
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "cprini.h"
#include "quadtree.h"

int main (int argc, char* argv[])
{
                                                                                

  double pi;
  pi = 4*atan(1.0);

  cprin_init("stdout","fort.13");  
  
  int npts = 5000;
  int nmax = 50;

  cprin_skipline(1);
  cprinf("npts = ", &npts, 1);

  //
  // generate some random points and build a tree on them
  //
  srand(0);
  double *xys = malloc(2*npts*sizeof(double));

  quadtree_rand(2*npts, xys);

  int i;
  double t, theta, x, y;
  theta = pi/4;
  for (i=0; i<npts; i++) {
    t = (xys[2*i]+xys[2*i+1])*pi;
    x = 2*cos(t);
    y = sin(t);
    xys[2*i] = x*cos(theta) - y*sin(theta);
    xys[2*i+1] = x*sin(theta) + y*cos(theta);
    //(xys[2*i]-.5)*2;
    //xys[2*i+1] = (xys[2*i+1]-.5)*2;
  }

  // add some points to test level-restricted-ness
  quadtree_rand(2*nmax, xys);
  for (i=0; i<nmax; i++) {
    xys[2*i] = -xys[2*i]*0.5;
    xys[2*i+1] = xys[2*i+1]*0.5;
  }

  int nmax2 = 1000;
  
  quadtree_rand(2*nmax2, &(xys[2*nmax]));
  for (i=0; i<nmax2; i++) {
    xys[2*nmax + 2*i] = xys[2*nmax + 2*i]*0.5;
    xys[2*nmax + 2*i+1] = xys[2*nmax + 2*i+1]*0.5;
  }

  
  

  cprind("xys = ", xys, 20);

  
  //
  // now build a tree on the sources and targets
  //
  int maxboxes = pow(2,20);
  struct quadtree_box *tree = malloc(maxboxes*sizeof(struct quadtree_box));
  struct quadtree_opts opts;

  //int nmax = 20;
  //if (nmax > (m+n)) opts.maxlev = 0;
  //opts.maxlev = (int) (log((double)((m+n)/2)/nmax)/log(4.0))+1;
  //opts.maxlev = opts.maxlev + 2;
  //opts.maxlev = 4 + nfac;

  opts.maxlev = 40;
  opts.maxboxes = maxboxes;  
  opts.nmax = nmax;
  opts.maxwidth = 10;

  cprinf("opts.maxlev = ", &opts.maxlev, 1);
  cprinf("opts.maxboxes = ", &opts.maxboxes, 1);
  cprinf("opts.nmax = ", &opts.nmax, 1);
  cprind("opts.maxwidth = ", &opts.maxwidth, 1);


  
  int nlev, nboxes;
  int *perm = malloc(npts*sizeof(int));
  double width, center[2];
  
  quadtree_get_extent(npts, xys, center, &width);
  quadtree_build(center, width, npts, xys, perm, opts, &nlev, &nboxes, tree);
  quadtree_print_tree(nboxes, tree);
  quadtree_plotboxes("2dtree", nboxes, tree, "Fully adaptive boxes");
  quadtree_plotleaves("2dtree_leaves", nboxes, tree, "Fully adaptive leaves");



  // now fix up the tree into one that is level-restricted
  quadtee_fixtree_lr(&nlev, &nboxes, tree);
  quadtree_print_tree(nboxes, tree);  
  quadtree_plotboxes("2dtree_lr", nboxes, tree, "LR boxes");
  quadtree_plotleaves("2dtree_leaves_lr", nboxes, tree, "LR leaves");

  // re-compute the colleagues since there are new boxes
  quadtree_gen_colleagues(nlev, nboxes, tree);



  // now check that the list generation is working, generate list1 for
  // a box
  int nlist1, ibox;
  ibox = 169;
  quadtree_plotcolleagues("coll1", nboxes, tree,
          &(tree[ibox]), "Colleagues for test box");

  quadtree_printbox("sample box info", &tree[ibox]);
  struct quadtree_box *list1[1000];
  quadtree_get_list1(&tree[ibox], &nlist1, list1);

//  quadtree_plot_list1("list1", nboxes, tree,
//          &(tree[ibox]), "List 1 for test box");



  exit(0);





  return 0;
}
