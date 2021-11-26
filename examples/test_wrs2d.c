
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
  
  int npts = 1000;
  
  cprin_skipline(1);
  cprinf("npts = ", &npts, 1);

  //
  // generate a curve and build a tree on them
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

  opts.maxlev = 30;
  opts.maxboxes = maxboxes;  
  opts.nmax = 20;
  opts.maxwidth = 10;

  cprinf("opts.maxlev = ", &opts.maxlev, 1);
  cprinf("opts.maxboxes = ", &opts.maxboxes, 1);
  cprinf("opts.nmax = ", &opts.nmax, 1);
  cprind("opts.maxwidth = ", &opts.maxwidth, 1);


  
  int nlev, nboxes;
  int *perm = malloc(npts*sizeof(int));
  double width, center[2];
  
  quadtree_get_extent(npts, xys, center, &width);
  quadtree_build_lr(center, width, npts, xys, perm, opts, &nlev, &nboxes, tree);
  quadtree_print_tree(nboxes, tree);
  quadtree_plotboxes("2dtree", nboxes, tree, "Level restricted tree");

  //void (*fkernel)() = &bf2d_kernel_h0;

  
  

  return 0;
}