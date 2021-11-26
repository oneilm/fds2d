
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
  
  int m, n;

  n = 4000;
  //n = 100;
  m = n;
  
  cprin_skipline(1);
  cprinf("m = ", &m, 1);
  cprinf("n = ", &n, 1);

  //
  // generate some random points and build a tree on them
  //
  srand(0);
  double *srcs = malloc(2*n*sizeof(double));
  double *targs = malloc(2*m*sizeof(double));

  quadtree_rand(2*n, srcs);
  quadtree_rand(2*m, targs);
  int i;
  double t, theta, x, y;
  theta = pi/4;
  for (i=0; i<n; i++) {
    t = (srcs[2*i]+srcs[2*i+1])*pi;
    x = 2*cos(t);
    y = sin(t);
    srcs[2*i] = x*cos(theta) - y*sin(theta);
    srcs[2*i+1] = x*sin(theta) + y*cos(theta);
    //(srcs[2*i]-.5)*2;
    //srcs[2*i+1] = (srcs[2*i+1]-.5)*2;
  }

  for (i=0; i<m; i++) {
    targs[2*i] = (targs[2*i]-.5)*2*2 + 8;
    targs[2*i+1] = (targs[2*i+1]-.5)*2*2;
    //targs[2*i] = (targs[2*i]-.5)*2;
    //targs[2*i+1] = (targs[2*i+1]-.5)*2;
  }

  cprind("sources = ", srcs, 10);
  cprind("targets = ", targs, 10);

  
  //
  // now build a tree on the sources and targets
  //
  int maxboxes = pow(2,20);
  struct quadtree_box *tree = malloc(maxboxes*sizeof(struct quadtree_box));
  struct quadtree_opts opts;

  int nmax = 20;
  if (nmax > (m+n)) opts.maxlev = 0;
  opts.maxlev = (int) (log((double)((m+n)/2)/nmax)/log(4.0))+1;
  opts.maxlev = opts.maxlev + 2;
  opts.maxlev = 30;
  //opts.maxlev = 4 + nfac;
  //cprinf("maxlev = ", &opts.maxlev, 1);

  
  opts.maxboxes = maxboxes;  
  opts.nmax = nmax;
  opts.maxwidth = 10;

  int nlev, nboxes;
  int *ps = malloc(n*sizeof(int));
  int *pt = malloc(m*sizeof(int));

  
  quadtree_build(m, targs, pt, n, srcs, ps, opts, &nlev, &nboxes, tree);
  //  quadtree_printtree(nboxes, tree);
  quadtree_plotboxes("2dtree", nboxes, tree, "the 2d boxes");


  return 0;
}
