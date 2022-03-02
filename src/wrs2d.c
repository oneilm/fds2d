//
// Feb 27, 2022 -- 
// Mike O'Neil
// moneil@flatironinstitute.org
// Center for Computational Mathematics
// Flatiron Institute
// New York, NY

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cprini.h"
#include "quadtree.h"
#include "wrs2d.h"
#include "kernels.h"
#include "cplot.h"

//
// These subroutines build a weak recursive skeletonization-based fast
// direct solver
//

void wrs2d_factor(int nboxes, struct quadtree_box *tree,
                  double diag, void (*fkernel)() ) {

  // This computes all the factors needed for a FDS based on the info
  // in tree and the kernel above. 
  //
  // input:
  //   nboxes - the total number of boxes in tree
  //   tree - the list of boxes
  //   diag - diagonal to add to the matrix, a constant for now
  //   kernel - the kernel to compress... more on this later
  //
  // output:
  //   factors - ??? now sure how to output things yet...
  //
  

  int i, j, maxlevel;
  int *processed = malloc(nboxes*sizeof(int));

  int isleaf, nproxy, npts;
  double center[2], width, radius;

  nproxy = 20;
  double *pxys = malloc(2*nproxy*sizeof(double));
  double *xys;
  

  maxlevel = 0;
  for (i=0; i<nboxes; i++) {
    processed[i] = 0;
    if (tree[i].level > maxlevel) maxlevel = tree[i].level;
  }

  int lev, ibox, ifproc, nlist1, ntarg, nsrc, ntmp, k, ijk;
  struct quadtree_box *box, *list1[1000];

  
  for (lev=maxlevel; lev>=0; lev--) {

    cprinf("processing level lev = ", &lev, 1);

    for (ibox=0; ibox<nboxes; ibox++) {

      if (tree[ibox].level != lev) continue;

      // process this box on level lev
      box = &tree[ibox];
      quadtree_isleaf(box, &isleaf);

      if ( (isleaf == 1) && (box->npts != 0) ) {
        //cprinf("isleaf = ", &isleaf, 1);

        center[0] = box -> center[0];
        center[1] = box -> center[1];
        cprind("box center = ", center, 2);
        width = box -> width;
        
        // generate the proxy surface to compress
        radius = 1.75*width/2;
        cprind("box width = ", &width, 1);
        cprind("radius for proxy points = ", &radius, 1);

        wrs2d_proxy(center, nproxy, radius, pxys);
        cprind("proxy points = ", pxys, 2*nproxy);

        // collect the points in the box
        npts = box->npts;
        xys = box->xys;
        cprind("points in the box = ", xys, 2*npts);

        // collect the points in list 1
        quadtree_get_list1(box, &nlist1, list1);
        quadtree_plot_list1("list1-wrs", nboxes, tree, box, nlist1, list1,
            "List 1 for box in WRS");
        cprinf("boxes in list 1 = ", &nlist1, 1);

        // first construct the incoming skeletons
        nsrc = nproxy;
        for (j=0; j<nlist1; j++) {
          if (box != list1[j]) {
            nsrc = nsrc + list1[j]->npts;
          }
        }

        double *srcs = malloc(2*nsrc*sizeof(double));
        ijk = 0;
        for (j=0; j<nlist1; j++) {
          if (box != list1[j]) {
            ntmp = list1[j]->npts;
            for (k=0; k<ntmp; k++) {
              srcs[2*ijk] = list1[j]->xys[2*k];
              srcs[2*ijk+1] = list1[j]->xys[2*k+1];
              ijk = ijk + 1;
            }
          }
        }

        for (k=0; k<nproxy; k++) {
          srcs[2*ijk] = pxys[2*k];
          srcs[2*ijk+1] = pxys[2*k+1];
          ijk = ijk + 1;
        }

        // plot these points to check
        cplot2(nsrc, srcs, "sources", npts, xys, "targets", "inc_skel");

        // construct the system matrix and skeletonize
        double *sysmat = malloc(npts*nsrc*sizeof(double));
        wrs2d_buildmat(npts, xys, nsrc, srcs, 


        free(srcs);

        exit(0);

        ntarg = ntarg + nproxy;
        double *targ = malloc(ntarg*sizeof(double));


        free(targ);


        exit(0);

        // for now, build interaction matrix of srcs/proxy points and
        // skeletonize

        exit(0);
      }
          
      
      // get box info
      
    }

    
    
  }


  
  
  return;
}






void wrs2d_proxy(double *center, int n, double radius, double *pxys) {

  // this routine generates a proxy surface centered at center and of
  // raidus radius with n points on it

  double one;
  one = 1.0;
  double pi = 4*atan(one);

  double h = 2*pi/n, t;
  int i;
  
  for (i = 0; i<n; i++) {
    t = h*i;
    pxys[2*i] = center[0] + radius*cos(t);
    pxys[2*i+1] = center[1] + radius*sin(t);
  }
  
  return;

}








