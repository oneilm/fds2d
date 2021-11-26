

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cprini.h"
#include "quadtree.h"
#include "wrs2d.h"
#include "kernels.h"


//
// These subroutines build a weak recursive skeletonization-based fast
// direct solver
//

void wrs2d_factor(int nboxes, struct quadtree_box *tree,
                  double diag, void (*fkernel)() ) {

  // this computes all the factors needed for a fds based on the info
  // in tree and the kernel above

  int i, j, maxlevel;
  int *processed = malloc(nboxes*sizeof(int));

  maxlevel = 0;
  for (i=0; i<nboxes; i++) {
    processed[i] = 0;
    if (tree[i].level > maxlevel) maxlevel = tree[i].level;
  }

  int lev, ibox, ifproc;
  struct quadtree_box *box;
  
  for (lev=0; lev<=maxlevel; lev++) {

    cprinf("pass through tree lev = ", &lev, 1);

    for (ibox=0; ibox<nboxes; ibox++) {

      // skip if already processed
      if (processed[ibox] != 0) continue;
      
      // if ibox is a leaf, process now
      box = &tree[ibox];
      ifproc = 0;


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








