

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cprini.h"
#include "quadtree.h"
#include "wrs2d.h"

//
// These subroutines build a weak recursive skeletonization-based fast
// direct solver
//

void wrs2d_factor() {



  
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








