

#include <math.h>


void l2d_single(double par0, double *src, double *targ, double *par1,
                double *par2, double *val) {

  // this routine computes just the green's function of the laplace
  // equation in 2D

  const double fac = 0.159154943091895335;
  
  double dx, dy, r;
  dx = src[0] - targ[0];
  dy = src[1] - targ[1];
  r = sqrt(dx*dx + dy*dy);

  *val = fac*log(r);

  return;
}

