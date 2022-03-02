//
// 
//

#include <string.h>
#include <stdio.h>
#include "cplot.h"

void cplot2(int n1, double *xys1, char *label1, 
    int n2, double *xys2, char *label2, char *filename) {
  //
  // This plots the two sets of points in R^2.
  //

  char buffer[100];
  FILE *fp;

  strcpy(buffer, filename);
  strcat(buffer, ".py");
  fp = fopen(buffer, "w");

  fprintf(fp, "import matplotlib.pyplot as plt\n");
  fprintf(fp, "import numpy as np\n");
  fprintf(fp, "\n");
  fprintf(fp, "fig, ax = plt.subplots()\n");
  fprintf(fp, "\n");

  // plot the first sequence
  fprintf(fp, "n1 = %d\n", n1);
  fprintf(fp, "x1 = np.zeros(n1)\n");
  fprintf(fp, "y1 = np.zeros(n1)\n");

  int i;  
  for (i=0; i<n1; i++) {
    fprintf(fp, "x1[%d] = %e\n", i, xys1[2*i]);
    fprintf(fp, "y1[%d] = %e\n", i, xys1[2*i+1]);
  }

  fprintf(fp, "plt.scatter(x1, y1, s=2, c='blue', alpha=0.5, label='%s')\n",
    label1);
  fprintf(fp, "\n");

  // plot the second sequence
  fprintf(fp, "n2 = %d\n", n2);
  fprintf(fp, "x2 = np.zeros(n2)\n");
  fprintf(fp, "y2 = np.zeros(n2)\n");

  for (i=0; i<n2; i++) {
    fprintf(fp, "x2[%d] = %e\n", i, xys2[2*i]);
    fprintf(fp, "y2[%d] = %e\n", i, xys2[2*i+1]);
  }

  fprintf(fp, "plt.scatter(x2, y2, s=2, c='red', alpha=0.5, label='%s')\n",
    label2);
  fprintf(fp, "\n");


  fprintf(fp, "plt.axis('equal')\n");
  fprintf(fp, "plt.axis('on')\n");
  fprintf(fp, "ax.legend()\n");
  fprintf(fp, "plt.tight_layout()\n");
  fprintf(fp, "plt.show()\n");

  fclose(fp);


  return;
}