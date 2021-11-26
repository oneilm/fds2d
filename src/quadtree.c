

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cprini.h"
#include "quadtree.h"





void quadtree_unpermute(int m, int n, int *p, double *xs) {
  //
  // unpermute the array xs in place
  //
  // Input:
  //   m,n - dimensions of xs
  //   p - the permutation
  //   xs - the data to permute
  //
  double *tmp = malloc(m*n*sizeof(double));

  cprinf("inside unpermute, m = ", &m, 1);
  cprinf("inside unpermute, n = ", &n, 1);


  //cprind("unsorted vals xs = ", xs, 10);
  

  int i,j;


  //for (i=0; i<n; ++i) {
  //  tmp[i] = xs[i];
  // }

  //for (i=0; i<n; ++i) {
  //  xs[p[i]] = tmp[i];
  // }

  //  return;


  for (i=0; i<m*n; ++i) {
    tmp[i] = xs[i];
  }

  for (j=0; j<n; ++j) {
    for (i=0; i<m; ++i) {
      xs[i+m*p[j]] = tmp[i+m*j];
    }
  }
  
  return;
}




void quadtree_rand(int n, double *rs) {

  int i;
  for (i=0; i<n; i++) {
    rs[i] = (double)rand() / RAND_MAX;
  }
  return;
}





void quadtree_getboxes_bylevel(int nboxes, struct quadtree_box *tree,
                             int lev, int *nb, int *boxes) {
  //
  // return all the boxes on a particular level
  //
  int i;
  *nb = 0;
  for (i=0; i<nboxes; ++i) {
    if (tree[i].level == lev) {
      boxes[*nb] = i;
      *nb = *nb + 1;
    }
  }

  return;
}





void quadtree_isleaf(struct quadtree_box *box, int *isleaf) {
  //
  // check whether box has any children
  //
  int i;
  for (i=0; i<4; ++i) {
    if (box->child[i] != NULL) {
      *isleaf = 0;
      return;
    }
  }

  *isleaf = 1;
  return;
}





void quadtree_getsiblings(struct quadtree_box *box, struct quadtree_box **pch) {
  //
  // return pointers to all the siblings of box, which will include box itself
  //

  struct quadtree_box *parent;

  parent = box->parent;


  int i;
  for (i=0; i<4; ++i) {
    pch[i] = parent->child[i];
  }
  
  return;
}




void quadtree_getleafs(int nboxes, struct quadtree_box *tree, int *nleafs,
                     int *leafs) {
  //
  // search for all the leafs
  //

  int i;
  *nleafs = 0;
  void *ptr0, *ptr1, *ptr2, *ptr3;
  for (i=0; i<nboxes; i++) {
    ptr0 = tree[i].child[0];
    ptr1 = tree[i].child[1];
    ptr2 = tree[i].child[2];
    ptr3 = tree[i].child[3];
    if (ptr0 != NULL) continue;
    if (ptr1 != NULL) continue;
    if (ptr2 != NULL) continue;
    if (ptr3 != NULL) continue;
    leafs[*nleafs] = tree[i].id;
    *nleafs = *nleafs +1;
  }

  
  return;
}




void quadtree_build(int m, double *targs, int *pt, int n, double *srcs, int *ps,
                  struct quadtree_opts opts, int *nlev,
                  int *nboxes, struct quadtree_box *tree) {
  //
  // this routine builds a uniform pruned tree on srcs and targs
  // according to the parameters in opts
  //
  //     input:
  // m - number of targets
  // targs - the targets in 2d
  // n - number of sources
  // srcs - the source points in 2d
  // opts - tree options, structure for futureproofing
  //
  //     output:
  // targs - on output, have been re-ordered
  // pt - permutation vector re-ordering targets to be contiguous
  // srcs - on output, have been re-ordered
  // ps - permutation vector re-ordering sources to be contiguous
  // nlev - the number of levels in the tree, 0 means just root
  // nboxes - total number of boxes created
  // tree - array of box structs, linked together by pointers
  //
  //
  
  // initialize the permutation vector
  int i;
  for (i=0; i<m; i++) {
    pt[i]=i;
  }

  for (i=0; i<n; i++) {
    ps[i]=i;
  }

  
  //
  // find the bounding box
  //
  double xmin, xmax, ymin, ymax, x, y;
  xmin = srcs[0];
  xmax = srcs[0];
  ymin = srcs[1];
  ymax = srcs[1];;
  
  for (i=0; i<n; i++) {
    x = srcs[2*i];
    y = srcs[2*i+1];
    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (x < xmin) xmin = x;
    if (y < ymin) ymin = y;
    x = targs[2*i];
    y = targs[2*i+1];
    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (x < xmin) xmin = x;
    if (y < ymin) ymin = y;
  }

  
  //
  // construct the inital box
  //
  *nboxes = 0;
  tree[0].id = *nboxes;
  tree[0].level = 0;
  tree[0].center[0] = (xmin+xmax)/2;
  tree[0].center[1] = (ymin+ymax)/2;

  double w = (xmax-xmin);
  if ((ymax-ymin) > w) w = (ymax-ymin);
  w = w*(1.0 + 1.0e-8);
  tree[0].width = w;

  tree[0].nsrcs = n;
  tree[0].srcs = srcs;
  tree[0].ps = ps;

  tree[0].ntargs = m;
  tree[0].targs = targs;
  tree[0].pt = pt;

  tree[0].parent = NULL;

  for (i=0; i<4; ++i) {
    tree[0].child[i] = NULL;
  }
  
  *nboxes = 1;
  
  //
  // now run the subdivisions
  //
  int lev, ifsplit, ntemp, ifdone, nnn, ifchild;
  //int ifsrc, iftarg;
  

  *nlev = 0;
  int maxboxes, nmax;
  maxboxes = opts.maxboxes;
  nmax = opts.nmax;

  
  for (lev=0; lev<opts.maxlev; lev++) {

    //cprinf("beginning outer loop, lev = ", &lev, 1);
    //cprinf("beginning outer loop, nboxes = ", nboxes, 1);

    ntemp = *nboxes;
    ifdone = 1;
    for (nnn=0; nnn<ntemp; nnn++) {

      //cprin_skipline(1);
      //cprinf("processing box: ", &nnn, 1);
      
      ifsplit = 0;
      ifchild = 1;
      for (i=0; i<4; ++i) {
        if (tree[nnn].child[i] != NULL) ifchild=0;
      }
      
      //
      // change these options to reflect the correct subdivision
      // criteria
      //
      if (ifchild == 1) {
        ifsplit = 1;
        //if ((tree[nnn].nsrcs+tree[nnn].ntargs) > opts.nmax) ifsplit=1;
        //if (tree[nnn].width > opts.maxwidth) ifsplit=1;
      }

      if (ifsplit == 1) {
        if (*nboxes+4 >maxboxes) {
          cprin_message("on split, maxboxes will be exceeded");
          cprinf("nboxes = ", nboxes, 1);
          exit(0);
        }
        ifdone = 0;
        quadtree_split(nnn, nboxes, tree);
      }
    }

    if (ifdone == 1) break;
    *nlev = *nlev+1;
  }

  return;
}







void quadtree_plotboxes(char *filename, int nboxes, struct quadtree_box *tree,
                      char *title) {
  //
  // plot all the boxes and points
  //
  char buffer[100];
  FILE *fp;

  strcpy(buffer, filename);
  strcat(buffer, ".py");
  fp = fopen(buffer, "w");

  
  fprintf(fp, "import matplotlib.pyplot as plt\n");
  fprintf(fp, "import numpy as np\n");
  fprintf(fp, "import matplotlib as mpl\n");
  fprintf(fp, "\n");
  fprintf(fp, "fig, ax = plt.subplots()\n");
  fprintf(fp, "patches = []\n");
  
  fprintf(fp, "print('. . . constructing boxes')\n");  

  int i;
  double xc, yc, w;
  for (i=0; i<nboxes; i++) {
  //for (i=0; i<5; i++) {
    w = tree[i].width;
    xc = tree[i].center[0] - w/2;
    yc = tree[i].center[1] - w/2;
    fprintf(fp, "\n");
    fprintf(fp, "rect = mpl.patches.Rectangle([%e,%e], %e, %e)\n", xc, yc, w, w);
    fprintf(fp, "patches.append(rect)\n");
  }

  fprintf(fp, "print('. . . plotting boxes')\n");  
  fprintf(fp, "collection = mpl.collections.PatchCollection(patches, alpha=0.2, facecolor='grey', edgecolor='black')\n");
  fprintf(fp, "ax.add_collection(collection)\n");


  // plot the sources now
  int nsrcs = tree[0].nsrcs;
  int ntargs = tree[0].ntargs;
  double *srcs = tree[0].srcs;
  double *targs = tree[0].targs;

  fprintf(fp, "print('. . . plotting sources')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", nsrcs);
  fprintf(fp, "ys = np.zeros(%d)\n", nsrcs);

  for (i=0; i<nsrcs; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, srcs[2*i]);
  }

  for (i=0; i<nsrcs; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, srcs[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='blue', alpha=0.4)\n");

  
  fprintf(fp, "print('. . . plotting targets')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", ntargs);
  fprintf(fp, "ys = np.zeros(%d)\n", ntargs);
  for (i=0; i<ntargs; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, targs[2*i]);
  }

  for (i=0; i<ntargs; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, targs[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='red', alpha=0.4)\n");
  
  fprintf(fp, "plt.axis('equal')\n");
  fprintf(fp, "plt.axis('on')\n");
  fprintf(fp, "plt.tight_layout()\n");
  fprintf(fp, "print('. . . showing plot')\n");  
  fprintf(fp, "plt.show()\n");

  fclose(fp);

  return;
}





void quadtree_split(int ibox, int *nboxes, struct quadtree_box *tree) {
  //
  // split a single box into four children based on the total number
  // of sources and targets in the box
  //
  // Input:
  //   ibox - box number to split
  //   nboxes - total number of boxes in tree
  //   tree - array of boxes
  //
  // Output:
  //   nboxes - new total number of boxes
  //   tree - updated array of boxes
  //
  
  int lev = tree[ibox].level;
  double xc = tree[ibox].center[0];
  double yc = tree[ibox].center[1];
  double w = tree[ibox].width;


  int nsrcs = tree[ibox].nsrcs;
  double *srcs = tree[ibox].srcs;
  int *ps = tree[ibox].ps;

  int ntargs = tree[ibox].ntargs;
  double *targs = tree[ibox].targs;
  int *pt = tree[ibox].pt;

  
  //
  // first collect the sources
  //
  int is[4], its[4], ns[4], nts[4];
  int i;
  for (i=0; i<4; ++i) {
    is[i] = 0;
    ns[i] = 0;
    its[i] = 0;
    nts[i] = 0;
  }


  //
  // do sorting, based on value of box, 0, 1, 2, 3
  // first count the zeros
  //
  int ikey;

  for (i=0; i<nsrcs; i++) {
    quadtree_key(tree[ibox].center, &srcs[2*i], &ikey);
    ns[ikey] = ns[ikey] + 1;
  }

  for (i=0; i<ntargs; i++) {
    quadtree_key(tree[ibox].center, &targs[2*i], &ikey);
    nts[ikey] = nts[ikey] + 1;
  }

  is[0] = 0;
  is[1] = ns[0];
  is[2] = ns[0]+ns[1];
  is[3] = ns[0]+ns[1]+ns[2];
  
  its[0] = 0;
  its[1] = nts[0];
  its[2] = nts[0]+nts[1];
  its[3] = nts[0]+nts[1]+nts[2];
  
  //
  // add the new boxes to the tree, setting up correct pointers now
  //
  int nb = *nboxes;
  int j;
  
  for (j=0; j<4; ++j) {
    
    if ((ns[j]+nts[j] > 0) ) {
      
      tree[nb].id = nb;
      tree[nb].level = lev+1;
      tree[nb].width = w/2;

      if (j==0) {
        tree[nb].center[0] = xc-w/4;
        tree[nb].center[1] = yc-w/4;
      } else if (j==1) {
        tree[nb].center[0] = xc+w/4;
        tree[nb].center[1] = yc-w/4;
      } else if (j==2) {
        tree[nb].center[0] = xc-w/4;
        tree[nb].center[1] = yc+w/4;
      } else if (j==3) {
        tree[nb].center[0] = xc+w/4;
        tree[nb].center[1] = yc+w/4;
      }

      tree[nb].nsrcs = ns[j];
      tree[nb].srcs = &(srcs[2*is[j]]);
      tree[nb].ps = &ps[is[j]];
      
      tree[nb].ntargs = nts[j];
      tree[nb].targs = &(targs[2*its[j]]);
      tree[nb].pt = &pt[its[j]];
      
      tree[nb].parent = &(tree[ibox]);

      for (i=0; i<4; ++i) {
        tree[nb].child[i] = NULL;
      }
      tree[ibox].child[j] = &(tree[nb]);
      
      nb = nb+1;
      *nboxes = nb;
    }

  }


  //
  // now collect all the sources that are in each box, permuting them
  // as need be
  //
  int ind;

  for (j=0; j<4; ++j) {
    ind = is[j];

    for (i=is[j]; i<nsrcs; ++i) {
      quadtree_key(tree[ibox].center, &srcs[2*i], &ikey);
      if (ikey == j) {
        quadtree_swapd(2, &srcs[2*i], &srcs[2*ind]);
        quadtree_swapi(1, &ps[i], &ps[ind]);
        ind=ind+1;
      }
    }
    
  }
  

  //
  // and collect all the targets that are in each box, permuting them
  // as need be
  //
  for (j=0; j<4; ++j) {
    ind = its[j];

    for (i=its[j]; i<ntargs; ++i) {
      quadtree_key(tree[ibox].center, &targs[2*i], &ikey);
      if (ikey == j) {
        quadtree_swapd(2, &targs[2*i], &targs[2*ind]);
        quadtree_swapi(1, &pt[i], &pt[ind]);
        ind=ind+1;
      }
    }
    
  }
  
  return;
}





void quadtree_key(double *center, double *src, int *ikey) {
  //
  // get the quadrant that src is in
  double x = src[0];
  double y = src[1];
  double xc = center[0];
  double yc = center[1];

  *ikey = -1;
  if ((x <= xc) && (y <= yc)) *ikey = 0;
  if ((x > xc) && (y <= yc)) *ikey = 1;
  if ((x <= xc) && (y > yc)) *ikey = 2;
  if ((x > xc) && (y > yc)) *ikey = 3;

  return;
}




void quadtree_printbox(char *name, struct quadtree_box *box) {
  //
  // print out details of this box
  //

  printf("\n");
  printf("- - - %s - - -\n", name);
  printf("self =        %p\n", box);
  printf("id =          %d\n", box->id);
  printf("level =       %d\n", box->level);
  printf("width =       %e\n", box->width);
  printf("center =      (%+e, %+e)\n", box->center[0], box->center[1]);
  printf("nsrcs =       %d\n", box->nsrcs);
  printf("srcs =        %p\n", box->srcs);
  printf("src perm =    %p\n", box->ps);
  printf("ntargs =      %d\n", box->ntargs);
  printf("targs =       %p\n", box->targs);
  printf("targ perm =   %p\n", box->pt);
  printf("parent =      %p\n", box->parent);
  printf("child[0] =    %p\n", box->child[0]);
  printf("child[1] =    %p\n", box->child[1]);
  printf("child[2] =    %p\n", box->child[2]);
  printf("child[3] =    %p\n", box->child[3]);
  printf("\n");
  


  return;
}





void quadtree_swapd(int n, double *xs, double *ys) {
  //
  // swap n doubles from arrays xs and ys
  //
  int i;
  double temp;
  for (i=0; i<n; i++) {
    temp = xs[i];
    xs[i] = ys[i];
    ys[i] = temp;
  }
  return;
}





void quadtree_swapi(int n, int *xs, int *ys) {
  //
  // swap n ints from arrays xs and ys
  //
  int i;
  int temp;
  for (i=0; i<n; i++) {
    temp = xs[i];
    xs[i] = ys[i];
    ys[i] = temp;
  }
  return;
}





void quadtree_print_tree(int nboxes, struct quadtree_box *tree) {
  //
  // print out details of this box
  //

  printf("\n");
  printf("- - - - - - - - - - - - - - - start of the tree - - - - - - - - - - - - - - -\n");

  int i = 0;
  int ip, ic0, ic1, ic2, ic3;

  printf("   id  level    width        center          nsrcs  ntargs parent   ch0   ch1   ch2   ch3\n");

  for (i=0; i<nboxes; i++) {

    ip = -1;
    ic0 = -1;
    ic1 = -1;
    ic2 = -1;
    ic3 = -1;
    if (tree[i].parent != NULL) ip = tree[i].parent->id;
    if (tree[i].child[0] != NULL) ic0 = tree[i].child[0]->id;
    if (tree[i].child[1] != NULL) ic1 = tree[i].child[1]->id;
    if (tree[i].child[2] != NULL) ic2 = tree[i].child[2]->id;
    if (tree[i].child[3] != NULL) ic3 = tree[i].child[3]->id;
    
    printf("%5d %5d  %8.2e (%+8.2e, %+8.2e) %6d %6d  %5d %5d %5d %5d %5d\n", tree[i].id,
           tree[i].level,
           tree[i].width, tree[i].center[0], tree[i].center[1], tree[i].nsrcs, tree[i].ntargs, 
           ip, ic0, ic1, ic2, ic3);

  }
  printf("- - - - - - - - - - - - - - - end of the tree - - - - - - - - - - - - - - -\n");
  

  return;
}
