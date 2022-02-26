// 
// This code is responsible for building and operating quadtree routines
// for use with fast direct solvers.
// 
// Jan 01, 2022 -- 
// Mike O'Neil
// moneil@flatironinstitute.org
//

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "cprini.h"
#include "quadtree.h"



void quadtree_get_list1(struct quadtree_box *box, int *nlist1,
    struct quadtree_box **list1) {
  //
  // This routine computes, on the fly, LIST1 for a box which is defined to be
  // the self box, plus all adjacent neighbors. IT IS ASSUMED that the tree
  // in which box resides is a level-restricted one, meaning that adjacent leaf
  // boxes are only up or down one level at most.
  //

  // initialize
  *nlist1 = 0;
  int isleaf;
  quadtree_isleaf(box, &isleaf);
  if (isleaf == 0) return;

  // get colleague info for box
  int i, j, ncoll, nl1;
  struct quadtree_box *coll, *child;
  double width;
  double dx, dy, dist;

  width = box->width;
  nl1 = 0;
  cprinf("finding list 1 for box = ", &(box->id), 1);
  cprind("box center = ", (box->center), 2);

  // First, check the colleagues: if any are leaves, they are in list1,
  // otherwise check their children
  ncoll = box->ncoll;
  cprinf("ncoll for box = ", &ncoll, 1);    
  for (i=0; i<ncoll; i++) {    
    coll = box->colls[i];
    quadtree_isleaf(coll, &isleaf);
    //cprinf("isleaf for colleague = ", &isleaf, 1);
    if (isleaf == 1) {
      list1[nl1] = coll;
      nl1 = nl1 + 1;
      cprinf("colleague id = ", &(coll->id), 1);
      cprind("colleague center = ", (coll->center), 2);
    }
    // check children of colleagues
    else {
      for (j=0; j<4; j++) {
        child = coll->child[j];
        // no need to check if it is a leaf, since by assumption the tree is
        // level restricted
        if (child == NULL) continue;
        dx = fabs(box->center[0] - child->center[0]);
        dy = fabs(box->center[1] - child->center[1]);
        dist = 0.75*width*1.000001;
        if ( (dx <= dist) && (dy <= dist) ) {
          list1[nl1] = child;
          nl1 = nl1 + 1;
        }


      }
    }
  }

  // lastly, check colleages of parent (using assumption that tree is
  // level-restricted)
  struct quadtree_box *parent;
  parent = box->parent;
  ncoll = parent->ncoll;
  cprinf("ncoll for parent = ", &ncoll, 1);    
  for (i=0; i<ncoll; i++) {    

    coll = parent->colls[i];
    quadtree_isleaf(coll, &isleaf);
    //cprinf("isleaf for parent colleague = ", &isleaf, 1);
    if (isleaf == 0) continue;

    // check the distance, if close enough it is in list 1
    dx = fabs(box->center[0] - coll->center[0]);
    dy = fabs(box->center[1] - coll->center[1]);
    dist = 1.5*width*1.000001;
    if ( (dx <= dist) && (dy <= dist) ) {
      list1[nl1] = coll;
      nl1 = nl1 + 1;
    }

  }

  *nlist1 = nl1;

  cprinf("found list 1 for box = ", &(box->id), 1);
  cprinf("nlist1 = ", nlist1, 1);

  return;
}





void quadtree_get_extent(int npts, double *xys, double *center, double *width) {

  // get the bounding box for the points in xys

  double xmin, xmax, ymin, ymax, x, y;
  xmin = xys[0];
  xmax = xys[0];
  ymin = xys[1];
  ymax = xys[1];;

  int i;
  for (i=0; i<npts; i++) {
    x = xys[2*i];
    y = xys[2*i+1];
    if (x > xmax) xmax = x;
    if (y > ymax) ymax = y;
    if (x < xmin) xmin = x;
    if (y < ymin) ymin = y;
  }

  center[0] = (xmin + xmax)/2;
  center[1] = (ymin + ymax)/2;

  double w = (xmax-xmin);
  if ((ymax-ymin) > w) {
    w = (ymax-ymin);
  }
  *width = w*(1.0 + 1.0e-8);

  return;
}
  



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





void quadtree_restrict1(int *nlev, int *nboxes, struct quadtree_box *tree) {
  
  // This routine makes a single pass through the tree (meaning the
  // leaves) and refines boxes whose neighbors are more than 1 level
  // down.

  int nleafs;
  int *leafs = malloc(*nboxes *sizeof(int));
  quadtree_getleafs(*nboxes, tree, &nleafs, leafs);

  //cprinf("number of leafs = ", &nleafs, 1);
  //cprinf("leafs = ", leafs, nleafs);

  // scan through the leafs
  int i, ilev, j, ifsplit;
  double dx, dy, center[2], width, sc;
  
  for (i=0; i<nleafs; i++) {

    ifsplit = 0;
    ilev = tree[leafs[i]].level;
    center[0] = tree[leafs[i]].center[0];
    center[1] = tree[leafs[i]].center[1];
    width = tree[leafs[i]].width;

    sc = .75*0.999999;
    for (j = 0; j<nleafs; j++) {

      if (tree[leafs[j]].level > ilev+1) {
        dx = fabs(tree[leafs[j]].center[0] - center[0]);
        dy = fabs(tree[leafs[j]].center[1] - center[1]);

        if ((dx < sc*width) && (dy < sc*width)) {
          ifsplit = 1;
        }
      }
      if (ifsplit == 1) break;
    }

    if (ifsplit == 1) {
      quadtree_split(tree[leafs[i]].id, nboxes, tree);
      return;
    }
    
  }

  return;
}





void quadtee_fixtree_lr(int *nlev, int *nboxes, struct quadtree_box *tree) {
  //
  // This routine takes an arbitrarily adaptive tree and fixes it so that it is
  // level-restricted.
  //

  int i, nprev, nnew;
  for (i=0; i<1000; i++) {
    cprinf("fixing tree, iteration i = ", &i, 1);
    nprev = *nboxes;
    //cprinf("before quadtree_restrict1, nboxes = ", &nprev, 1);
    quadtree_restrict1(nlev, nboxes, tree);
    //cprinf("after quadtree_restrict1, nboxes = ", &nboxes, 1);
    nnew = *nboxes - nprev;
    cprinf("new boxes created = ", &nnew, 1);
    if (nnew == 0) break;
    cprin_skipline(1);
  }

if (nnew != 0) {
    cprinf("level restricted tree not built, nnew = ", &nnew, 1);
    cprin_message(". . . exiting code");
    exit(0);
  }
  return;
}




void quadtree_build_lr(double *center, double width,
                    int npts, double *xys, int *perm, 
                  struct quadtree_opts opts, int *nlev,
                    int *nboxes, struct quadtree_box *tree) {
  //
  // this routine builds a level restricted tree by iteratively fixing
  // a fully adaptive one (not a fast algorithm)
  //

  quadtree_build(center, width, npts, xys, perm, opts, nlev, nboxes, tree);
  
  int i, nprev, nnew;
  for (i=0; i<1000; i++) {
    cprinf("fixing tree, iteration i = ", &i, 1);
    nprev = *nboxes;
    //cprinf("before quadtree_restrict1, nboxes = ", &nprev, 1);
    quadtree_restrict1(nlev, nboxes, tree);
    //cprinf("after quadtree_restrict1, nboxes = ", &nboxes, 1);
    nnew = *nboxes - nprev;
    cprinf("new boxes created = ", &nnew, 1);
    if (nnew == 0) break;
    cprin_skipline(1);
  }

  if (nnew != 0) {
    cprinf("level restricted tree not built, nnew = ", &nnew, 1);
    cprin_message(". . . exiting code");
    exit(0);
  }

  return;
}






void quadtree_build(double *center, double width,
                    int npts, double *xys, int *perm, 
                  struct quadtree_opts opts, int *nlev,
                    int *nboxes, struct quadtree_box *tree) {
  //
  // This routine builds an UNPRUNED space-filling fully adaptive quadtree on
  // the points xys, and will initially be used for a basic fast direct solver.
  //
  //     Input:
  //        center - 
  //        width - 
  //        npts - 
  //
  //     Output:
  //      nlev - max depth in the tree, 0, 1, ..., nlev
  //      nboxes - total number of boxes created
  //      tree - an array of boxes which is the tree
  //
  
  // initialize the permutation vector
  int i;
  for (i=0; i<npts; i++) {
    perm[i]=i;
  }
  
  //
  // construct the inital box
  //
  *nboxes = 0;
  *nlev = -1;
  if (npts <= 0) {
    cprin_message("exiting tree build, no points");
    return;
  }
  
  tree[0].id = *nboxes;
  tree[0].level = 0;
  tree[0].center[0] = center[0];
  tree[0].center[1] = center[1];
  tree[0].width = width;

  tree[0].npts = npts;
  tree[0].xys = xys;
  tree[0].perm = perm;

  tree[0].parent = NULL;

  for (i=0; i<4; ++i) {
    tree[0].child[i] = NULL;
  }

  *nboxes = 1;
  
  //
  // now run the subdivisions
  //
  int lev, ifsplit, ntemp, ifdone, nnn, ifchild;

  *nlev = 0;
  int maxboxes;
  maxboxes = opts.maxboxes;
  
  for (lev=0; lev<opts.maxlev; lev++) {

    ntemp = *nboxes;
    ifdone = 1;
    for (nnn=0; nnn<ntemp; nnn++) {

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
        if ((tree[nnn].npts) > opts.nmax) ifsplit=1;
        if (tree[nnn].width > opts.maxwidth) ifsplit=1;
      }

      if (ifsplit == 1) {
        if (*nboxes+4 > maxboxes) {
          cprin_message("on split, maxboxes will be exceeded");
          cprinf("nboxes = ", nboxes, 1);
          cprin_message("... exiting code");
          exit(0);
        }
        ifdone = 0;
        quadtree_split(nnn, nboxes, tree);        
      }
    }

    if (ifdone == 1) break;
    *nlev = *nlev+1;
  }

  //
  // now scan through all the boxes, level by level, and find the colleagues of
  // all possible boxes
  quadtree_gen_colleagues(*nlev, *nboxes, tree);

  int ifplot;
  ifplot = 0;

  // plot the colleagues of a particular box
  if (ifplot == 1) {  
    int ibox;
    char buffer[100], snum[16];
    for (ibox = 0; ibox<*nboxes; ibox ++){ 
      strcpy(buffer, "colltest_");
      strcpy(snum, "");
      sprintf(snum, "%d", ibox);
      strcat(buffer, snum);
      quadtree_plotcolleagues(buffer, *nboxes, tree,
          &(tree[ibox]), "Example colleague plot");
    }
  }

  return;
}





void quadtree_gen_colleagues(int nlev, int nboxes, struct quadtree_box *tree) {
  //
  // This routine initializes and generates all colleagues for all boxes. A box
  // A is a colleague of a box B if A and B are on the same level, and A and B
  // share a vertex.
  //
  // NOTE: This definition of colleagues means that A is a colleague of itself.
  //

  int lev, ibox;
  for (lev=0; lev < nlev+1; lev++) {
    for (ibox=0; ibox<nboxes; ibox++) {
      if (tree[ibox].level == lev) {
        quadtree_colleagues1( &tree[ibox] );
      }
    }
  }

  return;
}




void quadtree_colleagues1( struct quadtree_box *box ) {
  //
  // This routine finds the colleagues of a box by searching through children of
  // its parent and its parent's colleagues. The box struct is updated directly.
  //
  // Input:
  //  box - the box whose colleagues should be calculated
  //
  // Output:
  //  box - the structure is updated directly
  //

  int i;
  struct quadtree_box *parent;
  struct quadtree_box *pcoll, *pch;
  int ncoll, npcoll;
  int j;
  double center[2], dx, dy, width;

  // initialize all colleagues back to zero
  box->ncoll = 0;
  for (i=0; i<9; i++) {
    box->colls[i] = NULL;
  }

  // If the box is the root box, just set its colleagues to itself
  if (box->level == 0) {
    box->ncoll = 1;
    box->colls[0] = box;
    return;
  }

  // Otherwise search through children of the parent's colleages
  parent = box->parent;

  // scan through colleagues of parent and add children that are within a
  // certain distance

  ncoll = 0;
  npcoll = parent->ncoll;
  center[0] = box->center[0];
  center[1] = box->center[1];
  width = box->width;

  for (i=0; i<npcoll; i++) {

    pcoll = parent->colls[i];
    // scan through children of pcoll
    for (j=0; j<4; j++) {
      pch = pcoll->child[j];
      if (pch != NULL) {
        dx = pch->center[0] - center[0];
        dy = pch->center[1] - center[1];
        if ( (fabs(dx) < width*1.000001) && (fabs(dy) < width*1.000001) ) {
               // store pch as a colleague
               box->colls[ncoll] = pch;
               ncoll = ncoll + 1;
        }
      }

    }

  }

  box->ncoll = ncoll;

  // print out the colleague centers for the box
  int ifplot;
  ifplot = 0;
  if (ifplot == 1) {
    //printf("id = %d, ncoll = %d\n", box->id, ncoll);
    //cprinf("box id = ", &(box->id), 1);
    //cprinf("ncoll for this box = ", &ncoll, 1);
  }

  return;
}






void quadtree_plotcolleagues(char *filename, int nboxes,
    struct quadtree_box *tree, struct quadtree_box *box, char *title) {
  //
  // Plot all boxes and highlight box and its colleages.
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


  // plot the colleagues of box now
  fprintf(fp, "patches = []\n");

  int ncoll;
  ncoll = box->ncoll;
  for (i=0; i<ncoll; i++) {
    w = box->colls[i]->width;
    xc = box->colls[i]->center[0] - w/2;
    yc = box->colls[i]->center[1] - w/2;
    fprintf(fp, "\n");
    fprintf(fp, "rect = mpl.patches.Rectangle([%e,%e], %e, %e)\n", xc, yc, w, w);
    fprintf(fp, "patches.append(rect)\n");
  }
  fprintf(fp, "print('. . . highlighting colleagues')\n");  
  fprintf(fp, "collection = mpl.collections.PatchCollection(patches, alpha=0.2, facecolor='red', edgecolor='black')\n");
  fprintf(fp, "ax.add_collection(collection)\n");

  // highlight the self block
  fprintf(fp, "\n");
  fprintf(fp, "patches = []\n");
  w = box->width;
  xc = box->center[0] - w/2;
  yc = box->center[1] - w/2;
  fprintf(fp, "rect = mpl.patches.Rectangle([%e,%e], %e, %e)\n", xc, yc, w, w);
  fprintf(fp, "patches.append(rect)\n");
  fprintf(fp, "collection = mpl.collections.PatchCollection(patches, alpha=0.6, facecolor='red', edgecolor='black')\n");
  fprintf(fp, "ax.add_collection(collection)\n");



  // plot the sources now
  int npts = tree[0].npts;
  double *xys = tree[0].xys;

  fprintf(fp, "print('. . . plotting sources')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", npts);
  fprintf(fp, "ys = np.zeros(%d)\n", npts);

  for (i=0; i<npts; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, xys[2*i]);
  }

  for (i=0; i<npts; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, xys[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='blue', alpha=0.4)\n");

  fprintf(fp, "plt.axis('equal')\n");
  fprintf(fp, "plt.axis('on')\n");
  fprintf(fp, "plt.tight_layout()\n");
  fprintf(fp, "print('. . . showing plot')\n");  
  fprintf(fp, "plt.show()\n");

  fclose(fp);

  return;
}





void quadtree_plot_list1(char *filename, int nboxes,
    struct quadtree_box *tree, struct quadtree_box *box, 
    int nlist1, struct quadtree_box **list1, char *title) {
  //
  // Plot all boxes and highlight box and its colleages.
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


  // plot the list1 of box now
  fprintf(fp, "patches = []\n");

  for (i=0; i<nlist1; i++) {
    w = list1[i]->width;
    xc = list1[i]->center[0] - w/2;
    yc = list1[i]->center[1] - w/2;
    fprintf(fp, "\n");
    fprintf(fp, "rect = mpl.patches.Rectangle([%e,%e], %e, %e)\n", xc, yc, w, w);
    fprintf(fp, "patches.append(rect)\n");
  }
  fprintf(fp, "print('. . . highlighting list 1')\n");  
  fprintf(fp, "collection = mpl.collections.PatchCollection(patches, alpha=0.2, facecolor='green', edgecolor='black')\n");
  fprintf(fp, "ax.add_collection(collection)\n");

  // highlight the self block
  fprintf(fp, "\n");
  fprintf(fp, "patches = []\n");
  w = box->width;
  xc = box->center[0] - w/2;
  yc = box->center[1] - w/2;
  fprintf(fp, "rect = mpl.patches.Rectangle([%e,%e], %e, %e)\n", xc, yc, w, w);
  fprintf(fp, "patches.append(rect)\n");
  fprintf(fp, "collection = mpl.collections.PatchCollection(patches, alpha=0.6, facecolor='green', edgecolor='black')\n");
  fprintf(fp, "ax.add_collection(collection)\n");



  // plot the sources now
  int npts = tree[0].npts;
  double *xys = tree[0].xys;

  fprintf(fp, "print('. . . plotting sources')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", npts);
  fprintf(fp, "ys = np.zeros(%d)\n", npts);

  for (i=0; i<npts; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, xys[2*i]);
  }

  for (i=0; i<npts; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, xys[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='blue', alpha=0.4)\n");

  fprintf(fp, "plt.axis('equal')\n");
  fprintf(fp, "plt.axis('on')\n");
  fprintf(fp, "plt.tight_layout()\n");
  fprintf(fp, "print('. . . showing plot')\n");  
  fprintf(fp, "plt.show()\n");

  fclose(fp);

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
  int npts = tree[0].npts;
  double *xys = tree[0].xys;

  fprintf(fp, "print('. . . plotting sources')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", npts);
  fprintf(fp, "ys = np.zeros(%d)\n", npts);

  for (i=0; i<npts; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, xys[2*i]);
  }

  for (i=0; i<npts; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, xys[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='blue', alpha=0.4)\n");

  fprintf(fp, "plt.axis('equal')\n");
  fprintf(fp, "plt.axis('on')\n");
  fprintf(fp, "plt.tight_layout()\n");
  fprintf(fp, "print('. . . showing plot')\n");  
  fprintf(fp, "plt.show()\n");

  fclose(fp);

  return;
}




void quadtree_plotleaves(char *filename, int nboxes, struct quadtree_box *tree,
                      char *title) {
  //
  // Plot all the LEAF boxes and points.
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

  int i, isleaf;
  double xc, yc, w;
  for (i=0; i<nboxes; i++) {
    quadtree_isleaf(&tree[i], &isleaf);
    if (isleaf == 0) continue;
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
  int npts = tree[0].npts;
  double *xys = tree[0].xys;

  fprintf(fp, "print('. . . plotting sources')\n");
  fprintf(fp, "xs = np.zeros(%d)\n", npts);
  fprintf(fp, "ys = np.zeros(%d)\n", npts);

  for (i=0; i<npts; i++) {
    fprintf(fp, "xs[%d] = %e\n", i, xys[2*i]);
  }

  for (i=0; i<npts; i++) {
    fprintf(fp, "ys[%d] = %e\n", i, xys[2*i+1]);
  }

  fprintf(fp, "plt.scatter(xs, ys, s=5, c='blue', alpha=0.4)\n");

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
  // Split a single box into four children based on the total number
  // of sources and targets in the box. Note: This routine does NOT prune
  // empty boxes. The resulting tree is space-filling leafs, some of which
  // are empty.
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

  int npts = tree[ibox].npts;
  double *xys = tree[ibox].xys;
  int *perm = tree[ibox].perm;
  
  //
  // first collect the nodes
  //
  int is[4], ns[4];
  int i;
  for (i=0; i<4; ++i) {
    is[i] = 0;
    ns[i] = 0;
  }

  //
  // do sorting, based on value of box, 0, 1, 2, 3
  // first count the zeros
  //
  int ikey;

  for (i=0; i<npts; i++) {
    quadtree_key(tree[ibox].center, &xys[2*i], &ikey);
    ns[ikey] = ns[ikey] + 1;
  }

  is[0] = 0;
  is[1] = ns[0];
  is[2] = ns[0]+ns[1];
  is[3] = ns[0]+ns[1]+ns[2];
  
  //
  // add the new boxes to the tree, setting up correct pointers now
  //
  int nb = *nboxes;
  int j;
  
  for (j=0; j<4; ++j) {
    
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

      tree[nb].npts = ns[j];
      tree[nb].xys = &(xys[2*is[j]]);
      tree[nb].perm = &perm[is[j]];
      
      tree[nb].parent = &(tree[ibox]);

      for (i=0; i<4; ++i) {
        tree[nb].child[i] = NULL;
      }
      tree[ibox].child[j] = &(tree[nb]);

      nb = nb+1;
      *nboxes = nb;
  }

  //
  // now collect all the points that are in each box, permuting them
  // as need be
  //
  int ind;

  for (j=0; j<4; ++j) {
    ind = is[j];

    for (i=is[j]; i<npts; ++i) {
      quadtree_key(tree[ibox].center, &xys[2*i], &ikey);
      if (ikey == j) {
        quadtree_swapd(2, &xys[2*i], &xys[2*ind]);
        quadtree_swapi(1, &perm[i], &perm[ind]);
        ind=ind+1;
      }
    }
    
  }
  
  return;
}






void quadtree_key(double *center, double *xy, int *ikey) {
  //
  // get the quadrant that xy is in
  //
  double x = xy[0];
  double y = xy[1];
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
  printf("npts =        %d\n", box->npts);
  printf("xys =         %p\n", box->xys);
  printf("xys perm =    %p\n", box->perm);
  printf("parent =      %p\n", box->parent);
  printf("child[0] =    %p\n", box->child[0]);
  printf("child[1] =    %p\n", box->child[1]);
  printf("child[2] =    %p\n", box->child[2]);
  printf("ncoll =       %d\n", box->ncoll);
  printf("colls[0] =    %p\n", box->colls[0]);
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

  printf("   id  level    width        center          npts   parent   ch0   ch1   ch2   ch3\n");

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
    
    printf("%5d %5d  %8.2e (%+8.2e, %+8.2e) %6d %5d %5d %5d %5d %5d\n", tree[i].id,
           tree[i].level,
           tree[i].width, tree[i].center[0], tree[i].center[1], tree[i].npts,
           ip, ic0, ic1, ic2, ic3);

  }
  printf("- - - - - - - - - - - - - - - end of the tree - - - - - - - - - - - - - - -\n");
  

  return;
}
