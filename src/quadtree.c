

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cprini.h"
#include "quadtree.h"






void quadtree_findall_list1(int nboxes, struct quadtree_box *tree) {

  // This subroutine passes through the tree and 
  // generates lists of all boxes in list1. A box J is in list1 of a
  // box I if it is adjacent.
  //
  // 

  int i, ibox;
  //int *ifproc = malloc(nboxes*sizeof(int));

  //for (i = 0; i < nboxes; i++) {
  //  ifproc[i] = 0;
  //}

  // box 0 has no neighbors
  //tree[0].nnbors = 0;
  
  struct quadtree_box *parent;

  if (nboxes == 0) return;

  // initialize the parent box
  ibox = 0;
  //  tree[ibox].nlist1 = 0;

  // generate lists for all other boxes

  int j, nlp;
  
  for (ibox = 1; ibox<nboxes; ibox++) {

    parent = tree[ibox].parent;
    //    nlp = parent.nlist1;
    
    // scan through list1 of the parent
    
    ///for (j=0; j<

  }

  
  
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

    sc = .75*1.000001;

    for (j = 0; j<nleafs; j++) {

      if (tree[leafs[j]].level > ilev+1) {
        dx = abs(tree[leafs[j]].center[0] - center[0]);
        dy = abs(tree[leafs[j]].center[1] - center[1]);

        if ((dx <= sc*width) && (dy <= sc*width)) {
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





void quadtree_build_lr(double *center, double width,
                    int npts, double *xys, int *perm, 
                  struct quadtree_opts opts, int *nlev,
                    int *nboxes, struct quadtree_box *tree) {

  // this routine builds a level restricted tree by iteratively fixing
  // a fully adaptive one (not a fast algorithm)
  
  quadtree_build(center, width, npts, xys, perm, opts, nlev, nboxes, tree);
  
  int i, nprev;
  for (i=0; i<100; i++) {
    cprinf("fixing adaptive tree, iteration i = ", &i, 1);
    nprev = *nboxes;
    cprinf("before quadtree_restrict1, nboxes = ", &nprev, 1);
    quadtree_restrict1(nlev, nboxes, tree);
    cprinf("after quadtree_restrict1, nboxes = ", nboxes, 1);
    if (nprev == *nboxes) break;
    cprin_skipline(1);
  }

  return;
}





void quadtree_build(double *center, double width,
                    int npts, double *xys, int *perm, 
                  struct quadtree_opts opts, int *nlev,
                    int *nboxes, struct quadtree_box *tree) {

  // This routine builds a stadard quadtree on the points xys, and
  // will originally be used for a basic fast direct solver. The tree
  // is pruned and fully adaptive.
  //
  //     Input:
  //        center - 
  //        width - 
  //        npts - 
  //
  //     Output:
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

  // set list1 for the root
  tree[0].ncoll = 0;

  *nboxes = 1;
  
  //
  // now run the subdivisions
  //
  int lev, ifsplit, ntemp, ifdone, nnn, ifchild;

  *nlev = 0;
  int maxboxes, nboxes0, nnew;
  //int nmax;

  maxboxes = opts.maxboxes;
  //nmax = opts.nmax;


  
  
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
        if ((tree[nnn].npts) > opts.nmax) ifsplit=1;
        if (tree[nnn].width > opts.maxwidth) ifsplit=1;
      }

      if (ifsplit == 1) {
        if (*nboxes+4 > maxboxes) {
          cprin_message("on split, maxboxes will be exceeded");
          cprinf("nboxes = ", nboxes, 1);
          exit(0);
        }
        ifdone = 0;
        //nboxes0 = *nboxes;
        quadtree_split(nnn, nboxes, tree);
        //nnew = *nboxes - nboxes0;

        // for each new box, calculate colleagues
        //cprinf("num new boxes = ", &nnew, 1);

        //for (i=0; i<nnew; i++) {
        //  quadtree_colleagues1( &tree[nboxes0+i] );
        // }
        
      }
    }

    if (ifdone == 1) break;
    *nlev = *nlev+1;
  }


  // now scan through all the boxes, level by level, and find the colleagues of
  // all possible boxes
  int ibox;
  for (lev=0; lev<3; lev++) {

    cprinf("computing colleagues for level lev = ", &lev, 1);

    for (ibox=0; ibox<*nboxes; ibox++) {
      
      if (tree[ibox].level == lev) {
        quadtree_colleagues1( &tree[ibox] );

      }


    }

  }

  exit(0);

  return;

}





void quadtree_colleagues1( struct quadtree_box *box ) {

  // This routine finds the colleagues of box by searching through children of
  // its parent and its parent's colleagues. The box struct is updated directly.
  //
  // Input:
  //   box - the box whose colleagues should be calculated
  //

  // initialize all colleagues back to zero
  box->ncoll = 0;
  int i;
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
  struct quadtree_box *parent;
  parent = box->parent;

  // scan through colleagues of parent and add children that are within a
  // certain distance

  int ncoll, npcoll;
  ncoll = 0;
  npcoll = parent->ncoll;
  
  int j;
  double center[2], dx, dy, width;
  center[0] = box->center[0];
  center[1] = box->center[1];
  width = box->width;

  struct quadtree_box *pcoll, *pch;
  for (i=0; i<npcoll; i++) {

    pcoll = parent->colls[i];
    // scan through children of pcoll
    for (j=0; j<4; j++) {
      pch = pcoll->child[j];
      if (pch != NULL) {
        dx = pch->center[0] - center[0];
        dy = pch->center[1] - center[1];
        if ( (abs(dx) < width*1.000001) && 
             (abs(dy) < width*1.000001) ) {
               // store pch as a colleague
               box->colls[ncoll] = pch;
               ncoll = ncoll + 1;
        }
      }

    }

  }

  // print out the colleague centers for the box
  box->ncoll = ncoll;

  cprinf("ncoll for this box = ", &ncoll, 1);

  //exit(0);

  
  
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
    
    if ( ns[j] > 0 ) {
      
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



  // generate a list of colleagues for each new box
  // ... first add the siblings

  struct quadtree_box *parent, *newbox, *sibling;
  struct quadtree_box **pcolls, *child;
  int npcoll;
  int ncoll, k;
  double width, *center;
  
  parent = &(tree[ibox]);
  cprin_skipline(2);
  cprinf("parent index, ibox = ", &ibox, 1);
  cprind("parent center = ", parent->center, 2);

  for (j=0; j<4; j++) {

    newbox = parent->child[j];

    ncoll = 0;
    if (newbox != NULL) {
      for (i=0; i<4; i++) {
        sibling = parent->child[i];
        if ( (i != j) && (sibling != NULL) ) {
          newbox->colls[ncoll] = sibling;
          ncoll = ncoll + 1;
        }
      }

      cprin_skipline(2);
      cprinf("number of colleagues from siblins, ncoll = ", &ncoll, 1);

      // and now check children of the colleagues of the parent
      if (parent->parent != NULL) {

        npcoll = parent->ncoll;
        pcolls = parent->colls;
        cprinf("npcoll = ", &npcoll, 1);

        width = newbox->width;
        center = newbox->center;

        cprind("center of newbox = ", center, 2);

        ///// CHILD HASN"T BEEN CREATED YET!!!!
        cprind("center of parent colleague 0 = ",
               pcolls[0]->child[0]->center, 2);
        exit(0);
        ///// CHILD HASN"T BEEN CREATED YET!!!!
        
        for (i=0; i<npcoll; i++) {

          // check each child of colleague i of the parent
          for (k=0; k<4; k++){
            if (pcolls[i]->child[k] != NULL) {
              xc = pcolls[i]->child[k]->center[0];
              yc = pcolls[i]->child[k]->center[1];
              cprind("testing box with center = ",
                     pcolls[i]->child[k]->center, 2);
              if ( (abs(xc-center[0]) <= width*1.000001) &&
                   (abs(yc-center[1]) <= width*1.000001) ) {
                // store the colleague
                newbox->colls[ncoll] = pcolls[i]->child[k];
                ncoll = ncoll + 1;
              }
            }
          }

          
        } 
        // loop over i


        cprinf("for newbox, ncoll = ", &ncoll, 1);
      }
      // if parent isn't NULL...


      
      newbox->ncoll = ncoll;
      cprinf("for newbox, ncoll = ", &ncoll, 1);
    }
    
  }


  if (ibox == 1) exit(0);
  
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
