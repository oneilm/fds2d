
#include <complex.h>

//
// tree declariations for 2d butterfly code
//



struct quadtree_opts {

  int maxlev;
  int maxboxes;  
  int nmax;
  double maxwidth;
  //double rmax;
};





struct quadtree_mpoles {
  // struct to store multipole and local expansions for every node in
  // the tree
  int nterms;
  int nd;
  double sc;
  double complex *mpole;
  double complex *lpole;
};



struct quadtree_box {

  int id;
  int level;
  
  double center[2];
  double width;
  //double r;
  
  int nsrcs;
  double *srcs;
  int *ps;

  int ntargs;
  double *targs;
  int *pt;
  
  struct quadtree_box *parent;
  struct quadtree_box *child[4];

  struct quadtree_mpoles *fmm;
};





void quadtree_permute(int m, int n, int *p, double *xs);

void quadtree_unpermute(int m, int n, int *p, double *xs);

void quadtree_getsiblings(struct quadtree_box *box, struct quadtree_box **pch);
  
void quadtree_isleaf(struct quadtree_box *box, int *isleaf);

void quadtree_rand(int n, double *rs);

void quadtree_getboxes_bylevel(int nboxes, struct quadtree_box *tree,
                             int lev, int *nb, int *boxes);

void quadtree_getleafs(int nboxes, struct quadtree_box *tree, int *nleafs,
                     int *leafs);

void quadtree_build(int m, double *targs, int *pt, int n, double *srcs, int *ps,
                  struct quadtree_opts opts, int *nlev,
                  int *nboxes, struct quadtree_box *tree);

void quadtree_build_up(double *center, double width, 
                     int m, double *targs, int *pt, int n, double *srcs, int *ps,
                  struct quadtree_opts opts, int *nlev,
                  int *nboxes, struct quadtree_box *tree);

void quadtree_split(int n,  int *nboxes, struct quadtree_box *tree);

void quadtree_printbox(char *name, struct quadtree_box *box);

void quadtree_print_tree(int nboxes, struct quadtree_box *tree);

void quadtree_key(double *center, double *src, int *ikey);

void quadtree_plotboxes(char *filename, int nboxes, struct quadtree_box *tree,
                      char *title);

void quadtree_swapi(int n, int *xs, int *ys);

void quadtree_swapd(int n, double *xs, double *ys);

