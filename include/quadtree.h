
#include <complex.h>

//
// tree declarations for 2d fast direct solver code
//




struct quadtree_opts {

  int maxlev;
  int maxboxes;  
  int nmax;
  double maxwidth;
};





struct quadtree_box {

  int id;
  int level;
  double center[2];
  double width;
  
  int npts;
  double *xys;
  int *perm;

  struct quadtree_box *parent;
  struct quadtree_box *child[4];

  int ncoll;
  struct quadtree_box *colls[8];
  
};





void quadtree_findall_list1(int nboxes, struct quadtree_box *tree);

void quadtree_get_extent(int npts, double *xys, double *center, double *width);

void quadtree_permute(int m, int n, int *p, double *xs);

void quadtree_unpermute(int m, int n, int *p, double *xs);

void quadtree_getsiblings(struct quadtree_box *box, struct quadtree_box **pch);
  
void quadtree_isleaf(struct quadtree_box *box, int *isleaf);

void quadtree_rand(int n, double *rs);

void quadtree_getboxes_bylevel(int nboxes, struct quadtree_box *tree,
                             int lev, int *nb, int *boxes);

void quadtree_getleafs(int nboxes, struct quadtree_box *tree, int *nleafs,
                     int *leafs);

void quadtree_restrict1(int *nlev, int *nboxes, struct quadtree_box *tree);

void quadtree_build_lr(double *center, double width,
                    int npts, double *xys, int *perm, 
                  struct quadtree_opts opts, int *nlev,
                  int *nboxes, struct quadtree_box *tree);

void quadtree_build(double *center, double width,
                    int npts, double *xys, int *perm, 
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

