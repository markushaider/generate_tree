/* fast3tree.c
   A fast BSP tree implementation.
   (C) 2010, Peter Behroozi, Stanford University.
   Usage:
      #define FAST3TREE_TYPE struct mytype //mytype must have an array pos[]
      #define FAST3TREE_DIM X //Optional---the dimension of pos[], default 3
      #define FAST3TREE_POINTS_PER_LEAF X //Optional, num. points per leaf node
      #define FAST3TREE_PREFIX XYZ //Optional, if using multiple different trees
      #define FAST3TREE_FLOATTYPE float //Optional: or double, or long double
      #include "fast3tree.c"


   PUBLIC METHODS:
   Initialize a fast3tree from a list of points:
      struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p);

   Rebuild a fast3tree from a new (or the same) list of points:
      void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p);

   Rebuilds the tree boundaries, but keeps structure the same:
      void fast3tree_maxmin_rebuild(struct fast3tree *t);

   Frees the tree memory and sets tree pointer to NULL.
      void fast3tree_free(struct fast3tree **t);

   Initialize a fast3tree results structure:
      struct fast3tree_results *fast3tree_results_init(void);

   Find all points within a sphere centered at c[FD] with radius r:
   (FD = FAST3TREE_DIM, usually 3 unless you changed it).
      void fast3tree_find_sphere(struct fast3tree *t,
              struct fast3tree_results *res, float c[FD], float r);

   Find all points within a sphere centered at c[FD] with radius r,
   assuming periodic boundary conditions:
      int fast3tree_find_sphere_periodic(struct fast3tree *t, 
              struct fast3tree_results *res, float c[FD], float r);

   Reset memory stored in results structure:
      void fast3tree_results_clear(struct fast3tree_results *res);

   Free the results structure returned by fast3tree_find_sphere:
      void fast3tree_results_free(struct fast3tree_results *res);

   END PUBLIC METHODS
*/

#ifndef _FAST3TREE_C_
#define _FAST3TREE_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#ifndef FAST3TREE_PREFIX
#define FAST3TREE_PREFIX
#endif /* FAST3TREE_PREFIX */

#ifndef FAST3TREE_DIM
#define FAST3TREE_DIM 3
#endif /* FAST3TREE_DIM */

#ifndef POINTS_PER_LEAF
#define POINTS_PER_LEAF 40
#endif  /* POINTS_PER_LEAF */

#ifdef FAST3TREE_FLOATTYPE
#define float FAST3TREE_FLOATTYPE
#endif /* FAST3TREE_FLOATTYPE */

#ifndef FAST3TREE_TYPE
#error Usage:
#error #define FAST3TREE_TYPE point_structure
#error #include "fast3tree.c"
#define FAST3TREE_TYPE struct fast3tree_default_point 
#endif /* FAST3TREE_TYPE */

#define _F3TS(a,b)  a ## b
#define _F3TN(a,b) _F3TS(a, b)
#undef tree3_node
#define tree3_node _F3TN(FAST3TREE_PREFIX,tree3_node)
#undef fast3tree
#define fast3tree  _F3TN(FAST3TREE_PREFIX,fast3tree)
#undef fast3tree_results
#define fast3tree_results _F3TN(FAST3TREE_PREFIX,fast3tree_results)
#undef fast3tree_default_point
#define fast3tree_default_point _F3TN(FAST3TREE_PREFIX,fast3tree_default_point)

struct tree3_node;
struct fast3tree;
struct fast3tree_results;
struct fast3tree_default_point;

#undef fast3tree_init
#define fast3tree_init _F3TN(FAST3TREE_PREFIX,fast3tree_init)
struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p);

#undef fast3tree_rebuild
#define fast3tree_rebuild _F3TN(FAST3TREE_PREFIX,fast3tree_rebuild)
void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p);

#undef fast3tree_maxmin_rebuild
#define fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,fast3tree_maxmin_rebuild)
void fast3tree_maxmin_rebuild(struct fast3tree *t);

#undef fast3tree_results_init
#define fast3tree_results_init _F3TN(FAST3TREE_PREFIX,fast3tree_results_init)
struct fast3tree_results *fast3tree_results_init(void);

#undef fast3tree_find_sphere
#define fast3tree_find_sphere _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere)
static inline void fast3tree_find_sphere(struct fast3tree *t,
			   struct fast3tree_results *res, float c[3], float r);

#undef fast3tree_find_sphere_periodic
#define fast3tree_find_sphere_periodic _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere_periodic)
int fast3tree_find_sphere_periodic(struct fast3tree *t,
			   struct fast3tree_results *res, float c[3], float r);

#undef fast3tree_results_clear
#define fast3tree_results_clear _F3TN(FAST3TREE_PREFIX,fast3tree_results_clear)
void fast3tree_results_clear(struct fast3tree_results *res);

#undef fast3tree_results_free
#define fast3tree_results_free _F3TN(FAST3TREE_PREFIX,fast3tree_results_free)
void fast3tree_results_free(struct fast3tree_results *res);


#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

struct fast3tree_default_point { float pos[FAST3TREE_DIM]; };

struct tree3_node {
  float min[FAST3TREE_DIM], max[FAST3TREE_DIM];
  int64_t num_points;
  int32_t div_dim;
  struct tree3_node *left, *right, *parent;
#ifdef FAST3TREE_EXTRA_INFO
  FAST3TREE_EXTRA_INFO;
#endif /*FAST3TREE_EXTRA_INFO*/
  FAST3TREE_TYPE *points;
};

struct fast3tree {
  FAST3TREE_TYPE *points;
  int64_t num_points;
  struct tree3_node *root;
  int64_t num_nodes;
  int64_t allocated_nodes;
};

struct fast3tree_results {
  int64_t num_points;
  int64_t num_allocated_points;
  FAST3TREE_TYPE **points;
};


/* PRIVATE METHODS */

#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _F3TN(FAST3TREE_PREFIX,_fast3tree_check_realloc)
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason);

#undef _fast3tree_build
#define _fast3tree_build _F3TN(FAST3TREE_PREFIX,_fast3tree_build)
void _fast3tree_build(struct fast3tree *t);

#undef _fast3tree_maxmin_rebuild
#define _fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,_fast3tree_maxmin_rebuild)
void _fast3tree_maxmin_rebuild(struct tree3_node *n);



struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p) {
  struct fast3tree *new=NULL;
  new = _fast3tree_check_realloc(new,sizeof(struct fast3tree), "Allocating fast3tree.");
  if (!new) return 0;
  memset(new, 0, sizeof(struct fast3tree));
  fast3tree_rebuild(new, n, p);
  return new;
}

void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p) {
  t->points = p;
  t->num_points = n;
  _fast3tree_build(t);
}

void fast3tree_maxmin_rebuild(struct fast3tree *t) {
  _fast3tree_maxmin_rebuild(t->root);
}

#undef fast3tree_free
#define fast3tree_free _F3TN(FAST3TREE_PREFIX,fast3tree_free)
void fast3tree_free(struct fast3tree **t) {
  if (!t) return;
  struct fast3tree *u = *t;
  if (u) {
    free(u->root);
    free(u);
  }
  *t = NULL;
}

#undef _fast3tree_fast_box_not_intersect_sphere
#define _fast3tree_fast_box_not_intersect_sphere \
  _F3TN(FAST3TREE_PREFIX,_fast3tree_fast_box_not_intersect_sphere)
/* Some false negatives, but best compromise for speed and accuracy.  */
static inline int _fast3tree_fast_box_not_intersect_sphere(struct tree3_node *node, float c[FAST3TREE_DIM], float r) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if ((node->min[i]-c[i]) > r) return 1;
    if ((c[i]-node->max[i]) > r) return 1;
  }
  return 0;
}

/* Fast, accurate. */
#undef _fast3tree_box_inside_sphere
#define _fast3tree_box_inside_sphere _F3TN(FAST3TREE_PREFIX,_fast3tree_box_inside_sphere)
static inline int _fast3tree_box_inside_sphere(struct tree3_node *node, float c[FAST3TREE_DIM], float r) {
  int i;
  float dx, dx2, dist = 0, r2 = r*r;
  if (fabsf(c[0]-node->min[0]) > r) return 0; //Rapid short-circuit.
  for (i=0; i<FAST3TREE_DIM; i++) {
    dx = node->min[i] - c[i];
    dx *= dx;
    dx2 = c[i]-node->max[i];
    dx2 *= dx2;
    if (dx2 > dx) dx = dx2;
    dist += dx;
    if (dist > r2) return 0;
  }
  return 1;
}

#undef _fast3tree_sphere_inside_box
#define _fast3tree_sphere_inside_box _F3TN(FAST3TREE_PREFIX,_fast3tree_sphere_inside_box)
static inline int _fast3tree_sphere_inside_box(struct tree3_node *node, float c[FAST3TREE_DIM], float r) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (c[i]-r < node->min[i]) return 0;
    if (c[i]+r < node->max[i]) return 0;
  }
  return 1;
}

#undef _fast3tree_check_results_space
#define _fast3tree_check_results_space _F3TN(FAST3TREE_PREFIX,_fast3tree_check_results_space)
static inline void _fast3tree_check_results_space(struct tree3_node *n, struct fast3tree_results *res) {
  if (res->num_points + n->num_points > res->num_allocated_points) {
    res->num_allocated_points = res->num_points + n->num_points + 1000;
    res->points = _fast3tree_check_realloc(res->points, 
 res->num_allocated_points * sizeof(FAST3TREE_TYPE *), "Allocating fast3tree results");
  }
}

#undef _fast3tree_find_sphere
#define _fast3tree_find_sphere _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere)
void _fast3tree_find_sphere(struct tree3_node *n, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  int64_t i,j;
  float r2, dist, dx;

  if (_fast3tree_fast_box_not_intersect_sphere(n,c,r)) return;
  if (_fast3tree_box_inside_sphere(n,c,r)) { /* Entirely inside sphere */  
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++)
      res->points[res->num_points+i] = n->points+i;
    res->num_points += n->num_points;
    return;
  }

  r2 = r*r;
  if (n->div_dim < 0) { /* Leaf node */
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      j = dist = 0;
      for (; j<FAST3TREE_DIM; j++) {
	dx = c[j]-n->points[i].pos[j];
	dist += dx*dx;
      }
      if (dist < r2) {
	res->points[res->num_points] = n->points + i;
	res->num_points++;
      }
    }
    return;
  }
  /* Otherwise, search in descendant tree nodes. */
  //if ((n->left->min[n->div_dim] - c[n->div_dim]) < r && (c[n->div_dim] - n->left->max[n->div_dim]) < r)
  _fast3tree_find_sphere(n->left, res, c, r);
  _fast3tree_find_sphere(n->right, res, c, r);
}

struct fast3tree_results *fast3tree_results_init(void) {
  struct fast3tree_results *res = 
    _fast3tree_check_realloc(NULL, sizeof(struct fast3tree_results), "Allocating fast3tree results structure.");
  res->points = NULL;
  res->num_points = 0;
  res->num_allocated_points = 0;
  return res;
}

static inline void fast3tree_find_sphere(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  res->num_points = 0;
  _fast3tree_find_sphere(t->root, res, c, r);
}

#undef _fast3tree_find_sphere_periodic_dim
#define _fast3tree_find_sphere_periodic_dim _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere_periodic_dim)
void _fast3tree_find_sphere_periodic_dim(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r, float dims[FAST3TREE_DIM], int dim) {
  float c2[FAST3TREE_DIM];
  if (dim<0) {
    _fast3tree_find_sphere(t->root, res, c, r);
    return;
  }
  memcpy(c2, c, sizeof(float)*FAST3TREE_DIM);
  _fast3tree_find_sphere_periodic_dim(t, res, c2, r, dims, dim-1);
  if (c[dim]+r > t->root->max[dim]) {
    c2[dim] = c[dim]-dims[dim];
    _fast3tree_find_sphere_periodic_dim(t, res, c2, r, dims, dim-1);
  }
  if (c[dim]-r < t->root->min[dim]) {
    c2[dim] = c[dim]+dims[dim];
    _fast3tree_find_sphere_periodic_dim(t, res, c2, r, dims, dim-1);
  }
}

int fast3tree_find_sphere_periodic(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  float dims[FAST3TREE_DIM];
  int i;
  
  if (_fast3tree_sphere_inside_box(t->root, c, r)) {
    fast3tree_find_sphere(t, res, c, r);
    return 2;
  }

  for (i=0; i<FAST3TREE_DIM; i++) {
    dims[i] = t->root->max[i] - t->root->min[i];
    if (r*2.0 > dims[i]) return 0; //Avoid wraparound intersections.
  }

  res->num_points = 0;
  _fast3tree_find_sphere_periodic_dim(t, res, c, r, dims, FAST3TREE_DIM-1);
  return 1;
}


void fast3tree_results_clear(struct fast3tree_results *res) {
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
}


void fast3tree_results_free(struct fast3tree_results *res) {
  if (!res) return;
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
  free(res);
}


#undef _fast3tree_find_largest_dim
#define _fast3tree_find_largest_dim _F3TN(FAST3TREE_PREFIX,_fast3tree_find_largest_dim)
static inline int64_t _fast3tree_find_largest_dim(float * min, float * max) {
  int64_t i, dim = FAST3TREE_DIM-1;
  float d = max[FAST3TREE_DIM-1]-min[FAST3TREE_DIM-1], d2;
  for (i=0; i<(FAST3TREE_DIM-1); i++) {
    d2 = max[i] - min[i];
    if (d2 > d) { d=d2; dim = i; }
  }
  return dim;
}

#undef _fast3tree_sort_dim_pos
#define _fast3tree_sort_dim_pos _F3TN(FAST3TREE_PREFIX,_fast3tree_sort_dim_pos)
static inline int64_t _fast3tree_sort_dim_pos(struct tree3_node * node) {
  int64_t dim = node->div_dim = 
    _fast3tree_find_largest_dim(node->min, node->max);
  FAST3TREE_TYPE *p = node->points;
  int64_t i,j=node->num_points-1;
  FAST3TREE_TYPE temp;
  float lim = 0.5*(node->max[dim]+node->min[dim]);

  if (node->max[dim]==node->min[dim]) return node->num_points;

  for (i=0; i<j; i++) {
    if (p[i].pos[dim] > lim) {
      temp = p[j];
      p[j] = p[i];
      p[i] = temp;
      i--;
      j--;
    }
  }
  if ((i==j) && (p[i].pos[dim] <= lim)) i++;
  return i;
}

#undef _fast3tree_find_minmax
#define _fast3tree_find_minmax _F3TN(FAST3TREE_PREFIX,_fast3tree_find_minmax)
static inline void _fast3tree_find_minmax(struct tree3_node *node) {
  int64_t i, j;
  float x;
  FAST3TREE_TYPE * p = node->points;
  assert(node->num_points > 0);
  for (j=0; j<FAST3TREE_DIM; j++) node->min[j] = node->max[j] = p[0].pos[j];
  for (i=1; i<node->num_points; i++)  {
    for (j=0; j<FAST3TREE_DIM; j++) {
      x = p[i].pos[j];
      if (x<node->min[j]) node->min[j] = x;
      else if (x>node->max[j]) node->max[j] = x;
    }
  }
}

#undef _fast3tree_split_node
#define _fast3tree_split_node _F3TN(FAST3TREE_PREFIX,_fast3tree_split_node)
void _fast3tree_split_node(struct fast3tree *t, struct tree3_node *node) {
  int64_t num_left;
  struct tree3_node *null_ptr = NULL;
  struct tree3_node *left, *right;
  int64_t left_index, node_index;

  num_left = _fast3tree_sort_dim_pos(node);
  if (num_left == node->num_points || num_left == 0) 
  { //In case all node points are at same spot
    node->div_dim = -1;
    return;
  }

  node_index = node - t->root;
  if ((t->num_nodes+2) > t->allocated_nodes) {
    t->allocated_nodes += 1000;
    t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes");
    node = t->root + node_index;
  }

  left_index = t->num_nodes;
  t->num_nodes+=2;

  node->left = null_ptr + left_index;
  node->right = null_ptr + (left_index + 1);

  left = t->root + left_index;
  right = t->root + (left_index + 1);
  memset(left, 0, sizeof(struct tree3_node)*2);

  right->parent = left->parent = null_ptr + node_index;
  left->num_points = num_left;
  right->num_points = node->num_points - num_left;
  left->div_dim = right->div_dim = -1;
  left->points = node->points;
  right->points = node->points + num_left;

  _fast3tree_find_minmax(left);
  _fast3tree_find_minmax(right);

  if (left->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, left);

  right = t->root + (left_index + 1);
  if (right->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, right);
}

#undef _fast3tree_rebuild_pointers
#define _fast3tree_rebuild_pointers _F3TN(FAST3TREE_PREFIX,_fast3tree_rebuild_pointers)
void _fast3tree_rebuild_pointers(struct fast3tree *t) {
  int64_t i;
  struct tree3_node *nullptr = NULL;
  for (i=0; i<t->num_nodes; i++) {
#define REBUILD(x) x = t->root + (x - nullptr)
    REBUILD(t->root[i].left);
    REBUILD(t->root[i].right);
    REBUILD(t->root[i].parent);
#undef REBUILD
  }
}

#undef _fast3tree_build
#define _fast3tree_build _F3TN(FAST3TREE_PREFIX,_fast3tree_build)
void _fast3tree_build(struct fast3tree *t) {
  int64_t i, j;
  struct tree3_node *root;
  FAST3TREE_TYPE tmp;
  t->allocated_nodes = (3+t->num_points/(POINTS_PER_LEAF/2));
  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes"); //Estimate memory load
  t->num_nodes = 1;

  //Get rid of NaNs / infs
  for (i=0; i<t->num_points; i++) {
    for (j=0; j<FAST3TREE_DIM; j++) if (!isfinite(t->points[i].pos[j])) break;
    if (j<FAST3TREE_DIM) {
      tmp = t->points[i];
      t->num_points--;
      t->points[i] = t->points[t->num_points];
      t->points[t->num_points] = tmp;
      i--;
    }
  }

  root = t->root;
  memset(root, 0, sizeof(struct tree3_node));
  root->num_points = t->num_points;
  root->points = t->points;
  root->div_dim = -1;
  if (t->num_points) _fast3tree_find_minmax(root);
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->min[j]));
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->max[j]));

  if (root->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, root);

  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->num_nodes), "Tree nodes");
  t->allocated_nodes = t->num_nodes;
  _fast3tree_rebuild_pointers(t);
}

#undef _fast3tree_maxmin_rebuild
#define _fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,_fast3tree_maxmin_rebuild)
void _fast3tree_maxmin_rebuild(struct tree3_node *n) {
  int i;
  if (n->div_dim < 0 && n->num_points) {
    _fast3tree_find_minmax(n);
    return;
  }
  _fast3tree_maxmin_rebuild(n->left);
  _fast3tree_maxmin_rebuild(n->right);
  memcpy(n->min, n->left->min, sizeof(float)*FAST3TREE_DIM);
  memcpy(n->max, n->right->max, sizeof(float)*FAST3TREE_DIM);
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (n->min[i] > n->right->min[i]) n->min[i] = n->right->min[i];
    if (n->max[i] < n->left->max[i]) n->max[i] = n->left->max[i];
  }
}


#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _F3TN(FAST3TREE_PREFIX,_fast3tree_check_realloc)
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if (res == NULL) {
    fprintf(stderr, "[Error] Failed to allocate memory (%s)!\n", reason);
    exit(1);
  }
  return res;
}


#undef _fast3tree_set_minmax
#define _fast3tree_set_minmax _F3TN(FAST3TREE_PREFIX,fast3tree_set_minmax)
void _fast3tree_set_minmax(struct fast3tree *t, float min, float max) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    t->root->min[i] = min;
    t->root->max[i] = max;
  }
}

#undef fast3tree_find_next_closest_distance
#define fast3tree_find_next_closest_distance _F3TN(FAST3TREE_PREFIX,fast3tree_find_next_closest_distance)
float fast3tree_find_next_closest_distance(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM]) {
  int64_t i=0, j;
  float dist = 0, dx, min_dist = 0;
  struct tree3_node *nd = t->root;
  
  while (nd->div_dim >= 0) {
    if (c[nd->div_dim] <= (nd->left->max[nd->div_dim])) { nd = nd->left; }
    else { nd = nd->right; }
  }

  while (nd != t->root && (nd->min[nd->parent->div_dim] == nd->max[nd->parent->div_dim])) nd = nd->parent;

  min_dist = 0;
  for (i=0; i<nd->num_points; i++) {
    dist = 0;
    for (j=0; j<FAST3TREE_DIM; j++) {
      dx = c[j] - nd->points[i].pos[j];
      dist += dx*dx;
    }
    if (!min_dist || (dist && dist < min_dist)) min_dist = dist;
  }
  min_dist = sqrt(min_dist);

  fast3tree_find_sphere(t, res, c, min_dist*1.01);
  for (i=0; i<res->num_points; i++) {
    dist = 0;
    for (j=0; j<FAST3TREE_DIM; j++) { dx = res->points[i]->pos[j] - c[j]; dist += dx*dx; }
    if ((dist && dist < min_dist)) min_dist = dist;
  }
  //if (!min_dist) fprintf(stderr, "Huh???\n");
  return sqrt(min_dist);
}

#undef float

#endif /* _FAST3TREE_C_ */
