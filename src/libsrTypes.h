/*****************************************************************************
 * libsrTypes.h is the header file containing type declarations for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#ifndef SR_TYPES_H
#define SR_TYPES_H

#include"list.h"

#define SR_TYPES_VERSION_H "$Id: libsrTypes.h,v 1.29 2007/05/22 19:18:11 oph Exp $"

#define SR_MAX_STR_LEN 255

#ifdef SR_DEBUG_FLAG
#define SR_DEBUG 1
#else
#define SR_DEBUG 0
#endif

#define SR_FAILURE -1
#define SR_SUCCESS 1

#define SR_PI 3.14159265358979
#define SR_TOL 0.00000000001

/* stuffSR_ for bounds */
#define SR_BIG 1.0e32
#define SR_MED 1.0e16
#define SR_LITTLE -SR_BIG
#define SR_TINY 0.000000001

#define SR_X 0
#define SR_Y 1
#define SR_Z 2

#define SR_XMIN 0
#define SR_XMAX 1
#define SR_YMIN 2
#define SR_YMAX 3
#define SR_ZMIN 4
#define SR_ZMAX 5

/* simple int list data struct */
typedef struct {
  int val;
} intNode;

/* simple float list data struct */
typedef struct {
  float val;
} floatNode;

/* simple double list data struct */
typedef struct {
  double val;
} doubleNode;

/* boundary state */
enum BOUNDARYNESS {INTERIOR, BOUNDARY, BRANCHED_BOUNDARY};

/* structure for a vertex */
typedef struct {
  /* geometry */
  double x;
  double y;
  double z;

  /* if this vertex is selected */
  int selected;

  /* node number in its contour */
  int number;

  /* if this is a boundary vertex */
  int boundary;

  /* the label of this vertex */
  int label;
} vertex;

/* node for an edge in graph */
typedef struct {
  /* endpoints */
  vertex *v1;
  vertex *v2;

  /* length of the edge */
  double length;

  /* distance from the last vertical step in the edge-face graph */
  int lastVert;

  /* distance from the last horizontal step in the edge-face graph */
  int lastHorz;
} edge;

/* node for a face in graph */
typedef struct {
  /* left and right edges */
  edge *left;
  edge *right;

  /* vertex indices */
  int v1;
  int v2;
  int v3;

  /* "cost" of the face */
  double cost;
} face;

/* correspondence methods */
enum CORRESPONDENCE_METHOD {HISTOGRAM, KMEANS, DENSITY};
enum CORRESPONDENCE_SCOPE {GLOBAL, LOCAL};

/* histogram */
typedef struct {
  int numBins;
  double binSize;

  int *hist;

  double minval;
  double maxval;
} histogram;

/** contour graph **/

/** contour type **/
#define OPEN 0
#define CLOSED 1
#define BRANCHED_OPEN 2

/* for contour branches */
enum BRANCHTYPE {DISCONNECT, CONNECT};
enum BRANCHDIRECTION {FORWARD, BACKWARD};

typedef struct {
  /* 1D list of vertices making up the contour */
  list *vertices;

  /* whether this contour is closed or not */
  int closed;

  /* list of the contours on the next slice that this one is connected to */
  list *adjacentContours;
  list *adjacentBackwardContours;

  struct contour *origin;
} contour;

typedef struct {
  contour *c1;
  contour *c1Origin;

  contour *c2;
  contour *c2Origin;
} contourPair;

/* edge directions */
enum DIR {VERT,HORZ};

/* closure cases */
enum CLOSURE {BOTH_CLOSED, ONE_CLOSED, NONE_CLOSED};

/* tiling methods */
enum TILING_METHOD {OPTIMIZING};

/* capping methods */
enum CAPPING_METHOD {NO_CAPPING, TRIANGULATE, TILE_SKELETON, TILE_POINT};

/* node for a path through graph */
typedef struct {
  list *faceList; /* faces in path */
  double cost;    /* cost of the path with all terms */
  double baseCost;/* base cost of the path (usually surface area)
                   * keep this around so we can penalize based on
                   * surface area, not on previously applied penalties */
} path;

/* include the heap header here for its element declarations */
#include"heap.h"

/* struct for a graph */
typedef struct {
  int m,n;     /* size of the graph */

  /* closure cases */
  enum CLOSURE closureCase;

  edge ***v;    /* vertex matrix */
  double **wh;  /* horizontal edge weights */
  double **wv;  /* vertical edge weights */
  edge ***pi;   /* predecessor matrix */

  heap_it_handler **h; /* pointers to heap items for each edge */
} graph;

/* curvature modes, either vertex or face, depending on how curvature is stored */
enum CURVATURE_MODE {NO_CURVATURE, VERTEX_CURVATURE, FACE_CURVATURE};

/* surface file formats */
enum SURFACE_FORMAT {INVALID_SURFACE_FORMAT, OFF, OBJ, MGHSURF};

/* surface manifoldness state */
enum SURACE_MANIFOLDNESS_STATE {SURF_NON_MANIFOLD, SURF_MANIFOLD, SURF_UNKNOWN};

/* vertex manifoldness state */
enum VERTEX_MANIFOLDNESS_STATE {VERT_MANIFOLD,
                                VERT_MORE_THAN_TWO_TRIS_AT_EDGE_ERROR,
                                VERT_CONE_ERROR,
                                VERT_TWO_HOLE_ERROR,
                                VERT_UNKNOWN};

/** star edge for topo fixing **/
typedef struct {
  int edge[2];
  int neighborIndex[2];
  int face;
} staredge;

/** surface **/
typedef struct {
  /* vertex data and count */
  unsigned int numVertices;
  unsigned int maxVertices;
  double **vertices;
  double *vertexLabels;
  int *vertexBoundaryness;

  /* vertex normals */
  double **vertexNormals;
  int numVertexNormals;

  /* face data and count */
  unsigned int numFaces;
  unsigned int maxFaces;
  unsigned int **faces;

  /* face normals */
  double **faceNormals;
  unsigned int numFaceNormals;

  /* curvature data */
  int *curvature;
  int curvatureMode;

  /* topology data */
  int *NTaV;                // number of triangles at vertex
  int **LTaV;       // label of triangles at vertex
  int **v1TaV;        // first other vertex at triangle at vertex
  int **v2TaV;        // second other vertex at triangle at vertex
  int *NNoV;        // number of neighbors of vertex
  int maxNTaV;              // max number or triangles at any vertex
  int E;        // number of edges
  int **Neighbors;      // neighbors of vertex
  int **NTaE;       // number of triangles at each edge
  int **LT1aNE;       // label of the first triangle at neighbor edge
  int **LT2aNE;       // label of the second triangle at neighbor edge
  int manifoldness;     // manifoldness state
  int *vertexManifoldness;  // array of vertex manifoldness state
  int **OT;                 // opposite triangle
  int **OV;                 // opposite vertex
  list *CC;                 // list containing one vertex from each
  // connected component
} surface;

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/libsrTypes.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
