/*****************************************************************************
 * surfUtil.h is the header file for surface utility functions for
 * libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 * 
 *
 *****************************************************************************/

#ifndef SURF_UTIL_H
#define SURF_UTIL_H

#include"libsrTypes.h"
#include"list.h"
#include"tile.h"
#include"surfIO.h"
#include"surfTopology.h"
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_sort.h>
#include"surfUtil.extern"

#define DEFAULT_MAX_VERTICES 100000
//#define DEFAULT_MAX_VERTICES_EXPAND 10000 
#define DEFAULT_MAX_FACES 200000
//#define DEFAULT_MAX_FACES_EXPAND 20000

#define MAX_SKEL_V 2048
#define MAX_SKEL_N 1024

/** surface creation and deletion functions **/

/**
 * create a surface structure
 * returns a surface structure or null if creation flibvps
 */
surface *createSurfaceDefault();

/**
 * create a surface structure with a specified number of vertices and faces 
 * returns a surface structure or null if creation flibvps
 */
surface *createSurface(unsigned int nVerts, unsigned int nFaces);

/**
 * delete a surface and free all vertices and faces
 * returns success (1) or flibvpure (-1)
 */
int deleteSurface(surface *surf);

/**
 * expands the capacity of the vertex array of the surface to amount
 * returns success (1) or flibvpure (-1)
 */
int expandVertices(surface *surf, int amount);

/**
 * expands the capacity of the faces array of the surface to amount
 * returns success (1) or flibvpure (-1)
 */
int expandFaces(surface *surf, int amount);

/** surface conversion and manipulation tools **/

/**
 * add a vertex to a surface
 * this function just uses the x, y, and z fields of the vertex struct
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexStruct(surface *surf, vertex* v);

/**
 * add a vertex to a surface, without duplicating vertices this
 * function just the x, y, and z fields of the vertex struct to check
 * if the vertex is already in the surface. if so, the original index 
 * is returned.
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexStructNoDuplicate(surface *surf, vertex* v);

/**
 * add a vertex to a surface by coordinates
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexCoord(surface *surf, double x, double y, double z, double lab, int boundary);

/**
 * add a vertex to a surface by coordinates if it doesnt already exist 
 * in the surface
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexCoordNoDuplicate(surface *surf, double x, double y, double z, double lab, int boundary);

/**
 * add a vertex to a surface by an array of coords
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexArr(surface *surf, double *coords);

/**
 * add a vertex to a surface by an array of coords if coords not already
 * in the surface
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexArrNoDuplicate(surface *surf, double *coords);

/**
 * add a face to a surface
 * this function just uses the v1, v2, and v3 fields of the face struct
 * returns index of added face or flibvpure (-1)
 */
int addFaceStruct(surface *surf, face *f);

/**
 * add a face to a surface by indices
 * returns index of added face or flibvpure (-1)
 */
int addFaceInd(surface *surf, int v1, int v2, int v3);

/**
 * add a face to a surface by an array of indices
 * returns index of added face or flibvpure (-1)
 */
int addFaceArr(surface *surf, int *inds);

/**
 * add a facelist to a surface structure without duplicating vertices
 * already contained in a surface.
 */
void addFaceListToSurface(surface *surf, list *faceList);

/**
 * combine two existing surfaces
 * returns resulting surface
 */
surface *combineSurfaces(surface *surf1, surface *surf2);

/**
 * build a surface structure from a face list
 */
surface *buildSurfaceFromFaceList(list *faceList);

// contour capping methods

/**
 * caps a terminating surface
 * inputs are the contour to be capped and the z coordinate to add any vertices
 * returns a facelist, or NULL, if the contour cannot be capped
 */
list *capContour(contour *cont, enum CAPPING_METHOD cappingMethod, double z);

/**
 * caps a contour by converting its skeleton to a contour, then
 * triangulating them together
 * takes as input a contour and a z coordinate to place the skeleton
 * returns a face list, or NULL, if the method fails
 */
list *capContourSkeleton(contour *cont, double z);

/**
 * caps a contour by tiling it with a point at its centroid and the z
 * coordinate passed 
 * takes as input a contour and a z coordinate to place the point
 * returns a face list, or NULL, if the method fails
 */
list *capContourPoint(contour *cont, double z);

/**
 * determine if a contour is roughly circular
 */
int contourIsCircular(contour *cont);

/**
 * mark a vertex as on a boundary so we don't close a boundary we didnt cap
 */
void markAVertexAsBoundary(contour *cont);

/**
 * mark a vertex as on a boundary so we don't close a boundary we didnt cap
 * for all contour on a slice
 */
void markAVertexAsBoundarySlice(list *slice);

/**
 * removes a vertex from an adjacency list
 */
void removeAdjacentEntry(vertex **vlist, vertex *otherV);

/**
 * replace a vertex in an adjacency list with another
 */
void replaceAdjacentEntry(vertex **vlist, vertex *oldV, vertex *newV);

/**
 * appends a deleted nodes adjacency entries to the current list
 */
void appendAdjacentEntries(vertex ***vlist, vertex *v, vertex *otherV);

// branched hole filling methods

/**
 * fill holes
 */
void fillHoles(surface *surf, int numHolesLeft);

/**
 * fill holes created by the branching process. 
 */
void fillBranchedHoles(surface *surf);

/**
 * decides whether a contour pair is one side of a branched boundary hole
 */ 
int isBranchedBoundary(contourPair *cp);

/**
 * tests if two vertices from separate contours pairs were adjacent in the
 * original contour
 */
int verticesOriginallyAdjacent(contour *c1, contour *oc1, vertex *v1, vertex *ov1);

/**
 * given a surface and a vertex, determine which of two other vertices are
 * 1-neighbors of the passed vertex, if any.
 * RETURNS THE FIRST INSTANCE FOUND!
 */
vertex *getAdjacentVertex(surface *surf, vertex *v0, vertex *v10, vertex *v11);

/**
 * adds a three dimensional polygon to an existing surface.
 * input is an ordered set of points representing the boundary of the polygon.
 * WARNING! THIS IS A DIRTY HACK, IT WILL PRODUCE SELF_INTERSECTIONS!
 */
void add3DPolyToSurf(surface *surf, list *plist);

/**
 * recursively trace a poly boundary to recover an optimal tiling
 */ 
void tracePolyBoundary(surface *surf, int **lambda, int *vlist, int i, int k);

/**
 * determine the cost of adding this face to a polygonal tiling
 */
vector *getPolyTilingTriangleCost(surface *surf, int v1, int v2, int v3);

/** surface * surface query tools **/

/**
 * find a vertex by vertex struct
 * returns index of vertex or flibvpure (-1)
 */
int findVertexStruct(surface *surf, vertex *v);

/**
 * find a vertex by coords
 * returns index of vertex or flibvpure (-1)
 */
int findVertexCoords(surface *surf, double x, double y, double z);

/**
 * find a vertex by array
 * returns index of vertex or flibvpure (-1)
 */
int findVertexArr(surface *surf, double *coords);

/**
 * find a face by face struct
 * returns index of face or flibvpure (-1)
 */
int findFaceStruct(surface *surf, face *f);

/**
 * find a face by indices
 * returns index of face or flibvpure (-1)
 */
int findFaceInds(surface *surf, int v1, int v2, int v3);

/**
 * find a face by array
 * returns index of face or flibvpure (-1)
 */
int findFaceArr(surface *surf, int *inds);

/**
 * copy a contour
 */
void copyContour(contour *orig, contour *copy);

/**
 * use Triangle to get a delaunay triangulated contour
 * NOTE: uses the -qBPNEg flags when running Triangle
 */
list *triangulateContour(contour *c);

/**
 * copies an edge
 */
edge *copyEdge(edge *e);

/**
 * tests two vertices for equality
 */
int verticesEqual(vertex *v1, vertex *v2);


/**
 * tests two segments for intersection based on the 12 coordinates of
 * the end points of two three dimensional segments
 */
int segmentsIntersect(float x00, float y00, float z00, 
		      float x01, float y01, float z01, 
		      float x10, float y10, float z10, 
		      float x11, float y11, float z11);

/* /\** */
/*  * tests two segments for intersection based on the 8 coordinates of */
/*  * the end points of two two dimensional segments */
/*  *\/ */
/* int segmentsIntersect2d(float x00, float y00,  */
/* 			float x01, float y01,  */
/* 			float x10, float y10,  */
/* 			float x11, float y11); */

/**
 * tests two edges for intersection
 */
int edgesIntersect(edge *e1, edge *e2);

/**
 * test intersection of an edge and a face
 *
 * algorithm pieced from two places 
 * http://astronomy.swin.edu.au/~pbourke/geometry/planeline/
 * http://www.ecse.rpi.edu/Homepages/wrf/research/geom/pnpoly.html
 */
int edgeFaceIntersect(face *f, edge *e);

/**
 * center the surf at its center of mass
 */
void centerSurf(surface *surf);

/**
 * calculate the center of mass
 */
vector computeCenterOfMass(surface *surf);

/**
 * calculate the bounds in each dimension
 */
void computeSurfBounds(surface *surf, double *bounds);

/** 
 * transform the coordinates of the surface vertices by a linear transform
 */
int transformSurfaceVertices(surface *surf, float A[4][4]);

/**
 * compute the face normals for all faces
 */
void computeFaceNormals(surface *surf);

/**
 * compute the vertex normals for all vertices
 */
void computeVertexNormals(surface *surf);

/**
 * get the unit normal to a vertex
 */
void getVertexNormal(surface *surf, int vertexInd, vector *n);

/**
 * compute the face normals for all faces
 */
void computeFaceNormals(surface *surf);

/**
 * get the unit normal to a face
 */
void getFaceNormal(face *f, vector *n);

/**
 * get the unit normal to a face
 */
void getFaceNormalV(double *v1, double *v2, double *v3, vector *n);

/**
 * get the unit normal to a face
 */
void getFaceNormalVerts(vertex *v1, vertex *v2, vertex *v3, vector *n);

/**
 * cross procuct
 */
void cross(vector *v1, vector *v2, vector *c);

/**
 * normalize a vector
 */
void normalizeVector(vector *v);

/**
 * dumps a contour pair to a stream
 */
void dumpContourPair(FILE *fp, contourPair *cp);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/surfUtil.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
