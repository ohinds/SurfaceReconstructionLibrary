/*****************************************************************************
 * libsrUtil.h is the header file for utility functions for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#ifndef SR_UTIL_H
#define SR_UTIL_H

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<gsl/gsl_fit.h>

#include"libsrTypes.h"
#include"libvpUtil.h"
#include"surfUtil.h"
#include"list.h"

#include"libsrUtil.extern"

#define dist(a,b) sqrt(pow((a).x-(b).x,2)+pow((a).y-(b).y,2)+pow((a).z-(b).z,2))
#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)


/**
 * create a contour
 * returns a pointer to a new contour, or null, if flibvpure
 */
contour *createContour();

/**
 * clone a contour
 * returns a pointer to a new contour copied from the one passed,
 * or null, if flibvpure
 */
contour *cloneContour(contour *cont);

/**
 * deletes a contour
 */
void deleteContour(contour *cont);

/**
 * determines the area of the polygon bounded by a contour
 */
double contourArea(contour* cont);

/**
 * delete empty slices
 */
void deleteEmptySlices(list *slices);

/**
 * clones a pair of slices, maintaining adjacency
 */
void cloneSlices(list *slice1, list *slice2,
                 list **newSlice1, list **newSlice2);

/**
 * gets the average interslice distance
 */
double getAverageIntersliceDistance(list *slices);

/**
 * gets interslice distance
 */
double getIntersliceDistance(list *slice1, list *slice2);

/**
 * preprocesses slice contours before surface reconstruction
 */
int preprocessSliceContours(list *slices);

/**
 * straight up copy the forward to backward adjacency list
 */
void copyForwardToBackwardAdjacency(list *slice);

/**
 * builds a backward adjacency list between to slices
 */
void buildBackwardAdjacency(list *slice1, list *slice2);

/**
 * deletes a contour pair, and the contours themselves
 */
void deleteContourPairList(list *contourPairs);

/**
 * builds a backward adjacency list between two slices represented by
 * contour pairs. returns a list of contours representing the slice 2 contours
 */
list *buildBackwardAdjacencyFromContourPairs(list *contourPairs);

/**
 * tests to see if a slice has no vertices
 */
int sliceEmpty(list *slice);

/**
 * select all the vertices between two selected vertices
 */
int selectVerticesBetween(list *vertices);

/** some generic geometry utils */

/**
 * create a vertex structure with default values
 * returns a vertex structure or null if creation fails
 */
vertex *createVertex();

/**
 * copies a vertex
 * returns a copy of the passed vertex, or null, if cant allocate
 */
vertex *copyVertex(vertex *v);

/**
 * gets a new vertex that is the midpoint between two vertices
 */
vertex *getMidpoint(vertex *v1, vertex *v2);

/**
 * test for intersection of a contour and an edge
 */
int edgeContourIntersect(edge *e, contour *c);

/**
 * resample slice contours by integrated angle
 */
void resampleContoursByAngle(list *slices, double angle);

/**
 * resample one contour by integrated distance
 */
void resampleContourByDistance(contour *cont, double distance);

/**
 * resample slice contours by integrated distance
 */
void resampleContoursByDistance(list *slices, double distance);

/**
 * delete repeated vertices in a contour
 */
void deleteRepeatedContourVerts(list *slices);

/**
 * all vertices in the slices based on contour position, also assigns
 * boundaryness to the contour vertices
 */
void numberVertices(list *slices);

/**
 * label all vertices based on whether they are on a boundary or not
 */
void labelBoundaries(list *slices);

/**
 * delete repeated vertices in a contour
 */
void orientContours(list *slices);

/**
 * randomize contour indexing for closed contours
 */
void randomizeClosedContourIndexing(list *slices);

/* get the centroid of a contour */
void getContourCentroid(contour* c, vertex *cent);

/**
 * get the angle between two line segments
 */
double angleBetweenSegments(vertex *a, vertex *b, vertex *c);

/**
 * test a skeleton vertex for containment in the projection of a contour
 */
int inContourSkelVert(contour *c, vertex *v);

/**
 * get the skeleton of a closed contour
 * input is a closed contour
 * output is a list of vertices on the skeleton and a list of edges
 * connecting the vertices, both are null if skeleton cannot be found
 */
int getContourSkeleton(contour* c, list **vertices, list **edges);

/** Triangle interface utility functions **/

#define TRIANGLE_EXEC "triangle"

/**
 * write a Triangle .poly file from a contour
 * NOTE: this function assumes the contour vertices are connected in order
 */
void writeTrianglePolyFile(contour *c, char *filename);

/**
 * reads a node file output by Triangle into a list
 */
list *readTriangleNodeFile(char *filename);

/**
 * reads an edge file output by Triangle into a list
 */
list *readTriangleEdgeFile(char *filename);

/**
 * reads a node and corresponding edge file output by Triangle into a lists,
 * where the edge list has vertex pointers
 */
int readTriangleNodeAndEdgeFile(char *filebase, list **vertices, list **edges);

/**
 * runs Triangle on an existing .poly file, producing an off file
 * NOTE: uses the -qBPNEg flags when running Triangle
 */
int executeTriangle(char *polyInputFilename);

/**
 * runs Triangle on an existing .poly file with specified flags
 * returns the value returned by system()
 */
int executeTriangleWithFlags(char *polyInputFilename, char *flags);

/** math utils **/

/* phase unwrapping */
inline double unwrap(double a);

/**
 * wrapper for matrix multiply of vertex coords
 */
void matrixMult4byVert(float a[4][4], vertex b, vertex *c);

/**
 * perform a linear fit given a vector of observations and loadings via gsl
 */
double *linearFit(double *x, double *y, double *w, int n);

/* debugging funcs */
void printPath(FILE* fp, path *p);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/libsrUtil.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
