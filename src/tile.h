/*****************************************************************************
 * tile.h is the header file with functions to tile a surface between two
 * slices/contours for libsr 
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#ifndef TILE_H
#define TILE_H

#define TILE_VERSION_H "$Id: tile.h,v 1.7 2007/05/22 19:18:11 oph Exp $"

#include"tile.extern"

#include"libsrTypes.h"
#include"libsrUtil.h"
#include"list.h"
#include"surfUtil.h"
#include"branch.h"

/* some geometric macros */

/* length of an edge in the actual contour graph */
#define len(u,v) sqrt(pow(u->x-v->x,2) + pow(u->y-v->y,2) + pow(u->z-v->z,2))

/**
 * get a surface represented by contours arranged in slices
 */
surface *tileSlices(list *slices);

/**
 * branch and tile a corresponded slice pair
 */
void tileSliceByOptimizingMethod(surface *surf, list *slice1, list *slice2);

/**
 * tiles the contours on a slice by an optimizing method derived from the
 * tiling method of Fuchs, et al. (1977). The resulting surface is
 * nonintersecting and of minimum surface area.
 */
void tileSliceByOldOptimizingMethod(surface *surf, list *slice);

/**
 * tile a contour to a single vertex, return a facelist
 */
list* tileContourToSingleVertex(contour *cont, vertex* v);

/**
 * build a surface from a set of branched contours
 */
surface *getBranchedSurface(contour *preBranch, list *postBranch);

/**
 * add a surface for two contours to an existing surface
 */
void addContoursToSurf(surface *surf, contour *contour1, contour* contour2);

/**
 * get face list between two contours. just a wrapper around getMinPath
 */
list *getFaceList(contour *contour1, contour *contour2);

/**
 * get the minimum path for two contours
 */
path *getMinPath(contour *contour1, contour *contour2);

/**
 * build planar graph for fuch's algorithm to operate on
 */
graph *buildPlanarGraph(contour *contour1, contour *contour2);

/**
 * free the allocated contents of a planar graph
 */
void freePlanarGraph(graph *G);

/**
 * find all shortest paths from the m starting nodes to the n ending
 * nodes implementation of algorithm presented in fuchs (1977) for two
 * closed contours
 */
path **allPaths(graph *G);

/**
 * find all shortest paths from the m starting nodes to the n ending
 * nodes, but dont wrap to the begining of contour2. implementation of
 * algorithm presented in fuchs (1977) for first contour closed and
 * second open
 */
path **allPathsNoclosure(graph *G);

/**
 * find two paths, one from the first vertices of each contour, the
 * second from the last of the first contour and the first of the
 * second contour . implementation of algorithm presented in fuchs
 * (1977) for both contours open
 */
path **twoPathsNoclosure(graph *G);

/**
 * find the paths between two bounding paths
 */
void pathsBetween(graph *G, path **paths, int botPath, 
		  int topPath, int length, int close);

/**
 * run dijstra's single source shortest path algorithm 
 */
path *singlePath(graph *G, int source, int length, int close);

/**
 * relax function 
 */
double relax(graph *G, edge *u, double wu, edge *v, double wv, int dir);
  //listNode *relax(graph *G, listNode* u, listNode *v, int dir);

/**
 * build an actual edge list from a path in the planar graph
 */
list *makeEdgeList(path *p);

/**
 * builds a path from a graph
 */
path *buildPath(graph *G, int source, int length, int close);

/**
 * frees the contents of a facelist
 */
void freeFaceList(list *faceList);

/**
 * frees the contents of a path
 */
void freePath(path *p);

/**
 * finds the index of a node in a heap corresponding to the lower
 * neighbor of a given node
 */
int findAdj(graph *G, list *vertexHeap, edge *u, int dir);

/** 
 * tests that a contour is not refered to in the contour adjacency
 * list of any of a list of contours
 */
int hasBackwardAdjacentContours(list *contourList, contour *cont);

/**
 * gets the cost of a face made of three vertices
 */
double getFaceCost(vertex *a, vertex *b, vertex *c);

/**
 * gets the irregularity of the angles in the face, a measure that goes from
 * 0 to one
 */
double getFaceAngleIrregularity(vertex *a, vertex *b, vertex *c);

/**
 * gets the area of a face made of three vertices
 */
inline double getFaceArea(vertex *a, vertex *b, vertex *c);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/tile.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
