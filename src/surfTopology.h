/******************************************************************************
 * surfTopology contains functions that identify and try to correct
 * topological errors in surfaces.
 *
 * this code is modified from that contained in the SurfaceGeometry toolbox
 * 
 * Mukund Balasubramanian and Oliver Hinds <oph@bu.edu> 2006-08-18 
 * 
 *
 *
 *****************************************************************************/

/* Notes: */
/* V = number of vertices */
/* F = number of faces (i.e. polygons) */
/* E = number of edges */

/* NTaV[V]  = Number of Triangles at Vertex i */
/* LTaV[3F] = Labels of Triangles at Vertex i */
/*   The Triangle at Vertex i has three vertices: i, v1 and v2. */
/* v1TaV[3F] and v2TaV[3F] contain the vertices v1 and v2 */
/* Index[V] is used to index into LTaV */

#include "libsrTypes.h"
#include "surfUtil.h"


/**
 * prepare the surface for topological processing
 */
void prepareSurf(surface *surf);

/** 
 * unprepare a surface (just delete all topology check fields)
 */
void unprepareSurf(surface *surf);

/**
 * preprocess a surface
 */
void preprocessSurf(surface *surf);

/**
 * fix simple topological errors
 */
void topoFixer(surface *surf);

/**
 * get a list of boundaries, each of which is a list of vertices
 */
list *getBoundaries(surface *s);

/**
 * computes the number and label of triangles for each vertex in a surface
 */
void trianglesAtVertex(surface *surf);

/**
 * processes the edges of a surface for manifoldness and neighbor relations
 */
void edgeProcess(surface *surf);

/**
 *
 */
void consistentOrientation(surface *surf);

/**
 *
 */
void oppositeProcess(surface *surf);

/**
 *
 */
void consistentFaces(surface *surf);

/**
 *
 */
int findIndex(int *Array, int segmentLength,  int key);

/**
 * implementation of the gueziec (1998) local surface cutting method
 * for topology fixing
 */
void cutNonmanifoldVertex(surface *s, int v);

/**
 * build the star for a vertex, a la gueziec (1998)
 */
list *buildStar(surface *s, int v);

/**
 * build a map of the number of triangles at each vertex
 * note that this info is only needed for topofixing
 */
void buildNTaE(surface *s);

/**
 * map the connected components of the surface by building a list of
 * vertices in each connected component
 */
void buildConnectedComponents(surface *surf);

/**
 * compute the shortest number of edges from a siingle vert to all others
 * modified from the surface geometry toolbox
 */
int *neighborhood(surface *surf, int n0);

/**
 * fix a two hole topological error
 */
void fixTwoHoleError(surface *s, int v);

/**
 * fix a cone topological error
 */
void fixConeError(surface *s, int v);

/**
 * translated a given vertex towards its neighbors
 */
void translateVertexTowardNeighbors(surface *surf, int v, int *nList, int N);
