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

#include "surfTopology.h"

double vertexTranslateRatio = 0.1;

/**
 * prepare the surface for topological processing
 */
void prepareSurf(surface *surf) {
  if(surf == NULL) return;

  int i, V = surf->maxVertices;

  // delete any fields already allocated
  unprepareSurf(surf);

  // prepare the surface fields

  // NTaV
  surf->NTaV = (int*) malloc(V*sizeof(int));
  for(i = 0; i < V; i++) {
    surf->NTaV[i] = 0;
  }
  surf->maxNTaV = 0;

  // LTaV
  surf->LTaV = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->LTaV[i] = NULL;
  }

  // v1TaV
  surf->v1TaV = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->v1TaV[i] = NULL;
  }

  // v2TaV
  surf->v2TaV = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->v2TaV[i] = NULL;
  }

  // NNoV
  surf->NNoV = (int*) malloc(V*sizeof(int));

  for(i = 0; i < V; i++) {
    surf->NNoV[i] = 0;
  }

  // Neighbors
  surf->Neighbors = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->Neighbors[i] = NULL;
  }

  // number of triangles at edge
  surf->NTaE = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->NTaE[i] = NULL;
  }

  // LT1aNE
  surf->LT1aNE = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->LT1aNE[i] = NULL;
  }

  // LT2aNE
  surf->LT2aNE = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->LT2aNE[i] = NULL;
  }

  // vertex manifoldness
  surf->vertexManifoldness = (int*) malloc(V*sizeof(int));
  for(i = 0; i < V; i++) {
    surf->vertexManifoldness[i] = VERT_MANIFOLD;
  }

  // surface manifoldness
  surf->manifoldness = SURF_UNKNOWN;

  // opposite triangles
  surf->OT = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->OT[i] = NULL;
  }

  // opposite triangles
  surf->OV = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    surf->OV[i] = NULL;
  }

  // connected components
  surf->CC = newList(LIST);
}

/**
 * unprepare a surface (just delete all topology check fields)
 */
void unprepareSurf(surface *surf) {
  if(surf == NULL) return;

  int i, V = surf->maxVertices;

  // prepare the surface fields

  // NTaV
  if(surf->NTaV != NULL) {
    free(surf->NTaV);
    surf->NTaV = NULL;
  }

  // LTaV
  if(surf->LTaV != NULL) {
    for(i = 0; i < V; i++) {
      if(surf->LTaV[i] != NULL) {
	free(surf->LTaV[i]);
      }
    }
    free(surf->LTaV);
    surf->LTaV = NULL;
  }

  // v1TaV
  if(surf->v1TaV != NULL) {
    for(i = 0; i < V; i++) {
      if(surf->v1TaV[i] != NULL) {
	free(surf->v1TaV[i]);
      }
    }
    free(surf->v1TaV);
    surf->v1TaV = NULL;
  }

  // v2TaV
  if(surf->v2TaV != NULL) {
    for(i = 0; i < V; i++) {
      if(surf->v2TaV[i] != NULL) {
	free(surf->v2TaV[i]);
      }
    }
    free(surf->v2TaV);
    surf->v2TaV = NULL;
  }

  // NNoV
  if(surf->NNoV != NULL) {
    free(surf->NNoV);
    surf->NNoV = NULL;
  }

  // Neighbors
  if(surf->Neighbors != NULL) {
   for(i = 0; i < V; i++) {
      if(surf->Neighbors[i] != NULL) {
	free(surf->Neighbors[i]);
      }
    }
   free(surf->Neighbors);
   surf->Neighbors = NULL;
  }

  // NTaE
  if(surf->NTaE != NULL) {
   for(i = 0; i < V; i++) {
      if(surf->NTaE[i] != NULL) {
	free(surf->NTaE[i]);
      }
    }
   free(surf->NTaE);
   surf->NTaE = NULL;
  }

  // LT1aNE
  if(surf->LT1aNE != NULL) {
   for(i = 0; i < V; i++) {
      if(surf->LT1aNE[i] != NULL) {
	free(surf->LT1aNE[i]);
      }
    }
   free(surf->LT1aNE);
   surf->LT1aNE = NULL;
  }

  // LT2aNE
  if(surf->LT2aNE != NULL) {
   for(i = 0; i < V; i++) {
      if(surf->LT2aNE[i] != NULL) {
	free(surf->LT2aNE[i]);
      }
    }
   free(surf->LT2aNE);
   surf->LT2aNE = NULL;
  }

  // vertex manifoldness
  if(surf->vertexManifoldness != NULL) {
    free(surf->vertexManifoldness);
    surf->vertexManifoldness = NULL;
  }

  // surface manifoldness
  surf->manifoldness = SURF_UNKNOWN;

  // opposite triangles
  if(surf->OT != NULL) {
   for(i = 0; i < V; i++) {
      if(surf->OT[i] != NULL) {
	free(surf->OT[i]);
      }
    }
   free(surf->OT);
   surf->OT = NULL;
  }

  // opposite vertices
  if(surf->OV != NULL) {
    for(i = 0; i < V; i++) {
      if(surf->OV[i] != NULL) {
	free(surf->OV[i]);
      }
    }
    free(surf->OV);
    surf->OV = NULL;
  }

  // connected components
  if(surf->CC != NULL) {
    freeList(surf->CC);
    surf->CC = NULL;
  }
}

/**
 * preprocess surface
 */
void preprocessSurf(surface *surf) {
  if(surf == NULL) return;

  prepareSurf(surf);

  trianglesAtVertex(surf);
  edgeProcess(surf);

  // if manifold, orient
  if(surf->manifoldness == SURF_MANIFOLD) {
    oppositeProcess(surf);
    consistentOrientation(surf);
    consistentFaces(surf);
  }

  buildConnectedComponents(surf);
}

/**
 * fix simple topological errors
 */
void topoFixer(surface *surf) {
  if(surf == NULL || surf->manifoldness == SURF_UNKNOWN) {
    return;
  }

  int i, V = surf->numVertices;

//  // find and fix two hole errors by duplicating a vertex
//  for(i = 0; i < V; i++) {
//    if(surf->vertexManifoldness[i] == VERT_TWO_HOLE_ERROR) {
//      fixTwoHoleError(surf,i);
//    }
//  }
//
//  preprocessSurf(surf);
//  V = surf->numVertices;
//
//  // find and fix cone errors by duplicating a vertex
//  for(i = 0; i < V; i++) {
//    if(surf->vertexManifoldness[i] == VERT_CONE_ERROR) {
//      fixConeError(surf,i);
//    }
//  }


  // find errors and fix by duplicating a vertex
  buildNTaE(surf);
  for(i = 0; i < V; i++) {
    if(surf->vertexManifoldness[i] != VERT_MANIFOLD) {
      cutNonmanifoldVertex(surf,i);
    }
  }

  //preprocessSurf(surf);

}

/**
 * get a list of boundaries, each of which is a list of vertices
 */
list *getBoundaries(surface *s) {
  if(s == NULL || s->manifoldness != SURF_MANIFOLD) {
    return NULL;
  }

// NOTES on algorithm:
//
// 1) we assume that every vertex is either interior w/
// NNoV - NTaV = 0, or a boundary vertex with NNoV - NTaV = 1.
//
// 2) consistentOrientation means Neighbors is such that a vertex is
// circumnavigated by its neighbors in a consistent way (consistent
// across vertices).
//
// For a boundary vertex n1, let the first item in its Neighbor-list be
// n2. Then for n2, n1 must be the LAST item in n2's Neighbor-list for
// the orientation to be consistent.
//
// 3) consistentFaces makes each triangle [n0 n1 n2], [n0 n2 n3],
// etc. where n0 is the source, n1-n3 are its ordered neighbors.
//
// So if, in a flattening, each triangle is ordered in an
// anti-clockwise manner, then the ordered boundary vertices will go
// around in an anti-clockwise manner.



// PREPROCESSING NOTE: surfStruct must be a MANIFOLD. If the boundary
// vertices are non-manifold, it is not possible to traverse them in
// order (robustly). This function therefore assumes a single component
// manifold, and uses the fact that the triangulation is oriented
// consistently. However, it could be generalized to work on
// non-manifold meshes for which the non-manifoldness occurs away from
// the boundary.
//
//checkSurface(surfStruct, mfilename, inputname(1), 'ccE', 'manifE');

  listNode *ln;
  list *boundaries = newList(LIST);
  list *curBoundary;
  long i, n0, n1, n2, ni, found, tmp;
  int *indivMap;

  listNode *comp;
  // go through each connected component, find its boundaries
  for(comp = getListNode(s->CC,0); comp; comp = (listNode*) comp->next) {
    // bN contains the labels of all the boundary vertices
    list *bN = newList(LIST);

    // get the vertex indices in this connected component
    indivMap = neighborhood(s,(long)comp->data);
    for(i = 0; i < s->numVertices; i++) {
      if(indivMap[i] != -1 && s->NNoV[i] - s->NTaV[i] == 1) {
	enqueue(bN,(void*)i);
      }
    }
    free(indivMap);

    if(listSize(bN) == 0) {
      continue;
    }

    while(listSize(bN) > 0) {

      ln = pop(bN);
      n0 = (long) ln->data;
      free(ln);

      curBoundary = newList(LIST);
      enqueue(curBoundary,(void*)n0);
      enqueue(boundaries,curBoundary);

      // sort the boundary vertices into ordered_bN so that they
      // traverse the boundary of the surface
      n1 = n0;
      while(listSize(bN) > 0) {
	found = 0;
	for(i = 0; i < s->NNoV[n1]; i++) {
	  tmp = s->Neighbors[n1][i];
	  ni = findInListI(bN,(void*) tmp);

	  if(ni > -1 && (s->LT1aNE[n1][i] == -1 || s->LT2aNE[n1][i] == -1)) {
	    ln = removeListNode(bN,ni);
	    n2 = (long) ln->data;
	    // add n2 to the ordered list
	    enqueue(curBoundary,(void*)n2);
	    free(ln);
	    n1 = n2;
	    found = 1;
	    break;
	  }
	}

	if(!found) {
	  break;
	}

	//n2 = s->Neighbors[n1][0];
	// Above, we assume a consistently oriented triangulation.
	// Otherwise, then we would have to search through the neighbors of
	// n_1 for another vertex on the boundary (but other than the one
	// prior to n_1)

	//printf("n0=%d n1=%d n2=%d\n",n0,n1,n2);

//	if(n2 == n0) {
//
//	  break;
//	}
//	else {
//	  // add n2 to the ordered list
//	  enqueue(curBoundary,(void*)n2);
//	  // eliminate n2 from the original unordered list
//	  free(removeListNode(bN, findInListI(bN,(void*)n2)));
//	  n1 = n2;
//	}
      }
    }

    freeList(bN);
  }

  return boundaries;
}

/**
 * computes the number and label of triangles for each vertex in a surface
 */
void trianglesAtVertex(surface *surf) {
  if(surf == NULL) return;

  int i, j, *counter;
  unsigned int V = surf->numVertices, F = surf->numFaces;
  int vert, vert1, vert2, vert3;
  unsigned int **f = surf->faces;

  int *NTaV = surf->NTaV;
  int **LTaV = surf->LTaV;
  int **v1TaV = surf->v1TaV;
  int **v2TaV = surf->v2TaV;

  /* First Pass through f: calc #triangles at each vertex */
  surf->maxNTaV = 0;
  for(i = 0; i < F; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  vert = f[i][j]; /* next vertex from adjacency list f */
	  NTaV[vert]++; /* increment NTaV[vert] */

	  if( NTaV[vert] > surf->maxNTaV ) /* keep track of max NTaV[vertex] */
	    surf->maxNTaV = NTaV[vert];    /* in maxNTaV */
	}
    }

  // allocate room for each triangle label list and v1 and v2 list
  counter = (int*) malloc(V*sizeof(int));
  for(i = 0; i < V; i++) {
    LTaV[i] = (int*) malloc(NTaV[i]*sizeof(int));
    v1TaV[i] = (int*) malloc(NTaV[i]*sizeof(int));
    v2TaV[i] = (int*) malloc(NTaV[i]*sizeof(int));
    counter[i] = 0;
  }

  /* Second pass through f: we know the number of triangles at */
  /* each vertex. Now pass through f */
  /* a triangle at a time, putting the label of the triangle at the */
  /* appropriate place in LTaV[], and at corresponding locations, */
  /* put the labels of the other 2 vertices of that triangle in */
  /* v1TaV[] and v2TaV[] */


  for(i = 0; i < F; i++)
    {
      /* get vertices from adjacency list, f, */
      /* 3 at a time, i.e. a triangle at a time */
      vert1 = f[i][0]; /* vert1 = f[3*i]; */
      vert2 = f[i][1]; /* vert2 = f[3*i + 1]; */
      vert3 = f[i][2]; /* vert3 = f[3*i + 2]; */

      /* Put triangle in LTaV, at the appropriate location  */
      /* given by Index[vert]. Put other 2 vertices in the */
      /* corresponding locations in v1TaV and v2TaV */
      /*                 vertex 1                            */
      LTaV [vert1][counter[vert1]] = i; /* i = triangle label */
      v1TaV[vert1][counter[vert1]] = vert2;
      v2TaV[vert1][counter[vert1]] = vert3;
      /*                 vertex 2                            */
      LTaV [vert2][counter[vert2]] = i; /* i = triangle label */
      v1TaV[vert2][counter[vert2]] = vert1;
      v2TaV[vert2][counter[vert2]] = vert3;
      /*                 vertex 3                            */
      LTaV [vert3][counter[vert3]] = i; /* i = triangle label */
      v1TaV[vert3][counter[vert3]] = vert1;
      v2TaV[vert3][counter[vert3]] = vert2;

      /* Increment index to "point" to the next open slot in LTaV */
      /* for each vertex */
      counter[vert1]++;
      counter[vert2]++;
      counter[vert3]++;
    }

  free(counter);

  return;
}

/**
 * processes the edges of a surface for manifoldness and neighbor relations
 * NNoV[V] = Number of Neighbors of vertex i
 *         = Number of Edges at Vertex i
 * Neighbors[V][N] = for each vertex v, this array
 *           contains the labels of its neighbors
 * LT1aNE[V][N] = Label of Triangles at Neigbor-Edge: (a Neigbor-Edge
 * LT2aNE[V][N]   is the edge from a vertex to one of its neighbors)
 *              for each such edge, these 2 arrays contain
 *              the labels of the triangles on either side of
 *              the edge (LT2aNE[] = -1 on a boundary/hole edge)
 */
void edgeProcess(surface *surf) {
//		 /* Inputs */
//		 /* Outputs */
//		 int *NNoV, int *E,
//		 int **Neighbors, int **LT1aNE, int **LT2aNE,
//		 int *NON_MANIFOLD,
//		 int *vertexManifoldness)
//{

  int i, j;
  unsigned int V = surf->numVertices;
  int *vertexList, *triangList, *secondTriList;
  int key, item, endOfList, FOUND;
  int soloIndex, vert, tri1, tri2;
  int NON_MANIFOLD_VERTEX;
  int counter;

  vertexList = (int *) malloc( 2*(surf->maxNTaV)*sizeof(int) );
  triangList = (int *) malloc( 2*(surf->maxNTaV)*sizeof(int) );
  secondTriList = (int *) malloc( 2*(surf->maxNTaV)*sizeof(int) );

  int *NTaV = surf->NTaV;
  int **LTaV = surf->LTaV;
  int **v1TaV = surf->v1TaV;
  int **v2TaV = surf->v2TaV;

  surf->E = 0;

//  int *NNoV = surf->NNoV;
//  int **Neighbors = surf->Neighbors;
//  int **LT1aNE = surf->LT1aNE;
//  int **LT2aNE = surf->LT2aNE;
//
//  int *vertexManifoldness = surf->vertexManifoldness;
  surf->manifoldness = SURF_MANIFOLD;

  /* Loop: fill up Neighbors[], LT1aNE[] and LT2aNE[], and in */
  /* the process, calc E, NNoV[] and Index2[] */
  for(i = 0; i < V; i++)
    {

      surf->NNoV[i] = 0;
      counter = 0;
      NON_MANIFOLD_VERTEX = 0;
      surf->vertexManifoldness[i] = VERT_MANIFOLD;

      /* ------------------------------------------------------ */
      /* Fill up vertexList & triangList w/ NTaV[i] pairs */
      /* (i.e. 2 times the number of triangles at vertex i). */
      /* For vertexList, each pair is label of the 2 vertices */
      /* other than i (from v1TaV and v2TaV), and for triangList, */
      /* the pair is simply the triangle label twice (from LTaV) */
      for(j = 0; j < NTaV[i]; j++)
	{
	  vertexList[counter] = v1TaV[i][j];
	  triangList[counter] = LTaV[i][j];
	  counter++;

	  vertexList[counter] = v2TaV[i][j];
	  triangList[counter] =  LTaV[i][j];
	  counter++;
	}
      /* e.g. vertexList = E B  A E  B G */
      /*      triangList = 3 3  1 1  2 2 */
      /* Done filling up vertexList and triangList */

      /* (HOMEWORK POINT 1) */

      /* ------------------------------------------------------ */
      /* vertexList will contain 2 vertex labels for any edge */
      /* that is a part of 2 triangles. Next, search for the */
      /* repeat of a given vertex label, and delete this repeat */
      /* from vertexList, after storing the label of the triangle */
      /* the repeated vertex was a part of. */

      /* This is done as follows: take the first vertex label in */
      /* vertexList as the "key". Search vertexList "items" until */
      /* either 1) "item" is found that matches "key" or 2) the end */
      /* of the list is encountered. */

      /* If 1), then copy the triangle label at "item" (i.e. the */
      /* triangle label of the repeated vertex label) into */
      /* secondTriList[key]. Then delete "item", by copying */
      /* vertexList[endOfList] and triangList[endOfList] into */
      /* vertexList[item] and triangList[item] respectively, and */
      /* then decrementing endOfList */

      /* If 2), then put a -1 in secondTriList[key], indicating */
      /* that no second triangle was found for the vertex whose */
      /* label is in vertexList[key] */

      endOfList = 2*NTaV[i] - 1;

      for(j = 0; j <= endOfList ; j++)
	secondTriList[j] = -1;

      /* soloIndex will store the index of a vertex with no repeats,
         or stay equal to zero if all vertices have repeats */
      soloIndex = 0;

      key = 0; /* endOfList = 2*NTaV[i] - 1; */
      while(key < endOfList) /* W1 */
	{
	  item = key + 1;
	  FOUND = 0;
	  while(item <= endOfList) /* W2 */
	    {

	      if(vertexList[item] == vertexList[key])
		{
		  /* Copy the last item on the list to "item", the
		     current position, and decrement "endOfList",
		     the list length */
		  secondTriList[key] = triangList[item];
		  triangList[item] = triangList[endOfList];
		  vertexList[item] = vertexList[endOfList];
		  endOfList--;

		  /* Because "item" now contains a new item (which
		     used to be at the end of the list), we
		     shouldn't increment "item" if we want to see
		     this item. So we decrement "item" here, as it
		     will be incremented later */
		  item--;

		  if(FOUND == 0)
		    { /* This is the first time item matches key */
		      FOUND = 1;
		    }
		  else
		    { /* This means the vertex in vertexList[key] has
			 occurred 3 or more times in the list:
			 i.e.  this edge is in 3 or more triangles
			 i.e.  this is a non-manifold edge */

		      NON_MANIFOLD_VERTEX = 1;
		      surf->vertexManifoldness[i] = VERT_MORE_THAN_TWO_TRIS_AT_EDGE_ERROR;
		      surf->manifoldness = SURF_NON_MANIFOLD;

		    }

		}

	      item++;

	    }/* return to W2 */

	  if(FOUND == 0)
	    soloIndex = key;

	  key++;

	}/* return to W1 */

      /* For the example above, we'll have: (unless NON_MANIFOLD_VERTEX = 1)*/
      /* e.g. vertexList = E  B  A  G   \                  */
      /*      triangList = 3  3  1  2    > endOfList = 3   */
      /*   secondTriList = 1  2 -1 -1   /                  */
      /*                    */
      /* Done removing repeats from vertexList and corresponding */
      /* triangList entries, and filling up secondTriList */

      /* (HOMEWORK POINT 2) */

      if(NON_MANIFOLD_VERTEX == 0)
	{
	  /* --------------------------------------------------------- */
	  /* Sort vertexList, triangList and secondTriList, so that */
	  /* they are in the order you would encounter them as */
	  /* you circled around vertex i, starting from a boundary/hole */
	  /* (if a boundary/hole exists). */

	  /* First, swap 0 and soloIndex: this way, if there is only */
	  /* one boundary/hole, it occurs before vertexList[0]         */
	  /*                         and after vertexList[endOfList] */
	  vert = vertexList[0];
	  tri1 = triangList[0];
	  tri2 = secondTriList[0];

	  vertexList[0] = vertexList[soloIndex];    /* this puts */
	  triangList[0] = secondTriList[soloIndex]; /* the -1 in */
	  secondTriList[0] = triangList[soloIndex]; /* triangList */

	  vertexList[soloIndex] = vert;
	  triangList[soloIndex] = tri2;
	  secondTriList[soloIndex] = tri1;

	  /* Next, do the sort */
	  key = 0;
	  while(key < endOfList) /* W1 */
	    {
	      item = key + 1;
	      FOUND = 0;
	      while( (item <= endOfList) && !FOUND ) /* W2 */
		{
		  if(triangList[item] == secondTriList[key])
		    {
		      FOUND = 1;

		      vert = vertexList[key+1];
		      tri1 = triangList[key+1];
		      tri2 = secondTriList[key+1];

		      vertexList[key+1] = vertexList[item];
		      triangList[key+1] = triangList[item];
		      secondTriList[key+1] = secondTriList[item];

		      vertexList[item] = vert;
		      triangList[item] = tri1;
		      secondTriList[item] = tri2;
		    }

		  if(secondTriList[item] == secondTriList[key])
		    {
		      FOUND = 1;

		      vert = vertexList[key+1];
		      tri1 = triangList[key+1];
		      tri2 = secondTriList[key+1];

		      vertexList[key+1] = vertexList[item];
		      triangList[key+1] = secondTriList[item];
		      secondTriList[key+1] = triangList[item];

		      vertexList[item] = vert;
		      triangList[item] = tri2;
		      secondTriList[item] = tri1;
		    }

		  item++;
		}/* return to W2 */

	      if(!FOUND)
		{/* This means that we couldn't find a triangle in
		    triangList or secondTriList that matched
		    secondTriList[key]. This implies that i is a
		    (CONE-ERROR) non-manifold vertex */

		  NON_MANIFOLD_VERTEX = 2; /* Cone error */
		  surf->vertexManifoldness[i] = VERT_CONE_ERROR;
		  surf->manifoldness = SURF_NON_MANIFOLD;

		}

	      key++;

	    }/* return to W1 */

	  /* For the example above, we'll have:                */
	  /* e.g. vertexList =  G  B  E  A   \                  */
	  /*      triangList = -1  2  3  1    > endOfList = 3   */
	  /*   secondTriList =  2  3  1 -1   /                  */
	  /*                   */
	  /* Note the "/////" ordering:              */
	  /*      secondTriList[k] = triangList[k+1] */
	  /*                                         */
	  /* Done sorting vertexList, triangList and secondTriList */

	} /* Matches if NON_MANIFOLD_VERTEX == 0 */

      /* (HOMEWORK POINT 3) */

      /* ------------------------------------------------------*/
      surf->NNoV[i] = endOfList + 1;
      surf->E += surf->NNoV[i];

      if(NON_MANIFOLD_VERTEX == 0)
	{
	  if(surf->NNoV[i] - NTaV[i] >= 2)
	    {
	      NON_MANIFOLD_VERTEX = 3; /* 2-hole error */
	      surf->vertexManifoldness[i] = VERT_TWO_HOLE_ERROR;
	      surf->manifoldness = SURF_NON_MANIFOLD;
	    }
	}

      /* (HOMEWORK POINT 4) */

      /* Paste the sorted vertexList for vertex i into Neighbors */
      /* and triangList and secondTriList into LT1aNE and LT2aNE */
      /* respectively. Index2[i] tells you where the information */
      /* for vertex i begins, and NNoV[i] tells you the "length" */
      /* of the info */
      surf->Neighbors[i] = (int*) malloc(surf->NNoV[i]*sizeof(int));
      surf->NTaE[i] = (int*) malloc(surf->NNoV[i]*sizeof(int));
      surf->LT1aNE[i] = (int*) malloc(surf->NNoV[i]*sizeof(int));
      surf->LT2aNE[i] = (int*) malloc(surf->NNoV[i]*sizeof(int));
      for(j = 0; j <= endOfList; j++)
	{
	  surf->Neighbors[i][j] =  vertexList[j];
	  if(NON_MANIFOLD_VERTEX != 1)
	    {
	      surf->LT1aNE[i][j] =  triangList[j];
	      surf->LT2aNE[i][j] =  secondTriList[j];
	    }
	  else
	    {
	      surf->LT1aNE[i][j] =  -1;
	      surf->LT2aNE[i][j] =  -1;
	    }


	}

    }

  surf->E /= 2;

  free(vertexList);
  free(triangList);
  free(secondTriList);

  return;
}

/* Notes: */
/* V = number of vertices */
/* F = number of faces (i.e. polygons) */
/* E = number of edges */

/* NTaV[V]  = Number of Triangles at Vertex i */
/* LTaV[3F] = Labels of Triangles at Vertex i */
/*   The Triangle at Vertex i has three vertices: i, v1 and v2. */
/* v1TaV[3F] and v2TaV[3F] contain the vertices v1 and v2 */
/* Index[V] is used to index into LTaV, v1TaV, v2TaV */

/* NNoV[V] = Number of Neighbors of vertex i */
/*         = Number of Edges at Vertex i */
/* Neighbors[2E] = for each vertex v, this array */
/*           contains the labels of its neighbors */
/* LT1aNE[2E] = Label of Triangles at Neigbor-Edge: (a Neigbor-Edge */
/* LT2aNE[2E]   is the edge from a vertex to one of its neighbors) */
/*              for each such edge, these 2 arrays contain */
/*              the labels of the triangles on either side of */
/*              the edge (LT2aNE[] = -1 on a boundary/hole edge) */
/* Index2[V] is used to index into Neighbors, LT1aNE, LT2aNE */

/* For vertex = i, triangle = LT2aNE[ Index2[i] + j] */
/* OT[2E]: Opposite Triangle  = OT[ Index2[i] + j] */
/* OV[2E]: Opposite Vertex = OV[ Index2[i] + j] */
void oppositeProcess(surface *surf) {
		     /* Inputs */
  if(surf == NULL) return;
  int i, j, V = surf->numVertices;
  int t, vA, vB;
  int index2, indexPlus1, indexMinus1;

  int *NNoV = surf->NNoV;
  int **Neighbors = surf->Neighbors;
  int **LT1aNE = surf->LT1aNE;
  int **LT2aNE = surf->LT2aNE;
  int **OT = surf->OT;
  int **OV = surf->OV;

  /* Given Neighbors[], LT1aNE[], LT2aNE[] with Index2[] and NNoV[] to
     index into them, we can use this info to find opposite triangles
     and vertices and therefore fill in OT[] and OV[] as follows: */

  /* Given index = Index2[i] + j, t = LT2aNE[index] is the triangle
     for some vertex i. The other two vertices in t are given by vA =
     Neigbors[index] and vB = Neigbors[index+1]. */

  /* Go to Neighbors and find the edge (vA, vB), i.e. in Neighbors,
     starting at Index2[vA], search no more than NNoV[vA] entries, for
     vB. This is done by the function findIndex, which returns
     index2. */

  /* At the location of vB in Neighbors (index2), look at
     LT1aNE[index2]: this is one of the two triangles straddling the
     edge (vA,vB). LT1aNE[index2] will either be the initial triangle
     t, or the opposite triangle. */

  /* If LT1aNE[index2] == t, then tOpposite is in LT2aNE[index2], and
     as the order of traversal of the triangles is {t, tOpposite}, the
     order of traversal of the neighbors must be {vB, vOpposite}
     (remember that this version of OppositeProcess assumes that
     Neighbors, LT1aNE, and LT2aNE are sorted). --> vOpposite =
     Neighbors[ mod(index2 + 1) ]. */

  /* If LT1aNE[index2] != t, then tOpposite is in LT1aNE[index2], and
     as the order of traversal of the triangles is {tOpposite, t}, the
     order of traversal of the neighbors must be {vOpposite, vB}
     (remember that this version of OppositeProcess assumes that
     Neighbors, LT1aNE, and LT2aNE are sorted). --> vOpposite =
     Neighbors[ mod(index2 - 1) ]. */

  for(i = 0; i < V; i++)
    {
      OT[i] = (int*) malloc(NNoV[i]*sizeof(int));
      OV[i] = (int*) malloc(NNoV[i]*sizeof(int));
      for(j = 0; j < NNoV[i]; j++)
	{
	  t = LT2aNE[i][j];
	  if(t == -1)
	    {
	      OT[i][j] = -1;
	      OV[i][j] = -1;
	    }
	  else /* e1 */
	    {
	      vA = Neighbors[i][j];
	      /* if index is the last of i's entries, wrap */
	      /* around to the first of i's entries for vB */
	      if(j == NNoV[i] - 1)
		indexPlus1 = 0;
	      else
		indexPlus1 = j + 1;
	      vB = Neighbors[i][indexPlus1];

	      index2 = findIndex(Neighbors[vA], NNoV[vA], vB);

	      if(LT1aNE[vA][index2] == t)
		{
		  OT[i][j] = LT2aNE[vA][index2];

		  if(OT[i][j] == -1)
		    OV[i][j] = -1;
		  else
		    {
		      /* if index2 is the last of vA's entries, */
		      /* wrap around to the first of vA's entries */
		      if(index2 == NNoV[vA] - 1)
			indexPlus1 = 0;
		      else
			indexPlus1 = index2 + 1;

		      OV[i][j] = Neighbors[vA][indexPlus1];
		    }
		}
	      else /* i.e. LT1aNE[index2] != t */
		{
		  OT[i][j] = LT1aNE[vA][index2];

		  if(OT[i][j] == -1)
		    OV[i][j] = -1;
		  else
		    {
		      /* if index2 is the first of vA's entries, */
		      /* wrap around to the last of vA's entries */
		      if(index2 == 0)
			indexMinus1 = NNoV[vA] - 1;
		      else
			indexMinus1 = index2 - 1;

		      OV[i][j] = Neighbors[vA][indexMinus1];
		    }
		}/* Matches else for if tri == t */

	    }/* Matches e1: else for if t == -1 */

	}/* Matches for(j = 0; j < NNoV[i]; j++) */

    }/* Matches for(i = 0; i < V; i++) */

  return;

}


/*

  This reorders Neighbors[] so that 1-neighbors are traversed in a
  consistent (e.g. clockwise) order for each node.

  This requires also reordering LT1aNE[], LT2aNE[], OT[] and OV[]

  Inputs from matlab:
  Index2[V] - 1
  NNoV[V]
  NTaV[V]
  Neighbors[] - 1
  LT1aNE[] - 1
  LT2aNE[] - 1
  OT[] - 1
  OV[] - 1

  Outputs:
  Neighbors_out[]
  LT1aNE_out[]
  LT2aNE_out[]
  OT_out[]
  OV_out[]

  From matlab call with:
  [Neighbors, LT1aNE, LT2aNE, OT, OV] = ...
  consistentOrientation(Index2-1, NNoV, NTaV, ...
  Neighbors-1, LT1aNE-1, LT2aNE-1, OT-1, OV-1);

*/

void consistentOrientation(surface *surf) {
  if(surf == NULL) return;
  int i, j, V = surf->numVertices;
  int indexPlus1, index2, index3, index4;
  int nAft, nBef;
  int nRightOfnAft=-1, nRightOfnBef=-1;

  int *Q, Qstart, Qend, Qempty, *onQ;
  int *sense = NULL;

  int n0 = 0;

  /* Inputs */
  int *NNoV = surf->NNoV;
  int *NTaV = surf->NTaV;
  int **Neighbors = surf->Neighbors;
  int **LT1aNE = surf->LT1aNE;
  int **LT2aNE = surf->LT2aNE;
  int **OT = surf->OT;
  int **OV = surf->OV;
  /* Outputs */
  int **Neighbors_out;
  int **LT1aNE_out;
  int **LT2aNE_out;
  int **OT_out;
  int **OV_out;


  /* Initialize Q */
  Q = (int*) malloc( V*sizeof(int) );
  Qstart = 0; Qend = 0; Qempty = 1;

  onQ = (int*) malloc( V*sizeof(int) );
  for(i = 0; i < V; i++)
    onQ[i] = 0;
  onQ[ n0 ] = 1;

  Q[ Qend ] = n0;
  Qend++;
  Qempty = 0;

  sense = (int*) malloc( V*sizeof(int) );
  for(i = 0; i < V; i++)
    sense[i] = 0;
  sense[n0] = 1;

  while( !Qempty )
    {
      /* Take the information from the item at Qstart */
      nBef = Q[ Qstart ];

      /* Update Qstart */
      Qstart++;

      for(i = 0; i < NNoV[nBef]; i++)
	{
	  nAft = Neighbors[nBef][i];

	  /* if nAft:
	     1) is not on the Q, i.e. has no sense
	     then:
	     1) put it on the Q
	     2) figure out its sense
	  */

	  if( onQ[nAft] == 0 )
	    {

	      /* 1) Put nAft on the Q */
	      onQ[ nAft ] = 1;
	      Q[ Qend ] = nAft;
	      Qend++;

	      /* 2) Figure out the sense of nAft */

	      /* Go around nAft to find nBef, and
		 the nodes adjacent to it */
	      for(j = 0; j < NNoV[nAft]; j++)
		{
		  if(j < NNoV[nAft] - 1)
		    indexPlus1 = j + 1;
		  else
		    indexPlus1 = 0;

		  if( Neighbors[nAft][j] == nBef )
		    nRightOfnAft = Neighbors[nAft][indexPlus1];

		}

	      /* Go around nBef to find nAft, and
		 the nodes adjacent to it */
	      for(j = 0; j < NNoV[nBef]; j++)
		{
		  if(j < NNoV[nBef] - 1)
		    indexPlus1 = j + 1;
		  else
		    indexPlus1 = 0;

		  if( Neighbors[nBef][j] == nAft )
		    nRightOfnBef = Neighbors[nBef][indexPlus1];

		}

	      /* Figure out sense[nAft] */
	      if( nRightOfnAft != nRightOfnBef )
		sense[nAft] =  sense[nBef];
	      else
		sense[nAft] = -sense[nBef];

	      /* If we have a boundary node nBef propagating to another
		 boundary node nAft, we have to be careful: what we did
		 above to figure out the sense of nAft relative to the
		 sense of nBef is not guaranteed to work.

		 If nAft is neither the FIRST nor the LAST node on
		 nBef's neighbor list, do nothing. (i.e. sense is
		 correct). Otherwise:

		 The (boundary) nodes will have the same sense if: nAft
		 is the FIRST node on nBef's neighbor list AND nBef is
		 the LAST node on nAft's neighbor list.

		 Or if: nAft is the LAST node on nBef's neighbor list
		 AND nBef is the FIRST node on nAft's neighbor list.  */
	      if( NNoV[nBef] != NTaV[nBef] )
		if( NNoV[nAft] != NTaV[nAft] )
		  { /* Boundary node propagating to a boundary node. */

		    if( Neighbors[nBef][0] == nAft )
		      { /* nAft == FIRST on nBef's NeighborList */
			index2 = NNoV[nAft] - 1;
			if( Neighbors[nAft][index2] == nBef )
			  /* nBef == LAST on nAft's neighbor list. */
			  sense[nAft] =  sense[nBef];
			else
			  sense[nAft] = -sense[nBef];
		      }

		    if( Neighbors[nBef][NNoV[nBef] - 1] == nAft )
		      { /* nAft == LAST node on nBef's NeighborList */
			index2 = 0;
			if( Neighbors[nAft][index2] == nBef )
			  /* nBef == FIRST on nAft's neighbor list. */
			  sense[nAft] =  sense[nBef];
			else
			  sense[nAft] = -sense[nBef];
		      }

		  }/* matches Boundary node propagating to a boundary node. */

	    }/* matches if( onQ[nAft] == 0 ) */

	}/* matches for(i = 0; i < NNoV[nBef]; i++) */

      /* Check for Q empty condition */
      if( Qend == Qstart )
	Qempty = 1;

    }/* while( !Qempty ) */

  /* Now that we know the sense of all the nodes (relative to the n0 =
     node 0), reorder the list of Neighbors, LT1aNE, LT2aNE, OT, OV so
     that everything has a consistent orientation */
  V = surf->maxVertices;
  Neighbors_out = (int**) malloc(V*sizeof(int*));
  LT1aNE_out = (int**) malloc(V*sizeof(int*));
  LT2aNE_out = (int**) malloc(V*sizeof(int*));
  OT_out = (int**) malloc(V*sizeof(int*));
  OV_out = (int**) malloc(V*sizeof(int*));
  for(i = 0; i < V; i++) {
    Neighbors_out[i] = (int*) malloc(NNoV[i]*sizeof(int));
    LT1aNE_out[i] = (int*) malloc(NNoV[i]*sizeof(int));
    LT2aNE_out[i] = (int*) malloc(NNoV[i]*sizeof(int));
    OT_out[i] = (int*) malloc(NNoV[i]*sizeof(int));
    OV_out[i] = (int*) malloc(NNoV[i]*sizeof(int));
    if( sense[i] == -1 )
      for(j = 0; j < NNoV[i]; j++)
	{
	  index2 = NNoV[i] - 1 - j;
	  Neighbors_out[i][index2] = Neighbors[i][j];

	  if( j == 0 )
	    index3 = 0;
	  else
	    index3 = index2 + 1;
	  LT1aNE_out[i][index3] = LT1aNE[i][j];

	  if( j == NNoV[i]-1 )
	    index4 = j;
	  else
	    index4 = index2 - 1;
	  LT2aNE_out[i][index4] = LT2aNE[i][j];
	  OT_out[i][index4] = OT[i][j];
	  OV_out[i][index4] = OV[i][j];

	}
    else
      for(j = 0; j < NNoV[i]; j++)
	{
	  Neighbors_out[i][j] = Neighbors[i][j];
	  LT1aNE_out[i][j] = LT1aNE[i][j];
	  LT2aNE_out[i][j] = LT2aNE[i][j];
	  OT_out[i][j] = OT[i][j];
	  OV_out[i][j] = OV[i][j];

	}

    free(Neighbors[i]);
    free(LT1aNE[i]);
    free(LT2aNE[i]);
    free(OT[i]);
    free(OV[i]);
  }

  free(Neighbors);
  free(LT1aNE);
  free(LT2aNE);
  free(OT);
  free(OV);

  surf->Neighbors = Neighbors_out;
  surf->LT1aNE = LT1aNE_out;
  surf->LT2aNE = LT2aNE_out;
  surf->OT = OT_out;
  surf->OV = OV_out;

  free(Q);
  free(onQ);
  free(sense);

}

/*

  Given an ordered Neighbors[] and LT2aNE[], this produces an f
  (i.e. Fx3 face array) such that the vertices are ordered
  consistently.

  Inputs from matlab:
  NTaV[V]
  NNoV[V]
  Neighbors[] - 1
  LT2aNE[] - 1
  Index2[V] - 1
  F

  Outputs:
  f[]

  From matlab call with:
  f = consistentFaces(NTaV, NNoV, Neighbors-1, LT2aNE-1, Index2-1, F)

*/

void consistentFaces(surface *surf) {

  if(surf == NULL) return;
  int V = surf->numVertices;
  int n0, j, t, nA, nB;
  int indexPlus1;

  /* Inputs */
  int *NTaV = surf->NTaV;
  int *NNoV = surf->NNoV;
  int **Neighbors = surf->Neighbors;
  int **LT2aNE = surf->LT2aNE;

  for(n0 = 0; n0 < V; n0++)
    {

      for(j = 0; j < NTaV[n0]; j++)
	{
	  nA = Neighbors[n0][j];
	  /* if index is the last of n0's entries, wrap */
	  /* around to the first of n0's entries for nB */
	  if(j < NNoV[n0] - 1)
	    indexPlus1 = j + 1;
	  else
	    indexPlus1 = 0;
	  nB = Neighbors[n0][indexPlus1];

	  t  = LT2aNE[n0][j];
	  surf->faces[t][0] = (unsigned int) nA;
	  surf->faces[t][1] = (unsigned int) nB;
	  surf->faces[t][2] = (unsigned int) n0;
	}

    }

}

int findIndex(int *Array,
	      int segmentLength,
	      int key)
{
  int i;

  for(i = 0; i < segmentLength; i++)
    {
      if( (int) Array[i] == key )
	return i;
    }

  return(-1);

}

/**
 * implementation of the gueziec (1998) local surface cutting method
 * for topology fixing
 */
void cutNonmanifoldVertex(surface *s, int v) {
  if(s == NULL || v < 0 || v > s->numVertices) return;

  list *star = buildStar(s,v);

  list *compList = newList(LIST);
  list *curComp;
  listNode *ln, *ln2;
  staredge *se;

  int ends[2], i, si, ei, f, fi, vi;

  // decompose star into n connected components
  while(listSize(star) > 0) {
    ln = pop(star);
    se = ln->data;
    free(ln);

    ends[0] = se->edge[0];
    ends[1] = se->edge[1];

    curComp = newList(LIST);
    enqueue(compList,curComp);
    enqueue(curComp,se);

    // search the remaining edges for the current ends
    i = 0;
    ln = getListNode(star,0);
    while(ln) {
      se = (staredge*) ln->data;

      // try to match an end to an edge node
      if(se->edge[0] == ends[0]) {
	si = 0;
	ei = 1;
      }
      else if(se->edge[0] == ends[1]) {
	si = 0;
	ei = 0;
      }
      else if(se->edge[1] == ends[0]) {
	si = 1;
	ei = 1;
      }
      else if(se->edge[1] == ends[1]) {
	si = 1;
	ei = 0;
      }
      else {
	ln = (listNode*) ln->next;
	i++;
	continue;
      }

      // make sure the connecting edge is a manifold edge
      if(s->NTaE[v][se->neighborIndex[si]] < 3) {
	ends[si] = se->edge[ei];
	enqueue(curComp,se);
	removeListNode(star,i);

	if(ends[0] == ends[1]) {
	  ln = NULL;
	}
	else {
	  i = 0;
	  ln = getListNode(star,0);
	}
      }
      else {
	ends[si] = -1;
	ln = (listNode*) ln->next;
	i++;
      }

    }
  }

  // duplicate vertex n-1 times
  for(ln = (listNode*) getListNode(compList,0)->next; ln; ln = (listNode*) ln->next) {
    vi = s->numVertices;

    // duplicate the vertex and adjust the face references
    addVertexArr(s,s->vertices[v]);

    curComp = (list*) ln->data;
    for(ln2 = getListNode(curComp,0); ln2; ln2 = (listNode*) ln2->next) {
      f = ((staredge*) ln2->data)->face;
      for(fi = 0; fi < 3; fi++) {
	if(s->faces[f][fi] == v) {
	  s->faces[f][fi] = vi;
	}
      }
    }

    // update topo-fields in surf struct
    preprocessSurf(s);

    freeListAndData(curComp);
  }

  free(star);
  freeList(compList);

}

/**
 * build the star for a vertex, a la gueziec (1998)
 */
list *buildStar(surface *s, int v) {
  if(s == NULL || v < 0 || v > s->numVertices) return NULL;

  list *st = newList(LIST);
  listNode *ln, *ln2;
  staredge *se, *se2;
  int i,j,k,n,nn,found;

  // loop over neighbors, adding edges
  for(i = 0; i < s->NNoV[v]; i++) {
    n = s->Neighbors[v][i];
    // search the neighbors of this neighbor for neighbors of this vertex
    for(j = 0; j < s->NNoV[n]; j++) {
      nn = s->Neighbors[n][j];
      if(nn == v) {
	continue;
      }

      // find nn in v's Neighbors
      found = 0;
      for(k = 0; k < s->NNoV[v]; k++) {
	if(nn == s->Neighbors[v][k]) {
	  found = 1;
	  break;
	}
      }

      if(found) {
	// determine the face this edge is in with v

	// try to find n in v1TaV and nn in v2TaV or vice versa
	for(k = 0; k < s->NTaV[v]; k++) {
	  if((s->v1TaV[v][k] == n && s->v2TaV[v][k] == nn)
	     || (s->v2TaV[v][k] == n && s->v1TaV[v][k] == nn)) {

	    // add the edge, delete repeats later
	    se = (staredge*) malloc(sizeof(staredge));
	    se->edge[0] = n;
	    se->edge[1] = nn;
	    se->neighborIndex[0] = i;
	    se->neighborIndex[1] = k;
	    se->face = s->LTaV[v][k];

	    enqueue(st,se);
	    break;
	  }
	}
      }
    }
  }

  // delete repeated edges
  for(ln = getListNode(st,0); ln; ln = (listNode*) ln->next) {
    if(ln->delete) continue;
    se = (staredge*) ln->data;
    for(ln2 = (listNode*) ln->next; ln2; ln2 = (listNode*) ln2->next) {
      se2 = (staredge*) ln2->data;

      // delete if they represent the same edge
      if(((se->edge[0] == se2->edge[0]
	   && se->edge[1] == se2->edge[1])
	  ||
	  (se->edge[0] == se2->edge[1]
	   && se->edge[1] == se2->edge[0]))
	 && se->face == se2->face) {
	markForDeletion(ln2);
      }
    }
  }
  deleteMarkedNodes(st);

  return st;
}

/**
 * build a map of the number of triangles at each vertex
 * note that this info is only needed for topofixing
 */
void buildNTaE(surface *s) {
  if(s == NULL) return;

  int p[3][3] = {{0,1},{1,2},{0,2}};
  int f,i,v1,v2,n;
  for(f = 0; f < s->numFaces; f++) {
    for(i = 1; i < 3; i++) {
      v1 = s->faces[f][p[i][0]];
      v2 = s->faces[f][p[i][1]];

      // find v2 in v1's neighbor list
      for(n = 0; n < s->NNoV[v1]; n++) {
	if(s->Neighbors[v1][n] == v2) {
	  s->NTaE[v1][n]++;
	}
      }

      // find v1 in v2's neighbor list
      for(n = 0; n < s->NNoV[v2]; n++) {
	if(s->Neighbors[v2][n] == v1) {
	  s->NTaE[v2][n]++;
	}
      }
    }
  }
}

/** tools to build maps of connection components **/

/**
 * map the connected components of the surface by building a list of
 * vertices in each connected component
 */
void buildConnectedComponents(surface *surf) {
  if(surf == NULL || surf->numVertices < 1
     || surf->manifoldness == SURF_UNKNOWN) {
    return;
  }

  long i,lastFind,nextFind,done,V = surf->numVertices;
  int *foundCC = (int*) malloc(V*sizeof(int));
  int *indivMap;

  /* build a found vector to keep track of the vertices for which CC is known*/
  for(i = 0; i < V; i++) {
    foundCC[i] = 0;
  }

  /* get the map for successive unfound vertices until all are found */
  nextFind = 0;
  done = 0;
  while(!done) {
    indivMap = neighborhood(surf,nextFind);

    /* add this vertex to CC list */
    enqueue(surf->CC, (void*) nextFind);

    /* update foundCC */
    lastFind = nextFind;
    for(i = 0; i < V; i++) {
      if(foundCC[i] == 0 && indivMap[i] != -1) {
	foundCC[i] = 1;
      }
      else if(foundCC[i] == 0 && lastFind == nextFind) {
	nextFind = i;
      }
    }

    if(nextFind == lastFind) {
      done = 1;
    }

    free(indivMap);
  }

  free(foundCC);
}

/**
 * compute the shortest number of edges from a siingle vert to all others
 * modified from the surface geometry toolbox
 */
int *neighborhood(surface *surf, int n0) {
  if(surf == NULL || surf->numVertices < 1
     || surf->manifoldness == SURF_UNKNOWN) {
    return NULL;
  }

  int i, j, nB, V = surf->numVertices;
  int *Q, Qstart = 0, Qend = 0, Qempty = 1;

  /* Initiaize nbd */
  int *nbd = (int*) malloc(V*sizeof(int));
  for(i = 0; i < V; i++)
    nbd[i] = -1;
  nbd[n0] = 0;

  /* Initialize Q */
  Q = (int*) malloc( V*sizeof(int) );

  /* Put n0 on the Q */
  Q[ Qend ] = n0;
  Qend++;
  Qempty = 0;

  while( !Qempty )
    {
      /* pop vertex i from Q */
      i = Q[ Qstart ];
      Qstart++;

      /* Go all the way around i */
      for(j = 0; j < surf->NNoV[i]; j++)
	{
	  nB = surf->Neighbors[i][j];

	  /* if nB has not been hit */
	  if( nbd[nB] == -1 )
	    {
	      /* 1) set nbd of vertex nB */
	      nbd[nB] = nbd[i] + 1;

	      /* 2) put vertex nB on Q */
	      Q[ Qend ] = nB;
	      Qend++;
	      Qempty = 0;

	    } /* Matches: if nB has not been hit */

	} /* Gone all the way around i */

      /* Check for Q empty condition */
      if( Qend == Qstart )
	Qempty = 1;

    } /* Matches while( !Qempty ) */


  free(Q);

  return nbd;
}

/**
 * fix a two hole topological error
 */
void fixTwoHoleError(surface *s, int v) {
  if(s == NULL || v < 0 || v > s->numVertices) return;

  int *neighList = (int*) malloc(s->NNoV[v]*sizeof(int));
  int *triList = (int*) malloc(s->NNoV[v]*sizeof(int));
  int vi, f, fi, n, orign, done, count, curNeigh, curTri;
  for(n = 0; n < s->NNoV[v] && s->LT1aNE[v][n] != -1; n++);

  orign = n;
  done = 0;
  count = 0;
  while(!done) {
    curNeigh = curTri = 0;
    neighList[curNeigh++] = s->Neighbors[v][n];

    n++;
    if(n == s->NNoV[v]) {
      n = 0;
    }
    if(n == orign) {
      done = 1;
    }

    while(s->LT1aNE[v][n] != -1) {
      neighList[curNeigh++] = s->Neighbors[v][n];
      triList[curTri++] = s->LT1aNE[v][n];

      n++;
      if(n == s->NNoV[v]) {
	n = 0;
      }
      if(n == orign) {
	done = 1;
      }

    }

    vi = v;
    if(count > 0) {
      vi = s->numVertices;

      // duplicate the vertex and adjust the face references
      addVertexArr(s,s->vertices[v]);

      for(f = 0; f < curTri; f++) {
	for(fi = 0; fi < 3; fi++) {
	  if(s->faces[triList[f]][fi] == v) {
	    s->faces[triList[f]][fi] = vi;
	  }
	}
      }
    }


    // move the vertex epsilon toward these neighbors
    //translateVertexTowardNeighbors(s,vi,neighList,curNeigh);
    count++;
  }

  free(neighList);
  free(triList);
}

/**
 * fix a cone topological error
 */
void fixConeError(surface *s, int v) {
  if(s == NULL || v < 0 || v > s->numVertices) return;

  int *neighList = (int*) malloc(s->NNoV[v]*sizeof(int));
  int *triList = (int*) malloc(s->NNoV[v]*sizeof(int));
  int vi, f, fi, n, t1, done, count, curNeigh, curTri;

  n = 0;
  done = 0;
  count = 0;
  while(!done) {
    curNeigh = curTri = 0;

    t1 = s->LT1aNE[v][n];
    neighList[curNeigh++] = s->Neighbors[v][n];
    triList[curTri++] = s->LT2aNE[v][n];
    while(t1 != s->LT2aNE[v][n]) {
      n++;
      neighList[curNeigh++] = s->Neighbors[v][n];
      triList[curTri++] = s->LT2aNE[v][n];
    }

    vi = v;
    if(count > 0) {
      vi = s->numVertices;

      // duplicate the vertex and adjust the face references
      addVertexArr(s,s->vertices[v]);

      for(f = 0; f < curTri; f++) {
	for(fi = 0; fi < 3; fi++) {
	  if(s->faces[triList[f]][fi] == v) {
	    s->faces[triList[f]][fi] = vi;
	  }
	}
      }

      // move the vertex epsilon toward these neighbors
      //translateVertexTowardNeighbors(s,vi,neighList,curNeigh);
    }

    count++;

    n++;
    if(n == s->NNoV[v]) {
      done = 1;
    }
  }

}


/**
 * translated a given vertex towards its neighbors
 */
void translateVertexTowardNeighbors(surface *surf, int v, int *nList, int N) {
  int i,j;
  double mean[3];

  for(i = 0; i < 3; i++) {
    mean[i] = surf->vertices[nList[0]][i];
  }

  for(i = 1; i < N; i++) {
    for(j = 1; j < N; j++) {
      mean[j] += surf->vertices[nList[i]][j];
    }
  }

  for(i = 0; i < 3; i++) {
    mean[i]/=N;
    surf->vertices[v][i] += vertexTranslateRatio*(mean[i]-surf->vertices[v][i]);
  }
}
