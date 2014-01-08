/*****************************************************************************
 * surfUtil.c is the source file for surface utility functions for
 * libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#include"surfUtil.h"
#include"inhedron.h"

/* contour geometry magic numbers */
double thresholdForCircularity = 0.25;

/** surface creation and deletion functions **/

/**
 * create a surface structure
 * returns a surface structure or null if creation flibvps
 */
surface *createSurfaceDefault() {
  return createSurface(DEFAULT_MAX_VERTICES, DEFAULT_MAX_FACES);
}

/**
 * create a surface structure with room for a specified number of
 * vertices and faces
 * returns a surface structure or null if creation flibvps
 */
surface *createSurface(unsigned int maxVerts, unsigned int maxFaces) {
  int i;

  surface *surf = (surface*) malloc(sizeof(surface));
  /* validate */
  if(surf == NULL) return NULL;

  /* allocate vertices */
  surf->maxVertices = maxVerts;
  surf->numVertices = 0;
  surf->vertices = (double**) malloc(surf->maxVertices*sizeof(double*));
  /* validate */
  if(surf->vertices == NULL) {
    free(surf);
    return NULL;
  }

  surf->vertexLabels = (double*) malloc(surf->maxVertices*sizeof(double));
  /* validate */
  if(surf->vertexLabels == NULL) {
    free(surf);
    return NULL;
  }

  surf->vertexBoundaryness = (int*) malloc(surf->maxVertices*sizeof(int));

  /* allocate each vertex */
  for(i = 0; i < surf->maxVertices; i++) {
    surf->vertices[i] = (double*) malloc(3*sizeof(double));
    surf->vertexLabels[i] = -1;
  }
  surf->vertexNormals = NULL;
  surf->numVertexNormals = 0;

  /* allocate faces */
  surf->maxFaces = maxFaces;
  surf->numFaces = 0;
  surf->faces = (unsigned int**) malloc(surf->maxFaces*sizeof(unsigned int*));
  /* validate */
  if(surf->faces == NULL) {
    surf->maxFaces = 0;
    deleteSurface(surf);
    return NULL;
  }

  /* allocate each face */
  for(i = 0; i < surf->maxFaces; i++) {
    surf->faces[i] = (unsigned int*) malloc(3*sizeof(unsigned int));
  }

  surf->faceNormals = NULL;
  surf->numFaceNormals = 0;

  /* curvature stuff */
  surf->curvatureMode = NO_CURVATURE;
  surf->curvature = NULL;

  /* surface topology fields */
  surf->NTaV = NULL;
  surf->LTaV = NULL;
  surf->v1TaV = NULL;
  surf->v2TaV = NULL;
  surf->NNoV = NULL;
  surf->maxNTaV = 0;
  surf->E = 0;
  surf->Neighbors = NULL;
  surf->NTaE = NULL;
  surf->LT1aNE = NULL;
  surf->LT2aNE = NULL;
  surf->manifoldness = SURF_UNKNOWN;
  surf->vertexManifoldness = NULL;
  surf->OT = NULL;
  surf->OV = NULL;
  surf->CC = NULL;


  return surf;
}

/**
 * delete a surface and free all vertices and faces
 * returns success (1) or flibvpure (-1)
 */
int deleteSurface(surface *surf) {
  int i;

  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* unprepare (delete all topology fields) */
  unprepareSurf(surf);

  /* delete the vertices */
  for(i = 0; i < surf->maxVertices; i++) {
    free(surf->vertices[i]);
  }
  free(surf->vertices);

  /* delete the vertex normals */
  if(surf->numVertexNormals > 0) {
    for(i = 0; i < surf->numVertexNormals; i++) {
      free(surf->vertexNormals[i]);
    }
    free(surf->vertexNormals);
  }

  /* delete the vertex labels */
  if(surf->vertexLabels != NULL) {
    free(surf->vertexLabels);
  }

  /* delete the vertex boundaries */
  if(surf->vertexBoundaryness != NULL) {
    free(surf->vertexBoundaryness);
  }

  /* delete the faces */
  for(i = 0; i < surf->maxFaces; i++) {
    free(surf->faces[i]);
  }
  free(surf->faces);

  /* delete the face normals */
  if(surf->numFaceNormals > 0) {
    for(i = 0; i < surf->numFaceNormals; i++) {
      free(surf->faceNormals[i]);
    }
    free(surf->faceNormals);
  }

  /* free curvature, if necessary */
  if(surf->curvatureMode != NO_CURVATURE) {
    free(surf->curvature);
  }

  /* free the structure */
  free(surf);

  return SR_SUCCESS;
}

/**
 * expands the capacity of the vertex array of the surface to amount
 * returns success (1) or flibvpure (-1)
 */
int expandVertices(surface *surf, int amount) {
  int i;
  int maxOldVertices;
  double **oldVertices;
  double *oldLabels;
  int *oldBoundaryness;

  fprintf(stderr,"expanding a surface from %d to %d vertices.\n",
	  surf->maxVertices, amount);

  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(amount <= surf->maxVertices) return SR_FAILURE;

  /** re-preprocess **/
  // INEFFICIENT! REIMPLEMENT!
  unprepareSurf(surf);

  /* store the old vertices */
  maxOldVertices = surf->maxVertices;
  oldVertices = surf->vertices;
  oldLabels = surf->vertexLabels;
  oldBoundaryness = surf->vertexBoundaryness;

  /* allocate and copy to new vertex list */
  surf->maxVertices = amount;
  surf->vertices = (double**) malloc(surf->maxVertices*sizeof(double*));
  surf->vertexLabels = (double*) malloc(surf->maxVertices*sizeof(double));
  surf->vertexBoundaryness = (int*) malloc(surf->maxVertices*sizeof(int));
  for(i = 0; i < surf->maxVertices; i++) {
    /* if this vertex already exists, just copy it over */
    if(i < maxOldVertices) {
      surf->vertices[i] = oldVertices[i];
      surf->vertexLabels[i] = oldLabels[i];
      surf->vertexBoundaryness[i] = oldBoundaryness[i];
    }
    else { /* create it */
      surf->vertices[i] = (double*) malloc(3*sizeof(double));
      surf->vertexLabels[i] = -1;
      surf->vertexBoundaryness[i] = -1;
    }
  }

  /** re-preprocess **/
  // INEFFICIENT! REIMPLEMENT!
  preprocessSurf(surf);

  /* free the old vertices */
  free(oldVertices);
  free(oldLabels);

  return SR_SUCCESS;
}

/**
 * expands the capacity of the faces array of the surface to amount
 * returns success (1) or flibvpure (-1)
 */
int expandFaces(surface *surf, int amount) {
  int i;
  unsigned int maxOldFaces;
  unsigned int **oldFaces;

  fprintf(stderr,"expanding a surface from %d to %d faces.\n",
	  surf->maxFaces, amount);

  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(amount <= surf->maxFaces) return SR_FAILURE;

  /** re-preprocess **/
  // INEFFICIENT! REIMPLEMENT!
  unprepareSurf(surf);

  /* store the old faces */
  maxOldFaces = surf->maxFaces;
  oldFaces = surf->faces;

  /* allocate and copy to new face list */
  surf->maxFaces = amount;
  surf->faces = (unsigned int**) malloc(surf->maxFaces*sizeof(unsigned int*));
  for(i = 0; i < surf->maxFaces; i++) {
    /* if this face already exists, just copy it over */
    if(i < maxOldFaces) {
      surf->faces[i] = oldFaces[i];
    }
    else { /* create it */
      surf->faces[i] = (unsigned int*) malloc(3*sizeof(unsigned int));
    }
  }

  /** re-preprocess **/
  // INEFFICIENT! REIMPLEMENT!
  preprocessSurf(surf);

  /* free the old vertices */
  free(oldFaces);

  return SR_SUCCESS;
}

/** surface conversion and manipulation tools **/

/**
 * add a vertex to a surface
 * this function just uses the x, y, and z fields of the vertex struct
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexStruct(surface *surf, vertex* v) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(v == NULL) return SR_FAILURE;

  /* add coords */
  return addVertexCoord(surf,v->x,v->y,v->z,v->label,v->boundary);
}

/**
 * add a vertex to a surface, without duplicating vertices this
 * function just the x, y, and z fields of the vertex struct to check
 * if the vertex is already in the surface. if so, the original index
 * is returned.
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexStructNoDuplicate(surface *surf, vertex* v) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(v == NULL) return SR_FAILURE;

  /* add coords */
  return addVertexCoordNoDuplicate(surf,v->x,v->y,v->z,v->label,v->boundary);
}

/**
 * add a vertex to a surface by coordinates
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexCoord(surface *surf, double x, double y, double z, double lab, int boundary) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* expand if necessary */
  if(surf->numVertices >= surf->maxVertices) {
    if(!expandVertices(surf,2*surf->maxVertices)) {
      return SR_FAILURE;
    }
  }

  /* insert the vertex */
  surf->vertices[surf->numVertices][SR_X] = x;
  surf->vertices[surf->numVertices][SR_Y] = y;
  surf->vertices[surf->numVertices][SR_Z] = z;

  surf->vertexLabels[surf->numVertices] = lab;
  surf->vertexBoundaryness[surf->numVertices] = boundary;

  surf->numVertices++;

  return surf->numVertices-1;
}

/**
 * add a vertex to a surface by coordinates if it doesnt already exist
 * in the surface
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexCoordNoDuplicate(surface *surf, double x, double y, double z, double lab, int boundary) {
  int index = -1;

  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* expand if necessary */
  if(surf->numVertices >= surf->maxVertices) {
    if(!expandVertices(surf,2*surf->maxVertices)) {
      return SR_FAILURE;
    }
  }

  index = findVertexCoords(surf,x,y,z);
  if(index < 0) { /* didn't find it */
    /* insert the vertex */
    surf->vertices[surf->numVertices][SR_X] = x;
    surf->vertices[surf->numVertices][SR_Y] = y;
    surf->vertices[surf->numVertices][SR_Z] = z;

    surf->vertexLabels[surf->numVertices] = lab;
    surf->vertexBoundaryness[surf->numVertices] = boundary;

    index = surf->numVertices;
    surf->numVertices++;
  }
  else if(boundary == BRANCHED_BOUNDARY
	  &&
	  surf->vertexBoundaryness[index] != BRANCHED_BOUNDARY) {
    surf->vertexBoundaryness[index] = BRANCHED_BOUNDARY;
  }
  else if(boundary == BOUNDARY 
	  && surf->vertexBoundaryness[index] != BOUNDARY) {
    surf->vertexBoundaryness[index] = BOUNDARY;
  }

  return index;
}

/**
 * add a vertex to a surface by an array of coords
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexArr(surface *surf, double *coords) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(coords == NULL) return SR_FAILURE;

  /* add coords */
  return addVertexCoord(surf,coords[0],coords[1],coords[2],-1,-1);
}

/**
 * add a vertex to a surface by an array of coords if coords not already
 * in the surface
 * returns index of added vertex or flibvpure (-1)
 */
int addVertexArrNoDuplicate(surface *surf, double *coords) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(coords == NULL) return SR_FAILURE;

  /* add coords */
  return addVertexCoordNoDuplicate(surf,coords[0],coords[1],coords[2],-1,-1);
}

/**
 * add a face to a surface
 * this function just uses the v1, v2, and v3 fields of the face struct
 * returns index of added face or flibvpure (-1)
 */
int addFaceStruct(surface *surf, face *f) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(f == NULL) return SR_FAILURE;

  /* add face indices */
  return addFaceInd(surf,f->v1,f->v2,f->v3);
}

/**
 * add a face to a surface by indices
 * returns index of added face or flibvpure (-1)
 */
int addFaceInd(surface *surf, int v1, int v2, int v3) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* expand if necessary */
  if(surf->numFaces >= surf->maxFaces) {
    if(!expandFaces(surf,2*surf->maxFaces)) {
      return SR_FAILURE;
    }
  }

  /* insert the vertex */
  surf->faces[surf->numFaces][0] = v1;
  surf->faces[surf->numFaces][1] = v2;
  surf->faces[surf->numFaces][2] = v3;
  surf->numFaces++;

  return surf->numFaces-1;
}

/**
 * add a face to a surface by an array of indices
 * returns index of added face or flibvpure (-1)
 */
int addFaceArr(surface *surf, int *inds) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(inds == NULL) return SR_FAILURE;

  /* add face indices */
  return addFaceInd(surf,inds[0],inds[1],inds[2]);
}

/**
 * combine two existing surfaces
 * returns resulting surface, or NULL, if flibvpure
 */
surface *combineSurfaces(surface *surf1, surface *surf2) {
  surface *surf;
  int curV = 0,
    curV1 = 0,
    curV2 = 0,
    curF = 0,
    curF1 = 0,
    curF2 = 0;

  /* validate */
  if(surf1 == NULL || surf2 == NULL) return NULL;

  /* create the surface */
  surf = createSurface(surf1->numVertices
		       +surf2->numVertices
		       +DEFAULT_MAX_VERTICES,
		       surf1->numFaces
		       +surf2->numFaces
		       +DEFAULT_MAX_FACES);
  if(surf == NULL) return NULL;

  /* combine the vertices */
  while(surf->numVertices < surf1->numVertices + surf2->numVertices) {
    /* test if we should add more surf1 vertices */
    if(curV1 < surf1->numVertices) {
      addVertexArr(surf,surf1->vertices[curV1++]);
      curV++;
    }
    else if(curV2 < surf2->numVertices) {
      addVertexArr(surf,surf2->vertices[curV2++]);
      curV++;
    }
  }

  /* combine the faces */
  while(surf->numFaces < surf1->numFaces + surf2->numFaces) {
    /* test if we should add more surf1 faces */
    if(curF1 < surf1->numFaces) {
      addFaceInd(surf,
		 surf1->faces[curF1][0],
		 surf1->faces[curF1][1],
		 surf1->faces[curF1][2]);
      curF1++;
      curF++;
    }
    else if(curF2 < surf2->numFaces) {
      addFaceInd(surf,
		 surf2->faces[curF2][0]+surf1->numVertices,
		 surf2->faces[curF2][1]+surf1->numVertices,
		 surf2->faces[curF2][2]+surf1->numVertices);
      curF2++;
      curF++;
    }
  }

  return surf;
}

/**
 * build a surface structure from a face list
 */
surface *buildSurfaceFromFaceList(list *faceList) {
/*   int v1 = 0, v2 = 0, v3 = 0; */
/*   listNode *i; */
/*   face *f; */

  /* create a surface hopefully big enough to hold the verts and faces */
  surface * surf = createSurface(2*listSize(faceList),2*listSize(faceList));
  if(surf == NULL) return NULL;

  addFaceListToSurface(surf,faceList);

/*   for(i = getListNode(faceList,0); i; i = (listNode*) i->next) { */
/*     f = (face*) i->data; */

/*     /\* find or add each vertex *\/ */
/*     if(SR_FAILURE == (v1 = findVertexStruct(surf,f->left->v1))) { */
/*       v1 = addVertexStruct(surf,f->left->v1); */
/*     } */
/*     if(SR_FAILURE == (v2 = findVertexStruct(surf,f->left->v2))) { */
/*       v2 = addVertexStruct(surf,f->left->v2); */
/*     } */
/*     if(verticesEqual(f->left->v1,f->right->v1)) { */
/*       if(SR_FAILURE == (v3 = findVertexStruct(surf,f->right->v2))) { */
/* 	v3 = addVertexStruct(surf,f->right->v2); */
/*       } */
/*     } */
/*     else { */
/*       if(SR_FAILURE == (v3 = findVertexStruct(surf,f->right->v1))) { */
/* 	v3 = addVertexStruct(surf,f->right->v1); */
/*       } */
/*     } */

/*     /\* add the face *\/ */
/*     addFaceInd(surf,v1,v2,v3); */
/*   } */

  return surf;
}

/**
 * add a facelist to a surface structure without duplicating vertices
 * already contained in a surface.
 */
void addFaceListToSurface(surface *surf, list *faceList) {
  int v1 = 0, v2 = 0, v3 = 0;
  listNode *i;
  face *f;



  /* go through all faces in the list, adding each */
  for(i = getListNode(faceList,0); i; i = (listNode*) i->next) {
    f = (face*) i->data;

    /* find or add each vertex */
    if(SR_FAILURE == (v1 = findVertexStruct(surf,f->left->v1))) {
      v1 = addVertexStruct(surf,f->left->v1);
    }
    if(SR_FAILURE == (v2 = findVertexStruct(surf,f->left->v2))) {
      v2 = addVertexStruct(surf,f->left->v2);
    }
    if(verticesEqual(f->left->v1,f->right->v1)) {
      if(SR_FAILURE == (v3 = findVertexStruct(surf,f->right->v2))) {
	v3 = addVertexStruct(surf,f->right->v2);
      }
    }
    else {
      if(SR_FAILURE == (v3 = findVertexStruct(surf,f->right->v1))) {
	v3 = addVertexStruct(surf,f->right->v1);
      }
    }

    /* add the face */
    addFaceInd(surf,v1,v2,v3);
  }
}


/**
 * caps a terminating surface
 * inputs are the contour to be capped and the z coordinate to add any vertices
 * returns a facelist, or NULL, if the contour cannot be capped
 */
list *capContour(contour *cont, enum CAPPING_METHOD cappingMethod, double z) {
  list *retList = NULL;

  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 3
     || cont->closed == OPEN) {
    return NULL;
  }

  /* build a cap based on the method */
  switch(cappingMethod) {
  case TRIANGULATE:
    retList = triangulateContour(cont);
    break;
  case TILE_SKELETON:
    retList = capContourSkeleton(cont,z);
    break;
  case TILE_POINT:
    retList = capContourPoint(cont,z);
    break;
  default:
    markAVertexAsBoundary(cont);
    break;
  }

  return retList;
}

/**
 * caps a contour by converting its skeleton to a contour, then
 * triangulating them together
 * takes as input a contour and a z coordinate to place the skeleton
 * returns a face list, or NULL, if the method fails
 */
list *capContourSkeleton(contour *cont, double z) {
  list *faceList, *skelVerts = NULL, *skelEdges = NULL;
  listNode *ln, *ln2;
  contour *skelCont;
  vertex *firstV, *lastV, *firstLastV, *v, *otherV, *minV;
  vertex ***vlist;
  edge *e;
  double minAngle, curAngle, otherAngle, diffAngle;
  int i,numV,curV,tr,allLabel = 1;

  /* validate */
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) == 0
     || cont->closed == OPEN) {
    return NULL;
  }

  /* determine if all vertices on this contour are labeled, if so label all added vertices */
  for(ln = getListNode(cont->vertices,0); ln; ln = (listNode*) ln->next) {
    if(((vertex*)ln->data)->label == -1) {
      allLabel = 0;
      break;
    }
  }

  /* get the skeleton, tile to a point if necessary */
  if(contourIsCircular(cont)
     || SR_FAILURE == getContourSkeleton(cont,&skelVerts,&skelEdges)
     || listSize(skelVerts) < 2 || listSize(skelEdges) < 1) {
    return capContourPoint(cont,z);
  }

  vlist = (vertex***) malloc(sizeof(vertex**)*MAX_SKEL_V);

  /* assign the z coords as passed and labels as found */
  for(ln = getListNode(skelVerts,0), i = 0; ln; ln = (listNode*) ln->next, i++) {
    vlist[i] = (vertex**) malloc(sizeof(vertex*)*MAX_SKEL_N);
    ((vertex*)ln->data)->z = z;
    if(allLabel) {
      ((vertex*)ln->data)->label = 1;
    }
  }

  /* if there is only one vertex in skeleton, tile it */
  if(listSize(skelVerts) == 1) {
    faceList = tileContourToSingleVertex(cont,(vertex*)getListNode(skelVerts,0)->data);
    freeList(skelVerts);
    freeListAndData(skelEdges);
    return faceList;
  }

  /* build a vertex adjacency list */
  numV = listSize(skelVerts);
  for(ln = getListNode(skelVerts,0), i = 0; ln; ln = (listNode*) ln->next, i++) {
    v = (vertex*) ln->data;
    v->number = i;

    /* build an adjacency list for this vertex by searching the edges for it */
    curV = 0;
    for(ln2 = getListNode(skelEdges,0); ln2; ln2 = (listNode*) ln2->next) {
      e = (edge*) ln2->data;
      if(e->v1 == v && e->v2 != NULL) {
	vlist[i][curV] = e->v2;
	curV++;
      }
      else if(e->v2 == v && e->v1 != NULL) {
	vlist[i][curV] = e->v1;
	curV++;
      }

      if(curV > MAX_SKEL_N) {
	fprintf(stderr,"FOUND A SKEL VERTEX WITH TOO MANY NEIGHBORS!\n");
	curV--;
	break;
      }
    }
    vlist[i][curV] = NULL;
  }


  /* delete vertices that are too close together */
  if(!resample) {
    lastResampleDist = fabs(z-((vertex*)getListNode(cont->vertices,0)->data)->z);
  }

  /* resample skeleton to avoid nonmanifoldness */
  for(ln = getListNode(skelVerts,0); ln; ln = (listNode*) ln->next) {
    v = (vertex*) ln->data;
    for(i = 0; vlist[v->number][i]; i++) {
      otherV = vlist[v->number][i];
      if(dist(*v,*otherV) < lastResampleDist) {
	// found one to delete, shift the list to cover it
	removeAdjacentEntry(vlist[v->number],otherV);

	// copy the adjacent list for the deleted node here
	appendAdjacentEntries(vlist,v,otherV);

	// erase all adajcency entries for removed vertex
	vlist[otherV->number][0] = NULL;

	i--;
      }
    }
  }

  /* traverse the skeleton, creating a contour */
  skelCont = createContour();
  skelCont->closed = CLOSED;

  //firstV = v = ((edge*)getListNode(skelEdges,0)->data)->v2;
  //firstLastV = lastV = ((edge*)getListNode(skelEdges,0)->data)->v1;

  // find the first and last vertices
  v = firstV = firstLastV = lastV = NULL;
  for(i = 0; v == NULL && i < listSize(skelVerts); i++) {
    if(vlist[i][0] != NULL) {
      firstV = v = vlist[i][0];
      firstLastV = lastV = (vertex*) getListNode(skelVerts,i)->data;
    }
  }

  /* tile to a point if all vertices were deleted */
  if(v == NULL) {
    faceList = capContourPoint(cont,z);
  }
  else {

    // follow the adjacency list until we find the edge we just added
    do {
      /* add the current vertex as a vertex in the contour */
      if(-1 != findInListI(skelCont->vertices,v)) {
	enqueue(skelCont->vertices,v);
      }
      else {
	enqueue(skelCont->vertices,copyVertex(v));
      }

      /* search for the leftmost adjacent vertex */
      minV = v;
      minAngle = 0;//2*SR_PI;
      curAngle = unwrap(atan2(lastV->y-v->y,lastV->x-v->x));
      //edgeList = adjList[findInListI(skelVerts,v)];
      i = findInListI(skelVerts,v);

      /* if there is only one neighbor, make it v */
      //if(listSize(edgeList) == 1) {
      if(vlist[i][1] == NULL) {
	lastV = v;
	//v = (vertex*) getListNode(edgeList,0)->data;
	v = vlist[i][0];
	continue;
      }

      /* otherwise search neighbors */
      //for(ln = getListNode(edgeList,0); ln; ln = (listNode*) ln->next) {
      for(curV = 0; curV < MAX_SKEL_N && vlist[i][curV] != NULL; curV++) {
	//otherV = (vertex*) ln->data;
	otherV = vlist[i][curV];
	otherAngle = unwrap(atan2(otherV->y-v->y,otherV->x-v->x));

	diffAngle = unwrap(curAngle-otherAngle);
	if(lastV == otherV) {
	  diffAngle = 0;//2*SR_PI;
	}

	if(diffAngle > minAngle) {
	  minAngle = diffAngle;
	  minV = otherV;
	}
      }

      lastV = v;
      v = minV;

    } while(!(v == firstV && lastV == firstLastV));

    /* get the skeleton, tile to a point if necessary */
    if(listSize(skelCont->vertices) < 2) {
      faceList = capContourPoint(cont,z);
    }
    else {
      /* tile the contour with the skeleton contour */
      tr = testRepeat;
      //testRepeat = 1;
      faceList = getFaceList(cont,skelCont);
      testRepeat = tr;
    }
  }

  /* free the skeleton */
  for(i = 0; i < numV; i++) {
    free(vlist[i]);
  }
  free(vlist);
  freeListAndData(skelCont->vertices);
  freeList(skelCont->adjacentContours);
  freeList(skelCont->adjacentBackwardContours);
  free(skelCont);

  freeListAndData(skelEdges);
  freeListAndData(skelVerts);

  return faceList;
}

/**
 * caps a contour by tiling it with a point at its centroid and the z
 * coordinate passed
 * takes as input a contour and a z coordinate to place the point
 * returns a face list, or NULL, if the method fails
 */
list *capContourPoint(contour *cont, double z) {
  vertex v;
  int allLabel = 1;
  listNode *ln;

  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) == 0
     || cont->closed == OPEN) {
    return NULL;
  }

  /* validate */
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) == 0
     || cont->closed == OPEN) {
    return NULL;
  }

  /* determine if all vertices on this contour are labeled, if so label all added vertices */
  for(ln = getListNode(cont->vertices,0); ln; ln = (listNode*) ln->next) {
    if(((vertex*)ln->data)->label == -1) {
      allLabel = 0;
      break;
    }
  }

  /* get the centroid and assign the z coord */
  getContourCentroid(cont,&v);
  v.z = z;

  if(allLabel) {
    v.label = 1;
  }

  return tileContourToSingleVertex(cont,&v);
}

/**
 * determine if a contour is roughly circular
 */
int contourIsCircular(contour *cont) {
  listNode* ln;
  vertex centroid;
  double *dist;
  int cur = 0;
  double circularity;

  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 3) {
    return 0;
  }

  /* compute distance from each vertex to the centroid */
  getContourCentroid(cont,&centroid);
  dist = (double*) malloc(listSize(cont->vertices)*sizeof(double));
  for(ln = getListNode(cont->vertices,0); ln; ln = (listNode*) ln->next) {
    dist[cur++] = dist(centroid,*((vertex*) ln->data));
  }

  circularity = gsl_stats_sd(dist,1,listSize(cont->vertices))
    / gsl_stats_mean(dist,1,listSize(cont->vertices));

  free(dist);

  return circularity < thresholdForCircularity;
}

/**
 * mark a vertex as on a boundary so we don't close a boundary we didnt cap
 */
void markAVertexAsBoundary(contour *cont) {
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 1) {
    return;
  }

  ((vertex*) getListNode(cont->vertices,0)->data)->boundary = BOUNDARY;
}

/**
 * mark a vertex as on a boundary so we don't close a boundary we didnt cap
 * for all contour on a slice
 */
void markAVertexAsBoundarySlice(list *slice) {
  if(slice == NULL || listSize(slice) < 1) {
    return;
  }
  
  listNode* ln;
  for(ln = getListNode(slice,0); ln; ln = (listNode*) ln->next) {
    markAVertexAsBoundary((contour*) ln->data);
  }
}

/** tools for manipulating adject vertex arrays for the skeleton computation**/

/**
 * removes a vertex from an adjacency list
 */
void removeAdjacentEntry(vertex **vlist, vertex *otherV) {
  int i,j;

  for(i = 0; vlist[i]; i++) {
    if(vlist[i] == otherV) {
      for(j = i+1; vlist[j]; j++) {
	vlist[j-1] = vlist[j];
      }
      vlist[j-1] = NULL;
      i--;
    }
  }
}

/**
 * replace a vertex in an adjacency list with another
 */
void replaceAdjacentEntry(vertex **vlist, vertex *oldV, vertex *newV) {
  int i;

  for(i = 0; vlist[i]; i++) {
    if(vlist[i] == oldV) {
      vlist[i] = newV;
    }
  }
}

/**
 * appends a deleted nodes adjacency entries to the current list
 */
void appendAdjacentEntries(vertex ***vlist, vertex *v, vertex *otherV) {
  int i, cur;

  for(cur = 0; vlist[v->number][cur]; cur++);

  for(i = 0; vlist[otherV->number][i]; i++) {
    if(vlist[otherV->number][i] != v) {
      vlist[v->number][cur++] = vlist[otherV->number][i];
      replaceAdjacentEntry(vlist[vlist[otherV->number][i]->number],otherV,v);
    }
  }
  vlist[v->number][cur] = NULL;
}

// hole filling methods

/**
 * fill holes
 */
void fillHoles(surface *surf, int numHolesLeft) {
  if(surf == NULL) {
    return;
  }

  list *boundaries;
  listNode *ln;
  int i;

  preprocessSurf(surf);
  if(surf->manifoldness != SURF_MANIFOLD) {
    if(SR_VERBOSE) {
      fprintf(stdout,"surface is not manifold!\ntopoFixing...");
    }

    topoFixer(surf);
    preprocessSurf(surf);

    if(SR_VERBOSE) {
      if(surf->manifoldness != SURF_MANIFOLD) {
	fprintf(stdout,"failed\n");
      }
      else {
	fprintf(stdout,"succeded\n");
      }
    }

  }
  else if(SR_VERBOSE) {
    fprintf(stdout,"surface is manifold\n");
  }

  boundaries = getBoundaries(surf);

  if(SR_VERBOSE) {
    fprintf(stdout,"surface has %d boundaries, fixing...\n",
	    listSize(boundaries));
  }

  // sort the boundaries in decending order of number of vertices
  int *numVB = (int*) malloc(listSize(boundaries)*sizeof(int));
  size_t *sortedVB = (size_t*) malloc(listSize(boundaries)*sizeof(int));
  for(ln = getListNode(boundaries,0), i = 0; ln; ln = (listNode*) ln->next, i++) {
    numVB[i] = listSize((list*) ln->data);
  }
  gsl_sort_int_index(sortedVB,numVB,1,listSize(boundaries));
  
  // fill the desireded boundaries
  for(i = 0; i < listSize(boundaries)-numHolesLeft; i++) {
    if(SR_VERBOSE) {
      fprintf(stdout,"filling hole %d of %d\n",i+1,listSize(boundaries));
    }

    add3DPolyToSurf(surf,(list*)getListNode(boundaries,sortedVB[i])->data);
  }
  freeListAndData(boundaries);
  
  if(SR_VERBOSE) {
    boundaries = getBoundaries(surf);
    fprintf(stdout,"surface now has %d boundaries\n",
	    listSize(boundaries));
    freeListAndData(boundaries);
  }

  free(sortedVB);
}

/**
 * fill holes created by the branching process.
 */
void fillBranchedHoles(surface *surf) {
  if(surf == NULL) {
    return;
  }

  list *boundaries, *boundary;
  listNode *ln, *ln2;
  long v;

  preprocessSurf(surf);
  if(surf->manifoldness != SURF_MANIFOLD) {
    if(SR_VERBOSE) {
      fprintf(stdout,"surface is not manifold!\ntopoFixing...");
    }

    topoFixer(surf);
    preprocessSurf(surf);

    if(SR_VERBOSE) {
      if(surf->manifoldness != SURF_MANIFOLD) {
	fprintf(stdout,"failed\n");
	return;
      }
      else {
	fprintf(stdout,"succeded\n");
      }
    }

  }
  else if(SR_VERBOSE){
    fprintf(stdout,"surface is manifold\n");
  }

  boundaries = getBoundaries(surf);

  if(SR_VERBOSE) {
    fprintf(stdout,"surface has %d boundaries, fixing...\n",
	    listSize(boundaries));
  }

  // fill the desired boundaries
  for(ln = getListNode(boundaries,0); ln; ln = (listNode*) ln->next) {
    boundary = (list*) ln->data;

    // debugging
//    for(ln2 = getListNode(boundary,0), v = 0; ln2; ln2 = (listNode*) ln2->next, v++) {
//      printf("vert %d boundaryness == %d\n",v,surf->vertexBoundaryness[(int) ln2->data]);
//    }
//    printf("\n");

    // find a boundary that shouldnt be a boundary
    for(ln2 = getListNode(boundary,0); ln2; ln2 = (listNode*) ln2->next) {
      v = (long) ln2->data;
      if(surf->vertexBoundaryness[v] == BOUNDARY) {
	break;
      }
    }
    
    // test if we found none
    if(ln2 != NULL) {
      freeList(boundary);
      continue;
    }

    if(SR_VERBOSE) {
      fprintf(stdout,"filling hole with %d vertices\n",listSize(boundary));
    }

    add3DPolyToSurf(surf,boundary);
    freeList(boundary);
  }
  freeList(boundaries);
  
  preprocessSurf(surf);
  if(surf->manifoldness == SURF_MANIFOLD) {
    if(SR_VERBOSE) {
      boundaries = getBoundaries(surf);
      fprintf(stdout,"surface now has %d boundaries\n",
	      listSize(boundaries));
      freeListAndData(boundaries);
    }
  }
}

/**
 * fill holes created by the branching process.
 */
void fillBranchedHolesOld(list *contourPairs, surface *surf) {
  if(contourPairs == NULL || listSize(contourPairs) < 2 || surf == NULL) {
    return;
  }

  list *plist;
  listNode *i,*j;
  int c1n, c2n, oc1n, oc2n, foundSeat;
  vertex *c10, *c1e, *c20, *c2e, *oc10, *oc1e, *oc20, *oc2e;
  contourPair *cp, *ocp;
  contour *c1, *c2, *oc1, *oc2;

  // search the possible pair endpoints for surface branched boundary edges
  for(i = getListNode(contourPairs,0); i; i = (listNode*) i->next) {
    cp = (contourPair*) i->data;

    // test if both contours have branched open endpoints
    if(isBranchedBoundary(cp)) {

      // two cases:

      // one: the seat of the pants
      foundSeat = 0;
      for(j = (listNode*) i->next; j; j = (listNode*) j->next) {
	ocp = (contourPair*) j->data;

	if(isBranchedBoundary(ocp)) {
	  // try to match the branched contours

	  // decide which contours have a common origin, so would have
	  // adjacent contour verts
	  if(cp->c1->origin == ocp->c1->origin) {
	    c1 = cp->c1;
	    c2 = cp->c2;
	    oc1 = ocp->c1;
	    oc2 = ocp->c2;
	  }
	  else {
	    c1 = cp->c2;
	    c2 = cp->c1;
	    oc1 = ocp->c2;
	    oc2 = ocp->c1;
	  }

	  // get the first contour ordering and endpoints
	  c1n = listSize(c1->vertices);
          c2n = listSize(c2->vertices);

          c10 = (vertex*) getListNode(c1->vertices,0)->data;
          c1e = (vertex*) getListNode(c1->vertices,c1n-1)->data;
          c20 = getAdjacentVertex(surf, c10,
    			(vertex*) getListNode(c2->vertices,0)->data,
    			(vertex*) getListNode(c2->vertices,c2n-1)->data);
          c2e = (c20 == (vertex*) getListNode(c2->vertices,0)->data)
	    ? (vertex*) getListNode(c2->vertices,c2n-1)->data
	    : (vertex*) getListNode(c2->vertices,0)->data;


	  // get the other contours ordering and endpoints
	  oc1n = listSize(oc1->vertices);
	  oc2n = listSize(oc2->vertices);

	  oc10 = (vertex*) getListNode(oc1->vertices,0)->data;
	  oc1e = (vertex*) getListNode(oc1->vertices,oc1n-1)->data;
	  oc20 = getAdjacentVertex(surf, oc10,
		(vertex*) getListNode(oc2->vertices,0)->data,
		(vertex*) getListNode(oc2->vertices,oc2n-1)->data);
	  oc2e = (oc20 == (vertex*) getListNode(oc2->vertices,0)->data)
	    ? (vertex*) getListNode(oc2->vertices,oc2n-1)->data
	    : (vertex*) getListNode(oc2->vertices,0)->data;

	  // decide if there is a valid branched boundary. by testing if the
	  // contour endpoints were adjacent on the original contours if so,
	  // fill the boundary
	  if(verticesOriginallyAdjacent(c1,oc1,c10,oc10)
	     && verticesOriginallyAdjacent(c1,oc1,c1e,oc1e)) {

	    plist = newList(LIST);
	    enqueue(plist,c10);
	    enqueue(plist,c20);
	    enqueue(plist,c2e);
	    enqueue(plist,c1e);
	    enqueue(plist,oc1e);
	    enqueue(plist,oc2e);
	    enqueue(plist,oc20);
	    enqueue(plist,oc20);

	    add3DPolyToSurf(surf,plist);
	    freeList(plist);
	  }
	  else if(verticesOriginallyAdjacent(c1,oc1,c10,oc1e)
	     && verticesOriginallyAdjacent(c1,oc1,c1e,oc10)) {

	    plist = newList(LIST);
	    enqueue(plist,c10);
	    enqueue(plist,c20);
	    enqueue(plist,c2e);
	    enqueue(plist,c1e);
	    enqueue(plist,oc10);
	    enqueue(plist,oc20);
	    enqueue(plist,oc2e);
	    enqueue(plist,oc2e);

	    add3DPolyToSurf(surf,plist);
	    freeList(plist);
	  }
	}
      }

      // two: the thumb
      // TODO
    }
  }
}

/**
 * decides whether a contour pair is one side of a branched boundary hole
 */
int isBranchedBoundary(contourPair *cp) {
  if(cp == NULL || cp->c1 == NULL || cp->c2 == NULL
     || cp->c1->vertices == NULL || cp->c2->vertices == NULL) {
    return 0;
  }

  int c1n = listSize(cp->c1->vertices);
  int c2n = listSize(cp->c2->vertices);

  // test if both contours have branched open endpoints
  if(((vertex*)getListNode(cp->c1->vertices,0)->data)->boundary
     == BRANCHED_BOUNDARY
     &&
     ((vertex*)getListNode(cp->c1->vertices,c1n-1)->data)->boundary
     == BRANCHED_BOUNDARY
     &&
     ((vertex*)getListNode(cp->c2->vertices,0)->data)->boundary
     == BRANCHED_BOUNDARY
     &&
     ((vertex*)getListNode(cp->c2->vertices,c2n-1)->data)->boundary
     == BRANCHED_BOUNDARY
     ) {
    return 1;
  }
  return 0;
}

/**
 * tests if two vertices from separate contours pairs were adjacent in the
 * original contour
 */
int verticesOriginallyAdjacent(contour *c1, contour *oc1, vertex *v1, vertex *ov1) {
  if(c1 == NULL || oc1 == NULL || c1->origin != oc1->origin || v1 == NULL || ov1 == NULL) {
    return 0;
  }

  if(abs(v1->number-ov1->number) == 1
     || (v1->number == 0 && ov1->number
         == listSize(((contour*)oc1->origin)->vertices)-1)
     || (ov1->number == 0 && v1->number
         == listSize(((contour*)c1->origin)->vertices)-1)) {
    return 1;
  }

  return 0;
}

/**
 * given a surface and a vertex, determine which of two other vertices are
 * 1-neighbors of the passed vertex, if any.
 * RETURNS THE FIRST INSTANCE FOUND!
 */
vertex *getAdjacentVertex(surface *surf, vertex *v0, vertex *v10, vertex *v11) {
  if(surf == NULL || v0 == NULL || v10 == NULL || v11 == NULL) {
    return NULL;
  }

  vertex v;
  int i,
    v0i = -1,
    v10i = -1,
    v11i = -1;

  // search the vertices for v0, v10 and v11
  for(i = 0; i < surf->numVertices; i++) {
    v.x = surf->vertices[i][0];
    v.y = surf->vertices[i][1];
    v.z = surf->vertices[i][2];

    if(verticesEqual(v0,&v)) {
      v0i = i;
    }
    if(verticesEqual(v10,&v)) {
      v10i = i;
    }
    if(verticesEqual(v11,&v)) {
      v11i = i;
    }

    if(v0i != -1 && v10i != -1 && v11i != -1)
      break;
  }

  // search the faces for an instance of one edge
  for(i = 0; i < surf->numFaces; i++) {
    if(surf->faces[i][0] == v0i
       || surf->faces[i][1] == v0i
       || surf->faces[i][2] == v0i) {
      if(surf->faces[i][0] == v10i
	 || surf->faces[i][1] == v10i
	 || surf->faces[i][2] == v10i) {
	return v10;
      }
      if(surf->faces[i][0] == v11i
	 || surf->faces[i][1] == v11i
	 || surf->faces[i][2] == v11i) {
	return v11;
      }
    }
  }

  return NULL;
}

/**
 * adds a three dimensional polygon to an existing surface.
 * input is an ordered set of points representing the boundary of the polygon.
 * algorithm is a la liepa 2003 and barequet and sharir 1995
 */
void add3DPolyToSurf(surface *surf, list *plist) {
  if(surf == NULL || plist == NULL || listSize(plist) < 3) {
    return;
  }

  int n = listSize(plist);
  int *vlist = (int*) malloc(n*sizeof(int));
  vector **W = (vector**) malloc(n*sizeof(vector*));
  int **lambda = (int**) malloc(n*sizeof(int*));
  int i,j,k,m,allBoundaries = 1;
  long tmp;
  listNode *ln;
  vector *curCost;

  for(ln = getListNode(plist,0), i = 0; ln; ln = (listNode*) ln->next, i++) {
    tmp = (long) ln->data;
    vlist[i] = (int) tmp;
    W[i] = (vector*) malloc(n*sizeof(vector));
    lambda[i] = (int*) malloc(n*sizeof(int));
    if(surf->vertexLabels[vlist[i]] == -1) {
      allBoundaries = 0;
    }
  }

  // initialize W
  for(i = 0; i < n-1; i++) {
    if(i < n-2) {
      curCost = getPolyTilingTriangleCost(surf,vlist[i],vlist[i+1],vlist[i+2]);
      W[i][i+2] = *curCost;
      free(curCost);
    }
    W[i][i+1].x = 0;
    W[i][i+1].y = 0;
  }

  // iterate the solution
  for(j = 3; j < n; j++) {
    for(i = 0; i < n-j; i++) {
      k = i+j;

      // find the minimum cost triangle at this step
      lambda[i][k] = i+1;
      curCost = getPolyTilingTriangleCost(surf,vlist[i],vlist[i+1],vlist[k]);
      W[i][k].x = W[i][i+1].x + W[i+1][k].x + curCost->x;
      W[i][k].y = W[i][i+1].y + W[i+1][k].y + curCost->y;
      //fprintf(stdout,"%d %d %d %g\n", i, k, i+1, W[i][k].x);
      free(curCost);

      for(m = i+2; m < k; m++) {
	curCost = getPolyTilingTriangleCost(surf,vlist[i],vlist[m],vlist[k]);
	curCost->x += W[i][m].x + W[m][k].x;
	curCost->y += W[i][m].y + W[m][k].y;
	//fprintf(stdout,"%d %d %d %g\n", i, k, m, curCost->x);

	if(curCost->x < W[i][k].x 
	   || (fabs(curCost->x - W[i][k].x) < SR_TOL 
	       && curCost->y < W[i][k].y)) {
	  W[i][k] = *curCost;
	  lambda[i][k] = m;
	  
	}
	free(curCost);
      }

    }
  }

  tracePolyBoundary(surf,lambda,vlist,0,n-1);

  // free stuff
  for(i = 0; i < n; i++) {
    free(W[i]);
    free(lambda[i]);
  }
  free(W);
  free(lambda);
  free(vlist);
}

/**
 * recursively trace a poly boundary to recover an optimal tiling
 */ 
void tracePolyBoundary(surface *surf, int **lambda, int *vlist, int i, int k) {
  int o;
  
  if(i+2 == k) {
    addFaceInd(surf,vlist[i],vlist[i+1],vlist[k]);
    preprocessSurf(surf);
  }
  else {
    o = lambda[i][k];
    if(o != i+1) tracePolyBoundary(surf,lambda, vlist, i, o);
    addFaceInd(surf,vlist[i],vlist[o],vlist[k]);
    preprocessSurf(surf);
    if(o != k-1) tracePolyBoundary(surf,lambda, vlist, o, k);
  }
}

/**
 * determine the cost of adding this face to a polygonal tiling
 */
vector *getPolyTilingTriangleCost(surface *surf, int v1, int v2, int v3) {
  vector *cost = (vector*) malloc(sizeof(vector));
  if(surf == NULL || v1 < 0 || v1 > surf->numVertices
     || v2 < 0 || v2 > surf->numVertices
     || v3 < 0 || v3 > surf->numVertices) {
    cost->x = cost->y = SR_BIG;
    return cost;
  }

  unsigned int n, *f;
  vector norm1, norm2;
  vertex a,b,c;

  a.x = surf->vertices[v1][0];
  a.y = surf->vertices[v1][1];
  a.z = surf->vertices[v1][2];
  b.x = surf->vertices[v2][0];
  b.y = surf->vertices[v2][1];
  b.z = surf->vertices[v2][2];
  c.x = surf->vertices[v3][0];
  c.y = surf->vertices[v3][1];
  c.z = surf->vertices[v3][2];

  cost->y = getFaceArea(&a,&b,&c);

  // find max dihedral angle (max additive inverseof the dot product) in
  // neighboring triangles
  getFaceNormalV(surf->vertices[v1],
		 surf->vertices[v2],
		 surf->vertices[v3],
		 &norm1);
  n = findIndex(surf->Neighbors[v1],surf->NNoV[v1],v2);
  if(n > -1) {
    f = surf->faces[surf->LT1aNE[v1][n] == -1 
		    ? surf->LT2aNE[v1][n] : surf->LT1aNE[v1][n]];
    getFaceNormalV(surf->vertices[f[0]],
		   surf->vertices[f[1]],
		   surf->vertices[f[2]],
		   &norm2);
    if(-dot(norm1,norm2) > cost->x) {      
      cost->x = -dot(norm1,norm2);
    }
  }

  n = findIndex(surf->Neighbors[v1],surf->NNoV[v1],v3);
  if(n > -1) {
    f = surf->faces[surf->LT1aNE[v1][n] == -1 
		    ? surf->LT2aNE[v1][n] : surf->LT1aNE[v1][n]];
    getFaceNormalV(surf->vertices[f[0]],
		   surf->vertices[f[1]],
		   surf->vertices[f[2]],
		   &norm2);
    if(-dot(norm1,norm2) > cost->x) {      
      cost->x = -dot(norm1,norm2);
    }
  }

  n = findIndex(surf->Neighbors[v2],surf->NNoV[v2],v3);
  if(n > -1) {
    f = surf->faces[surf->LT1aNE[v2][n] == -1 
		    ? surf->LT2aNE[v2][n] : surf->LT1aNE[v2][n]];
    getFaceNormalV(surf->vertices[f[0]],
		   surf->vertices[f[1]],
		   surf->vertices[f[2]],
		   &norm2);
    if(-dot(norm1,norm2) > cost->x) {      
      cost->x = -dot(norm1,norm2);
    }
  }

  return cost;
}

/**
 * adds a three dimensional polygon to an existing surface.
 * input is an ordered set of points representing the boundary of the polygon.
 * WARNING! THIS IS A DIRTY HACK, IT WILL PRODUCE SELF_INTERSECTIONS!
 */
void add3DPolyToSurfOld(surface *surf, list *plist) {
  if(surf == NULL || plist == NULL || listSize(plist) < 3) {
    return;
  }

  list *flist = newList(LIST);
  face *f;
  listNode *first, *cur, *next;

  // loop through, connect the first vertex to each, then each to next.
  first = getListNode(plist,0);
  for(cur = getListNode(plist,1), next = getListNode(plist,2); next;
      cur = next, next = (listNode*) next->next) {
    f = (face*) malloc(sizeof(face));
    f->left = (edge*) malloc(sizeof(edge));
    f->right = (edge*) malloc(sizeof(edge));

    f->left->v1 = (vertex*) first->data;
    f->right->v1 = (vertex*) first->data;

    f->left->v2 = (vertex*) cur->data;
    f->right->v2 = (vertex*) next->data;
  }

  addFaceListToSurface(surf,flist);

  freeFaceList(flist);
}

/** surface query tools **/

/**
 * find a vertex by vertex struct
 * returns index of vertex or flibvpure (-1)
 */
int findVertexStruct(surface *surf, vertex *v) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(v == NULL) return SR_FAILURE;

  return findVertexCoords(surf,v->x,v->y,v->z);
}

/**
 * find a vertex by coords
 * returns index of vertex or flibvpure (-1)
 */
int findVertexCoords(surface *surf, double x, double y, double z) {
  int i;

  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* iterate over vertices, test if within tolreance */
  for(i = surf->numVertices-1; i >= 0; i--) {
    if(fabs(surf->vertices[i][SR_X] - x) < SR_TOL &&
       fabs(surf->vertices[i][SR_Y] - y) < SR_TOL &&
       fabs(surf->vertices[i][SR_Z] - z) < SR_TOL) {
      return i;
    }
  }

  return SR_FAILURE;
}

/**
 * find a vertex by array
 * returns index of vertex or flibvpure (-1)
 */
int findVertexArr(surface *surf, double *coords) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(coords == NULL) return SR_FAILURE;

  return findVertexCoords(surf,coords[SR_X],coords[SR_Y],coords[SR_Z]);
}

/**
 * find a face by face struct
 * returns index of face or flibvpure (-1)
 */
int findFaceStruct(surface *surf, face *f) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(f == NULL) return SR_FAILURE;

  return findFaceInds(surf,f->v1,f->v2,f->v3);
}

/**
 * find a face by indices
 * returns index of face or flibvpure (-1)
 */
int findFaceInds(surface *surf, int v1, int v2, int v3) {
  int i;

  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* iterate over faces looking for vertex indices */
  for(i = 0; i < surf->numFaces; i++) {
    if(surf->faces[i][0] == v1 &&
       surf->faces[i][1] == v2 &&
       surf->faces[i][2] == v3) {
      return i;
    }
  }

  return SR_FAILURE;
}

/**
 * find a face by array
 * returns index of face or flibvpure (-1)
 */
int findFaceArr(surface *surf, int *inds) {
  /* validate */
  if(surf == NULL) return SR_FAILURE;
  if(inds == NULL) return SR_FAILURE;

  return findFaceInds(surf,inds[0],inds[1],inds[2]);
}

/**
 * tests two vertices for equality
 */
int verticesEqual(vertex *v1, vertex *v2) {
  if(v1 == NULL || v2 == NULL) return 0;

  return fabs(v1->x-v2->x) < SR_TOL
    && fabs(v1->y-v2->y) < SR_TOL
    && fabs(v1->z-v2->z) < SR_TOL;
}

/**
 * copies an edge
 */
edge *copyEdge(edge *e) {
  edge *newe = (edge*) malloc(sizeof(edge));
  newe->v1 = createVertex();
  newe->v2 = createVertex();

  *(newe->v1) = *(e->v1);
  *(newe->v2) = *(e->v2);

  newe->length = e->length;

  return newe;
}

/**
 * copy a contour
 */
void copyContour(contour *orig, contour *copy) {
  /* copy each field */
  copy->closed = orig->closed;

  copy->vertices = cloneList(orig->vertices);
  copy->adjacentContours = cloneList(orig->adjacentContours);
}

/**
 * use Triangle to get a delaunay triangulated contour, returns a face list
 */
list *triangulateContour(contour *c) {
  char polyFile[SR_MAX_STR_LEN];
  char offFile[SR_MAX_STR_LEN];
  surface *surf;
  list *faceList;
  face *f;
  vertex *v;
  int i;
  double z;

  

  /* create the filenames */
  static int calls = 0;
  sprintf(polyFile,"/tmp/surfRecon_triangle%d.poly",calls);
  sprintf(offFile,"/tmp/surfRecon_triangle%d.1.off",calls);
  calls++;

  /* validate the contour vertex list */
  if(c == NULL || c->vertices == NULL || listSize(c->vertices) == 0) {
    return NULL;
  }

  /* determine the z coordinate of this contour */
  z = ((vertex*) getListNode(c->vertices,0)->data)->z;

  /* write an input file for triangle */
  writeTrianglePolyFile(c, polyFile);

  /* execute Triangle */
  if(executeTriangle(polyFile) != 0) {
    fprintf(stderr, "error: Triangle returned an error, failed in triangulation. of %s\n", polyFile);
    return NULL;
  }

  surf = readOFF(offFile);

  if(surf == NULL) {
    return NULL;
  }

  /* build the face list from the surface */
  faceList = newList(LIST);
  for(i = 0; i < surf->numFaces; i++) {
    f = (face*) malloc(sizeof(face));
    f->left = (edge*) malloc(sizeof(edge));
    f->right = (edge*) malloc(sizeof(edge));

    /* allocate a vertex */
    v = createVertex();

    /* get the current vertex */
    v->x = surf->vertices[surf->faces[i][0]][0];
    v->y = surf->vertices[surf->faces[i][0]][1];
    v->z = z;

    v->boundary = FALSE;

    /* add the vertex to the face */
    f->left->v1 = v;

    /* allocate a vertex */
    v = createVertex();

    /* get the current vertex */
    v->x = surf->vertices[surf->faces[i][1]][0];
    v->y = surf->vertices[surf->faces[i][1]][1];
    v->z = z;

    v->boundary = FALSE;

    /* add the vertex to the face */
    f->left->v2 = v;
    f->right->v2 = v;

    /* allocate a vertex */
    v = createVertex();

    /* get the current vertex */
    v->x = surf->vertices[surf->faces[i][2]][0];
    v->y = surf->vertices[surf->faces[i][2]][1];
    v->z = z;

    v->boundary = FALSE;

    /* add the vertex to the face */
    f->right->v1 = v;

    /* add the face to the list */
    enqueue(faceList,f);
  }

  deleteSurface(surf);
  return faceList;
}

/**
 * tests two edges for intersection
 */
int edgesIntersect(edge *e1, edge *e2) {
  return segmentsIntersect(e1->v1->x, e1->v1->y, e1->v1->z,
			   e1->v2->x, e1->v2->y, e1->v2->z,
			   e2->v1->x, e2->v1->y, e2->v1->z,
			   e2->v2->x, e2->v2->y, e2->v2->z);
}

/**
 * tests two segments for intersection based on the 12 coordinates of
 * the end points of two three dimensional segments
 */
int segmentsIntersect(float x00, float y00, float z00,
		      float x01, float y01, float z01,
		      float x10, float y10, float z10,
		      float x11, float y11, float z11) {
  vector u,v,w,w2;
  float dist,tmp,s,t,t0,t1;

  u.x = x01 - x00;
  u.y = y01 - y00;
  u.z = z01 - z00;

  v.x = x11 - x10;
  v.y = y11 - y10;
  v.z = z11 - z10;

  w.x = x00 - x10;
  w.y = y00 - y10;
  w.z = z00 - z10;

  dist = perp(u,v);

  /* parallel? */
  if(fabs(dist) < SR_TOL) {
    /* colinear ? */
    if(perp(u,w) != 0 || perp(v,w) != 0) {
      return FALSE;
    }

    /* collinear */
    w2.x = x01 - x10;
    w2.y = y01 - y10;
    w2.z = z01 - z10;
    if(v.x != 0) {
      t0 = w.x / v.x;
      t1 = w2.x / v.x;
    }
    else {
      t0 = w.y / v.y;
      t1 = w2.y / v.y;
    }
    if(t0 > t1) {
      tmp = t0;
      t0 = t1;
      t1 = tmp;
    }
    if(t0 > 1 || t1 < 0) {
      return FALSE;
    }
    t0 = t0<0 ? 0 : t0;
    t1 = t1>1 ? 1 : t1;
    if(t0 == t1) {
      return TRUE;
    }

    return TRUE;
  }

  /* non parallel */

  /* test e1 */
  s = perp(v,w) / dist;
  if(s < 0 || s > 1) {
    return FALSE;
  }

  /* test e2 */
  t = perp(u,w) / dist;
  if(t < 0 || t > 1) {
    return FALSE;
  }

  return TRUE;
}

/**
 * tests if an edge intersects a face
 */
int edgeFaceIntersect(face *f, edge *e) {

  /* build the datastructures */
  tPointi T;
  tPointd q, r, p;
  char code;

  Vertices[0][X] = f->left->v1->x;
  Vertices[0][Y] = f->left->v1->y;
  Vertices[0][Z] = f->left->v1->z;

  Vertices[1][X] = f->left->v2->x;
  Vertices[1][Y] = f->left->v2->y;
  Vertices[1][Z] = f->left->v2->z;

  if(verticesEqual(f->left->v1,f->right->v1)) {
    Vertices[2][X] = f->right->v2->x;
    Vertices[2][Y] = f->right->v2->y;
    Vertices[2][Z] = f->right->v2->z;
  }
  else {
    Vertices[2][X] = f->right->v1->x;
    Vertices[2][Y] = f->right->v1->y;
    Vertices[2][Z] = f->right->v1->z;
  }

  T[0] = 0;
  T[1] = 1;
  T[2] = 2;

  q[0] = e->v1->x;
  q[1] = e->v1->y;
  q[2] = e->v1->z;

  r[0] = e->v2->x;
  r[1] = e->v2->y;
  r[2] = e->v2->z;

  code = SegTriInt(T,q,r,p);

  if(code == 'f') {
    return 1;
  }

  return 0;
}

/* /\** */
/*  * tests two segments for intersection based on the 8 coordinates of */
/*  * the end points of two two dimensional segments */
/*  *\/ */
/* int segmentsIntersect2d(float x00, float y00,  */
/* 			float x01, float y01,  */
/* 			float x10, float y10,  */
/* 			float x11, float y11) { */

/*   Vector    u = S1.P1 - S1.P0; */
/*   Vector    v = S2.P1 - S2.P0; */
/*   Vector    w = S1.P0 - S2.P0; */
/*   float     D = perp(u,v); */

/*   // test if they are parallel (includes either being a point) */
/*   if (fabs(D) < SMALL_NUM) {          // S1 and S2 are parallel */
/*     if (perp(u,w) != 0 || perp(v,w) != 0) { */
/*       return 0;                   // they are NOT collinear */
/*     } */
/*     // they are collinear or degenerate */
/*     // check if they are degenerate points */
/*     float du = dot(u,u); */
/*     float dv = dot(v,v); */
/*     if (du==0 && dv==0) {           // both segments are points */
/*       if (S1.P0 != S2.P0)         // they are distinct points */
/* 	return 0; */
/*       *I0 = S1.P0;                // they are the same point */
/*       return 1; */
/*     } */
/*     if (du==0) {                    // S1 is a single point */
/*       if (inSegment(S1.P0, S2) == 0)  // but is not in S2 */
/* 	return 0; */
/*       *I0 = S1.P0; */
/*       return 1; */
/*     } */
/*     if (dv==0) {                    // S2 a single point */
/*       if (inSegment(S2.P0, S1) == 0)  // but is not in S1 */
/* 	return 0; */
/*       *I0 = S2.P0; */
/*       return 1; */
/*     } */
/*     // they are collinear segments - get overlap (or not) */
/*     float t0, t1;                   // endpoints of S1 in eqn for S2 */
/*     Vector w2 = S1.P1 - S2.P0; */
/*     if (v.x != 0) { */
/*       t0 = w.x / v.x; */
/*       t1 = w2.x / v.x; */
/*     } */
/*     else { */
/*       t0 = w.y / v.y; */
/*       t1 = w2.y / v.y; */
/*     } */
/*     if (t0 > t1) {                  // must have t0 smaller than t1 */
/*       float t=t0; t0=t1; t1=t;    // swap if not */
/*     } */
/*     if (t0 > 1 || t1 < 0) { */
/*       return 0;     // NO overlap */
/*     } */
/*     t0 = t0<0? 0 : t0;              // clip to min 0 */
/*     t1 = t1>1? 1 : t1;              // clip to max 1 */
/*     if (t0 == t1) {                 // intersect is a point */
/*       *I0 = S2.P0 + t0 * v; */
/*       return 1; */
/*     } */

/*     // they overlap in a valid subsegment */
/*     *I0 = S2.P0 + t0 * v; */
/*     *I1 = S2.P0 + t1 * v; */
/*     return 2; */
/*   } */

/*   // the segments are skew and may intersect in a point */
/*   // get the intersect parameter for S1 */
/*   float     sI = perp(v,w) / D; */
/*   if (sI < 0 || sI > 1)               // no intersect with S1 */
/*     return 0; */

/*   // get the intersect parameter for S2 */
/*   float     tI = perp(u,w) / D; */
/*   if (tI < 0 || tI > 1)               // no intersect with S2 */
/*     return 0; */

/*   *I0 = S1.P0 + sI * u;               // compute S1 intersect point */
/*   return 1; */
/* } */

/**
 * test intersection of an edge and a face
 *
 * algorithm pieced from two places
 * http://astronomy.swin.edu.au/~pbourke/geometry/planeline/
 * http://www.ecse.rpi.edu/Homepages/wrf/research/geom/pnpoly.html
 */
//int edgeFaceIntersect(edge *e, face *f) {
//  vector faceNorm,dist1,dist2;
//  vertex onePt,*fp0,*fp1,*fp2,intersectionPt;
//  float num, den, r;
//  int in;
//
//  /* test for zero length edge */
//  if(fabs(e->v1->x-e->v2->x) < SR_TOL
//     && fabs(e->v1->y-e->v2->y) < SR_TOL
//     && fabs(e->v1->z-e->v2->z) < SR_TOL) {
//    return FALSE;
//  }
//
////  printf("%lf %lf %lf :: %lf %lf %lf\n",
////	 e->v1->x,e->v1->y,e->v1->z,e->v2->x,e->v2->y,e->v2->z);
////
//
////  if(fabs(e->v2->x-(-25.970018)) < 0.00001
////     && fabs(e->v2->y-49.276383) < 0.00001
////     && fabs(e->v1->x-(-25.78949)) < 0.00001
////     && fabs(e->v1->y-49.225258) < 0.00001) {
////    printf("found one\n");
////  }
//
//  /* get face normal*/
//  getFaceNormal(f,&faceNorm);
//
//  /* store one point for reference */
//  onePt.x = f->left->v1->x;
//  onePt.y = f->left->v1->y;
//  onePt.z = f->left->v1->z;
//
//  /* make some distance vectors */
//  dist1.x = onePt.x - e->v2->x;
//  dist1.y = onePt.y - e->v2->y;
//  dist1.z = onePt.z - e->v2->z;
//
//  dist2.x = e->v1->x - e->v2->x;
//  dist2.y = e->v1->y - e->v2->y;
//  dist2.z = e->v1->z - e->v2->z;
//
//  /* test segment for intersection with plane */
//  num = dot(faceNorm,dist1);
//  den = dot(faceNorm,dist2);
//
//  /* test for parallel segment */
//  if(fabs(den) < SR_TOL) {
//    return FALSE;
//  }
//
//  r = num/den;
//  /* do intersection with plane test */
//  if(r < 0 || r > 1) {
//    return FALSE;
//  }
//
//  /* find intersection point */
//  intersectionPt.x = e->v2->x + (dist2.x)*r;
//  intersectionPt.y = e->v2->y + (dist2.y)*r;
//  intersectionPt.z = e->v2->z + (dist2.z)*r;
//
//  /* assign the points of the face */
//  fp0 = f->left->v1;
//  fp1 = f->left->v2;
//
//  if(verticesEqual(fp1,f->right->v1)) {
//    fp2 = f->right->v2;
//  }
//  else {
//    fp2 = f->right->v1;
//  }
//
//  /* see if the intersection is a vertex */
//  if(verticesEqual(&intersectionPt,fp0)
//     || verticesEqual(&intersectionPt,fp1)
//     || verticesEqual(&intersectionPt,fp2)) {
//    return FALSE;
//  }
//
//  /* test intersection point for inclusion in polygon */
//  in = FALSE;
//
//  /* test points 0 and 1 */
//  if((((fp0->z<=intersectionPt.z) && (intersectionPt.z<fp1->z)) ||
//      ((fp1->z<=intersectionPt.z) && (intersectionPt.z<fp0->z))) &&
//     (intersectionPt.x < (fp1->x - fp0->x) *
//      (intersectionPt.z - fp0->z)
//      / (fp1->z - fp0->z) + fp0->x)) {
//    in = !in;
//  }
//  /* test points 1 and 2 */
//  if((((fp1->z<=intersectionPt.z) && (intersectionPt.z<fp2->z)) ||
//      ((fp2->z<=intersectionPt.z) && (intersectionPt.z<fp1->z))) &&
//     (intersectionPt.x < (fp2->x - fp1->x) *
//      (intersectionPt.z - fp1->z)
//      / (fp2->z - fp1->z) + fp1->x)) {
//    in = !in;
//  }
//
//  return in;
//}

/**
 * center the surf at its center of mass
 */
void centerSurf(surface *surf) {
  vector c;
  int i;

  /* validate */
  if(surf == NULL) return;

  c = computeCenterOfMass(surf);
  /* subtract off the center */
  for(i = 0; i < surf->numVertices; i++) {
    surf->vertices[i][0] -= c.x;
    surf->vertices[i][1] -= c.y;
    surf->vertices[i][2] -= c.z;
  }
}

/**
 * calculate the center of mass
 */
vector computeCenterOfMass(surface *surf) {
  int i;
  vector c;
  c.x = c.y = c.z = 0;

  /* validate */
  if(surf == NULL) return c;

  for(i = 0; i < surf->numVertices; i++) {
    c.x += surf->vertices[i][0];
    c.y += surf->vertices[i][1];
    c.z += surf->vertices[i][2];
  }

  c.x/=surf->numVertices;
  c.y/=surf->numVertices;
  c.z/=surf->numVertices;

  return c;
}

/**
 * calculate the bounds in each dimension
 */
void computeSurfBounds(surface *surf, double *bounds) {
  unsigned int i;

  if(surf == NULL || bounds == NULL) {
    return;
  }

  /* initialize bounds */
  bounds[SR_XMIN] = SR_BIG;
  bounds[SR_XMAX] = SR_LITTLE;
  bounds[SR_YMIN] = SR_BIG;
  bounds[SR_YMAX] = SR_LITTLE;
  bounds[SR_ZMIN] = SR_BIG;
  bounds[SR_ZMAX] = SR_LITTLE;

  /* search through vert coords */
  for(i = 0; i < surf->numVertices; i++) {
    if(surf->vertices[i][SR_X] < bounds[SR_XMIN])
      bounds[SR_XMIN] = surf->vertices[i][SR_X];
    if(surf->vertices[i][SR_X] > bounds[SR_XMAX])
      bounds[SR_XMAX] = surf->vertices[i][SR_X];
    if(surf->vertices[i][SR_Y] < bounds[SR_YMIN])
      bounds[SR_YMIN] = surf->vertices[i][SR_Y];
    if(surf->vertices[i][SR_Y] > bounds[SR_YMAX])
      bounds[SR_YMAX] = surf->vertices[i][SR_Y];
    if(surf->vertices[i][SR_Z] < bounds[SR_ZMIN])
      bounds[SR_ZMIN] = surf->vertices[i][SR_Z];
    if(surf->vertices[i][SR_Z] > bounds[SR_ZMAX])
      bounds[SR_ZMAX] = surf->vertices[i][SR_Z];
  }
}

/**
 * transform the coordinates of the surface vertices by a linear transform
 */
int transformSurfaceVertices(surface *surf, float A[4][4]) {
  int i;
  float u[4],v[4];

  /* validate */
  if(surf == NULL) return SR_FAILURE;

  /* iterate over vertices, transforming each */
  for(i = 0; i < surf->numVertices; i++) {
    v[0] = (float) surf->vertices[i][0];
    v[1] = (float) surf->vertices[i][1];
    v[2] = (float) surf->vertices[i][2];
    v[3] = 1;

    matrixMult4by1(A,v,u);

    surf->vertices[i][0] = (double) u[0];
    surf->vertices[i][1] = (double) u[1];
    surf->vertices[i][2] = (double) u[2];
  }

  return SR_SUCCESS;
}

/**
 * compute the vertex normals for all vertices
 */
void computeVertexNormals(surface *surf) {
  int i;
  vector *v,sum;
  list **faceNormalsAtVertex;
  list *faceNormalList;
  listNode *ln2;

  /* validate */
  if(surf == NULL) return;

  /* check if the normals have already been computed */
  if(surf->vertexNormals != NULL) {
    for(i = 0; i < surf->numVertexNormals; i++) {
      free(surf->vertexNormals[i]);
    }
    free(surf->vertexNormals);
  }

  faceNormalsAtVertex = (list**) malloc(surf->numVertices*sizeof(list*));

  /* create a list for each vertex */
  for(i = 0; i < surf->numVertices; i++) {
    faceNormalsAtVertex[i] = newList(LIST);
  }

  /* iterate over faces, adding its normal to the list for each vertex */
  for(i = 0; i < surf->numFaces; i++) {
    v = (vector*) malloc(sizeof(vector));

    /* get this face normal */
    getFaceNormalV(surf->vertices[surf->faces[i][0]],
		   surf->vertices[surf->faces[i][1]],
		   surf->vertices[surf->faces[i][2]], v);

    /* add the face normal to the list of faces for the 3 vertices */
    enqueue(faceNormalsAtVertex[surf->faces[i][0]],v);
    enqueue(faceNormalsAtVertex[surf->faces[i][1]],v);
    enqueue(faceNormalsAtVertex[surf->faces[i][2]],v);
  }

  surf->vertexNormals = (double**) malloc(surf->numVertices*sizeof(double*));
  /* compute the normal for each vertex */
  for(i = 0; i < surf->numVertices; i++) {
    faceNormalList = faceNormalsAtVertex[i];

    /* iterate over adjacent faces, adding the normals */
    sum.x = 0.0;
    sum.y = 0.0;
    sum.z = 0.0;
    for(ln2 = getListNode(faceNormalList,0); ln2; ln2 = (listNode*)ln2->next) {
      v = (vector*) ln2->data;
      sum.x += v->x;
      sum.y += v->y;
      sum.z += v->z;
    }

    /* normalize the normal */
    normalizeVector(&sum);

    /* allocate and assign */
    surf->vertexNormals[i] = (double*) malloc(3*sizeof(double));
    surf->vertexNormals[i][0] = sum.x;
    surf->vertexNormals[i][1] = sum.y;
    surf->vertexNormals[i][2] = sum.z;
  }

  /* clean up the temporary lists */
  //  for(ln = getListNode(faceNormalsAtVertex,0); ln; ln = (listNode*) ln->next) {
  //    for(ln2 = getListNode((list*)ln->data,0); ln2; ln2=(listNode*)ln2->next) {
  //      free(ln2->data);
  //    }
  //    freeList((list*) ln->data);
  //  }
  //  freeList(faceNormalsAtVertex);
  //
  surf->numVertexNormals = surf->numVertices;
}

/**
 * compute the face normals for all faces
 */
void computeFaceNormals(surface *surf) {
  int i;
  vector v;

  /* validate */
  if(surf == NULL) return;

  /* check if the normals have already been computed */
  if(surf->faceNormals != NULL) {
    for(i = 0; i < surf->numFaceNormals; i++) {
      free(surf->faceNormals[i]);
    }
    free(surf->faceNormals);
  }

  surf->faceNormals = (double**) malloc(surf->numFaces*sizeof(double*));
  /* compute the normal for each face */
  for(i = 0; i < surf->numFaces; i++) {
    getFaceNormalV(surf->vertices[surf->faces[i][0]],
		   surf->vertices[surf->faces[i][1]],
		   surf->vertices[surf->faces[i][2]], &v);

    /* allocate and assign */
    surf->faceNormals[i] = (double*) malloc(3*sizeof(double));
    surf->faceNormals[i][0] = v.x;
    surf->faceNormals[i][1] = v.y;
    surf->faceNormals[i][2] = v.z;
  }

  surf->numFaceNormals = surf->numFaces;
}

/**
 * get the unit normal to a face
 */
void getFaceNormal(face *f, vector *n) {
  vector d1,d2;

  /* assign the coordinates for the difference vectors */
  d1.x = f->left->v1->x - f->left->v2->x;
  d1.y = f->left->v1->y - f->left->v2->y;
  d1.z = f->left->v1->z - f->left->v2->z;

  d2.x = f->right->v1->x - f->right->v2->x;
  d2.y = f->right->v1->y - f->right->v2->y;
  d2.z = f->right->v1->z - f->right->v2->z;

  cross(&d1,&d2,n);
  normalizeVector(n);
  n->y *= -1;
}

/**
 * get the unit normal to a face
 */
void getFaceNormalV(double *v1, double *v2, double *v3, vector *n) {
  vector d1,d2;

  /* assign the coordinates for the difference vectors */
  d1.x = v1[0]-v2[0];
  d1.y = v1[1]-v2[1];
  d1.z = v1[2]-v2[2];

  d2.x = v1[0]-v3[0];
  d2.y = v1[1]-v3[1];
  d2.z = v1[2]-v3[2];

  cross(&d1,&d2,n);
  normalizeVector(n);
  n->y *= -1;
}

/**
 * get the unit normal to a face
 */
void getFaceNormalVerts(vertex *v1, vertex *v2, vertex *v3, vector *n) {
  vector d1,d2;

  /* assign the coordinates for the difference vectors */
  d1.x = v1->x-v2->x;
  d1.y = v1->y-v2->y;
  d1.z = v1->z-v2->z;

  d2.x = v1->x-v3->x;
  d2.y = v1->y-v3->y;
  d2.z = v1->z-v3->z;

  cross(&d1,&d2,n);
  normalizeVector(n);
  n->y *= -1;
}

/**
 * cross procuct
 */
void cross(vector *v1, vector *v2, vector *c) {
  c->x = v1->y*v2->z - v1->z*v2->y;
  c->y = v1->x*v2->z - v1->z*v2->x;
  c->z = v1->x*v2->y - v1->y*v2->x;
}

/**
 * normalize a vector
 */
void normalizeVector(vector *v) {
  float len = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
  v->x/=len;
  v->y/=len;
  v->z/=len;
}

/**
 * dumps a contour pair to a stream
 */
void dumpContourPair(FILE *fp, contourPair *cp) {
  vertex *c10, *c1N, *c20, *c2N;

  /* validate */
  if(fp == NULL || cp == NULL || cp->c1 == NULL || cp->c2 == NULL) {
    return;
  }

  c10 = (vertex*) getListNode(cp->c1->vertices,0)->data;
  c1N = (vertex*) getListNode(cp->c1->vertices,listSize(cp->c1->vertices)-1)->data;
  c20 = (vertex*) getListNode(cp->c2->vertices,0)->data;
  c2N = (vertex*) getListNode(cp->c2->vertices,listSize(cp->c2->vertices)-1)->data;

  fprintf(fp, "%p & %p: \n%d: (%3.1f,%3.1f) to (%3.1f,%3.1f) & %d: (%3.1f,%3.1f) to (%3.1f,%3.1f)\n",
	  cp->c1, cp->c2,
	  listSize(cp->c1->vertices),
	  c10->x, c10->y,
	  c1N->x, c1N->y,
	  listSize(cp->c2->vertices),
	  c20->x, c20->y,
	  c2N->x, c2N->y
	  );
}

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/surfUtil.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
