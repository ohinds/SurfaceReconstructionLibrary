/*****************************************************************************
 * tile.c is the source file with functions to tile a surface between two
 * slices/contours for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#define TILE_VERSION_C "$Id: tile.c,v 1.56 2007/05/22 19:18:11 oph Exp $"

#include"tile.h"

/* min and max slice indices to reconstruct */
int minSliceInd = 0;
int maxSliceInd = INT_MAX;
int savePaths = 0;

int printNow = 0;

/* some constants */
const double INF = 8.9885e+307;

/* algorithm methods parameters */
enum TILING_METHOD tilingMethod = OPTIMIZING;
enum CAPPING_METHOD cappingMethod = TILE_SKELETON;
//enum CAPPING_METHOD cappingMethod = TILE_POINT;
//enum CAPPING_METHOD cappingMethod = -1;

/** cost function penalty terms **/
int testAngles = TRUE;
double anglePenalty = 1.0;

int testEdgeCrossing = TRUE;
double edgeCrossingPenalty = 0.1;

int testRepeat = FALSE;
double repeatedEdgePenalty = 0.2;

int testManhattan = TRUE;
double manhattanPenalty = 0.3;

/**
 * for debugging, output a matlab script to display the graph for a path
 */
void printGraph(graph *g, int source, int close);

/**
 * get a surface represented by contours
 */
surface *tileSlices(list *slices) {
  listNode *i, *j, *k;
  int indi, numVerts;
  surface *surf = NULL;
  list *slice, *faceList;
  contour *cont, *adjcont;
  double sliceDist,z;

  /* assign the boundary vertices of the contours */
  //assignBoundaries(ds);

  if(skipSlices < 1) skipSlices = 1;

  /* count the vertices so we can estimate the surface size */
  numVerts = 0;
  for(i = getListNode(slices,0); i; i = (listNode*) i->next) {
    slice = (list*) i->data;
    for(j = getListNode(slice,0); j; j = (listNode*) j->next) {
      cont = (contour*) j->data;
      if(cont == NULL || cont->vertices == NULL
	 || listSize(cont->vertices) == 0) continue;

      numVerts+=listSize(cont->vertices);
    }
  }

  /* allocate surface */
  surf = createSurface((int)(1.3*numVerts),(int)(2.6*numVerts));

  /* iterate over slices */
  if(SR_DEBUG) {
    fprintf(stdout,"performing surface tiling on %d slices...\n",
	    listSize(slices));
  }

  for(i = getListNode(slices,0), indi=0; i; i = (listNode*) i->next, indi++) {

    slice = (list*) i->data;

    /* check index bounds */
    if(indi < minSliceInd) {
      continue;
    }
    else if(indi == minSliceInd && indi == maxSliceInd) {
      markAVertexAsBoundarySlice(slice);
      markAVertexAsBoundarySlice((list*)((listNode*) i->next)->data);
    }
    else if(indi == minSliceInd && indi != 0) {
      markAVertexAsBoundarySlice(slice);
    }
    else if(indi == maxSliceInd && indi != listSize(slices)-1) {
      markAVertexAsBoundarySlice((list*)((listNode*) i->next)->data);
    }
    else if(indi > maxSliceInd) {
      fprintf(stdout,"done\n");
      return surf;
    }

    if(SR_VERBOSE) {
      fprintf(stdout, "tiling slice %d of %d...\n",indi,listSize(slices)-1);
    }

    /* check for empty slice contour list */
    if(listSize(slice) > 0 && listSize(((contour*) getListNode(slice,0)->data)->vertices) < 1) {
      continue;
    }

    /* cap closed contours on first slice */
    if(indi == 0) {
      /* get the interslice distance */
      sliceDist = getIntersliceDistance(slice,(list*)((listNode*)i->next)->data);

      for(j = getListNode(slice,0); j; j = (listNode*) j->next) {
	cont = (contour*) j->data;
	if(cont == NULL || cont->vertices == NULL || cont->closed == OPEN
	   || listSize(cont->vertices) < 3) continue;
	z = ((vertex*) getListNode(cont->vertices,0)->data)->z;
	faceList = NULL;
	faceList = capContour(cont,cappingMethod,z-sliceDist/2.0);

	if(faceList != NULL) {
	  addFaceListToSurface(surf,faceList);
	  freeFaceList(faceList);
	}
	else {
	  markAVertexAsBoundary(cont);
	}

      }
    }

    // prevent boundary filling if the last contour is closed and not capped
    if(indi == listSize(slices)-2 && cappingMethod == NO_CAPPING) {
      for(j = getListNode(slice,0); j; j = (listNode*) j->next) {
	cont = (contour*) j->data;
	if(cont == NULL || cont->vertices == NULL || cont->closed == OPEN
	   || listSize(cont->vertices) < 3) continue;
	for(k = getListNode(cont->adjacentContours,0); k; k = (listNode*) k->next) {
	  adjcont = (contour*) k->data;
	  if(adjcont == NULL || adjcont->vertices == NULL || adjcont->closed == OPEN
	     || listSize(adjcont->vertices) < 3) continue;

	  markAVertexAsBoundary(adjcont);
	}
      }
    }

    /* cap closed contours on last slice */
    if(indi == listSize(slices)-1) {
      /* get the interslice distance */
      sliceDist 
	= getIntersliceDistance(slice,(list*)getListNode(slices,indi-1)->data);

      for(j = getListNode(slice,0); j; j = (listNode*) j->next) {
	cont = (contour*) j->data;
	if(cont == NULL || cont->vertices == NULL || cont->closed == OPEN
	   || listSize(cont->vertices) < 3) continue;

	z = ((vertex*) getListNode(cont->vertices,0)->data)->z;
	faceList = NULL;
	faceList = capContour(cont,cappingMethod,z+sliceDist/2.0);

	if(faceList != NULL) {
	  addFaceListToSurface(surf,faceList);
	  freeFaceList(faceList);
	}
	else {
	  markAVertexAsBoundary(cont);
	}
      }
    }


    /* tile the current slice */
    switch(tilingMethod) {
    case OPTIMIZING:
      if(i->next == NULL) {
	continue;
      }

      tileSliceByOptimizingMethod(surf,slice,(list*)((listNode*)i->next)->data);
      break;
    default:
      break;
    }

  }

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  //freeList(slices);

  return surf;
}

/**
 * branch and tile a corresponded slice pair
 */
void tileSliceByOptimizingMethod(surface *surf, list *slice1, list *slice2) {
  listNode *i,*j;
  list *contourPairs, *faceList;
  contourPair *cp;//, *otherCp;
  contour *cont;//, *c1, *c2;
  int found;//, c10, c11, c20, c21;
  double sliceDist,z;
  
  contourPairs = branchSlices(slice1,slice2);

  // check that all contours were connected to something. if not, and if
  // closed, cap the contour

  // slice 1
  sliceDist = getIntersliceDistance(slice2,slice1);
  for(i = getListNode(slice1,0); i; i = (listNode*) i->next) {
    cont = (contour*) i->data;
    found = 0;
    for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
      cp = (contourPair*) j->data;
      if((contour*) cp->c1->origin == (contour*) cont->origin) {
	found = 1;
	break;
      }
    }

    if(!found && listSize(cont->vertices) > 0) { // cap
      z = ((vertex*) getListNode(cont->vertices,0)->data)->z + sliceDist/2.0;
      faceList = capContour((contour*) i->data,cappingMethod,z);
      if(faceList != NULL) {
	addFaceListToSurface(surf,faceList);
	freeFaceList(faceList);
      }
    }
  }

  // slice 2
  sliceDist*=-1;
  for(i = getListNode(slice2,0); i; i = (listNode*) i->next) {
    cont = (contour*) i->data;
    found = 0;
    for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
      cp = (contourPair*) j->data;
      if((contour*) cp->c2->origin == (contour*) cont->origin) {
	found = 1;
	break;
      }
    }

    if(!found && listSize(cont->vertices) > 0) { // cap
      z = ((vertex*) getListNode(cont->vertices,0)->data)->z + sliceDist/2.0;
      faceList = capContour((contour*) i->data,cappingMethod,z);
      if(faceList != NULL) {
	addFaceListToSurface(surf,faceList);
	freeFaceList(faceList);
      }
    }
  }

  /* connect all the contour pairs */
  for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
    cp = (contourPair*) j->data;

    /* fix the contour pairs that have whole parents so they aren't null */
    if(cp->c1 == NULL) cp->c1 = cp->c1Origin;
    if(cp->c2 == NULL) cp->c2 = cp->c2Origin;

    /* build a surface between the two contours, add it to existing surface*/
    addContoursToSurf(surf,cp->c1,cp->c2);
  }


  deleteContourPairList(contourPairs);
}

/**
 * tiles the contours on a slice by an optimizing method derived from the
 * tiling method of Fuchs, et al. (1977). The resulting surface is
 * nonintersecting and of minimum surface area.
 */
void tileSliceByOldOptimizingMethod(surface *surf, list *slice) {
  listNode *j,*k;
  int indj;
  list *contourPairs, *newContourPairs, *branchContourPairs, *handledPairs;
  contourPair *cp, *otherCp;
  contour *cont, *adjCont;

  printNow = 0;

  /* iterate over contours, adding all adjacent contours as pairs */
  contourPairs = newList(LIST);
  for(j = getListNode(slice,0),indj=0; j; j = (listNode*) j->next,indj++) {
    cont = (contour*) j->data;

    /* check for empty contour */
    if(cont == NULL || cont->vertices == NULL
       || listSize(cont->vertices) <= 0) {
      continue;
    }

    /* add pairs for all adjacent contours */
    for(k = getListNode(cont->adjacentContours,0); k; k = (listNode*) k->next) {
      adjCont = (contour*) k->data;

      /* check for non-empty adjacent contour */
      if(adjCont == NULL || adjCont->vertices == NULL
	 || listSize(adjCont->vertices) <= 0) {
	continue;
      }

      /* allocate a contour pair, and keep the contours NULL, for now */
      cp = (contourPair*) malloc(sizeof(contourPair));
      cp->c1 = cp->c2 = NULL;

      /* assign the original contours for each */
      cp->c1Origin = cont;
      cp->c2Origin = adjCont;

      /* add the pair to the list */
      enqueue(contourPairs,cp);
    }
  }

  /* first resolve the many to one branches */
  newContourPairs = newList(LIST);
  handledPairs = newList(LIST);
  for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
    cp = (contourPair*) j->data;

    /* check if we've already handled this pair*/
    if(listContains(handledPairs,cp)) continue;

    branchContourPairs = newList(LIST);

    /* search the rest of the pairs for other contours of the same c2orgin */
    for(k = getListNode(contourPairs,0); k; k = (listNode*) k->next) {
      otherCp = (contourPair*)k->data;
      if(otherCp->c2Origin == cp->c2Origin) { /* found a many to one */
	/* store it */
	enqueue(branchContourPairs,otherCp);
	enqueue(handledPairs,otherCp);
      }
    }

    /* break up the many to one branch, if there is one */
    if(listSize(branchContourPairs) > 1) {
      splitContoursManyToOne(branchContourPairs);
      appendList(newContourPairs,cloneList(branchContourPairs));
    }
    else {
      enqueue(newContourPairs,cp);
    }

    freeList(branchContourPairs);
  }

  /* store the modified list */
  freeList(contourPairs);
  freeList(handledPairs);
  contourPairs = newContourPairs;

  /* next solve the one to many branches */
  newContourPairs = newList(LIST);
  handledPairs = newList(LIST);
  for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
    cp = (contourPair*) j->data;

    /* check if we've already handled this pair*/
    if(listContains(handledPairs,cp)) continue;

    branchContourPairs = newList(LIST);

    /* search the rest of the pairs for other contours of the same c1orgin */
    for(k = getListNode(contourPairs,0); k; k = (listNode*) k->next) {
      otherCp = (contourPair*)k->data;
      if(otherCp->c1Origin == cp->c1Origin) { /* found a one to many */
	/* store it */
	enqueue(branchContourPairs,otherCp);
	enqueue(handledPairs,otherCp);
      }
    }

    /* break up the many to one branch, if there is one */
    if(listSize(branchContourPairs) > 1) {
      splitContoursOneToMany(branchContourPairs);
      appendList(newContourPairs,cloneList(branchContourPairs));
    }
    else {
      enqueue(newContourPairs,cp);
    }

    freeList(branchContourPairs);
  }

  /* store the modified list */
  freeList(contourPairs);
  freeList(handledPairs);
  contourPairs = newContourPairs;

  /* connect all the contour pairs */
  for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
    cp = (contourPair*) j->data;

    /* fix the contour pairs that have whole parents so they aren't null */
    if(cp->c1 == NULL) cp->c1 = cp->c1Origin;
    if(cp->c2 == NULL) cp->c2 = cp->c2Origin;

    /* build a surface between the two contours, add it to existing surface*/
    addContoursToSurf(surf,cp->c1,cp->c2);
  }

  freeList(contourPairs);

}

/**
 * tile a contour to a single vertex, return a facelist
 */
list* tileContourToSingleVertex(contour *cont, vertex* v) {
  list *faceList;
  listNode *ln;
  vertex *v1, *v2;
  face *f;
  edge *e;
  int i;


  /*validate*/
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 3 
     || v == NULL) {
    return NULL;
  }

  faceList = newList(LIST);
  for(ln = getListNode(cont->vertices,0), i = 0; 
      ln->next; ln = (listNode*) ln->next, i++) {
    v1 = (vertex*) ln->data;
    v2 = (vertex*) ((listNode*)ln->next)->data;

    f = (face*) malloc(sizeof(face));

    /* assign the left edge */
    e = (edge*) malloc(sizeof(edge));
    e->v1 = v1;
    e->v2 = v;
    f->left = copyEdge(e);
    f->left->v1->number = i;
    f->left->v2->number = 0;
    free(e);

    /* assign the right edge */
    e = (edge*) malloc(sizeof(edge));
    e->v1 = v2;
    e->v2 = v;
    f->right = copyEdge(e);
    f->right->v1->number = i+1;
    f->right->v2->number = 0;
    free(e);

    /* assign the vertices from the edges */
    f->v1 = f->left->v1->number;
    f->v2 = f->left->v2->number;
    f->v3 = f->right->v1->number;

    /* add the edge */
    push(faceList,f);
  }

  /* if closed, add wrapping face */
  if(cont->closed && listSize(cont->vertices) > 2) {
    v1 = (vertex*) ((listNode*)cont->vertices->tail)->data;
    v2 = (vertex*) ((listNode*)cont->vertices->head)->data;

    f = (face*) malloc(sizeof(face));

    /* assign the left edge */
    e = (edge*) malloc(sizeof(edge));
    e->v1 = v1;
    e->v2 = v;
    f->left = copyEdge(e);
    f->left->v1->number = i;
    f->left->v2->number = 0;
    free(e);

    /* assign the right edge */
    e = (edge*) malloc(sizeof(edge));
    e->v1 = v2;
    e->v2 = v;
    f->right = copyEdge(e);
    f->right->v1->number = i+1;
    f->right->v2->number = 0;
    free(e);

    /* assign the vertices from the edges */
    f->v1 = f->left->v1->number;
    f->v2 = f->left->v2->number;
    f->v3 = f->right->v1->number;

    /* add the edge */
    push(faceList,f);    
  }

  return faceList;
}

/**
 * build a surface from a set of branched contours
 */
surface *getBranchedSurface(contour *preBranch, list *postBranch) {
  enum BRANCHTYPE branchType = CONNECT;
  listNode *curPostCont, *lastPostCont;
  edge e;

  /* validate */
  if(preBranch == NULL || postBranch == NULL || listSize(postBranch) == 0) {
    return NULL;
  }

  /* aloocate the edge vertices */
  e.v1 = createVertex();
  e.v2 = createVertex();

  /* test the edges connecting centroids of the post branch contours
     for intersection with the prebranch contour to decide the branch
     type */
  lastPostCont = getListNode(postBranch,0);
  for(curPostCont = getListNode(postBranch,1); curPostCont;
      lastPostCont = curPostCont,
	curPostCont = (listNode*) curPostCont->next) {

    /* get the two centroids */
    getContourCentroid((contour*) lastPostCont->data, e.v1);
    getContourCentroid((contour*) curPostCont->data, e.v2);

    /* test the edge for interesection with the prebranch contour */
    if(edgeContourIntersect(&e,preBranch)) {
      branchType = DISCONNECT;
    }
    else {
      branchType = CONNECT;
    }
  }

  /* free the vertices */
  free(e.v1);
  free(e.v2);

  /* handle the branch differently based on branchtype */
  printf("found a branch we will %s\n", branchType == DISCONNECT ? "disconnect pre" : "connect post");

  return NULL;
}

/**
 * add a surface for two contours to an existing surface
 */
void addContoursToSurf(surface *surf, contour *contour1, contour* contour2) {
  path *p;
  list *f;
  face *f1, *f2;
  vector n, n1, n2;

  /* validate input */
  if(surf == NULL || contour1 == NULL || contour2 == NULL) {
    return;
  }

  /* orient the contours */
  if(contourArea(contour1) < 0) {
    reverseList(contour1->vertices);
  }
  if(contourArea(contour2) < 0) {
    reverseList(contour2->vertices);
  }

  /* test for the contours to have more than one vertex each */
  if(listSize(contour1->vertices) < 2 && listSize(contour2->vertices) < 2) { 
    /* can't tile */
    return;
  }
  else if(listSize(contour1->vertices) == 2 && listSize(contour2->vertices) == 2) { 
    // both have two, special case
    f = newList(LIST);
    f1 = (face*) malloc(sizeof(face));
    f1->left = (edge*) malloc(sizeof(edge));
    f1->right = (edge*) malloc(sizeof(edge));
    enqueue(f,f1);

    f2 = (face*) malloc(sizeof(face));
    f2->left = (edge*) malloc(sizeof(edge));
    f2->right = (edge*) malloc(sizeof(edge));
    enqueue(f,f2);

    // make one face
    f1->left->v1 = copyVertex((vertex*) getListNode(contour1->vertices,0)->data);
    f1->left->v2 = copyVertex((vertex*) getListNode(contour2->vertices,0)->data);
    f1->right->v1 = copyVertex((vertex*) getListNode(contour1->vertices,1)->data);
    f1->right->v2 = copyVertex((vertex*) getListNode(contour2->vertices,0)->data);
    // make the other face

    // decide which edge should complete the face via minimizing curvature
    // at edge == maximizing dot between candidate normals 
    getFaceNormalVerts((vertex*) getListNode(contour1->vertices,0)->data,
		       (vertex*) getListNode(contour1->vertices,1)->data,
		       (vertex*) getListNode(contour2->vertices,0)->data,
		       &n);
    getFaceNormalVerts((vertex*) getListNode(contour1->vertices,1)->data,
		       (vertex*) getListNode(contour2->vertices,1)->data,
		       (vertex*) getListNode(contour2->vertices,0)->data,
		       &n1);
    getFaceNormalVerts((vertex*) getListNode(contour1->vertices,0)->data,
		       (vertex*) getListNode(contour2->vertices,0)->data,
		       (vertex*) getListNode(contour2->vertices,1)->data,
		       &n2);

    if(dot(n1,n) < dot(n2,n)) {
      f2->left->v1 = copyVertex((vertex*) getListNode(contour1->vertices,0)->data);
      f2->left->v2 = copyVertex((vertex*) getListNode(contour2->vertices,0)->data);
      f2->right->v1 = copyVertex((vertex*) getListNode(contour1->vertices,0)->data);
      f2->right->v2 = copyVertex((vertex*) getListNode(contour2->vertices,1)->data);
    }
    else {
      f2->left->v1 = copyVertex((vertex*) getListNode(contour1->vertices,1)->data);
      f2->left->v2 = copyVertex((vertex*) getListNode(contour2->vertices,0)->data);
      f2->right->v1 = copyVertex((vertex*) getListNode(contour1->vertices,1)->data);
      f2->right->v2 = copyVertex((vertex*) getListNode(contour2->vertices,1)->data);

    }

    addFaceListToSurface(surf,f);
    freeFaceList(f);
  }
  else if(listSize(contour1->vertices) < 2) { 
    /* tile if contour1 has one vertex */
    if(listSize(contour1->vertices) < 1) {
      return;
    }

    /* tile a one vertex surface */
    f = tileContourToSingleVertex(contour2,(vertex*)getListNode(contour1->vertices,0)->data);
    addFaceListToSurface(surf,f);
    freeFaceList(f);
  }
  else if(listSize(contour2->vertices) < 2) { 
    /* tile if contour2 has one vertex */
    if(listSize(contour2->vertices) < 1) {
      return;
    }

    /* tile a one vertex surface */
    f = tileContourToSingleVertex(contour1,(vertex*)getListNode(contour2->vertices,0)->data);
    addFaceListToSurface(surf,f);
    freeFaceList(f);
  }
  else {
    /* get the path and add the faces to the surface, then free the path */
    p = getMinPath(contour1,contour2);

/*     listNode *ln; */
/*     face *fc; */
/*     for(ln = getListNode(p->faceList,0); ln; ln = (listNode*) ln->next) { */
/*       fc = (face*) ln->data; */
/*       fprintf(stdout,"%d(%f) %d(%f) %d(%f) %d(%f)\n", */
/* 	      fc->left->v1->label,fc->left->v1->z, */
/* 	      fc->left->v2->label,fc->left->v2->z, */
/* 	      fc->right->v1->label,fc->right->v1->z, */
/* 	      fc->right->v2->label,fc->right->v2->z); */
/*     } */

    addFaceListToSurface(surf,p->faceList);
    freePath(p);
  }
}

/**
 * get face list between two contours. just a wrapper around getMinPath
 */
list *getFaceList(contour *contour1, contour *contour2) {
  path *p;
  list *faceList;

  /* validate input */
  if(contour1 == NULL || contour2 == NULL
     || contour1->vertices == NULL || contour2->vertices == NULL
     || listSize(contour1->vertices) == 0
     || listSize(contour2->vertices) == 0) {
    return NULL;
  }

  p = getMinPath(contour1,contour2);
  if(p == NULL) return NULL;

  faceList = p->faceList;
  free(p);

  return faceList;
}

/**
 * get a the minimum path representing a surface between two contours
 */
path *getMinPath(contour *contour1, contour *contour2) {
  graph *G = buildPlanarGraph(contour1, contour2);
  int numPaths = 0, i,minPathInd = 0;
  path **paths = NULL, **forwPaths = NULL, **backPaths = NULL;

  /* get the paths based on the closure case */
  if(G->closureCase == BOTH_CLOSED) {
    /* get the shortest paths from each starting vertex */
    paths = allPaths(G);
    numPaths = G->m+1;
  }
  else if(G->closureCase == ONE_CLOSED) {
    /* get the shortest paths from each starting vertex, leave one out
       and run both directions */
    numPaths = 2*(G->m+1);
    paths = (path**) malloc(numPaths*sizeof(path));

    forwPaths = allPathsNoclosure(G);

    /* reverse the open contour  and get the second set of paths */
    if(contour1->closed == OPEN) {
      contour1 = cloneContour(contour1);
      reverseList(contour1->vertices);

      /* free the last graph and build a new one */
      freePlanarGraph(G);
      G = buildPlanarGraph(contour1, contour2);
      deleteContour(contour1);
    }
    else {
      contour2 = cloneContour(contour2);
      reverseList(contour2->vertices);

      /* free the last graph and build a new one */
      freePlanarGraph(G);
      G = buildPlanarGraph(contour1, contour2);
      deleteContour(contour2);
    }

    /* get the second set of paths */
    backPaths = allPathsNoclosure(G);

    /* copy the back paths over */
    for(i = 0; i < numPaths; i++) {
      if(i < G->m+1) {
	paths[i] = forwPaths[i];
      }
      else {
	paths[i] = backPaths[i-(G->m+1)];
      }
    }

    /* free path arrays */
    free(forwPaths);
    free(backPaths);
  }
  else {
    /* get the shortest paths from two possible starting vertices */
    numPaths = 2;
    paths = (path**) malloc(numPaths*sizeof(path));
    paths[0] = singlePath(G,0,G->n,OPEN);

    /* reverse contour2 and get the second path */
    contour2 = cloneContour(contour2);
    reverseList(contour2->vertices);

    /* free the last graph and build a new one */
    freePlanarGraph(G);
    G = buildPlanarGraph(contour1, contour2);
    deleteContour(contour2);

    /* get the second path */
    paths[1] = singlePath(G,0,G->n,OPEN);
  }

  /* search the generated paths for the shortest path, free path if needed */
  path *minPath = paths[0], *tmpPath;

//  if(SR_DEBUG) {
//    fprintf(stderr,"path %d cost = %lf\n",0,paths[0]->cost);
//  }

  static int curPath = 0;
  char name[100]; 

  if(savePaths) {
    //for(i = 0; i < numPaths; i++) {
    sprintf(name,"path%d.obj",curPath++);
    writeOBJ(buildSurfaceFromFaceList(paths[0]->faceList),name);
    //fprintf(stderr,"path %d cost = %lf\n",i,paths[i]->cost);
    //}
  }

  for(i = 1; i < numPaths; i++) {
    if(savePaths) {
      //for(i = 0; i < numPaths; i++) {
	sprintf(name,"path%d.obj",curPath++);
	writeOBJ(buildSurfaceFromFaceList(paths[i]->faceList),name);
	//fprintf(stderr,"path %d cost = %lf\n",i,paths[i]->cost);
	//}
    }

//    if(SR_DEBUG) {
//      fprintf(stderr,"path %d cost = %lf\n",i,paths[i]->cost);
//    }

    if(paths[i]->cost < minPath->cost) {
      tmpPath = minPath;
      minPath = paths[i];
      minPathInd = i;
      freePath(tmpPath);
    }
    else {
      freePath(paths[i]);
    }
  }

  if(printNow) {
    printGraph(G,minPathInd,contour1->closed);
  }

  /* free the path array and graph */
  free(paths);
  freePlanarGraph(G);

  if(SR_DEBUG) {
    fprintf(stderr,"min path cost = %lf\n",minPath->cost);
  }

  // check number of edge crossings
//  int numCrossings = 0;
//  listNode *ln, *ln2;
//  face *f;
//  for(ln = getListNode(minPath->faceList,0); ln; ln = (listNode*) ln->next) {
//    for(ln2 = getListNode(minPath->faceList,0); ln2; ln2 = (listNode*) ln2->next) {
//      f = (face*) ln2->data;
//
//      /* check each edge for self intersecting surface */
//      if(edgeFaceIntersect((face*) ln->data, f->left)
//	 || edgeFaceIntersect((face*) ln->data, f->right)) {
//	numCrossings++;
//      }
//    }
//  }
//  printf("chose path with %d crossings.\n",numCrossings);


  return minPath;
}

/**
 * build planar graph for fuchs' algorithm to operate on this graph is
 * a 2m+1 x n+1 replicated structure for easy path searching without
 * cycles. for most closure cases the entire graph will not be used,
 * but we build it anyhow.
 */
graph *buildPlanarGraph(contour *contour1, contour *contour2) {
  graph *G = (graph*) malloc(sizeof(graph));
  edge *e;
  vertex *v1,*v2;
  contour *tmpContour;
  listNode *c1v = NULL, *c2v = NULL; 

  /* determine the closure case */
  if(contour1->closed && contour2->closed) {
    G->closureCase = BOTH_CLOSED;
  }
  else if(contour1->closed) {
    G->closureCase = ONE_CLOSED;
  }
  else if(contour2->closed) {
    /* swap them so the closed contour is contour1 */
    tmpContour = contour2;
    contour2 = contour1;
    contour1 = tmpContour;

    G->closureCase = ONE_CLOSED;
  }
  else {
    G->closureCase = NONE_CLOSED;
  }

  /* size of the graph */
  int m = listSize(contour1->vertices);
  int n = listSize(contour2->vertices);
  int i, j;

  /* assign the sizes */
  G->m = m;
  G->n = n;

  /* allocate the graph */
  G->v = (edge***) malloc((2*m+1)*sizeof(edge**));
  G->pi = (edge***) malloc((2*m+1)*sizeof(edge**));
  G->wh = (double**) malloc((2*m+1)*sizeof(double*));
  G->wv = (double**) malloc((2*m+1)*sizeof(double*));
  G->h = (heap_it_handler**) malloc((2*m+1)*sizeof(heap_it_handler*));

  /* iterate over nodes, first allocate, then assign vertices and edges */
  for(i = 0; i < 2*m+1; i++) {
    G->v[i] = (edge**) malloc((n+1)*sizeof(edge*));
    G->pi[i] = (edge**) malloc((n+1)*sizeof(edge*));
    G->wh[i] = (double*) malloc((n+1)*sizeof(double));
    G->wv[i] = (double*) malloc((n+1)*sizeof(double));
    G->h[i] = (heap_it_handler*) malloc((n+1)*sizeof(heap_it_handler));

    if(!c1v) {
      c1v = getListNode(contour1->vertices,0);
    }
    v1 = (vertex*) c1v->data;
    c2v = NULL;

    /* assign the edges by topology and weights from geometry  */
    for(j = 0; j < n+1; j++) {
      /* set the predecessor to nil */
      G->pi[i][j] = NULL;

      /* allocate this edge */
      e = (edge*) malloc(sizeof(edge));

      /** assign the vertices for this edge **/

      /* get and copy vertex 1 */
      //v1 = (vertex*) getListNode(contour1->vertices,i%m)->data;
      e->v1 = createVertex();
      *(e->v1) = *v1; /* copy data */
      e->v1->number = i; /* modify the vertex number */

      /* get and copy vertex 2 */
      if(!c2v) {
	c2v = getListNode(contour2->vertices,0);
      }

      v2 = (vertex*) c2v->data;
      c2v = (listNode*) c2v->next;

      e->v2 = createVertex();
      *(e->v2) = *v2;
      e->v2->number = j;

      /* calculate the length */
      e->length = len(e->v1,e->v2);

      G->v[i][j] = e;

      /* calculate the edges for this face and the cost */
      if(i == 0 || j == 0) { /* not enough info to calc weights here */
	continue;
      }
      else if(i > 0 && i < m && j > 0 && j < n) {
	/* weight of face from i-1,j to i,j */
	G->wv[i][j]=getFaceCost(G->v[i-1][j]->v1,
				G->v[i][j]->v2,
				G->v[i][j]->v1);

	/* weight of face from i,j-1 to i,j */
	G->wh[i][j]=getFaceCost(G->v[i][j-1]->v2,
				G->v[i][j]->v2,
				G->v[i][j]->v1);

      }
      else if(i == m || j == n) { /* assign the 0th index nodes */
	/* weight of face from i-1,j to i,j */
	G->wv[i%m][j] = 
	G->wv[i][j%n] = 
	G->wv[i%m][j%n] = 
	G->wv[i][j] =
	  getFaceCost(G->v[i-1][j]->v1, G->v[i][j]->v2, G->v[i][j]->v1);

	/* weight of face from i,j-1 to i,j */
	G->wh[i%m][j] = 
	G->wh[i][j%n] = 
	G->wh[i%m][j%n] = 
	G->wh[i][j] =
	  getFaceCost(G->v[i][j-1]->v2, G->v[i][j]->v2, G->v[i][j]->v1);
      }
      else { /* we've already calculated these */
	G->wv[i][j] = G->wv[i%m][j%n];
	G->wh[i][j] = G->wh[i%m][j%n];
      }

      if(SR_DEBUG) {
	if(isnan(G->wv[i][j])) {
	  G->wv[i][j] = SR_MED;
	  fprintf(stdout,
		  "vert w nan: i=%d j=%d == [%0.20f %0.20f %0.20f] [%0.20f %0.20f %0.20f] [%0.20f %0.20f %0.20f]\n",
		  i,j,
		  G->v[i-1][j]->v1->x, G->v[i-1][j]->v1->y, G->v[i-1][j]->v1->z,
		  G->v[i][j]->v2->x, G->v[i][j]->v2->y, G->v[i][j]->v2->z,
		  G->v[i][j]->v1->x, G->v[i][j]->v1->y, G->v[i][j]->v1->z
		  );
	}
	if(isnan(G->wh[i][j])) {
	  G->wh[i][j] = SR_MED;
	  fprintf(stdout,
		  "horz w nan: i=%d j=%d == [%0.20f %0.20f %0.20f] [%0.20f %0.20f %0.20f] [%0.20f %0.20f %0.20f]\n",
		  i,j,
		  G->v[i-1][j]->v2->x, G->v[i-1][j]->v2->y, G->v[i-1][j]->v2->z,
		  G->v[i][j]->v2->x, G->v[i][j]->v2->y, G->v[i][j]->v2->z,
		  G->v[i][j]->v1->x, G->v[i][j]->v1->y, G->v[i][j]->v1->z
		  );
	}
      }

      G->h[i][j] = NULL;
    }

    // increment contour 1 vertex pointer
    c1v = (listNode*) c1v->next;  
  }

  return G;
}

/**
 * free the allocated contents of a planar graph
 */
void freePlanarGraph(graph *G) {
  int i,j;

  /* iterate over nodes, freeing each */
  for(i = 0; i < 2*G->m+1; i++) {
    for(j = 0; j < G->n+1; j++) {
      free(G->v[i][j]->v1);
      free(G->v[i][j]->v2);
      free(G->v[i][j]);
    }
    free(G->v[i]);
    free(G->wv[i]);
    free(G->wh[i]);
    free(G->pi[i]);
    free(G->h[i]);
  }
  free(G->v);
  free(G->wv);
  free(G->wh);
  free(G->pi);
  free(G->h);
  free(G);
}

/**
 * find all shortest paths from the m starting nodes to the n ending
 * nodes implementation of algorithm presented in fuchs (1977) for two
 * closed contours
 */
path **allPaths(graph *G) {
  path **paths = (path**) malloc((G->m+1)*sizeof(path*));

  /* get the paths */
  paths[0] = singlePath(G,0,G->n+1,CLOSED);
  paths[G->m] = singlePath(G,G->m,G->n+1,CLOSED);
  pathsBetween(G,paths,0,G->m,G->n+1,CLOSED);

  return paths;
}

/**
 * find all shortest paths from the m starting nodes to the n ending
 * nodes, but dont wrap to the begining of contour2. implementation of
 * algorithm presented in fuchs (1977) for first contour closed and
 * second open
 */
path **allPathsNoclosure(graph *G) {
  path **paths = (path**) malloc((G->m+1)*sizeof(path*));

  /* get the paths */
  paths[0] = singlePath(G,0,G->n,CLOSED);
  paths[G->m] = singlePath(G,G->m,G->n,CLOSED);
  pathsBetween(G,paths,0,G->m,G->n,CLOSED);

  return paths;
}

/**
 * find two paths, one from the first vertices of each contour, the
 * second from the last of the first contour and the first of the
 * second contour . implementation of algorithm presented in fuchs
 * (1977) for both contours open
 */
path **twoPathsNoclosure(graph *G) {
  path **paths = (path**) malloc(2*sizeof(path*));

  /* get the paths */
  paths[0] = singlePath(G,0,G->n,OPEN);
  paths[1] = singlePath(G,G->m-1,G->n,OPEN);

  return paths;
}

/**
 * find the paths between two bounding paths,
 * implementation of algorithm presented in fuchs (1977)
 */
void pathsBetween(graph *G, path **paths, int botPath,
		  int topPath, int length, int close) {
  int midPath = (topPath+botPath)/2;

  if(botPath < midPath) {
    paths[midPath] = singlePath(G,midPath,length,close);
    pathsBetween(G,paths,botPath,midPath,length,close);
    pathsBetween(G,paths,midPath,topPath,length,close);
  }
}

/**
 * run dijstra's single source shortest path algorithm
 * implemntation is specific to the fuchs (1977) planar graph
 */
path *singlePath(graph *G, int source, int length, int close) {
  int i,j;
  edge *u,*v;
  //listNode *ln;
  //listNode *hn,*vn;
  double wu,wv,wvNew;
  heap_it_handler adjEl;
  path *p;

  //int node;
  int ind;

#ifdef MEX
  //mexPrintf("finding a single path from source %d\n", source);
#else
    //fprintf(stdout,"finding a single path from source %d\n", source);
#endif

/*   if(SR_DEBUG) { */
/*     fprintf(stderr,"source = %d\n",source); */
/*   } */

  /* heap for the vertex list */
  //list *vertexHeap = newList(MINHEAP);
  heap_t vertexHeap;
  vertexHeap.heap = NULL;

  heap_init(&vertexHeap, length*(source+G->m+1));

  /* construct the vertex heap, with weights */
  for(i = source; i < source+G->m+1; i++) {
    for(j = 0; j < length; j++) {
      /* reset step counts */
      G->v[i][j]->lastVert = 0;
      G->v[i][j]->lastHorz = 0;

      /* reset predecessor */
      G->pi[i][j] = NULL;

      if(i == source && j == 0) { /* insert the source node as zero distance */
	//insertHeapNode(vertexHeap,0,G->v[i][j]);
	G->h[i][j] = heap_insert(&vertexHeap,G->v[i][j],0.0);
      }
      else { /* assume infinite distance to begin with */
	//insertHeapNode(vertexHeap,INF,G->v[i][j]);
	G->h[i][j] = heap_insert(&vertexHeap,G->v[i][j],INF);
      }
    }
  }

  /* find the shortest path to each node */
  ind = 0;
  while(!heap_isempty(&vertexHeap)) {
    ind++;

    //ln = getHeapTop(vertexHeap);
    //u = (edge*) ln->data;
    //w = ln->value;
    wu = (double) heap_topkey(&vertexHeap);
    u = (edge*) *((edge**)heap_removetop(&vertexHeap));

/*     if(SR_DEBUG) fprintf(stderr,"-----------\nexamining v%d%d=%g\n", */
/* 		      u->v1->number,u->v2->number,w); */

    //fprintf(stderr,"v1=%d v2=%d lv=%d lh=%d w=%lf\n",u->v1->number,u->v2->number,u->lastVert,u->lastHorz,wu);

    /* test for end of path */
    if(u->v1->number == source + G->m - (close==CLOSED?0:1)
       && u->v2->number == length-1) {
      //free(ln);
      //freeList(vertexHeap);
      heap_destroy(&vertexHeap);

      p = buildPath(G,source,length,close);

      if(SR_DEBUG) {
	fprintf(stderr,"found path from %d of length %d and cost %lf\n",
		source,length,p->cost);
      }

      return p;
    }

    /* remove the adjacent vertices and relax them */
    /* NOTE: needs to be optimized */

    /* relax the vertical edge */
    /* test for within range */
    if(u->v1->number < source + G->m) {
      //node = findAdj(G,vertexHeap,u,VERT);
      adjEl = G->h[u->v1->number+1][u->v2->number];
      
/*       if(SR_DEBUG && node > 0) { */
/* 	listNode *tmpln = getListNode(vertexHeap,node); */
/* 	fprintf(stderr,"relaxing vert neigh (%g) v%d%d=%g to ", */
/* 		G->wv[((edge*)tmpln->data)->v1->number] */
/* 		     [((edge*)tmpln->data)->v2->number], */
/* 		((edge*)tmpln->data)->v1->number, */
/* 		((edge*)tmpln->data)->v2->number, */
/* 		tmpln->value); */
/*       } */

      
      wv = adjEl->key;
      v = adjEl->user_data;

      //vn = relax(G,ln,removeListNode(vertexHeap,node),VERT);
      wvNew = relax(G,u,wu,v,wv,VERT);

      /* test for this edge being the one that completes a straight path,
	 thus causing a edge repeat if it is, penalize it heavily */
      if(testManhattan && u->lastHorz >= (G->m - (close==CLOSED?0:1))-1) {
	wvNew += manhattanPenalty*wvNew;
      }

      if(wv - wvNew > 0) {
	//insertHeapNode(vertexHeap, vn->value, vn->data);
	adjEl->key = wvNew;
	v->lastVert = 0;
	v->lastHorz = u->lastHorz+1;
	heapfloat_heapify(&vertexHeap, adjEl->p_ctrl.heapdata);
      }

/*       if(SR_DEBUG) { */
/* 	fprintf(stderr, "%g\n", vn == NULL ? -1.0: vn->value); */
/*       } */
    }

    /* relax the horizontal edge */
    if(u->v2->number < length-1) { /* test for within range */
      //node = findAdj(G,vertexHeap,u,HORZ);
      adjEl = G->h[u->v1->number][u->v2->number+1];

/*       if(SR_DEBUG && node > 0) { */
/* 	listNode *tmpln = getListNode(vertexHeap,node); */
/* 	fprintf(stderr,"relaxing horz neigh (%g) v%d%d=%g to ", */
/* 		G->wh[((edge*)tmpln->data)->v1->number] */
/* 		     [((edge*)tmpln->data)->v2->number], */
/* 		((edge*)tmpln->data)->v1->number, */
/* 		((edge*)tmpln->data)->v2->number, */
/* 		tmpln->value); */
/*       } */

      wv = adjEl->key;
      v = adjEl->user_data;
      
      //hn = relax(G,ln,removeListNode(vertexHeap,node),HORZ);
      wvNew = relax(G,u,wu,v,wv,HORZ);


      /* test for this edge being the one that makes this a true manhattan path
	 if it is, penalize it heavily */
      if(testManhattan && u->lastVert >= length-2) {
	wvNew += manhattanPenalty*wvNew;
      }

      if(wv - wvNew > 0) {
	//insertHeapNode(vertexHeap, hn->value, hn->data);
	adjEl->key = wvNew;
	v->lastVert = u->lastVert+1;
	v->lastHorz = 0;
	heapfloat_heapify(&vertexHeap, adjEl->p_ctrl.heapdata);
      }

/*       if(SR_DEBUG) { */
/* 	fprintf(stderr, "%g\n", hn == NULL ? -1.0 : hn->value); */
/*       } */

    }


    /* free list node */
    //free(ln);
  }

  //freeList(vertexHeap);
  heap_destroy(&vertexHeap);

  return NULL;
}

/**
 * relax function
 */
double relax(graph *G, edge *u, double wu, edge *v, double wv, int dir) {
  if(u == NULL || v == NULL) return wv;

  int vi = v->v1->number;
  int vj = v->v2->number;
  int i,j,tmp,numRepeats=0;

  /* shortest dist to this node through u */
  double vdist = wu + (dir == HORZ ? G->wh[vi][vj] : G->wv[vi][vj]);

  /* check if this path is shortest so far*/
  if(wv > vdist) {
    /* adjustment the weight if this edge already appeared */
    if(testRepeat) {
      // run backward through the current path, looking for this edge */
      i = u->v1->number;
      j = u->v2->number;
      while(G->pi[i][j]) {
	if(verticesEqual(G->pi[i][j]->v1,v->v1) 
	   && verticesEqual(G->pi[i][j]->v2,v->v2) ) {
	  numRepeats++;
	}
	
	tmp = G->pi[i][j]->v1->number;
	j = G->pi[i][j]->v2->number;
	i = tmp;
      }

      vdist += numRepeats*repeatedEdgePenalty*vdist;
    }

    if(wv > vdist) {
      wv = vdist;
      G->pi[vi][vj] = u;
    }
  }

  return wv;
}

/**
 * relax function
 */
listNode *relaxLN(graph *G, listNode *u, listNode *v, int dir) {
  if(u == NULL || v == NULL) return NULL;

  /* get u indices */
  edge *ue = (edge*) u->data;

  /* get v indices */
  edge *ve = (edge*) v->data;
  int vi = ve->v1->number;
  int vj = ve->v2->number;

  /* shortest dist to this node through u */
  double vdist = u->value + (dir == HORZ ? G->wh[vi][vj] : G->wv[vi][vj]);

  /* check if this path is shortest so far*/
  if(v->value > vdist) {
    v->value = vdist;
    G->pi[vi][vj] = ue;
  }

  return v;
}

/**
 * build an actual edge list from a path in the planar graph
 */
list *makeEdgeList(path *p) {
  list *edgeList = newList(LIST);
  listNode *i;
  face *f;

  /* iterate over faces, adding novel edges to list */
  for(i = getListNode(p->faceList,0); i; i = (listNode*) i->next) {
    f = (face*) i->data;

    /* add the face edges, if they don't exist */
    if(f->left && !listContains(edgeList, (void*) f->left)) {
      enqueue(edgeList,f->left);
    }
    if(f->right && !listContains(edgeList, (void*) f->right)) {
      enqueue(edgeList,f->right);
    }
  }

  return edgeList;
}

/**
 * builds a path from a graph
 */
path *buildPath(graph *G, int source, int length, int close) {
  path *p = (path*) malloc(sizeof(path));
  int i = source + G->m - (close==CLOSED?0:1);
  int j = length-1;
  face *f;
  int tmp;

  /* create the face list */
  p->faceList = newList(STACK);
  p->cost = p->baseCost = 0;

  /* follow the path backward, adding each face by its predecessor */
  while(G->pi[i][j]) {
    /* create a face */
    f = (face*) malloc(sizeof(face));

    /* assign the left edge */
    f->left = copyEdge(G->pi[i][j]);
    f->left->v1->number %= G->m;
    f->left->v2->number %= G->n;

    /* assign the right edge */
    f->right = copyEdge(G->v[i][j]);
    f->right->v1->number %= G->m;
    f->right->v2->number %= G->n;

    /* assign the vertices from the edges */
    f->v1 = f->left->v1->number;
    f->v2 = f->left->v2->number;

    /* decide the three face vertices */
    if(verticesEqual(f->left->v1,f->right->v1)) {
      f->v3 = f->right->v2->number;

      /* get the cost */
      f->cost = getFaceCost(f->left->v1,f->left->v2,f->right->v2);
    }
    else {
      f->v3 = f->right->v1->number;

      /* get the cost */
      f->cost = getFaceCost(f->left->v1,f->left->v2,f->right->v1);
    }

    p->baseCost += f->cost;

    /* add the face cost to the path cost */
/*     if(SR_DEBUG) { */
/*       fprintf(stderr,"adding face cost %lf for a total cost of %lf\n",f->cost,p->cost); */
/*     } */

    /* add the edge */
    push(p->faceList,f);

    /* update the indices */
    tmp = G->pi[i][j]->v1->number;
    j = G->pi[i][j]->v2->number;
    i = tmp;
  }
  p->cost = p->baseCost;

  if(testEdgeCrossing) {
    /* test that each face doesnt intersect any others */
    listNode* ln, *ln2;
    int numCrossings = 0;
    for(ln = getListNode(p->faceList,0); ln; ln = (listNode*) ln->next) {
      for(ln2 = getListNode(p->faceList,0); ln2; ln2 = (listNode*) ln2->next) {
	f = (face*) ln2->data;

	/* check each edge for self intersecting surface */
	if(edgeFaceIntersect((face*) ln->data, f->left)
	   || edgeFaceIntersect((face*) ln->data, f->right)) {
//	  	if(SR_DEBUG) {
//	  	  fprintf(stderr,"face-face intersection found!\nsource=%d,i=%d,j=%d\n",
//	  		  source,i,j);
//	  	}
	  numCrossings++;
	}
      }
    }
    //    printf("nc=%d cost=%lf\n",numCrossings,p->cost);
    p->cost += edgeCrossingPenalty*numCrossings*p->baseCost;
  }

  return p;
}

/**
 * frees the contents of a path
 */
void freePath(path *p) {
  /* free the face list */
  freeFaceList(p->faceList);
  free(p);
}

/**
 * frees the contents of a facelist
 */
void freeFaceList(list *faceList) {
  listNode *ln;
  face *f;

  for(ln = getListNode(faceList,0); ln; ln = (listNode*) ln->next) {
    f = (face*) ln->data;

    /* free all contents */
    free(f->left->v1);
    free(f->left->v2);
    free(f->right->v1);
    free(f->right->v2);

    free(f->left);
    free(f->right);

    free(f);
  }

  freeList(faceList);
}

/**
 * finds the index of a node in a heap corresponding to the lower
 * neighbor of a given node
 */
int findAdj(graph *G, list *vertexHeap, edge *u, int dir) {
  int i = dir == VERT ? (u->v1->number+1) : u->v1->number;
  int j = dir == VERT ? u->v2->number : (u->v2->number+1);

  listNode *k;
  int ind;
  edge *v;

  /* cycle through edges in heap till we find the desired edge */
  for(k = getListNode(vertexHeap,0), ind = 0; k; k = (listNode*) k->next, ind++) {
    v = (edge *) k->data;

    /* test for the desired node */
    if(v->v1->number == i && v->v2->number == j) {
      return ind;
    }
  }

  /* return flibvpure */
  return -1;
}

/**
 * tests if a contour is refered to in the contour adjacency
 * list of any of a list of contours
 */
int hasBackwardAdjacentContours(list *contourList, contour *cont) {
  listNode *ln, *ln2;
  contour *curContour;

  for(ln = getListNode(contourList,0); ln; ln = (listNode*) ln->next) {
    curContour = (contour*) ln->data;

    for(ln2 = getListNode(curContour->adjacentContours,0); ln2; ln2 = (listNode*) ln2->next) {
      if(ln2->data == cont) return TRUE;
    }
  }

  return FALSE;
}

/**
 * gets the vertices in a contour that are closest to one contour of a set
 */
void getVerticesForContour(contour *splitCont, contour *refCont,
			   list *contours, list *verts) {
  vertex **centroids = (vertex**) malloc((1+listSize(contours))*sizeof(vertex*));
  int i, numCont = 1+listSize(contours);
  listNode *ln;

  /* build the centroid list */

  /* ref contour centroid */
  centroids[0] = createVertex();
  getContourCentroid(refCont,centroids[0]);

  /* build from the other contours */
  for(ln = getListNode(contours,0), i = 1; ln; ln = (listNode*) ln->next, i++){
    centroids[i] = createVertex();
    getContourCentroid((contour*) ln->data,centroids[i]);
  }

  /* iterate over points in the split contour, finding the closest
     centroid to it, then adding vertices closest to the reference
     contour*/
  for(ln = getListNode(splitCont->vertices,0); ln; ln = (listNode*) ln->next) {

  }

  /* free the centroids */
  for(i = 0; i < numCont; i++) {
    free(centroids[i]);
  }
  free(centroids);
}

/**
 * gets the irregularity of the angles in the face, a measure that goes from
 * 0 to one
 */
double getFaceAngleIrregularity(vertex *a, vertex *b, vertex *c) {
  double bestAngle = 60.0*SR_PI/180.0;
  double maxAngleDiff = sqrt(2*pow(0.0-60.0,2)+2*pow(180.0-60.0,2))*SR_PI/180.0;

  return sqrt(
	      pow(angleBetweenSegments(a,b,c)-bestAngle,2) +
	      pow(angleBetweenSegments(b,c,a)-bestAngle,2) +
	      pow(angleBetweenSegments(c,a,b)-bestAngle,2))
    / maxAngleDiff;
}

/**
 * gets the cost of a face made of three vertices
 */
double getFaceCost(vertex *a, vertex *b, vertex *c) {
  if(testAngles == FALSE) {
    return getFaceArea(a,b,c) ;
  }

  return (1 + anglePenalty*getFaceAngleIrregularity(a,b,c)) 
    * getFaceArea(a,b,c) ;
}

/**
 * gets the area of a face made of three vertices
 * i know you hate me
 */
inline double getFaceArea(vertex *a, vertex *b, vertex *c) {
  return 0.5*sqrt(-2*a->x*a->x*b->z*c->z-2*a->x*b->y*b->y*c->x-2*a->x*c->y*c->y*b->x-2*a->y*a->y*b->z*c->z-2*a->x*a->x*b->y*c->y+2*b->x*a->y*c->x*b->y+2*b->x*c->y*c->x*a->y-2*b->x*c->y*c->x*b->y-2*a->y*b->z*b->y*a->z+2*a->y*b->z*b->y*c->z+2*a->y*b->z*c->y*a->z+2*a->y*c->z*b->y*a->z-2*a->y*c->z*c->y*a->z+2*a->y*c->z*c->y*b->z+2*b->y*a->z*c->y*b->z-2*b->y*c->z*c->y*b->z-2*a->x*b->z*b->x*a->z+2*a->x*b->z*b->x*c->z+2*a->x*b->z*c->x*a->z+2*a->x*c->z*b->x*a->z-2*a->x*c->z*c->x*a->z+2*a->x*c->z*c->x*b->z+2*b->x*a->z*c->x*b->z+2*b->x*c->z*c->x*a->z-2*b->x*c->z*c->x*b->z-2*b->x*b->x*a->y*c->y-2*b->x*a->y*a->y*c->x-2*c->x*c->x*a->y*b->y-2*a->y*b->z*b->z*c->y-2*a->y*c->z*c->z*b->y-2*b->y*b->y*a->z*c->z-2*b->y*a->z*a->z*c->y+2*a->x*c->y*c->x*b->y-2*a->x*c->y*c->x*a->y+2*a->x*c->y*b->x*a->y+2*a->x*b->y*c->x*a->y+2*a->x*b->y*b->x*c->y-2*a->x*b->y*b->x*a->y+a->x*a->x*b->y*b->y+a->x*a->x*c->y*c->y+b->x*b->x*a->y*a->y+b->x*b->x*c->y*c->y+c->x*c->x*a->y*a->y+c->x*c->x*b->y*b->y+a->y*a->y*b->z*b->z+a->y*a->y*c->z*c->z+b->y*b->y*a->z*a->z+b->y*b->y*c->z*c->z+c->y*c->y*a->z*a->z+c->y*c->y*b->z*b->z+a->x*a->x*b->z*b->z+a->x*a->x*c->z*c->z+b->x*b->x*a->z*a->z+b->x*b->x*c->z*c->z+c->x*c->x*a->z*a->z+c->x*c->x*b->z*b->z-2*c->y*c->y*a->z*b->z+2*b->y*c->z*c->y*a->z-2*a->x*b->z*b->z*c->x-2*a->x*c->z*c->z*b->x-2*b->x*b->x*a->z*c->z-2*b->x*a->z*a->z*c->x-2*c->x*c->x*a->z*b->z);
}

/** 
 * get the angle between two segments
 */

/**
 * for debugging, output the graph nodes for a path
 */
void printGraph(graph *g, int source, int close) {
  int length = g->n-(close==CLOSED?0:1);
  int i = source + g->m - (close==CLOSED?0:1);
  int j = length-1;
  int tmp;

  /* print the path */
  while(g->pi[i][j]) {
    fprintf(stdout,"%d %d\n",g->v[i][j]->v1->number,g->v[i][j]->v2->number);
    tmp = g->pi[i][j]->v1->number;
    j = g->pi[i][j]->v2->number;
    i = tmp;
  }
  fprintf(stdout,"\n");

  /* print the vertical weights */
  for(i = 0; i < 2*g->m+1; i++) {
    for(j = 0; j < g->n+1; j++) {
      fprintf(stdout,"%f ",g->wv[i][j]);
    }
    fprintf(stdout,"\n");
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"\n");

  /* print the horizontal weights */
  for(i = 0; i < 2*g->m+1; i++) {
    for(j = 0; j < g->n+1; j++) {
      fprintf(stdout,"%f ",g->wh[i][j]);
    }
    fprintf(stdout,"\n");
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"\n");

  /* print the vertices asscoiated with the graph nodes */
  for(i = 0; i < g->m; i++) {
    fprintf(stdout,"%lf %lf %lf\n",
	    g->v[i][0]->v1->x,g->v[i][0]->v1->y,g->v[i][0]->v1->z);
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"\n");

  for(i = 0; i < g->n; i++) {
    fprintf(stdout,"%lf %lf %lf\n",
	    g->v[0][i]->v2->x,g->v[0][i]->v2->y,g->v[0][i]->v2->z);
  }

}

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/tile.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
