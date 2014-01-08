/*****************************************************************************
 * libsrUtil.c is the source file for utility functions for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#include"libsrUtil.h"

int SR_VERBOSE = 0;

/* whether to resample contours or not */
int resample = FALSE;

/* distance to use for contour resampling */
double resampleDistance = 2.0; /* percent of thinnest linear dimension */
double lastResampleDist = 2.0; /* last distance used */
int resampleAbs = FALSE; /* whether to interpret resample distance as absolute */
int resampleEqual = FALSE; /* whether to resmaple the contours to equal vertex spacing same as between slices */

// number of slices to increment between recons
int skipSlices = 1;

/** contour utils **/

/**
 * create a contour
 * returns a pointer to a new contour, or null, if failure
 */
contour *createContour() {
  contour *cont = (contour*) malloc(sizeof(contour));

  if(cont == NULL) return NULL;

  /* allocate the vertices and contour connectivity list */
  cont->vertices = newList(LIST);
  cont->adjacentContours = newList(LIST);
  cont->adjacentBackwardContours = newList(LIST);

  /* assume the contour is open by default */
  cont->closed = OPEN;

  cont->origin = NULL;

  return cont;
}

/**
 * clone a contour
 * returns a pointer to a new contour copied from the one passed,
 * or null, if failure
 */
contour *cloneContour(contour *cont) {
  listNode *i;
  vertex *v, *newv;

  if(cont == NULL) return NULL;

  contour *newCont = (contour*) malloc(sizeof(contour));

  if(newCont == NULL) return NULL;

  /* copy the vertices */
  newCont->vertices = newList(LIST);
  for(i = getListNode(cont->vertices,0); i; i = (listNode*) i->next) {
    v = (vertex*) i->data;
    newv = copyVertex(v);
    enqueue(newCont->vertices, newv);
  }

  /* copy the adjacent contours */
  newCont->adjacentContours = newList(LIST);
  for(i = getListNode(cont->adjacentContours,0); i; i = (listNode*) i->next) {
    enqueue(newCont->adjacentContours, i->data);
  }

  /* copy the adjacent backward contours */
  newCont->adjacentBackwardContours = newList(LIST);
  for(i = getListNode(cont->adjacentBackwardContours,0);
      i; i = (listNode*) i->next) {
    enqueue(newCont->adjacentBackwardContours, i->data);
  }

  /* copy the open/closed state */
  newCont->closed = cont->closed;

  newCont->origin = (struct contour*) cont->origin;

  return newCont;
}

/**
 * deletes a contour
 */
void deleteContour(contour *cont) {
  listNode *i;

  if(cont == NULL) return;

  /* delete the vertices */
  for(i = getListNode(cont->vertices,0); i; i = (listNode*) i->next) {
    free(i->data);
  }
  freeList(cont->vertices);

  /* free the adjacent contour pointer list */
  freeList(cont->adjacentContours);

  /* copy the adjacent backward contour pointer list */
  freeList(cont->adjacentBackwardContours);

  free(cont);
}

/**
 * determines the area of the polygon bounded by a contour
 */
double contourArea(contour* cont) {
  listNode *i;
  vertex *v,*u;
  double area = 0;

  /* validate */
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 3) {
    return 0;
  }

  u =  (vertex*) getListNode(cont->vertices,listSize(cont->vertices)-1)->data;
  for(i = getListNode(cont->vertices,0); i; i = (listNode*) i->next) {
    v = (vertex*) i->data;

    area += u->x*v->y - v->x*u->y;

    u = v;
  }

  return area/2;
}

/**
 * get the min and max extents of the contours in a slice contour set
 */
void computeSliceContourBounds(list *slices, vector *minbounds, vector *maxbounds) {
  listNode *i,*j,*k;
  list *sl, *cv;
  vertex *v;

  minbounds->x = SR_BIG;
  minbounds->y = SR_BIG;
  maxbounds->x = -SR_BIG;
  maxbounds->y = -SR_BIG;
  for(i = getListNode(slices,0); i; i = (listNode*) i->next) {
    sl = (list*) i->data;

    for(j = getListNode(sl,0); j; j = (listNode*) j->next) {
      cv = ((contour*) j->data)->vertices;

      for(k = getListNode(cv,0); k; k = (listNode*) k->next) {
	v = (vertex*) k->data;

	/* check bounds */
	if(v->x < minbounds->x) minbounds->x = v->x;
	if(v->x > maxbounds->x) maxbounds->x = v->x;
	if(v->y < minbounds->y) minbounds->y = v->y;
	if(v->y > maxbounds->y) maxbounds->y = v->y;
      }
    }
  }

}

/**
 * clones a pair of slices, maintaining adjacency
 */
void cloneSlices(list *slice1, list *slice2,
		 list **newSlice1, list **newSlice2) {

  listNode *i,*j;
  int ind;
  contour *newCont = NULL;

  if(slice1 == NULL || slice2 == NULL) return;

  /* allocate */
  *newSlice1 = newList(LIST);
  *newSlice2 = newList(LIST);

/* clone slice 2 first, so we can maintain adjacency */
  for(i = getListNode(slice2,0); i; i = (listNode*) i->next) {
    newCont = cloneContour((contour*) i->data);
    enqueue(*newSlice2,newCont);
  }

  /* clone slice 1 */
  for(i = getListNode(slice1,0); i; i = (listNode*) i->next) {
    newCont = cloneContour((contour*) i->data);
    enqueue(*newSlice1,newCont);

    /* fix adjacency */
    for(j = getListNode(newCont->adjacentContours,0); j; j = (listNode*) j->next) {
      ind = findInListI(slice2,j->data);
      j->data = getListNode(*newSlice2,ind)->data;
    }
  }
}

/**
 * delete empty slices
 */
void deleteEmptySlices(list *slices) {
  listNode *slice, *cont;
  list *verts;
  int found;

  /* delete empty slices  */
  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    found = 0;
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      verts = ((contour*)cont->data)->vertices;
      if(listSize(verts) > 0) {
	found = 1;
	break;
      }
    }

    // delete all contours if empty
    if(!found) {
      for(cont = getListNode((list*)slice->data,0); cont;
	  cont = (listNode*) cont->next) {
	deleteContour((contour*) cont->data);
      }
      freeList((list*)slice->data);
      markForDeletion(slice);
    }

  }

  deleteMarkedNodes(slices);
}

/**
 * gets the average interslice distance
 */
double getAverageIntersliceDistance(list *slices) {
  listNode *cur, *last;
  double curDist,sum = 0.0;
  int n = 0;

  for(last = getListNode(slices,0), cur = getListNode(slices,1); last && cur;
      last = cur, cur = (listNode*) cur->next) {
    curDist = getIntersliceDistance((list*) cur->data, (list*) last->data);
    if(fabs(curDist) > SR_TOL) {
      n++;
      sum += curDist;
    }
  }

  return sum/n;
}

/**
 * gets interslice distance as z(curSlice) - z(lastSlice)
 */
double getIntersliceDistance(list *curSlice, list *lastSlice) {
  listNode *i;
  contour *lastFirst = NULL, *curFirst = NULL;

  if(curSlice == NULL || listSize(curSlice) < 1
     || lastSlice == NULL || listSize(lastSlice) < 1) {
    return 0.0;
  }

  for(i = getListNode(lastSlice,0); i; i = (listNode*) i->next) {
    lastFirst = (contour*) i->data;

    if(lastFirst != NULL && listSize(lastFirst->vertices) > 0) {
      break;
    }
  }

  for(i = getListNode(curSlice,0); i; i = (listNode*) i->next) {
    curFirst = (contour*) i->data;

    if(curFirst != NULL && listSize(curFirst->vertices) > 0) {
      break;
    }
  }

  if((lastFirst == NULL || listSize(lastFirst->vertices) < 1)
     || (curFirst == NULL || listSize(curFirst->vertices) < 1)) {
    return 0.0;
  }

  return fabs(((vertex*) getListNode(curFirst->vertices,0)->data)->z
	      - ((vertex*) getListNode(lastFirst->vertices,0)->data)->z);
}

/**
 * preprocesses slice contours before surface reconstruction
 *
 * this function deletes repreated vertices and resamples contours
 */
int preprocessSliceContours(list *slices) {
  /* validate */
  if(slices == NULL || listSize(slices) < 1) {
    return SR_FAILURE;
  }

  if((minSliceInd > 0 || maxSliceInd < listSize(slices)-1) && skipSlices > 1) {
    fprintf(stderr,"WARNING!!!: skipping slices (%d) with minSliceInd (%d) and/or maxSliceInd (%d) overridden!!\n",
	    skipSlices,minSliceInd,maxSliceInd);

    /* fix minSlice and maxSlice  for the new skip */
    //minSliceInd = floor(minSliceInd/skipSlices);
    //maxSliceInd = floor(maxSliceInd/skipSlices);
  }

  /* skip slices */
  takeEachNthElement(slices,skipSlices);

  if(SR_VERBOSE) {
    fprintf(stdout,"deleting empty slices\n");
  }
  deleteEmptySlices(slices);

  if(SR_VERBOSE) {
    fprintf(stdout,"deleting repeated vertices\n");
  }

  deleteRepeatedContourVerts(slices);

  if(resample) {
    if(resampleEqual) {
      lastResampleDist = getAverageIntersliceDistance(slices);
    }
    else if(!resampleAbs) {
      lastResampleDist = resampleDistance*getAverageIntersliceDistance(slices);
    }
    else {
      lastResampleDist = resampleDistance;
    }

    if(SR_VERBOSE) {
      fprintf(stdout,"resampling contours with vertex distance %lf\n",
	      lastResampleDist);
    }

    resampleContoursByDistance(slices,lastResampleDist);

    //resampleContoursByAngle(slices,resampleAngle);
  }

  if(SR_VERBOSE) {
    fprintf(stdout,"orienting contours\n");
  }

  orientContours(slices);

  //randomizeClosedContourIndexing(slices);

  numberVertices(slices);

  labelBoundaries(slices);

  return SR_SUCCESS;
}

/**
 * straight up copy the forward to backward adjacency list
 */
void copyForwardToBackwardAdjacency(list *slice) {
  listNode *i, *j;
  contour *curCont, *otherSliceCont;

  /* validate */
  if(slice == NULL || listSize(slice) < 1) {
    return;
  }

  /* initialize slice2 backward connectivity list */
  for(i = getListNode(slice,0); i; i = (listNode*) i->next) {
    curCont = (contour*) i->data;
    freeList(curCont->adjacentBackwardContours);
    curCont->adjacentBackwardContours = newList(LIST);
  }

  /* iterate over adjacent contours in slice, adding to slice back
     adjacent */
  for(i = getListNode(slice,0); i; i = (listNode*) i->next) {
    curCont = (contour*) i->data;

    /* insert contour as an adjacent one to each contour in list */
    for(j = getListNode(curCont->adjacentContours,0); j;
	j = (listNode*) j->next) {
      otherSliceCont = (contour*) j->data;
      enqueue(curCont->adjacentBackwardContours,otherSliceCont);
    }
  }
}

/**
 * builds a backward adjacency list between two slices
 */
void buildBackwardAdjacency(list *slice1, list *slice2) {
  listNode *i, *j;
  contour *curCont, *otherSliceCont;

  /* validate */
  if(slice1 == NULL || slice2 == NULL || listSize(slice1) < 1 || listSize(slice2) < 1) {
    return;
  }

  /* initialize slice2 backward connectivity list */
  for(i = getListNode(slice2,0); i; i = (listNode*) i->next) {
    curCont = (contour*) i->data;
    freeList(curCont->adjacentBackwardContours);
    curCont->adjacentBackwardContours = newList(LIST);
  }

  /* iterate over contours in slice1, adding to slice 2 contour back
     adjacent */
  for(i = getListNode(slice1,0); i; i = (listNode*) i->next) {
    curCont = (contour*) i->data;

    /* insert this contour as an adjacent one to each contour in list */
    for(j = getListNode(curCont->adjacentContours,0); j;
	j = (listNode*) j->next) {
      otherSliceCont = (contour*) j->data;
      enqueue(otherSliceCont->adjacentBackwardContours,curCont);
    }
  }
}

/**
 * builds a backward adjacency list between two slices represented by
 * contour pairs. returns a list of contours representing the slice 2 contours
 */
list *buildBackwardAdjacencyFromContourPairs(list *contourPairs) {
  list *slice2;
  listNode *i;
  contourPair *cp;

  slice2 = newList(LIST);

  /* go through all the c2 contours, add them as contours in the slice */
  for(i = getListNode(contourPairs,0); i; i = (listNode*) i->next) {
    cp = (contourPair*) i->data;

    /* see if the contour already exists */
    if(NULL == findInListLN(slice2,cp->c2)) { /* add it */
      enqueue(slice2,cp->c2);

      if(cp->c2->adjacentBackwardContours) {
	freeList(cp->c2->adjacentBackwardContours);
      }

      cp->c2->adjacentBackwardContours = newList(LIST);
    }

    /* add the correspondence */
    enqueue(cp->c2->adjacentBackwardContours,cp->c1);
  }

  return slice2;
}

/**
 * deletes a contour pair, and the contours themselves
 */
void deleteContourPairList(list *contourPairs) {
  contourPair *cp;
  listNode *ln;

  if(contourPairs == NULL) return;

  for(ln = getListNode(contourPairs,0); ln; ln = (listNode*) ln->next) {
    cp = (contourPair*) ln->data;
    deleteContour(cp->c1);
    deleteContour(cp->c2);
  }
  freeListAndData(contourPairs);
}

/**
 * Tests to see if a slice has no vertices
 */
int sliceEmpty(list *slice) {
  listNode *ln;
  contour *cont;

  if(slice == NULL || listSize(slice) < 1) return TRUE;

  for(ln = getListNode(slice,0); ln; ln = (listNode*) ln->next) {
    cont = (contour*) ln->data;
    if(cont != NULL && listSize(cont->vertices) > 0) {
      return FALSE;
    }
  }

  return TRUE;
}

/**
 * select all the vertices between two selected vertices
 */
int selectVerticesBetween(list *vertices) {
  listNode *ln, *first = NULL, *last = NULL;
  vertex *cur;

  /*validate */
  if(vertices == NULL || listSize(vertices) == 0) {
    return SR_FAILURE;
  }

  /* loop over vertices, finding the first and last */
  for(ln = getListNode(vertices,0);
      ln && last == NULL; ln = (listNode*) ln->next) {
    cur = (vertex*) ln->data;
    if(cur->selected) {
      if(first == NULL) {
	first = ln;
	continue;
      }
      else if(last == NULL) {
	last = ln;
	continue;
      }
    }
  }

  /* test for fitst and last */
  if(first == NULL || last == NULL) {
    return SR_FAILURE;
  }

  /* fill with selected */
  for(ln = first; ln && ln != last; ln = (listNode*) ln->next) {
    ((vertex*) ln->data)->selected = TRUE;
  }

  if(ln != last) {
    return SR_FAILURE;
  }

  return SR_SUCCESS;
}

/** some generic geometry utils */

/**
 * create a vertex structure with default values
 * returns a vertex structure or null if creation fails
 */
vertex *createVertex() {
  vertex *v = (vertex*) malloc(sizeof(vertex));

  if(v == NULL) return NULL;

  v->x = 0.0f;
  v->y = 0.0f;
  v->z = 0.0f;

  v->selected = FALSE;

  v->number = -1;

  v->boundary = -1;

  v->label = -1;

  return v;
}

/**
 * copies a vertex
 * returns a copy of the passed vertex, or null, if cant allocate
 */
vertex *copyVertex(vertex *v) {
  /* allocate vertex */
  vertex *newV = createVertex();
  if(newV == NULL) {
    return NULL;
  }

  /* copy individual fields */
  newV->x = v->x;
  newV->y = v->y;
  newV->z = v->z;

  newV->selected = v->selected;
  newV->number = v->number;
  newV->boundary = v->boundary;
  newV->label = v->label;

  //printf("%d\n",newV->label);

  return newV;
}

/**
 * gets a new vertex that is the midpoint between two vertices
 */
vertex *getMidpoint(vertex *v1, vertex *v2) {
  if(v1 == NULL || v2 == NULL) return NULL;

  vertex *m = copyVertex(v1);
  m->x+=0.5*(v2->x-v1->x);
  m->y+=0.5*(v2->y-v1->y);
  m->z+=0.5*(v2->z-v1->z);

  return m;
}

/**
 * test for intersection of a contour and an edge
 */
int edgeContourIntersect(edge *e, contour *c) {
  edge contEdge;
  listNode *ln;

  /* iterate over contour edges, testing each */
  contEdge.v1 = (vertex*) getListNode(c->vertices,0)->data;
  for(ln = getListNode(c->vertices,1); ln; ln = (listNode*) ln->next) {
    contEdge.v2 = (vertex*) ln->data;

    /* test for intersection */
    if(edgesIntersect(&contEdge,e)) {
      return TRUE;
    }

    /* update the edge vertices */
    contEdge.v1 = contEdge.v2;
  }

  return FALSE;

}

/**
 * resample one contour by integrated distance
 */
void resampleContourByDistance(contour *cont, double distance) {
  if(distance <= 0 || cont == NULL || cont->vertices == NULL || listSize(cont->vertices) < 2) {
    return;
  }

  double curd = 0.0, d;
  list *verts = cont->vertices;
  list *newVerts = newList(LIST);
  listNode *v;
  vertex *vert,*prev,*newVert;

  newVert = prev = (vertex*) getListNode(verts,0)->data;
  enqueue(newVerts,copyVertex(prev));

  /* integrate distance travelled, delete middle vertices with insufficient
     distance */
  for(v = getListNode(verts,1); v; v = (listNode*) v->next) {
    vert = (vertex*) v->data;
    d=dist((*vert), (*prev));
    curd+=d;
    while(curd >= distance) {
      newVert = copyVertex(vert);
      newVert->x = prev->x + ((distance-(curd-d)) / d) * (vert->x-prev->x);
      newVert->y = prev->y + ((distance-(curd-d)) / d) * (vert->y-prev->y);
      enqueue(newVerts,newVert);
      curd-=distance;
    }
    prev = vert;
  }

  if(cont->closed == OPEN) {
    /* handle the last vertex */
    if(curd > distance/2) {
      newVert = copyVertex(prev);
      enqueue(newVerts,newVert);
    }
    else {
      newVert->x = prev->x;
      newVert->y = prev->y;
    }
  }
  else {
    vert = (vertex*) getListNode(verts,0)->data;
    d=dist((*vert), (*prev));
    curd+=d;
    while(curd >= distance) {
      newVert = copyVertex(vert);
      newVert->x = prev->x + ((distance-(curd-d)) / d) * (vert->x-prev->x);
      newVert->y = prev->y + ((distance-(curd-d)) / d) * (vert->y-prev->y);
      enqueue(newVerts,newVert);
      curd-=distance;
    }
  }

  cont->vertices = newVerts;
  freeListAndData(verts);
}

/**
 * resample slice contours by integrated angle
 * ADD LABEL PRESERVATION
 */
void resampleContoursByAngle(list *slices, double angle) {
  list *verts, *newVerts;
  listNode *s,*c,*v;
  contour *cont;
  int i,j;
  double cura;
  vertex *vert,*prev,*prevprev;

  for(s = getListNode(slices,0),j=0; s; s = (listNode*) s->next,j++) {
    for(c = getListNode((list*)s->data,0); c; c = (listNode*) c->next) {
      cont = (contour*) c->data;

      cura = 0.0;
      verts = cont->vertices;
      newVerts = newList(LIST);

      prevprev = (vertex*) getListNode(verts,0)->data;
      prev = (vertex*) getListNode(verts,1)->data;
      enqueue(newVerts,copyVertex(prevprev));

      /* integrate angle swept, delete middle vertices w/ insufficient angle */
      for(v = getListNode(verts,2), i = 2; v; v = (listNode*) v->next, i++) {
	vert = (vertex*) v->data;
	cura+=fabs(SR_PI-angleBetweenSegments(prevprev,prev,vert));
	if(cura >= angle) {
	  enqueue(newVerts,copyVertex(vert));
	  prevprev = prev;
	  prev = vert;
	  cura = 0.0;
	}
      }

      cont->vertices = newVerts;
      freeListAndData(verts);
    }
  }
}

/**
 * resample slice contours by integrated distance
 */
void resampleContoursByDistance(list *slices, double distance) {
  listNode *s,*c;
  contour *cont;

  for(s = getListNode(slices,0); s; s = (listNode*) s->next) {
    for(c = getListNode((list*)s->data,0); c; c = (listNode*) c->next) {
      cont = (contour*) c->data;
      resampleContourByDistance(cont,distance);
    }
  }
}

/**
 * delete repeated vertices in a contour
 */
void deleteRepeatedContourVerts(list *slices) {
  listNode *slice, *cont;
  list *verts;
  listNode *vert;
  vertex *prev = NULL;
  int i;

  /* delete repeated vertices  */
  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      prev = NULL;
      verts = ((contour*)cont->data)->vertices;
      if(listSize(verts) < 2) continue;

      for(vert = getListNode(verts,0), i = 0; vert;
	  vert = (listNode*) vert->next, i++) {

	if(verticesEqual((vertex*)vert->data,prev)) {
	  markForDeletion(vert);
	}
	else {
	  prev = (vertex*)vert->data;
	}
      }

      /* test last and first */
      if(verticesEqual((vertex*)getListNode(verts,0)->data,prev)) {
	markForDeletion(getListNode(verts,0));
      }

      deleteMarkedNodes(verts);
    }
  }

}

/**
 * orient contours
 */
void orientContours(list *slices) {
  listNode *slice, *cont;
  list *verts;
  double a;

  /* orient contours by flipping contours with negative area */
  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      verts = ((contour*)cont->data)->vertices;
      a = contourArea((contour*)cont->data);
      if(a < 0) {
	reverseList(verts);
      }
    }
  }
}

/**
 * randomize contour indexing for closed contours
 */
void randomizeClosedContourIndexing(list *slices) {
  listNode *slice, *cont;
  contour *c;
  list *verts, *split;
  int ind;
  srand(time(NULL));

  /* randomize indexing by spliting and joining */
  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      c = (contour*) cont->data;

      if(c == NULL || c->vertices == NULL || listSize(c->vertices) < 2
	 || c->closed != CLOSED) {
	continue;
      }

      verts = c->vertices;
      ind = rand()/(LONG_MAX/listSize(verts)-1);
      split = splitList(verts,ind);
      appendList(split,verts);
      c->vertices = split;
    }
  }
}

/**
 * all vertices in the slices based on contour position
 */
void numberVertices(list *slices) {
  listNode *slice, *cont, *vert;
  list *verts;
  int cur;

  /* number vertices sequentially in each contour */
  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      verts = ((contour*)cont->data)->vertices;
      cur = 0;
      for(vert = getListNode(verts,0); vert; vert = (listNode*) vert->next) {
	((vertex*)vert->data)->number = cur++;
      }

    }
  }
}

/**
 * label all vertices based on whether they are on a boundary or not
 */
void labelBoundaries(list *slices) {
  listNode *slice, *cont, *vert;
  list *verts;

  for(slice = getListNode(slices,0); slice; slice = (listNode*) slice->next) {
    for(cont = getListNode((list*)slice->data,0); cont;
	cont = (listNode*) cont->next) {
      verts = ((contour*)cont->data)->vertices;

      if(listSize(verts) < 1) continue;

//      if(((contour*)cont->data)->closed == CLOSED && (!first && !last) &&

      ((vertex*) getListNode(verts,0)->data)->boundary
	= ((contour*)cont->data)->closed == CLOSED ? INTERIOR : BOUNDARY;

      for(vert = getListNode(verts,1); vert && vert->next;
	  vert = (listNode*) vert->next) {
	((vertex*)vert->data)->boundary = INTERIOR;
      }

      ((vertex*) getListNode(verts,listSize(verts)-1)->data)->boundary
	= ((contour*)cont->data)->closed == CLOSED ? INTERIOR : BOUNDARY;
    }
  }
}

/**
 * get the angle between two line segments
 */
double angleBetweenSegments(vertex *a, vertex *b, vertex *c) {
  vertex p;
  double theta;

  p.x = (c->x - b->x) * (a->x - b->x) + (c->y - b->y) * (a->y - b->y);
  p.y = (c->x - b->x) * (a->y - b->y) - (c->y - b->y) * (a->x - b->x);

  if(p.x == 0.0 && p.y == 0.0) {
    return 0.0;
  }

  theta = atan2(p.y, p.x);

  while(theta < 0.0) {
    theta = theta + 2.0 * SR_PI;
  }

  while(theta > SR_PI) {
    theta = 2.0 * SR_PI - theta;
  }

  return theta;
}

/* get the centroid of the vertices of a contour */
void getContourCentroid(contour* c, vertex *cent) {
  listNode *ln;

  /*initialize the centroid coordinates */
  cent->x = cent->y = cent->z = 0;

  /* sum the coordinates */
  for(ln = getListNode(c->vertices,0); ln; ln = (listNode*) ln->next) {
    cent->x += ((vertex*)ln->data)->x;
    cent->y += ((vertex*)ln->data)->y;
    cent->z += ((vertex*)ln->data)->z;
  }

  /* divide by the number of vertices */
  cent->x/=listSize(c->vertices);
  cent->y/=listSize(c->vertices);
  cent->z/=listSize(c->vertices);
}

/**
 * test a skeleton vertex for containment in the projection of a contour
 */
int inContourSkelVert(contour *c, vertex *v) {
  int cn = 0;
  listNode *i;
  vertex *v1, *v2;
  double vt;

  // loop through all edges
  for(i = getListNode(c->vertices,0); i; i = (listNode*) i->next) {
    v1 = (vertex*) i->data;
    if(i->next) {
      v2 = (vertex*) ((listNode*)i->next)->data;
    }
    else {
      v2 = (vertex*) getListNode(c->vertices,0)->data;
    }

    if(((v1->y <= v->y) && (v2->y > v->y))
       || ((v1->y > v->y) && (v2->y <= v->y))) {
      vt = (v->y - v1->y) / (v2->y - v1->y);
      if (v->x < v1->x + vt * (v2->x - v1->x)) {
	cn++;
      }
    }
  }
  return (cn&1);

}

/**
 * get the skeleton of a closed contour
 * input is a closed contour
 * output is a list of vertices on the skeleton and a list of edges
 * connecting the vertices, both are null if skeleton cannot be found
 * returns SR_SUCCESS or SR_FAILURE
 */
int getContourSkeleton(contour* c, list **vertices, list **edges) {
  char *trifnameBase = "/tmp/skel_tri";
  char trifnamePoly[SR_MAX_STR_LEN];
  char trifnameVor[SR_MAX_STR_LEN];
  listNode *ln, *ln2;
  vertex *v, *replaceV;
  edge *e;
  int lab = -1;

  /* validate input */
  if(c == NULL || c->vertices == NULL || listSize(c->vertices) == 0
     || c->closed == OPEN) {
    *vertices = *edges = NULL;
    return SR_FAILURE;
  }

  /* iterate over vertices, testing if they all belong to one label */
  lab = ((vertex*) getListNode(c->vertices,0)->data)->label;
  for(ln = getListNode(c->vertices,0); ln; ln = (listNode*) ln->next) {
    v = (vertex*) ln->data;
    if(v->label != lab) {
      lab = -1;
      break;
    }
  }

  /* if there are 3 or less vertices in the contour, just make skel==centroid*/
  if(listSize(c->vertices) <= 3) {
    v = createVertex();
    v->label = lab;

    getContourCentroid(c,v);
    *vertices = newList(LIST);
    enqueue(*vertices,v);
    *edges = newList(LIST);
    return SR_SUCCESS;
  }

  /** compute the voronoi diagram */

  /* build filenames */
  sprintf(trifnamePoly,"%s.poly",trifnameBase);
  sprintf(trifnameVor,"%s.1.v",trifnameBase);

  /* write a triangle poly file of the contour */
  writeTrianglePolyFile(c, trifnamePoly);

  /* execute triangle */
  if(0 != executeTriangleWithFlags(trifnamePoly,"-BPNEgQv -s")) {
    fprintf(stderr,"error invoking the triangle program. make sure it is available and in the path.\nto get Triangle, visit http://www.cs.cmu.edu/~quake/triangle.html");
    return SR_FAILURE;
  }

  /* read the result */
  readTriangleNodeAndEdgeFile(trifnameVor,vertices,edges);

  /* go through the edges to remove the external edges */
  for(ln = getListNode(*edges,0); ln; ln = (listNode*) ln->next) {
    e = (edge*) ln->data;
    if(e->v1 == NULL || e->v2 == NULL) {
      //       || !inContourSkelVert(c,e->v1) || !inContourSkelVert(c,e->v2)) {
      markForDeletion(ln);
    }
  }
  deleteMarkedNodes(*edges);

  /* find and remove vertices outside the contour */
  for(ln = getListNode(*vertices,0); ln; ln = (listNode*) ln->next) {
    v = (vertex*) ln->data;
    if(!inContourSkelVert(c,v)) {
      // find instances in edges and remove them
      replaceV = NULL;
      for(ln2 = getListNode(*edges,0); ln2; ln2 = (listNode*) ln2->next) {
	e = (edge*) ln2->data;
	if(e->v1 == v) {
	  if(replaceV == NULL) {
	    replaceV = e->v2;
	    markForDeletion(ln2);
	  }
	  else {
	    e->v1 = replaceV;
	  }
	}
	else if(e->v2 == v) {
	  if(replaceV == NULL) {
	    replaceV = e->v1;
	    markForDeletion(ln2);
	  }
	  else {
	    e->v2 = replaceV;
	  }
	}
      }
    }
  }
  deleteMarkedNodes(*edges);

  /* assign the z coords */
  for(ln = getListNode(*vertices,0); ln; ln = (listNode*) ln->next) {
    v = (vertex*) ln->data;
    v->label = lab;
    v->z = ((vertex*)getListNode(c->vertices,0)->data)->z;
  }

  return SR_SUCCESS;
}

/**
 * write a Triangle .poly file from a contour
 * NOTE: this function assumes the contour vertices are connected in order
 */
void writeTrianglePolyFile(contour *c, char *filename) {
  FILE *fp;
  list *vertices = c->vertices;
  listNode *ln;
  vertex *vert;
  int i;

  /* open the file */
  fp = fopen(filename,"w");
  if(fp == NULL) {
    fprintf(stderr,"error: failed in opening file %s for writing poly file.\n",
	    filename);
    return;
  }

  /* write the header */
  fprintf(fp, "%d 2 0 0\n", listSize(vertices));

  /* write the vertices */
  for(ln = getListNode(vertices,0), i = 0; ln; ln = (listNode*) ln->next, i++) {
    vert = (vertex*) ln->data;
    fprintf(fp,"%d %0.30lf %0.30lf\n", i+1, vert->x, vert->y);
  }

  /* write the edge header */
  fprintf(fp,"%d 0\n",
	  (c->closed == CLOSED) ? listSize(vertices) : listSize(vertices) - 1);

  /* write the edges */
  for(i = 0; i < listSize(vertices) - 1; i++) {
    fprintf(fp,"%d %d %d\n", i+1, i+1, i+2);
  }

  /* close, if a closed contour */
  if(c->closed) {
    fprintf(fp,"%d %d %d\n", i+1, i+1, 1);
  }

  /* print no holes */
  fprintf(fp,"0\n");

  fclose(fp);
}

/**
 * reads a node file output by Triangle into a list
 */
list *readTriangleNodeFile(char *filename) {
  list *nodes;
  vertex *v;
  int node, numNodes, trash;
  FILE *fp;

  /* open and validate */
  if(filename == NULL || !(fp = fopen(filename,"r"))) {
    return NULL;
  }

  nodes = newList(LIST);
  if(nodes == NULL) {
    return NULL;
  }

  /* read the number of nodes */
  skipComments(fp,'#');
  fscanf(fp,"%d %d 0 0",&numNodes,&trash);

  /* read each node */
  for(node = 0; node < numNodes; node++) {
    v = createVertex();
    fscanf(fp,"%d %lf %lf",&trash,&v->x,&v->y);
    v->z = 0;
    enqueue(nodes,v);
  }

  return nodes;
}

/**
 * reads an edge file output by Triangle into a list
 */
list *readTriangleEdgeFile(char *filename) {
  list *edges;
  int *e, edge, numEdges, trash;
  double dtrash;
  FILE *fp;

  /* open and validate */
  if(filename == NULL || !(fp = fopen(filename,"r"))) {
    return NULL;
  }

  edges = newList(LIST);
  if(edges == NULL) {
    return NULL;
  }

  /* read the number of edges */
  skipComments(fp,'#');
  fscanf(fp,"%d %d",&numEdges,&trash);

  /* read each edge */
  for(edge = 0; edge < numEdges; edge++) {
    e = (int*) malloc(2*sizeof(int));
    fscanf(fp,"%d %d %d",&trash,e,e+1);
    e[0]--;
    e[1]--;

    /* if this was an edge with a divergent endpoint, read the normal and
       trash it */
    if(e[1] == -2) {
      fscanf(fp,"%lf %lf",&dtrash,&dtrash);
    }

    enqueue(edges,e);
  }

  return edges;
}

/**
 * reads a node and corresponding edge file output by Triangle into a lists,
 * where the edge list has vertex pointers
 */
int readTriangleNodeAndEdgeFile(char *filebase, list **vertices, list**edges) {
  edge *e;
  int edg, numEdges, trash, v1, v2;
  double dtrash;
  FILE *fp;
  char nodefilename[SR_MAX_STR_LEN];
  char edgefilename[SR_MAX_STR_LEN];
  listNode *ln;

  /* validate */
  if(filebase == NULL) {
    return SR_FAILURE;
  }

  /* read the nodes as normal */
  sprintf(nodefilename,"%s.node",filebase);
  *vertices = readTriangleNodeFile(nodefilename);
  if(*vertices == NULL) {
    return SR_FAILURE;
  }

  /* open the edge file */
  sprintf(edgefilename,"%s.edge",filebase);
  if(!(fp = fopen(edgefilename,"r"))) {
    freeListAndData(*vertices);
    return SR_FAILURE;
  }

  /* read the edge file */
  *edges = newList(LIST);

  if(*edges == NULL) {
    return SR_FAILURE;
  }

  /* read the number of edges */
  skipComments(fp,'#');
  fscanf(fp,"%d %d",&numEdges,&trash);

  /* read each node */
  for(edg = 0; edg < numEdges; edg++) {
    e = (edge*) malloc(sizeof(edge));

    fscanf(fp,"%d %d %d",&trash,&v1,&v2);

    /* if this was an edge with a divergent endpoint, read the normal and
       trash it */
    if(v2 == -1) {
      fscanf(fp,"%lf %lf",&dtrash,&dtrash);
    }

    /* get the vertex pointers */
    ln = getListNode(*vertices,v1-1);
    e->v1 = ln == NULL ? NULL : (vertex*) ln->data;

    ln = getListNode(*vertices,v2-1);
    e->v2 = ln == NULL ? NULL : (vertex*) ln->data;

    enqueue(*edges,e);
  }

  return SR_SUCCESS;
}

/**
 * runs Triangle on an existing .poly file, producing an off file
 * NOTE: uses the -qBPNEg flags when running Triangle
 * returns the value returned by system()
 */
int executeTriangle(char *polyInputFilename) {
  char *triangleFlags = "-BPNEgQ -q 30 -s";
  char execStr[100] = "";

  /* build a string fopr the command to be executed */
  sprintf(execStr,"%s %s %s",TRIANGLE_EXEC,triangleFlags,polyInputFilename);

  /* execute Triangle */
  return system(execStr);
}

/**
 * runs Triangle on an existing .poly file with specified flags
 * returns the value returned by system()
 */
int executeTriangleWithFlags(char *polyInputFilename, char *flags) {
  char execStr[100] = "";

  /* build a string for the command to be executed */
  sprintf(execStr,"%s %s %s",TRIANGLE_EXEC,flags,polyInputFilename);

  /* execute Triangle */
  return system(execStr);
}

/* phase unwrapping, takes (-pi, pi] to [0, 2pi)*/
inline double unwrap(double a) {
  if(a < 0) return 2*SR_PI+a;
  return a;
}

/**
 * wrapper for matrix multiply of vertex coords
 */
void matrixMult4byVert(float a[4][4], vertex b, vertex *c) {
  float arr[4], res[4];
  arr[0] = b.x;
  arr[1] = b.y;
  arr[2] = b.z;
  arr[3] = 1.0f;

  matrixMult4by1(a,arr,res);

  c->x = res[0];
  c->y = res[1];
  c->z = res[2];
}

/**
 * perform a linear fit given a vector of observations and loadings via gsl
 */
double *linearFit(double *x, double *y, double *w, int n) {
  double *b = malloc(2*sizeof(double));
  double cov[3];
  double chi2;

  gsl_fit_wlinear(x,1,w,1,y,1,n,&b[0],&b[1],&cov[0],&cov[1],&cov[2],&chi2);

  return b;
}

/* debugging funcs */
void printPath(FILE* fp, path *p) {
  listNode *ln;
  face *f;

  if(p == NULL || p->faceList == NULL || listSize(p->faceList) == 0) {
    return;
  }

  fprintf(fp,"printing path with %d faces\n", listSize(p->faceList));
  for(ln = getListNode(p->faceList,0); ln; ln = (listNode*) ln->next) {
    f = (face*) ln->data;
    fprintf(fp,"face: %d %d %d\n",f->v1,f->v2,f->v3);
  }
  fprintf(fp,"done\n\n");
}

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/libsrUtil.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
