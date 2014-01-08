/*****************************************************************************
 * branch.c is the source file for an implementation of a solution to the
 * branching problem for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#include"branch.h"

int cliqueRadius = 4; /* number of neighbors on each side that get votes */

/**
 * take two sets of contours for adjacent slices that already have
 * correspondences identified. For each contour on slice one that is
 * connected to multiple contours on slice two, break it up. returns a list
 * of contour pairs to be tiled.
 *
 * 1) for each vertex on prebranch contour, assign closest postbranch
 * contour
 *
 * 2) use a maximum likelihood estimator to find boundaries between
 * branches.  basically, at any vertex of the prebranch contour, find the
 * maximum number of votes for any postbranch contour among vertices within
 * a clique.
 *
 * 3) split the prebranch contours into new contour pairs based on regions
 *    of constant assignment.
 */
list *branchSlices(list *slice1, list *slice2) {
		   //list *unmatchedSlice1, list *unmatchedList2) {
  list *newSlice,
    *newSlice2,
    *forwardContourPairs, 
    *backContourPairs, 
    *tmpContourPairs,
    *contourPairs;
  listNode *i,*j;
  contour *curCont;
  contour **closestContours;
  contourPair *cp, *newcp;
  int multi = 1;

  // DEBUGGING
  //int ind;
  //int j;
  //edge e;
  //fprintf(stdout,"branching slices with %d and %d contours\n",listSize(slice1),listSize(slice2));


  /* validate */
  if(slice1 == NULL || slice2 == NULL || listSize(slice1) < 1 || listSize(slice2) < 1) {
    return NULL;
  }

  // make copies of everything to avoid messing with original data
  cloneSlices(slice1,slice2,&newSlice,&newSlice2);

  /* label the origins of all contours in old and new slices */
  for(i = getListNode(newSlice,0), j = getListNode(slice1,0); i; 
      i = (listNode*) i->next, j = (listNode*) j->next) {
    ((contour*) i->data)->origin = (struct contour*) i->data;
    ((contour*) j->data)->origin = (struct contour*) i->data;
  }
  for(i = getListNode(newSlice2,0), j = getListNode(slice2,0); i; 
      i = (listNode*) i->next, j = (listNode*) j->next) {
    ((contour*) i->data)->origin = (struct contour*) i->data;
    ((contour*) j->data)->origin = (struct contour*) i->data;
  }

  slice1 = newSlice;
  slice2 = newSlice2;

  // copy the forward to backward contours (confusing, but it must me done
  // so that the loop can just look at backward contours)
  copyForwardToBackwardAdjacency(newSlice);
  

  // this loop iteratively splits contours while there are still multiply
  // connected contours. each time the contours are iterated over they are
  // tested for multiple connectedness
  while(multi) {
    multi = 0;

    forwardContourPairs = newList(LIST);

    /* find prebranch contours in slice1 */
    for(i = getListNode(newSlice,0); i; i = (listNode*) i->next) {
      curCont = (contour*) i->data;

      if(listSize(curCont->adjacentBackwardContours) == 1) {
	/* just add this as a one to one correspondence */
	cp = (contourPair*) malloc(sizeof(contourPair));

	cp->c1 = cp->c1Origin = curCont;
	cp->c2 = cp->c2Origin = 
	  (contour*) getListNode(curCont->adjacentBackwardContours,0)->data;

	enqueue(forwardContourPairs,cp);
      }
      else if(listSize(curCont->adjacentBackwardContours) < 1) {
	continue;
      }
      else {
	multi = 1;

	/* found one, split it up */
	splitContoursMLE(curCont,curCont->adjacentBackwardContours,
			 &closestContours);

	// DEBUGGING
	//      for(ind = 0; ind < listSize(curCont->vertices); ind++) {
	//      	fprintf(stdout,"%d %p\n", ind, closestContours[ind]);
	//      }

	/* add contour pairs between branches */
	tmpContourPairs = newList(LIST);
	addContourPairsFromClosest(curCont,closestContours,tmpContourPairs);
	appendList(forwardContourPairs,tmpContourPairs);
	free(tmpContourPairs);

	if((contour*) curCont->origin != curCont) {
	  free(curCont);
	}
	free(closestContours);
      }
    }

    /* print for debugging */
    if(SR_DEBUG) {
      fprintf(stderr,"------------------------\n");
      for(i = getListNode(forwardContourPairs,0); i; i = (listNode*) i->next) {
	dumpContourPair(stdout,(contourPair*) i->data);
      }
      fprintf(stderr,"------------------------\n");
    }

    /* find prebranch contours in slice2 */
    newSlice = buildBackwardAdjacencyFromContourPairs(forwardContourPairs);
    freeListAndData(forwardContourPairs);

    backContourPairs = newList(LIST);
    for(i = getListNode(newSlice,0); i; i = (listNode*) i->next) {
      curCont = (contour*) i->data;

      if(listSize(curCont->adjacentBackwardContours) == 1) {
	/* just add this as a one to one correspondence */
	cp = (contourPair*) malloc(sizeof(contourPair));
	cp->c1 = cp->c1Origin =curCont;
	cp->c2 = cp->c2Origin =
	  (contour*) getListNode(curCont->adjacentBackwardContours,0)->data;
	enqueue(backContourPairs,cp);
      }
      else if(listSize(curCont->adjacentBackwardContours) < 1) {
	continue;
      }
      else {
	multi = 1;

	/* found one, split it up */

	splitContoursMLE(curCont,curCont->adjacentBackwardContours,
			 &closestContours);

	// DEBUGGING
	//      for(ind = 0; ind < listSize(curCont->vertices); ind++) {
	//      	fprintf(stdout,"%d %p\n", ind, closestContours[ind]);
	//      }
	//for(j = 0; j < listSize(curCont->vertices); j++) {
	//	fprintf(stdout,"%0.2lf %0.2lf %d %p\n", 
	//		((vertex*) getListNode(curCont->vertices,j)->data)->x, 
	//		((vertex*) getListNode(curCont->vertices,j)->data)->y, 
	//		j, closestContours[j]);
	//}

	/* add contour pairs between branches */
	tmpContourPairs = newList(LIST);
	addContourPairsFromClosest(curCont,closestContours,tmpContourPairs);
	appendList(backContourPairs,tmpContourPairs);
	free(tmpContourPairs);

	if((contour*) curCont->origin != curCont) {
	  free(curCont);
	}
	free(closestContours);
      }
    }
    freeList(newSlice);

    /* print for debugging */
    if(SR_DEBUG) {
      fprintf(stderr,"------------------------\n");
      for(i = getListNode(backContourPairs,0); i; i = (listNode*) i->next) {
	dumpContourPair(stdout,(contourPair*) i->data);
      }
      fprintf(stderr,"------------------------\n");
    }

    if(multi) {
      newSlice = buildBackwardAdjacencyFromContourPairs(backContourPairs);
      freeListAndData(backContourPairs);
    }
  }

//  forwardContourPairs = newList(LIST);
//  /* find prebranch contours in newSlice again */
//  for(i = getListNode(newSlice,0); i; i = (listNode*) i->next) {
//    curCont = (contour*) i->data;
//
//    if(listSize(curCont->adjacentBackwardContours) == 1) {
//      /* just add this as a one to one correspondence */
//      cp = (contourPair*) malloc(sizeof(contourPair));
//      cp->c1 = cp->c1Origin =curCont;
//      cp->c2 = cp->c2Origin =(contour*) getListNode(curCont->adjacentBackwardContours,0)->data;
//      enqueue(forwardContourPairs,cp);
//    }
//    else if(listSize(curCont->adjacentBackwardContours) < 1) {
//      continue;
//    }
//    else {
//      /* found one, split it up */
//      splitContoursMLE(curCont,curCont->adjacentBackwardContours,&closestContours);
//
//// DEBUGGING
////      for(ind = 0; ind < listSize(curCont->vertices); ind++) {
////      	fprintf(stdout,"%d %p\n", ind, closestContours[ind]);
////      }
//
//      /* add contour pairs between branches */
//      tmpContourPairs = newList(LIST);
//      addContourPairsFromClosest(curCont,closestContours,tmpContourPairs);
//      appendList(forwardContourPairs,tmpContourPairs);
//      free(tmpContourPairs);
//
//      if((contour*) curCont->origin != curCont) {
//	free(curCont);
//      }
//      free(closestContours);
//    }
//  }
//  freeList(newSlice);
//
//  /* print for debugging */
//  if(SR_DEBUG) {
//    fprintf(stderr,"------------------------\n");
//    for(i = getListNode(forwardContourPairs,0); i; i = (listNode*) i->next) {
//      dumpContourPair(stdout,(contourPair*) i->data);
//    }
//    fprintf(stderr,"------------------------\n");
//  }
//
//  newSlice = buildBackwardAdjacencyFromContourPairs(forwardContourPairs);
//  freeListAndData(forwardContourPairs);
//
//  backContourPairs = newList(LIST);
//  // do backward contours again
//  for(i = getListNode(newSlice,0); i; i = (listNode*) i->next) {
//    curCont = (contour*) i->data;
//
//    if(listSize(curCont->adjacentBackwardContours) == 1) {
//      /* just add this as a one to one correspondence */
//      cp = (contourPair*) malloc(sizeof(contourPair));
//      cp->c1 = cp->c1Origin =curCont;
//      cp->c2 = cp->c2Origin =(contour*) getListNode(curCont->adjacentBackwardContours,0)->data;
//      enqueue(backContourPairs,cp);
//    }
//    else if(listSize(curCont->adjacentBackwardContours) < 1) {
//      continue;
//    }
//    else {
//      /* found one, split it up */
//
//      splitContoursMLE(curCont,curCont->adjacentBackwardContours,
//		       &closestContours);
//
//// DEBUGGING
////      for(ind = 0; ind < listSize(curCont->vertices); ind++) {
////      	fprintf(stdout,"%d %p\n", ind, closestContours[ind]);
////      }
//      //for(j = 0; j < listSize(curCont->vertices); j++) {
//      //	fprintf(stdout,"%0.2lf %0.2lf %d %p\n", 
//      //		((vertex*) getListNode(curCont->vertices,j)->data)->x, 
//      //		((vertex*) getListNode(curCont->vertices,j)->data)->y, 
//      //		j, closestContours[j]);
//      //}
//
//      /* add contour pairs between branches */
//      tmpContourPairs = newList(LIST);
//      addContourPairsFromClosest(curCont,closestContours,tmpContourPairs);
//      appendList(backContourPairs,tmpContourPairs);
//      free(tmpContourPairs);
//
//      if((contour*) curCont->origin != curCont) {
//	free(curCont);
//      }
//      free(closestContours);
//    }
//  }
//  freeList(newSlice);


  /* print for debugging */
  if(SR_DEBUG) {
    fprintf(stderr,"------------------------\n");
    for(i = getListNode(backContourPairs,0); i; i = (listNode*) i->next) {
      dumpContourPair(stdout,(contourPair*) i->data);
    }
    fprintf(stderr,"------------------------\n");
  }

  /* reverse and clone the contour pair for tiling */
  contourPairs = newList(LIST);
  for(i = getListNode(backContourPairs,0); i; i = (listNode*) i->next) {
    cp = (contourPair*) i->data;
    newcp = (contourPair*) malloc(sizeof(contourPair));
    newcp->c1 = newcp->c1Origin = cloneContour(cp->c2);
    newcp->c2 = newcp->c2Origin = cloneContour(cp->c1);
    enqueue(contourPairs,newcp);
  }

  // delete all contours 
  for(i = getListNode(backContourPairs,0); i; i = (listNode*) i->next) {
    cp = (contourPair*) i->data;
    if((contour*) cp->c1->origin != cp->c1) {
      deleteContour(cp->c1);
    }
    if((contour*) cp->c2->origin != cp->c2) {
      deleteContour(cp->c2);
    }
  }
  freeListAndData(backContourPairs);  

  /* close contours that should be */
  for(i = getListNode(contourPairs,0); i; i = (listNode*) i->next) {
    cp = (contourPair*) i->data;
    
    if(cp->c1->closed == CLOSED && cp->c2->closed == CLOSED) continue;

    if(cp->c1->closed == BRANCHED_OPEN 
       && listSize(((contour*)cp->c1->origin)->vertices) 
       == listSize(cp->c1->vertices)) {
      cp->c1->closed = CLOSED;
    }

    if(cp->c2->closed == BRANCHED_OPEN 
       && listSize(((contour*)cp->c2->origin)->vertices) 
       == listSize(cp->c2->vertices)) {
      cp->c2->closed = CLOSED;
    }
  }


  /* free all the stuff we created */
  for(i = getListNode(slice1,0); i; i = (listNode*) i->next) {
    deleteContour((contour*) i->data);
  }
  freeList(slice1);

  for(i = getListNode(slice2,0); i; i = (listNode*) i->next) {
    deleteContour((contour*) i->data);
  }
  freeList(slice2);

  return contourPairs;
}

// REMOVED FUNCTIONALITY
/**
 * branch slices using vertex index arrays instead of vertex pointers and lists 
 */
//list *branchSlicesVertIndex(list *slice1, list *slice2) {
//  list *contourPairs = newList(LIST);
//  list *vertIndices1 = newList(LIST);
//  list *vertIndices2 = newList(LIST);
//  list *verts;
//
//  int *indArr;
//
//  listNode *i,*j;
//
//  // build the vert arrays for slice 1 contours
//  for(i = getListNode(slice1,0); i; i = (listNode*) i->next) {
//    verts = ((contour*) i->data)->vertices;
//    indArr = (int*) malloc(listSize(verts)*sizeof(int));
//    enqueue(vertIndices1,indArr);
//  }
//  
//  // build the vert arrays for slice 2 contours
//  for(i = getListNode(slice2,0); i; i = (listNode*) i->next) {
//    verts = ((contour*) i->data)->vertices;
//    indArr = (int*) malloc(listSize(verts)*sizeof(int));
//    enqueue(vertIndices2,indArr);
//  }
//  
//  // 
//}

/**
 * split contours via maximum likelihood estimation of boundaries
 */
void splitContoursMLE(contour *cont, list *adjacent,
		      contour ***closestContours) {
  listNode *i;
  int indi, indj, maxVotes, maxCont;
  int *votes;


  /* get the closest contour for each vertex */
  (*closestContours)
    = (contour**) malloc(listSize(cont->vertices)*sizeof(contour*));
  for(i = getListNode(cont->vertices,0), indi = 0; i;
      i = (listNode*) i->next, indi++) {
    (*closestContours)[indi] = getClosestContour((vertex*) i->data,adjacent);
    //fprintf(stderr,"%d %p\n",indi,(*closestContours)[indi]);
  }


  /* perform boundary finding */
  for(indi = 0; indi < listSize(cont->vertices); indi++) {
    /* initialize the votes for this vertex */
    votes = (int*) malloc(listSize(adjacent)*sizeof(int));
    for(indj = 0; indj < listSize(adjacent); indj++) {
      votes[indj] = 0;
    }

    /* go over clique, accumulating votes for each contour */
    for(indj = indi-cliqueRadius; indj <= indi+cliqueRadius; indj++) {
      if(indj < 0 || indj >= listSize(cont->vertices)) continue;
      votes[findInListI(adjacent,(*closestContours)[indj])]++;
    }

    //    fprintf(stdout,"%d %d %d %d\n", indi, votes[0], votes[1], votes[2]);

    /* tally votes for this vertex, reassign closest contour */
    maxVotes = votes[0];
    maxCont = 0;
    for(indj = 1; indj < listSize(adjacent); indj++) {
      if(votes[indj] > maxVotes) {
	maxVotes = votes[indj];
	maxCont = indj;
      }
    }

    (*closestContours)[indi] = (contour*) getListNode(adjacent,maxCont)->data;

    free(votes);
  }
}

/**
 * create contour pairs and add them to a list from a set of closest contours
 */
void addContourPairsFromClosest(contour *curCont, contour **closestContours,
				list *contourPairs) {
  contourPair *cp;//, *otherCp;
  listNode *i;//, *j;
  vertex *v;
  int indi;//, lastMatchInd, thisMatchInd;
  //int found;//, 
  //int n = listSize(curCont->vertices);

  /* set up and add the first contour pair */
  cp = (contourPair*) malloc(sizeof(contourPair));
  cp->c1 = createContour();
  cp->c1->origin = curCont->origin;
  cp->c1Origin = curCont;

  cp->c2Origin = closestContours[0];
  cp->c2 = closestContours[0];

  enqueue(contourPairs,cp);


  /* closure */
  if(curCont->closed == CLOSED) {
    cp->c1->closed = BRANCHED_OPEN;
  }
  else {
    cp->c1->closed = OPEN;
  }

  if(cp->c2->closed == CLOSED) {
    cp->c2->closed = BRANCHED_OPEN;
  }
  else {
    cp->c2->closed = OPEN;
  }


  /* add the first vertex */
  enqueue(cp->c1->vertices,
	  copyVertex((vertex*)getListNode(curCont->vertices,0)->data));

// REMOVED FUNCTIONALITY
//  // get the vertex closest match
//  lastMatchInd = getClosestVertexInContourInd(cp->c2,
//		     (vertex*)getListNode(curCont->vertices,0)->data);
//

// REMOVED FUNCTIONALITY
  /* test for first and last not the same so we can label a branched boundary*/
//  if(curCont->closed != OPEN 
//     && closestContours[0] != closestContours[n-1]) {
//    v = (vertex*) getListNode(cp->c1->vertices,0)->data;
//    v->boundary = BRANCHED_BOUNDARY;
//    v = (vertex*) getListNode(cp->c1->vertices,listSize(cp->c1->vertices)-1)->data;
//    v->boundary = BRANCHED_BOUNDARY;
//  }


  /* find branches */
  for(i = getListNode(curCont->vertices,1), indi = 1; i;
      i = (listNode*) i->next, indi++) {
// REMOVED FUNCTIONALITY
//    thisMatchInd = getClosestVertexInContourInd(cp->c2,(vertex*) i->data);

    /* look for branch as either a discontinutity in the closest contour or
       a half contour discontinutity in the closest vertex index HEURISTIC */
    if(closestContours[indi-1] != closestContours[indi] 
// REMOVED FUNCTIONALITY
//       || abs(thisMatchInd-lastMatchInd) > listSize(cp->c2->vertices)/2
       ) {

// REMOVED FUNCTIONALITY
      /* add a copy of this vertex to the last contour */
      //if(closestContours[indi]->closed == OPEN) {
      //enqueue(cp->c1->vertices,copyVertex((vertex*)i->data));
	//}

      // label back side of the branch as on a branched boundary
      //v = (vertex*) getListNode(cp->c1->vertices,
      //				listSize(cp->c1->vertices)-1)->data;
      //v->boundary = BRANCHED_BOUNDARY;

// REMOVED FUNCTIONALITY
//      /* see if this contour already has a contour pair */
//     found = 0;
//     for(j = getListNode(contourPairs,0); j; j = (listNode*) j->next) {
//	otherCp = (contourPair*) j->data;
//	if(otherCp->c2Origin == closestContours[indi]) {
//	  found = 1;
//	  cp = otherCp;
//	  break;
//	}
//     }      

//      if(!found) {
	/* setup a new contour pair */
      cp = (contourPair*) malloc(sizeof(contourPair));
      cp->c1 = createContour();
      cp->c1->origin = curCont->origin;
      cp->c1Origin = curCont;
	
      cp->c2Origin = closestContours[indi];
      cp->c2 = closestContours[indi];
	
      /* closure */
      if(curCont->closed == CLOSED) {
	cp->c1->closed = BRANCHED_OPEN;
      }
      else {
	cp->c1->closed = OPEN;
      }

      if(cp->c2->closed == CLOSED) {
	cp->c2->closed = BRANCHED_OPEN;
      }
      else {
	cp->c2->closed = OPEN;
      }
      enqueue(contourPairs,cp);
//      }

      // label front side of the branch as on a branched boundary
      v = copyVertex((vertex*)i->data);
      //v->boundary = BRANCHED_BOUNDARY;
      enqueue(cp->c1->vertices,v);

// REMOVED FUNCTIONALITY
//      thisMatchInd = getClosestVertexInContourInd(cp->c2,(vertex*) i->data); 
    }
    else {
      /* just add the vertex */
      enqueue(cp->c1->vertices,copyVertex((vertex*)i->data));
    }


// REMOVED FUNCTIONALITY
//    lastMatchInd = thisMatchInd;
  }


// REMOVED FUNCTIONALITY
//  /* orient the contours */
//  for(i = getListNode(contourPairs,0); i; i = (listNode*) i->next) {
//    cp = (contourPair*) i->data;
//    if(contourArea(cp->c1) < 0) {
//      reverseList(cp->c1->vertices);
//    }
//    if(contourArea(cp->c2) < 0) {
//      reverseList(cp->c2->vertices);
//    }
//  }
}

/**
 * take in a list of contour pairs that are connected to the split
 * contour and split the target contour into many smaller ones, each
 * of which corresponds to one of the contours in the list of contour
 * pairs passed.  returns a vector of contour pair index assignments
 * for each vertex in the target contour
 *
 * NOTE! POTENTIAL ERROR: IF A POSTBRANCH CONTOUR HAS NO VERTICES IN
 * THE PREBRANCH CONTOUR THAT ARE CLOSEEST TO IT, SEGFAULT WILL HAPPEN!
 *
 * step by step:
 * 1) assign a closest contour to each vertex in the prebranch contour
 * 2) find possible transitions
 * 3) assign each region of constant contour a closest contour, respecting
 *    the number of contours, by iterating over permutations of closest
 *    contour and iterating over possible transition regions within that
 * 4) choose the configuration of transitions and closest contour assignments
 *    with the fewest misclassified vertices
 */
void splitContoursTransitions(contour *splitContour, list *contours,
			      int *contourAssignments) {
  contour **closestContours;
  vertex *vert;
  list *possibleTransitions = newList(LIST);
  listNode *ln;
  int i, numTransitions, numMisclassified, bestMisclassified, curTransition,
    numVertices = listSize(splitContour->vertices),
    numSections = listSize(contours);
  intNode *in;
  gsl_permutation *order, *bestOrder;
  gsl_combination *transitions = NULL, *bestTransitions;

  /* allocate space to keep the list of closest contours */
  closestContours = (contour**) malloc(numVertices * sizeof(contour*));

  /* iterate over vertices, finding the closest contour for each */
  for(ln = getListNode(splitContour->vertices,0), i = 0; ln;
      ln = (listNode*) ln->next, i++){
    vert = (vertex*) ln->data;
    closestContours[i] = getClosestContour(vert, contours);
  }

  /* identify possible transitions */
  for(i = 0; i < numVertices-1; i++) {
    /* test for transition */
    if(closestContours[i] != closestContours[i+1]) {
      in = (intNode*) malloc(sizeof(intNode));
      in->val = i;
      enqueue(possibleTransitions,in);
    }
  }

  /** test permutations for smallest number of misclassified vertices **/

  /* initialize perm and comb engines */
  numTransitions = listSize(possibleTransitions);
  order = gsl_permutation_calloc(numSections);

  /* initialize the best */
  bestOrder = gsl_permutation_calloc(numSections);
  bestTransitions = gsl_combination_calloc(numTransitions,numSections-1);
  bestMisclassified = getNumMisclassifiedVertices(contours,
						  numVertices,
						  closestContours,
						  bestOrder,
						  possibleTransitions,
						  bestTransitions);

  /* run the perms to find the best section order */
  do {
    transitions = gsl_combination_calloc(numTransitions,numSections-1);
    /* run the combinations to find the best transition location */
    do {
      /* test this order and combination */
      numMisclassified = getNumMisclassifiedVertices(contours,
						     numVertices,
						     closestContours,
						     order,
						     possibleTransitions,
						     transitions);
      /* see if its best */
      if(numMisclassified < bestMisclassified) {
	bestMisclassified = numMisclassified;
	gsl_permutation_memcpy(bestOrder,order);
	gsl_combination_memcpy(bestTransitions,transitions);
      }
    } while(gsl_combination_next(transitions) == GSL_SUCCESS);

    gsl_combination_free(transitions);
  } while(gsl_permutation_next(order) == GSL_SUCCESS);

  /* build the contour assignment vector based on the best order */
  curTransition = 0;
  for(i = 0; i < numVertices; i++) {
    /* assign this contour */
    contourAssignments[i] = gsl_permutation_get(bestOrder,curTransition);

    /* transition to the next contour, if necessasry */
    if(curTransition < gsl_combination_k(bestTransitions)
       && i == ((intNode*)
		getListNode(possibleTransitions,
			    gsl_combination_get(bestTransitions,
						curTransition))->data)->val) {
      curTransition++;
    }
  }

  gsl_combination_free(bestTransitions);
}

/**
 * take in a list of contour pairs that all have the same c2Origin and
 * split the target contour into many smaller ones, each of which
 * corresponds to one of the c1 contours in the list of contour pairs
 * passed.
 *
 */
void splitContoursManyToOne(list *contourPairs) {
  contour *splitContour,*addContour;
  int *contourAssignments;
  list *contours = newList(LIST);
  listNode *ln;

  /* get the contour to be split */
  splitContour = ((contourPair*) getListNode(contourPairs,0)->data)->c2Origin;

  /* build the list of contours */
  for(ln = getListNode(contourPairs,0); ln; ln = (listNode*) ln->next) {
    /* add the modified contour, if possible, default to original one */
    addContour = ((contourPair*)ln->data)->c1;
    if(addContour == NULL) {
      addContour = ((contourPair*)ln->data)->c1Origin;
    }

    enqueue(contours,addContour);
  }

  /* allocate array of contour assignment */
  contourAssignments =
    (int*) malloc(listSize(splitContour->vertices)*sizeof(int));

  /* split the contour */
  splitContoursTransitions(splitContour, contours, contourAssignments);

  /* create new contour pairs from the best segmentation */
  splitC2Origin(contourPairs,contourAssignments);

  /** free stuff **/
  freeList(contours);
  free(contourAssignments);
}

/**
 * take in a list of contour pairs that all have the same c1Origin and
 * split the target contour into many smaller ones, each of which
 * corresponds to one of the c2 contours in the list of contour pairs
 * passed
 */
void splitContoursOneToMany(list *contourPairs) {
  contour *splitContour, *addContour;
  int *contourAssignments;
  list *contours = newList(LIST);
  listNode *ln;

  /* get the contour to be split  */
  splitContour = ((contourPair*) getListNode(contourPairs,0)->data)->c1Origin;

  /* build the list of contours */
  for(ln = getListNode(contourPairs,0); ln; ln = (listNode*) ln->next) {
    /* get the contour to be added from the handled contours, if possible */
    addContour = ((contourPair*)ln->data)->c2;
    if(addContour == NULL) {
      addContour = ((contourPair*)ln->data)->c2Origin;
    }

    enqueue(contours,addContour);
  }

  /* allocate array of contour assignment */
  contourAssignments =
    (int*) malloc(listSize(splitContour->vertices)*sizeof(int));

  /* split the contour */
  splitContoursTransitions(splitContour, contours, contourAssignments);

  /* create new contour pairs from the best segmentation */
  splitC1Origin(contourPairs,contourAssignments);

  /** free stuff **/
  freeList(contours);
  free(contourAssignments);
}

/**
 * use a contour assignment vector to create actual contour pairs
 * corresponding to a branch where all original c2 contours were the
 * same. so break up the c2 origin and create new contours
 */
void splitC2Origin(list *contourPairs, int *contourAssignments) {
  list *origContourVerts = ((contourPair*) getListNode(contourPairs,0)->data)->c2Origin->vertices;
  listNode *ln;
  contourPair *curContourPair, *lastContourPair;
  int curVertexInd;

  /* assign proper contour vertices */
  curContourPair = (contourPair*) getListNode(contourPairs,0)->data;
  if(curContourPair->c2 == NULL) {
    origContourVerts = curContourPair->c2Origin->vertices;
  }
  else {
    origContourVerts = curContourPair->c2->vertices;
  }

  /* allocate new contours */
  for(ln = getListNode(contourPairs,0); ln; ln = (listNode*) ln->next) {
    curContourPair = (contourPair*) ln->data;

    /* copy over the original vertex list if this one is null */
    if(curContourPair->c1 == NULL) {
      curContourPair->c1 = cloneContour(curContourPair->c1Origin);
    }
    curContourPair->c2 = (contour*) malloc(sizeof(contour));
    curContourPair->c2->closed = curContourPair->c2Origin->closed;
    curContourPair->c2->adjacentContours = NULL;
    curContourPair->c2->vertices = newList(LIST);
  }

  /* iterate over vertices, adding each to its contour */
  lastContourPair = (contourPair*)
    getListNode(contourPairs, contourAssignments[0])->data;
  for(ln = getListNode(origContourVerts,0), curVertexInd = 0;
      ln; ln = (listNode*) ln->next, curVertexInd++) {
    curContourPair =
      (contourPair*) getListNode(contourPairs,
				 contourAssignments[curVertexInd])->data;

    /* at transitions duplicate end vertex to eliminate discontinuity in surf*/
    if(lastContourPair != curContourPair) {
      enqueue(lastContourPair->c2->vertices, (vertex*) ln->data);
    }

    /* add the vertex to the c2 */
    enqueue(curContourPair->c2->vertices, (vertex*) ln->data);
    lastContourPair = curContourPair;
  }
}

/**
 * use a contour assignment vector to create actual contour pairs
 * corresponding to a branch where all original c1 contours were the
 * same. so break up the c1 origin and create new contours
 */
void splitC1Origin(list *contourPairs, int *contourAssignments) {
  list *origContourVerts;
  listNode *ln;
  contourPair *curContourPair, *lastContourPair;
  int curVertexInd;

  /* assign proper contour vertices */
  curContourPair = (contourPair*) getListNode(contourPairs,0)->data;
  if(curContourPair->c1 == NULL) {
    origContourVerts = curContourPair->c1Origin->vertices;
  }
  else {
    origContourVerts = curContourPair->c1->vertices;
  }

  /* allocate new contours */
  for(ln = getListNode(contourPairs,0); ln; ln = (listNode*) ln->next) {
    curContourPair = (contourPair*) ln->data;

    /* copy over the original vertex list if this one is null */
    if(curContourPair->c2 == NULL) {
      curContourPair->c2 = cloneContour(curContourPair->c2Origin);
    }
    curContourPair->c1 = (contour*) malloc(sizeof(contour));
    curContourPair->c1->closed = curContourPair->c1Origin->closed;
    curContourPair->c1->adjacentContours = NULL;
    curContourPair->c1->vertices = newList(LIST);
  }

  /* iterate over vertices, adding each to its contour. at
     transitions, add to both */
  lastContourPair = (contourPair*)
    getListNode(contourPairs, contourAssignments[0])->data;
  for(ln = getListNode(origContourVerts,0), curVertexInd = 0;
      ln; ln = (listNode*) ln->next, curVertexInd++) {
    curContourPair =
      (contourPair*) getListNode(contourPairs,
				 contourAssignments[curVertexInd])->data;

    /* at transitions duplicate end vertex to eliminate discontinuity in surf*/
    if(lastContourPair != curContourPair) {
      enqueue(lastContourPair->c1->vertices, (vertex*) ln->data);
    }

    /* add the vertex to the c1 */
    enqueue(curContourPair->c1->vertices, (vertex*) ln->data);
    lastContourPair = curContourPair;
  }
}

/**
 * get the closest vertex in a to a contour to a given vertex
 */
vertex *getClosestVertexInContour(contour *cont, vertex *v) {
  listNode *ln;
  vertex *closestV;
  double minDist,curDist;

  /* validate */
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) == 0
     || v == NULL) {
    return NULL;
  }

  /* set the min to the first */
  minDist = dist(*v,*((vertex*)getListNode(cont->vertices,0)->data));
  closestV = (vertex*)getListNode(cont->vertices,0)->data;

  /* iterate over vertices, testing the distance for each */
  for(ln = getListNode(cont->vertices,1); ln; ln = (listNode*) ln->next) {
    curDist = dist(*v, *(vertex*)ln->data);
    if(curDist < minDist) {
      minDist = curDist;
      closestV = (vertex*) ln->data;
    }
  }

  return closestV;
}

/**
 * get the closest vertex in a to a contour to a given vertex, return index
 */
int getClosestVertexInContourInd(contour *cont, vertex *v) {
  listNode *ln;
  int i, minind = 0;
  double minDist,curDist;

  /* validate */
  if(cont == NULL || cont->vertices == NULL || listSize(cont->vertices) == 0
     || v == NULL) {
    return -1;
  }

  /* set the min to the first */
  minDist = dist(*v,*((vertex*)getListNode(cont->vertices,0)->data));

  /* iterate over vertices, testing the distance for each */
  for(ln = getListNode(cont->vertices,1), i = 0; ln; i++, 
	ln = (listNode*) ln->next) {
    curDist = dist(*v, *(vertex*)ln->data);
    if(curDist < minDist) {
      minDist = curDist;
      minind = i;
    }
  }

  return minind;
}

/**
 * get the closest contour to a vertex
 */
contour *getClosestContour(vertex *v, list *contours) {
  listNode *ln;
  float minDist,curDist;
  contour *minContour;

  /* validate */
  if(v == NULL || contours == NULL || listSize(contours) == 0) return NULL;

  /* set the first contour as the min */
  ln = getListNode(contours,0);
  minDist = getMinDistance(v, ((contour*) ln->data)->vertices);
  minContour = (contour*) ln->data;

  /* iterate over contours, finding one with smallest min distance */
  for(ln = getListNode(contours,1); ln; ln = (listNode*) ln->next) {
    curDist = getMinDistance(v, ((contour*) ln->data)->vertices);
    if(curDist < minDist) {
      minDist = curDist;
      minContour = (contour*) ln->data;
    }
  }

  return minContour;
}

/**
 * get the minimum distance from a vertex to a set of vertices
 */
float getMinDistance(vertex *v, list *vertices) {
  if(v == NULL || vertices == NULL || listSize(vertices) < 1) {
    return SR_BIG;
  }

  float minDist;
  float curDist;
  listNode *ln;

  /* set the min to the first */
  minDist = dist(*v,*(vertex*)getListNode(vertices,0)->data);

  /* iterate over vertices, testing the distance for each */
  for(ln = getListNode(vertices,1); ln; ln = (listNode*) ln->next) {
    curDist = dist(*v,*(vertex*)ln->data);
    if(curDist < minDist) minDist = curDist;
  }

  return minDist;
}

/**
 * evaluate a particular order and location of constant conoutr assignment
 */
int getNumMisclassifiedVertices(list *contours,
				int numVertices,
				contour **closestContours,
				gsl_permutation *order,
				list *possibleTransitions,
				gsl_combination *transitions) {
  int misclassified = 0;
  contour *curContour = (contour*) getListNode(contours,
					       gsl_permutation_get(order,
								   0))->data;
  int curVertexInd, curTransitionInd = 0;

  /* iterate over vertices, testing each for misclassification */
  for(curVertexInd = 0; curVertexInd < numVertices; curVertexInd++) {
    /* test for misclassification */
    if(curContour != closestContours[curVertexInd]) {
      misclassified++;
    }

    /* transition, if needed */
    if(curTransitionInd < gsl_combination_k(transitions)
       && curVertexInd == ((intNode*) getListNode(possibleTransitions,gsl_combination_get(transitions,curTransitionInd))->data)->val) {
      curTransitionInd++;
      curContour = (contour*) getListNode(contours,gsl_permutation_get(order,curTransitionInd))->data;
    }
  }

  return misclassified;
}

/**
 * searches two branched contours for vertices that were adjacent on the
 * original contour
 */
int matchBranchedContours(contour *c1, contour *c2, 
			  int *c10, int *c11, int *c20, int *c21) {
  if(c1 == NULL || c1->vertices == NULL || listSize(c1->vertices) < 1 
     || c2 == NULL || c2->vertices == NULL || listSize(c2->vertices) < 1 
     || c2->origin != c1->origin) {
    return 0;
  }

  int c1ind, c2ind,c1OrigInd,c2OrigInd;
  listNode *i, *j;

  // find a discontinuity in c1
    // find the adjacent entry in c2
  c1OrigInd = ((vertex*) getListNode(c1->vertices,listSize(c1->vertices)-1)->data)->number;
  for(i = getListNode(c1->vertices,0), c1ind = 0; i; i = (listNode*) i->next, c1ind++) {
    if(abs(((vertex*) i->data)->number-c1OrigInd) > 1) {
      *c10 = c1OrigInd;
      *c11 = ((vertex*) i->data)->number;
      *c20 = -1;
      *c21 = -1;
      c2OrigInd = ((vertex*) getListNode(c2->vertices,listSize(c2->vertices)-2)->data)->number;
      for(j = getListNode(c2->vertices,0),c2ind = 0; j; j = (listNode*) j->next, c2ind++) {
	if((abs(c2OrigInd-*c10) == 1 && 
	    abs(((vertex*) j->data)->number-*c11) == 1)
	   ||
	   (abs(c2OrigInd-*c11) == 1 && 
	    abs(((vertex*) j->data)->number-*c10) == 1)) {	   

	  // found it, return
	  *c10 = (c1ind-1 < 0) ? listSize(c1->vertices)-1 : c1ind-1;
	  *c11 = c1ind;
	  *c20 = (c2ind-1 < 0) ? listSize(c2->vertices)-1 : c2ind-1;
	  *c21 = c2ind;
	  return 1;
	}
	c2OrigInd = ((vertex*) j->data)->number;
      }
    }
    c1OrigInd = ((vertex*) i->data)->number;
  }

  return 0;
}


/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/branch.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
