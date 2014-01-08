/*****************************************************************************
 * correspondence.c is the source file with functions to take a guess at
 * a solution to the correspondence problem for the brain surface
 * reconstruction problem for libsr.
 * Oliver Hinds <oph@bu.edu> 2005-10-18
 *
 *
 *
 *****************************************************************************/

#include"correspondence.h"

int correspondenceNumMeans = 2;

/* method to use for correspondence */
enum CORRESPONDENCE_METHOD correspondenceMethod = DENSITY;
enum CORRESPONDENCE_SCOPE correspondenceScope = GLOBAL;

/* min and max percentiles to consider in global tiling */
float minPercentile = 5.0;
float maxPercentile = 85.0;

int saveContourDistances = 0;
char contourDistancesFilename[SR_MAX_STR_LEN] = "contour_distances.txt";

/**
 * builds a guess at contour correspondence
 */
int buildCorrespondenceGuess(list *slices) {
  switch(correspondenceScope) {
  case LOCAL:
    return buildCorrespondenceGuessLocal(slices);
    break;
  case GLOBAL:
    return buildCorrespondenceGuessGlobal(slices);
    break;
  default:
    return buildCorrespondenceGuessGlobal(slices);
    break;
  }
}

/**
 * builds a guess at contour correspondence slice by slice
 */
int buildCorrespondenceGuessLocal(list *slices) {
  listNode *i;
  list *thisSlice, *nextSlice;
  double **dm;
  int indi;
  long m, n;
  contour *c1,*c2;

  /* validate */
  if(listSize(slices) == 0) {
    return SR_FAILURE;
  }

  clearAllCorrespondence(slices);

  /* for each pair of slices, correspond the contour lists */
  for(i = getListNode(slices,0), indi = 0; i->next;
      i = (listNode*) i->next, indi++) {
    thisSlice = (list*) i->data;
    nextSlice = (list*) ((listNode*)i->next)->data;

    /* if there are not contours on either slice, continue */
    if(listSize(thisSlice) == 0 || listSize(nextSlice) == 0) {
      continue;
    }

    if(SR_DEBUG) {
      fprintf(stderr,"corresponding contours on slices %d and %d\n",indi,indi+1);
    }

    /* CAUTION! if both slices contain only one contour, connect them by
       default, this could be wrong!!!*/
    if(listSize(thisSlice) == 1 && listSize(nextSlice) == 1) {
      c1 = (contour*) getListNode(thisSlice,0)->data;
      c2 = (contour*) getListNode(nextSlice,0)->data;

      enqueue(c1->adjacentContours,c2);

      if(SR_DEBUG) {
	fprintf(stderr,"only one contour on each slice, connecting them\n");
      }

      continue;
    }

    /* get the distance matrix for the two slices */
    dm = getContourDistanceMatrixLocal(thisSlice, nextSlice, &m, &n);

    /* make the correspondence */
    switch(correspondenceMethod) {
    case HISTOGRAM:
    case DENSITY:
      correspondContoursHistogram(thisSlice,nextSlice,dm,m,n);
      break;
    case KMEANS:
      correspondContoursKmeans(thisSlice,nextSlice,dm,m,n);
      break;
    default:
      correspondContoursKmeans(thisSlice,nextSlice,dm,m,n);
      break;
    }

    freeDistanceMatrix(dm,m,n);
  }

  return SR_SUCCESS;
}

/**
 * builds a guess at contour correspondence from all slices simultaneously
 */
int buildCorrespondenceGuessGlobal(list *slices) {
  listNode *i,*dmln,*mln,*nln;
  list *thisSlice, *nextSlice,
    *dms = newList(LIST),
    *ms = newList(LIST),
    *ns = newList(LIST);


  long m,n;
  double connectionThreshold = SR_BIG, **dm;

  /* validate */
  if(listSize(slices) == 0) {
    return SR_FAILURE;
  }

  clearAllCorrespondence(slices);

  if(SR_VERBOSE) {
    fprintf(stdout, "computing contour distance matrix...");
    fflush(stdout);
  }

  /* for each pair of slices, store the dms for analysis */
  for(i = getListNode(slices,0); i->next; i = (listNode*) i->next) {
    thisSlice = (list*) i->data;
    nextSlice = (list*) ((listNode*)i->next)->data;

    /* if there are not contours on a slice, continue */
    if(listSize(thisSlice) == 0 || listSize(nextSlice) == 0) {
      continue;
    }

    /* store the distance matrix for the two slices */
    enqueue(dms,getContourDistanceMatrix(thisSlice, nextSlice, &m, &n));

    enqueue(ms,(int*)m); /* WATCH OUT, THE INTEGER VALUES ARE */
    enqueue(ns,(int*)n); /* USED AS POINTERS! */
  }

  if(SR_VERBOSE) {
    fprintf(stdout, "done\n");
  }

  /* make the global correspondence */
  switch(correspondenceMethod) {
  case DENSITY:
    connectionThreshold = getGlobalConnectionThresholdDensity(dms,ms,ns);
    break;
  case HISTOGRAM:
    connectionThreshold = getGlobalConnectionThresholdHistogram(dms,ms,ns);
    break;
  case KMEANS:
    connectionThreshold = getGlobalConnectionThresholdKmeans(dms,ms,ns);
    break;
  default:
    connectionThreshold = getGlobalConnectionThresholdKmeans(dms,ms,ns);
    break;
  }

  if(SR_VERBOSE) {
    fprintf(stdout, "assigning correspondence with threshold %lf...",
	    connectionThreshold);
  }

  /* assign the adjacency */
  setAdjacentContours(slices,connectionThreshold,dms,ms,ns);

  /* free dm lists */
  for(dmln = getListNode(dms,0),
	mln  = getListNode(ms,0),
	nln  = getListNode(ns,0) ;
      dmln;
      dmln = (listNode*) dmln->next,
	mln  = (listNode*) mln->next,
	nln  = (listNode*) nln->next
      ) {

    dm = (double**) dmln->data;
#ifdef BITS64
    m = (long) mln->data; /* TREATING POINTERS AS INTEGERS ON PURPOSE */
    n = (long) nln->data;
#else
    m = (int) mln->data; /* TREATING POINTERS AS INTEGERS ON PURPOSE */
    n = (int) nln->data;
#endif
    freeDistanceMatrix(dm,m,n);
  }

  freeList(dms);
  freeList(ms);
  freeList(ns);

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  return SR_SUCCESS;
}

/**
 * set adjacent contours by a correspondence threshold
 */
void setAdjacentContours(list *slices, double correspondenceThreshold,
			 list* dms, list* ms, list* ns) {
  listNode *i,*j,*k,*dmln,*mln,*nln;
  list *thisSlice, *nextSlice;
  contour *c1,*c2;
  int indi,indj,indk;
  double **dm;

  /* for each pair of slices, connect the contours with below threshold dist */
  for(indi = 0,
	i = getListNode(slices,0),
	dmln = getListNode(dms,0),
	mln  = getListNode(ms,0),
	nln  = getListNode(ns,0) ;
      i->next;
      indi++,
	i = (listNode*) i->next,
	dmln = (listNode*) dmln->next,
	mln  = (listNode*) mln->next,
	nln  = (listNode*) nln->next
      ) {
    thisSlice = (list*) i->data;
    nextSlice = (list*) ((listNode*)i->next)->data;

    /* CAUTION! if both slices contain only one contour, connect them by
       default, this could be wrong!!!*/
    if(listSize(thisSlice) == 1 && listSize(nextSlice) == 1) {
      c1 = (contour*) getListNode(thisSlice,0)->data;
      c2 = (contour*) getListNode(nextSlice,0)->data;

      enqueue(c1->adjacentContours,c2);

      if(SR_DEBUG) {
	fprintf(stderr,"only one contour on each slice, connecting them\n");
      }

      continue;
    }

    dm = (double**) dmln->data;

    /* connect all contours with below threshold distances */
    for(j = getListNode(thisSlice,0), indj=0; j; j=(listNode*) j->next, indj++) {
      for(k = getListNode(nextSlice,0), indk=0; k; k=(listNode*)k->next, indk++) {
	//fprintf(stderr,"s=%d, d%d%d=%0.3lf, t=%0.3lf\n",
	//	indi,indj,indk,dm[indj][indk],correspondenceThreshold);

	if(dm[indj][indk] < correspondenceThreshold) {
	  c1 = (contour*) j->data;
	  c2 = (contour*) k->data;
	  enqueue(c1->adjacentContours,c2);
	}
      }
    }
  }
}

/**
 * set the adjacency manually
 */
int setAdjacentContoursManual(list *slices, double connectionThreshold) {
  listNode *i,*dmln,*mln,*nln;
  list *thisSlice, *nextSlice,
    *dms = newList(LIST),
    *ms = newList(LIST),
    *ns = newList(LIST);
  long m,n;
  double **dm;

  /* validate */
  if(listSize(slices) == 0) {
    return SR_FAILURE;
  }

  clearAllCorrespondence(slices);

  /* for each pair of slices, store the dms for analysis */
  for(i = getListNode(slices,0); i->next; i = (listNode*) i->next) {
    thisSlice = (list*) i->data;
    nextSlice = (list*) ((listNode*)i->next)->data;

    /* if there are not contours on a slice, continue */
    if(listSize(thisSlice) == 0 || listSize(nextSlice) == 0) {
      continue;
    }

    /* store the distance matrix for the two slices */
    enqueue(dms,getContourDistanceMatrix(thisSlice, nextSlice, &m, &n));
    enqueue(ms,(long*)m); /* WATCH OUT, THE INTEGER VALUES ARE */
    enqueue(ns,(long*)n); /* USED AS POINTERS! */
  }

  setAdjacentContours(slices,connectionThreshold,dms,ms,ns);

  /* free dm lists */
  for(dmln = getListNode(dms,0),
	mln  = getListNode(ms,0),
	nln  = getListNode(ns,0) ;
      dmln->next;
      dmln = (listNode*) dmln->next,
	mln  = (listNode*) mln->next,
	nln  = (listNode*) nln->next
      ) {

    dm = (double**) dmln->data;
    m = (long) mln->data; /* TREATING POINTERS AS INTEGERS ON PURPOSE */
    n = (long) nln->data;

    freeDistanceMatrix(dm,m,n);
  }

  freeList(dms);
  freeList(ms);
  freeList(ns);

  return SR_SUCCESS;
}

/**
 * get contour distance matrix
 *
 * CAUTION!  this function adds a row to the dm when slice1 has only one
 * contour, and a column to the dm when slice2 has only one contour. this
 * 'dummy' contour aids in identifying terminating contours and many to one
 * splits in isolation.
 */
double **getContourDistanceMatrixLocal(list *slice1, list *slice2, long *m, long *n) {
  listNode *i, *j;
  int ind, c1ind, c2ind;
  double **dm, i2j, j2i, maxdist = 0.0, curmaxdist, dummydist;

  /* validate */
  if(listSize(slice1) == 0 || listSize(slice2) == 0) {
    return NULL;
  }

  /* assign dims */
  if(listSize(slice1) > 1) {
    *m = (long) listSize(slice1);
  }
  else {
    *m = 2;
  }

  if(listSize(slice2) > 1) {
    *n = (long) listSize(slice2);
  }
  else {
    *n = 2;
  }

  /* allocate the distance matrix */
  dm = (double**) malloc(*m*sizeof(double*));
  for(ind = 0; ind < *m; ind++) {
    dm[ind] = (double*) malloc(*n*sizeof(double));
  }

  /* match each contour across slice */
  for(i = getListNode(slice1,0), c1ind = 0; i; i = (listNode*) i->next, c1ind++) {
    for(j = getListNode(slice2,0), c2ind = 0; j; j = (listNode*) j->next, c2ind++) {
      i2j = getContoursDist((contour*) i->data, (contour*) j->data,
			    &curmaxdist);
      if(maxdist < curmaxdist) {
	maxdist = curmaxdist;
      }

      j2i = getContoursDist((contour*) j->data, (contour*) i->data,
			    &curmaxdist);
      if(maxdist < curmaxdist) {
	maxdist = curmaxdist;
      }

      dm[c1ind][c2ind] = min(i2j,j2i);
    }
  }

  /* caluculate dummy distance arbitrarily as double the max min dist */
  dummydist = 2*maxdist;


  /* add dummy row if needed */
  if(*m > listSize(slice1)) {
    for(c2ind = 0; c2ind < *n; c2ind++) {
      dm[*m-1][c2ind] = dummydist;
    }
  }

  /* add dummy col if needed */
  if(*n > listSize(slice2)) {
    for(c1ind = 0; c1ind < *m; c1ind++) {
      dm[c1ind][*n-1] = dummydist;
    }
  }

  return dm;
}

/**
 * get contour distance matrix
 */
double **getContourDistanceMatrix(list *slice1, list *slice2, long *m, long *n) {
  listNode *i, *j;
  int ind, c1ind, c2ind;
  double **dm, i2j, j2i, mindist = SR_BIG, curmindist;

  /* validate */
  if(listSize(slice1) == 0 || listSize(slice2) == 0) {
    return NULL;
  }

  /* assign dims */
  *m = (long) listSize(slice1);
  *n = (long) listSize(slice2);

  /* allocate the distance matrix */
  dm = (double**) malloc(*m*sizeof(double*));
  if(dm == (double**) 0x0814f0d8) {
    printf("here\n");
  }

  for(ind = 0; ind < *m; ind++) {
    dm[ind] = (double*) malloc(*n*sizeof(double));
  }

  /* match each contour across slice */
  for(i = getListNode(slice1,0), c1ind = 0; i; i = (listNode*) i->next, c1ind++) {
    for(j = getListNode(slice2,0), c2ind = 0; j; j = (listNode*) j->next, c2ind++) {
      i2j = getContoursDist((contour*) i->data, (contour*) j->data,
			    &curmindist);
      if(isnan(mindist) && mindist > curmindist) {
	mindist = curmindist;
      }

      j2i = getContoursDist((contour*) j->data, (contour*) i->data,
			    &curmindist);
      if(isnan(mindist) && mindist > curmindist) {
	mindist = curmindist;
      }

      dm[c1ind][c2ind] = min(i2j,j2i);
    }
  }

  for(c1ind = 0; c1ind < *m; c1ind++) {
    for(c2ind = 0; c2ind < *n; c2ind++) {
      if(isnan(dm[c1ind][c2ind])) {
	dm[c1ind][c2ind] = mindist;
      }
    }
  }

  return dm;
}

#define gaussian(x,y) 1/(sqrt(2*SR_PI)*y)*exp(-(x*x)/(2*y*y))
/**
 * get a correspondence threshold via estimating the pdf of the contour
 * distances, then finding the first mode of the pdf.
 */
double getGlobalConnectionThresholdDensity(list *dms, list *ms, list *ns) {
  double connectionThreshold = -1.0, sigma, d, bw;

  int numDatapoints = 0, i, j, n = 50;
  double *data = NULL;
  double *x = NULL;
  double *pdf = NULL;
  double *slope = NULL;
  double *fit = NULL;

  double *nx = NULL;
  double *ny = NULL;
  double *nw = NULL;
  int numNeigh = 0;

  if(SR_VERBOSE) {
    fprintf(stdout,"contour correspondence via global density estimation.\n");
  }

  /* build the data and mask beyond percentile thresholds */
  buildDataArray(dms,ms,ns,&data,&numDatapoints);
  maskPercentiles(data,&numDatapoints,minPercentile,maxPercentile);

  if(numDatapoints < 1) {
    return 0/0.0;
  }

  /* estimate the width of the gaussian as 5% of the range */
  sigma = 0.05*(data[numDatapoints-1]-data[0]);

  /* build the locations to estimate the pdf at */
  x = malloc(n*sizeof(double));
  for(i = 0; i < n; i++) {
    x[i] = data[0]+(data[numDatapoints-1]-data[0])*i/(double)(n-1);
  }

  /* estimate the pdf */
  pdf = malloc(n*sizeof(double));
  for(i = 0; i < n; i++) {
    pdf[i] = 0;
    for(j = 0; j < numDatapoints; j++) {
      pdf[i] += gaussian(fabs(x[i]-data[j]),sigma);
    }
    pdf[i]=pdf[i]/(n*sigma);
  }

  /* estimate the band width of the slope estimator as 10% of the range */
  bw = 0.10*(data[numDatapoints-1]-data[0]);

  /* estimate the slope of the pdf */
  slope = malloc(n*sizeof(double));

  /* set up the neighborhood data arrays */
  nx = malloc(n*sizeof(double));
  ny = malloc(n*sizeof(double));
  nw = malloc(n*sizeof(double));

  /* estimate slope at each data point */
  for(i = 0; i < n; i++) {

    /* build the neighbor arrays */
    numNeigh = 0;
    for(j = 0; j < n; j++) {
      d = fabs(x[i] - x[j]);
      if(d < bw) {
	nx[numNeigh] = x[j];
	ny[numNeigh] = pdf[j];
	nw[numNeigh] = 1-(d/bw)*(d/bw);
	numNeigh++;
      }
    }

    /* estimate the slope via linear least squares */
    if(numNeigh == 0) {
      slope[i] = 0/0.0;
      continue;
    }

    fit = linearFit(nx,ny,nw,numNeigh);
    slope[i] = fit[1];
    free(fit);
  }

  /* find the first rising zero crossing as the threshold */
  for(i = 1; i < n; i++) {
    if(!isnan(slope[i-1]) && !isnan(slope[i])
	      && slope[i-1] < 0 && slope[i] >= 0) {
      connectionThreshold = x[i-1]
	- slope[i-1]*(x[i-1]-x[i])/(slope[i-1]-slope[i]);
      //      connectionThreshold *= 1.5; // DELETE ME
      break;
    }
  }

  //if(SR_DEBUG) {
  //fprintf(stderr,"connection threshold = %lf\n", connectionThreshold);
  //}

  /* free stuff*/
  free(data);
  free(x);
  free(pdf);
  free(slope);
  free(nx);
  free(ny);
  free(nw);

  return connectionThreshold;
}

/**
 * determine the optimal connection threshold for contours via kmeans
 */
double getGlobalConnectionThresholdKmeans(list *dms, list *ms, list *ns) {
  double connectionThreshold = -1.0, eps = 0.0000001;

  int numDatapoints = 0, indi;
  double *data = NULL;
  double **kdata = NULL;
  double **means = NULL;
  int *clusters = NULL;

  int minMeanInd;

  if(SR_VERBOSE) {
    fprintf(stdout,"contour correspondence via global kmeans clustering.\n");
  }

  /* build the data and mask beyond percentile thresholds */
  buildDataArray(dms,ms,ns,&data,&numDatapoints);
  maskPercentiles(data,&numDatapoints,minPercentile,maxPercentile);

  /* convert to kmeans data format */
  kdata = (double**) malloc(numDatapoints*sizeof(double*));
  for(indi = 0; indi < numDatapoints; indi++) {
    kdata[indi] = (double*) malloc(sizeof(double));
    kdata[indi][0] = data[indi];
  }

  /* kmeans cluster */
  kmeans(kdata, numDatapoints, 1, correspondenceNumMeans, &clusters, &means);

  /* print debugging stuff */
  if(SR_DEBUG) {
    fprintf(stderr,"data:\t\tcluster:\n-----\t\t-------\n");
    for(indi = 0; indi < numDatapoints; indi++) {
      fprintf(stderr, "%2d: %5.5g\t%d\n", indi, data[indi],clusters[indi]);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"means:\n------\n");
    for(indi = 0; indi < correspondenceNumMeans; indi++) {
      fprintf(stderr, "%2d: %5.5g\n", indi, means[indi][0]);
    }
    fprintf(stderr,"\n");
  }

  /* find the lowest mean  */
  minMeanInd = 0;
  for(indi = 1; indi < correspondenceNumMeans; indi++) {
    if(means[indi][0] < means[minMeanInd][0]) {
      minMeanInd = indi;
    }
  }

  /* determine threshold */
  for(indi = 0; indi < numDatapoints; indi++) {
    if(clusters[indi] != minMeanInd) {
      continue;
    }

    if(data[indi] > connectionThreshold) {
      connectionThreshold = data[indi] + eps;
    }
  }

  /* free the data array */
  for(indi = 0; indi < numDatapoints; indi++) {
    free(kdata[indi]);
  }
  free(kdata);
  free(data);

  /* free clusters and means */
  for(indi = 0; indi < correspondenceNumMeans; indi++) {
    free(means[indi]);
  }
  free(means);
  free(clusters);

  return connectionThreshold;
}

/**
 * connect contours based on kmeans clustering of the log minimum distance
 * to each other.
 */
int correspondContoursKmeans(list *slice1, list* slice2, double **dm,
			     int m, int n) {
  int indi, indj;
  listNode *i, *j;
  contour *c1, *c2;

  int numDatapoints;
  double **data;
  double **means = NULL;
  int *clusters = NULL;

  int minMeanInd;

  /* validate */
  if(dm == NULL || listSize(slice1) == 0 || listSize(slice2) == 0) {
    return SR_FAILURE;
  }

  /* if there is only one slice on each contour, connect them and return */
  if(m == 1 && n == 1) {
    c1 = (contour*) getListNode(slice1,0)->data;
    c2 = (contour*) getListNode(slice2,0)->data;

    enqueue(c1->adjacentContours,c2);
    return SR_SUCCESS;
  }

  if(SR_VERBOSE) {
    fprintf(stdout,"contour correspondence via local kmeans clustering.\n");
  }

  /* build data vector */
  numDatapoints = m*n;
  data = (double**) malloc(numDatapoints*sizeof(double*));
  for(indi = 0; indi < numDatapoints; indi++) {
    data[indi] = (double*) malloc(sizeof(double));
    data[indi][0] = dm[indi/n][indi%n];
  }

  kmeans(data, numDatapoints, 1, correspondenceNumMeans, &clusters, &means);

  /* print debugging stuff */
  if(SR_DEBUG) {
    fprintf(stderr,"data:\t\tcluster:\n-----\t\t-------\n");
    for(indi = 0; indi < numDatapoints; indi++) {
      fprintf(stderr, "%2d: %5.5g\t%d\n", indi, data[indi][0],clusters[indi]);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"means:\n------\n");
    for(indi = 0; indi < correspondenceNumMeans; indi++) {
      fprintf(stderr, "%2d: %5.5g\n", indi, means[indi][0]);
    }
    fprintf(stderr,"\n");
  }

  /* find the lowest mean, call contours belonging to that cluster connected */
  minMeanInd = 0;
  for(indi = 1; indi < correspondenceNumMeans; indi++) {
    if(means[indi][0] < means[minMeanInd][0]) {
      minMeanInd = indi;
    }
  }

  /* connect contours in the minMean cluster */
  for(i = getListNode(slice1,0), indi = 0; i; i = (listNode*) i->next, indi++) {
    for(j = getListNode(slice2,0), indj = 0; j; j = (listNode*) j->next, indj++) {
      if(clusters[n*indi+indj] == minMeanInd) {
	c1 = (contour*) i->data;
	c2 = (contour*) j->data;

	enqueue(c1->adjacentContours,c2);
      }
    }
  }

  /* free data array */
  for(indi = 0; indi < numDatapoints; indi++) {
    free(data[indi]);
  }
  free(data);

  /* free clusters and means */
  for(indi = 0; indi < correspondenceNumMeans; indi++) {
    free(means[indi]);
  }
  free(means);
  free(clusters);

  return SR_SUCCESS;
}

/**
 * connect contours based on histogramming minimum distance to each
 * other. if there are multiple low distance contours, connect them all.
 */
int correspondContoursHistogram(list *slice1, list* slice2, double **dm,
				int m, int n) {
  listNode *i, *j;
  int indi, indj, numBins = 10;
  histogram *hist = NULL;
  contour *c1, *c2;
  int foundNonZero = FALSE;
  int foundThresh = FALSE;
  int maxBins = 1000;
  int binIncr = 20;
  double defaultMultiplier = 10.0;
  double connectionThreshold = 0.0; // dynamically determined

  /* validate */
  if(dm == NULL || listSize(slice1) == 0 || listSize(slice2) == 0) {
    return SR_FAILURE;
  }

  /* if there is only one slice on each contour, connect them and return */
  if(m == 1 && n == 1) {
    c1 = (contour*) getListNode(slice1,0)->data;
    c2 = (contour*) getListNode(slice2,0)->data;

    enqueue(c1->adjacentContours,c2);
    return SR_SUCCESS;
  }

  if(SR_VERBOSE) {
    fprintf(stdout,"contour correspondence via local histograming.\n");
  }

  /* find the first zero count in the histogram */
  while(!foundThresh) {
    if(hist != NULL) {
      deleteHistogram(hist);
      numBins += binIncr;
    }

    /* stopping condition, just in case */
    if(numBins > maxBins) {
      foundThresh = TRUE;
      connectionThreshold = defaultMultiplier*hist->minval;
    }

    /* build a histogram with the current number of bins */
    hist = buildHistogramDM(dm,m,n,numBins);

    for(indi = 0; indi < hist->numBins; indi++) {
      if(!foundNonZero && hist->hist[indi] == 0) {
	continue;
      }
      else if(hist->hist[indi] == 0) {
	foundThresh = TRUE;
	connectionThreshold = hist->minval + indi*hist->binSize + hist->binSize/2.0;
	break;
      }
      else {
	foundNonZero = TRUE;
      }
    }
  }

  /* print debugging stuff */
  if(SR_DEBUG) {
    fprintf(stderr,"threshold=%f\ndata:\t\tcluster:\n-----\t\t-------\n",connectionThreshold);
    for(i = getListNode(slice1,0), indi = 0; i; i = (listNode*) i->next, indi++) {
      for(j = getListNode(slice2,0), indj = 0; j; j = (listNode*) j->next, indj++) {
	fprintf(stderr, "%2d: %5.5g\t%d\n", indi, dm[indi][indj], dm[indi][indj] >= connectionThreshold);
      }
    }
    fprintf(stderr,"\n");
  }

  /* find the contours with distances below threshold */
  for(i = getListNode(slice1,0), indi = 0; i; i = (listNode*) i->next, indi++) {
    for(j = getListNode(slice2,0), indj = 0; j; j = (listNode*) j->next, indj++) {
      if(dm[indi][indj] < connectionThreshold) {
	c1 = (contour*) i->data;
	c2 = (contour*) j->data;

	enqueue(c1->adjacentContours,c2);
      }
    }
  }

  deleteHistogram(hist);
  return SR_SUCCESS;
}

/**
 * determine the optimal connection threshold for contours via histogram
 */
double getGlobalConnectionThresholdHistogram(list *dms, list *ms, list *ns) {
  int numBins = 10;
  histogram *hist = NULL;
  int numDatapoints = 0, indi;
  double *data = NULL;
  int foundNonZero = FALSE;
  int foundThresh = FALSE;
  int maxBins = 1000;
  int binIncr = 20;
  double defaultMultiplier = 10.0;
  double connectionThreshold = 0.0; // dynamically determined

  if(SR_VERBOSE) {
    fprintf(stdout,"contour correspondence via global histograming.\n");
  }

  buildDataArray(dms,ms,ns,&data,&numDatapoints);
  maskPercentiles(data,&numDatapoints,minPercentile,maxPercentile);

  /* find the first zero count in the histogram */
  while(!foundThresh) {
    if(hist != NULL) {
      deleteHistogram(hist);
      numBins += binIncr;
    }

    /* stopping condition, just in case */
    if(numBins > maxBins) {
      foundThresh = TRUE;
      connectionThreshold = defaultMultiplier*hist->minval;
    }

    /* build a histogram with the current number of bins */
    hist = buildHistogramMultiDM(data,numDatapoints,numBins);

    for(indi = 0; indi < hist->numBins; indi++) {
      if(!foundNonZero && hist->hist[indi] == 0) {
	continue;
      }
      else if(hist->hist[indi] == 0) {
	foundThresh = TRUE;
	connectionThreshold = hist->minval + indi*hist->binSize + hist->binSize/2.0;
	break;
      }
      else {
	foundNonZero = TRUE;
      }
    }
  }

  deleteHistogram(hist);
  free(data);

  return connectionThreshold;
}

/**
 * gets distance to the nearest vertex on a contour from a point
 */
double getDistanceToNearestVertex(contour *c, vertex *v) {
  listNode *i;
  double curDist, minDist;
  vertex *u;

  /* validate */
  if(c == NULL || c->vertices == NULL || v == NULL) {
    return SR_BIG;
  }

  minDist = SR_BIG;
  for(i = getListNode(c->vertices,0); i; i = (listNode*) i->next) {
    u = (vertex*) i->data;

    curDist = dist(*u,*v);

    if(curDist < minDist) {
      minDist = curDist;
    }
  }

  return minDist;
}

/**
 * gets the minimum distance to the contour from a vertex
 */
double getDistanceToContour(contour *c, vertex *v) {
  listNode *begPt, *endPt;
  vertex *v2;
  double curDist, minDist = SR_BIG;

  /* validate */
  if(c == NULL || c->vertices == NULL || listSize(c->vertices) == 0 || v == NULL) {
    return SR_BIG;
  }

  /* check for only one vertex */
  if(listSize(c->vertices) == 1) {
    // return the distance between the vertices
    v2 = (vertex*)getListNode(c->vertices,0)->data;
    return dist((*v2),(*v));
  }

  /* iterate over line segments, finding the min distance to any */
  for(begPt = getListNode(c->vertices,0), endPt = getListNode(c->vertices,1);
      endPt; begPt = endPt, endPt = (listNode*) endPt->next) {
    curDist = getDistanceToSegment((vertex*) begPt->data,
				   (vertex*) endPt->data,
				   v);

    if(curDist < minDist) {
      minDist = curDist;
    }
  }

  return minDist;
}

/**
 * get the distance from a point to the nearest point on a line segment
 */
double getDistanceToSegment(vertex *p1, vertex *p2, vertex *v) {
  double dv1, /* length from p2 to p1  */
         dv2, /* length from v to p1 */
         dv3; /* length from v to p2 */

  double theta; /* angle between the segments p2,p1 and p1,v */
  double r;     /* length down the segment for the perpenticular */
  double d;     /* min dist */

  /* compute neccesary vals */
  dv1 = dist((*p2),(*p1));
  dv2 = dist((*v),(*p1));
  dv3 = dist((*v),(*p2));

  theta = angleBetweenSegments(p2,p1,v);
  d = dv2*sin(theta);
  r = sqrt(dv2*dv2-d*d)/dv1;

  if(theta > SR_PI/2.0 || r > 1.0) {
    return min(dv2,dv3);
  }

  return d;
}

/**
 * get the median distance between two contours in feature space.  this
 * function computes the following: the distance between two contours is the
 * minimum of the sum of the meidan minimum distance of each contour vertex
 * to the nearest point on the other contour computed for each contour with
 * respect to the other.
 *
 * RETURNS LOG(DIST)
 */
double getContoursDist(contour *cont1, contour *cont2, double *mindist) {
  listNode *i;
  vertex *u;
  //double minDistSum = 0.0;
  double curDist;
  double *distArr;
  int ind;

  *mindist = SR_BIG;

  /* validate */
  if(cont1 == NULL || cont2 == NULL
     || listSize(cont1->vertices) == 0 || listSize(cont1->vertices) == 0) {
    return (*mindist = SR_BIG);
  }

  distArr = (double*) malloc(listSize(cont1->vertices)*sizeof(double));

  /* iterate over vertices of the first contour, finding the minimum
     distance to the other contour
     */
  for(i = getListNode(cont1->vertices,0), ind = 0; i;
      i = (listNode*) i->next, ind++) {
    u = (vertex*) i->data;
    distArr[ind] = getDistanceToContour(cont2,u);
    if(distArr[ind] < *mindist) *mindist = distArr[ind];
    //minDistSum += curdist;
  }

  /* sort to find the median */
  gsl_sort(distArr,1,ind);

  if(ind%2 == 0) {
    curDist = distArr[ind/2];
  }
  else {
    curDist = 0.5*(distArr[ind/2] + distArr[ind/2+1]);
  }

  free(distArr);

  return log(curDist);
}

/**
 * clear all correspondence relations
 */
int clearAllCorrespondence(list *slices) {
  list *slice;
  listNode *i, *j;
  contour *c;

  /* validate */
  if(listSize(slices) == 0) {
    return SR_FAILURE;
  }

  for(i = getListNode(slices,0); i; i = (listNode*) i->next) {
    slice = (list*) i->data;
    for(j = getListNode(slice,0); j; j = (listNode*) j->next) {
      c = (contour*) j->data;
      freeList(c->adjacentContours);
      c->adjacentContours = newList(LIST);
    }
  }

  return SR_SUCCESS;
}

/**
 * free a created distance matrix
 */
int freeDistanceMatrix(double **dm, int m, int n) {
  int i;

  /* validate */
  if(dm == NULL) {
    return SR_FAILURE;
  }

  /* free each row */
  for(i = 0; i < m; i++) {
    free(dm[i]);
  }

  free(dm);

  return SR_SUCCESS;
}

/**
 * build data array from multiple distance matrices
 */
void buildDataArray(list *dms, list* ms, list *ns, double **data, int *numDP) {
  listNode *dmln, *mln, *nln;
  double **dm;
  int curDP = 0, m, n, indi;

  FILE *saveFP = NULL;
  if(saveContourDistances) {
    saveFP = fopen(contourDistancesFilename,"w");
  }

  /* count the number of datapoints */
  for(mln = getListNode(ms,0), nln = getListNode(ns,0); mln;
      mln = (listNode*) mln->next, nln = (listNode*) nln->next) {
    (*numDP) += (long)mln->data * (long)nln->data;
  }

  /* get the data array */
  (*data) = (double*) malloc(*numDP*sizeof(double));
  for(dmln = getListNode(dms,0),
	mln = getListNode(ms,0),
	nln = getListNode(ns,0);
      dmln;
      dmln = (listNode*) dmln->next,
	mln = (listNode*) mln->next,
	nln = (listNode*) nln->next) {
    dm = (double**) dmln->data;
    m  = (long) mln->data;
    n  = (long) nln->data;

    for(indi = 0; indi < m*n; indi++) {
      (*data)[curDP] = dm[indi/n][indi%n];

      if(saveFP) {
	fprintf(saveFP,"%0.9f ", (*data)[curDP]);
      }

      curDP++;
    }
  }

  if(saveFP) {
    fclose(saveFP);
  }
}

/**
 * mask a data array, discarding data outside of the given percentiles
 */
void maskPercentiles(double *data, int *numDatapoints,
		     double minPercentile, double maxPercentile) {
  int i, minInd, maxInd, beginSkip = 0;

  /* remove nans */
  for(i = 0; i < (*numDatapoints); i++) {
    if(isnan(data[i])) {
      data[i] = SR_BIG;
    }
    if(isinf(data[i]) == -1) {
      beginSkip++;
    }
  }

  /* sort the data */
  gsl_sort(data,1,*numDatapoints);
  minInd = rint(minPercentile/100.0*(*numDatapoints)) + beginSkip;
  maxInd = rint(maxPercentile/100.0*(*numDatapoints));

  /* shift the data */
  (*numDatapoints) = maxInd-minInd;
  for(i = 0; i < (*numDatapoints); i++) {
    data[i] = data[i+minInd];
  }

  if(*numDatapoints < 1) {
    (*numDatapoints) = 0;
  }
}

/**
 * build a histogram of multiple matrices entries
 */
histogram *buildHistogramMultiDM(double *data, int numDP, int numBins) {
  histogram *hist;
  int i;

  /*validate */
  if(data == NULL || numDP < 1 || numBins < 1) {
    return NULL;
  }

  hist = (histogram*) malloc(sizeof(histogram));
  if(hist == NULL) return NULL;

  hist->numBins = numBins;
  hist->hist = (int*) malloc(hist->numBins*sizeof(int));
  if(hist->hist == NULL) return NULL;

  /* zero the bins */
  for(i = 0; i < hist->numBins; i++) {
    hist->hist[i] = 0;
  }

  /* find the minimum and maximum values */
  hist->maxval = hist->minval = data[0];
  for(i = 0; i < numDP; i++) {
    if(data[i] < hist->minval) {
      hist->minval = data[i];
    }
    if(data[i] > hist->maxval) {
      hist->maxval = data[i];
    }
  }

  /* get the bin size */
  hist->binSize = (hist->maxval-hist->minval)/hist->numBins;

  /* accumulate */
  for(i = 0; i < numDP; i++) {
    hist->hist[(int)floor((data[i]-hist->minval)/hist->binSize)]++;
  }

  return hist;
}

/**
 * build a histogram of a matrix's entries
 */
histogram *buildHistogramDM(double **A, int m, int n, int numBins) {
  histogram *hist;
  int i,j;

  /*validate */
  if(A == NULL || n < 1 || m < 1 || numBins < 1) {
    return NULL;
  }

  hist = (histogram*) malloc(sizeof(histogram));
  if(hist == NULL) return NULL;

  hist->numBins = numBins;
  hist->hist = (int*) malloc(hist->numBins*sizeof(int));
  if(hist->hist == NULL) return NULL;

  /* zero the bins */
  for(i = 0; i < hist->numBins; i++) {
    hist->hist[i] = 0;
  }

  /* find the minimum and maximum values */
  hist->maxval = hist->minval = A[0][0];
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      if(A[i][j] < hist->minval) {
	hist->minval = A[i][j];
      }
      if(A[i][j] > hist->maxval) {
	hist->maxval = A[i][j];
      }
    }
  }

  /* get the bin size */
  hist->binSize = (hist->maxval-hist->minval)/hist->numBins;

  /* accumulate */
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      hist->hist[(int)floor((A[i][j]-hist->minval)/hist->binSize)]++;
    }
  }

  return hist;
}

/**
 * frees a histogram
 */
int deleteHistogram(histogram *hist) {
  if(hist == NULL || hist->numBins < 1){
    return SR_FAILURE;
  }

  free(hist->hist);
  free(hist);

  return SR_SUCCESS;
}

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/correspondence.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
