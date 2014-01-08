/******************************************************************************
 * clustering.c is the source file for functions to perform cluster analysis
 *
 * Oliver Hinds <oph@bu.edu> 2005-12-12
 *
 *
 *
 *****************************************************************************/

#include"clustering.h"

/**
 * perform a k means cluster analysis
 */
int kmeans(double **data, int numDatapoints, int dataDimensions, int numMeans,
	   int **clusters, double ***means) {
  int i;
  int numFlips = 1;

  /* validate */
  if(data == NULL || numDatapoints == 0 || dataDimensions == 0
     || numMeans == 0 || numDatapoints < numMeans) {
    return SR_FAILURE;
  }

  /* set up the cluster assignments and means */
  (*clusters) = (int*) malloc(numDatapoints*sizeof(int));

  (*means) = (double**) malloc(numMeans*sizeof(double*));
  for(i = 0; i < numMeans; i++) {
    (*means)[i] = (double*) malloc(dataDimensions*sizeof(double));
  }

  /* intialize clusters randomly */
  srand(time(NULL));
  for(i = 0; i < numDatapoints; i++) {
    (*clusters)[i] = rand()/(LONG_MAX/numMeans);
  }

  /* run the iterative mean finding/cluster assignment */
  i = 0;
  while(numFlips > 0) {
    updateMeans(data,numDatapoints,dataDimensions,
		(*clusters),numMeans,(*means));
    assignClustersByClosestMean(data,numDatapoints,dataDimensions,numMeans,
				(*means),(*clusters),&numFlips);

  }

  return SR_SUCCESS;
}

/**
 * assign cluster numbers to data points via closest mean calculation
 */
int assignClustersByClosestMean(double **data, int numDatapoints,
				int dataDimensions, int numMeans,
				double **means, int *clusters, int *flips) {
  int i,el;

  *flips = 0;

  /* validate */
  if(data == NULL || numDatapoints == 0 || dataDimensions == 0
     || numMeans == 0 || numDatapoints < numMeans || means == NULL) {
    return SR_FAILURE;
  }

  /* iterate over datapoints, assigning each */
  for(i = 0; i < numDatapoints; i++) {
    el = findClosestArrayElement(means,numMeans,data[i],dataDimensions);

    /* check for flips */
    if(el != clusters[i]) {
      (*flips)++;
      clusters[i] = el;
    }
  }

  return SR_SUCCESS;
}

/**
 * update the locations of n mean clusters
 */
int updateMeans(double **data, int numDatapoints, int dataDimensions,
		int *clusters, int numMeans, double **means) {
  int i,j;
  double **sums;
  int *counts;

  /* validate */
  if(data == NULL || numDatapoints < 1 || dataDimensions < 1 
     || clusters == NULL || numMeans < 1 || means == NULL) {
    return SR_FAILURE;
  }

  /* intialize sums and counts */
  sums = (double**) malloc(numMeans*sizeof(double*));
  counts = (int*) malloc(numMeans*sizeof(int));
  for(i = 0; i < numMeans; i++) {
    sums[i] = (double*) malloc(dataDimensions*sizeof(double));
    counts[i] = 0;
    for(j = 0; j < dataDimensions; j++) {
      sums[i][j] = 0.0;
    }
  }

  /* iterate over datapoints, accumulating sums for each cluster */
  for(i = 0; i < numDatapoints; i++) {
    vectorSum(sums[clusters[i]],data[i],dataDimensions,sums[clusters[i]]);
    counts[clusters[i]]++;
  }

  /* make the sums means */
  for(i = 0; i < numMeans; i++) {
    for(j = 0; j < dataDimensions; j++) {
      means[i][j] = sums[i][j]/counts[i];
    }
  }

  /* free the sums and counts */
  /* make the sums means */
  for(i = 0; i < numMeans; i++) {
    free(sums[i]);
  }
  free(sums);
  free(counts);

  return SR_SUCCESS;
}


/**
 * finds the closest array element to a given element.  assuming that the
 * second dimension of the first array is the same as the same first
 * dimension of the second
 * uses euclidean distance
 */
int findClosestArrayElement(double **arr, int m, double *el, int n) {
  int i, minind;
  double curdist,mindist;

  /* validate */
  if(arr == NULL || m < 1 || el == NULL || n < 1) {
    return -1;
  }

  /* search through arr for closest to el */
  minind = 0;
  mindist = nd_dist(arr[0],el,n);
  for(i = 1; i < m; i++) {
    if((curdist = nd_dist(arr[i],el,n)) < mindist) {
      mindist = curdist;
      minind = i;
    }
  }

  return minind;
}

/**
 * finds the distance between two n-dimensional vectors
 */
double nd_dist(double *a, double *b, int n) {
  int i;
  double sum = 0.0;

  /* validate */
  if(a == NULL || b == NULL || n < 1) {
    return 0.0;
  }

  /* sum square differences */
  for(i = 0; i < n; i++) {
    sum += pow(a[i]-b[i],2);
  }

  return sqrt(sum);
}

/**
 * vector sum
 */
int vectorSum(double *a, double *b, int n, double *c) {
  int i;

  /* validate */
  if(a == NULL || b == NULL || n < 1 || c == NULL) {
    return SR_FAILURE;
  }

  /* iterate over dims */
  for(i = 0; i < n; i++) {
    c[i] = a[i]+b[i];
  }

  return SR_SUCCESS;
}

/**
 * factorial
 */
//int factorial(int a) {
//  if(a < 2) return 1;
//  return a*factorial(a-1);
//}

/************************************************************************
*** $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/clustering.c,v $
*** Local Variables:
*** mode: c
*** fill-column: 76
*** End:
************************************************************************/
