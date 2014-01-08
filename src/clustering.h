/******************************************************************************
 * clustering.h is the header file for functions to perform cluster analysis
 *
 * Oliver Hinds <oph@bu.edu> 2005-12-12
 *
 *
 *
 *****************************************************************************/

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits.h>
#include<time.h>
#include<string.h>
#include <gsl/gsl_combination.h>

#include"libsrTypes.h"

/**
 * perform a k means cluster analysis
 */
int kmeans(double **data, int numDatapoints, int dataDimension, int numMeans,
           int **clusters, double ***means);

/**
 * assign cluster numbers to data points via closest mean calculation
 */
int assignClustersByClosestMean(double **data, int numDatapoints,
                                int dataDimensions, int numMeans,
                                double **means, int *clusters, int *flips);

/**
 * update the locations of n mean clusters
 */
int updateMeans(double **data, int numDataPoints, int dataDimensions,
                int *clusters, int numMeans, double **means);

/**
 * finds the closest array element to a given element.  assuming that the
 * second dimension of the first array is the same as the same first
 * dimension of the second
 * uses euclidean distance
 */
int findClosestArrayElement(double **arr, int m, double *el, int n);

/**
 * finds the distance between two n-dimensional vectors
 */
double nd_dist(double *a, double *b, int n);

/**
 * vector sum
 */
int vectorSum(double *a, double *b, int n, double *c);

/**
 * factorial
 */
int factorial(int a);

#endif

/************************************************************************
 *** $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/clustering.h,v $
 *** Local Variables:
 *** mode: c
 *** fill-column: 76
 *** End:
 ************************************************************************/
