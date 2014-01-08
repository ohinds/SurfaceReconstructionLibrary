/*****************************************************************************
 * correspondence.h is the header file with functions to take a guess at 
 * a solution to the correspondence problem for the brain surface 
 * reconstruction problem for libsr. 
 * Oliver Hinds <oph@bu.edu> 2005-10-18
 *
 *
 *
 *****************************************************************************/

#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H

#define CORRESPONDENCE_VERSION_H "$Id: correspondence.h,v 1.15 2007/05/22 19:18:11 oph Exp $"

#include"gsl/gsl_sort.h"

#include"libsrTypes.h"
#include"libsrUtil.h"
#include"list.h"
#include"surfUtil.h"
#include"branch.h"
#include"clustering.h"

#include"correspondence.extern"

/**
 * builds a guess at contour correspondence
 */
int buildCorrespondenceGuess(list *slice);

/**
 * builds a guess at contour correspondence slice by slice
 */
int buildCorrespondenceGuessLocal(list *slices);

/**
 * builds a guess at contour correspondence by examining all slices
 * simultaneuosly
 */
int buildCorrespondenceGuessGlobal(list *slices);

/**
 * set adjacent contours by a correspondence threshold
 */
void setAdjacentContours(list *slices, double correspondenceThreshold,
			 list* dms, list* ms, list* ns);

/**
 * set the adjacency manually
 */
int setAdjacentContoursManual(list *slices, double connectionThreshold);

/**
 * get contour distance matrix
 */
double **getContourDistanceMatrix(list *slice1, list *slice2, long *m, long *n);

/**
 * get contour distance matrix for local correspondence
 */
double **getContourDistanceMatrixLocal(list *slice1, list *slice2, long *m, long *n);

/**
 * get a correspondence threshold via estimating the pdf of the contour
 * distances, then finding the first mode of the pdf.
 */
double getGlobalConnectionThresholdDensity(list *dms, list *ms, list *ns);

/**
 * connect contours based on histogramming minimum distance to each
 * other. if there are multiple low distance contours, connect them all.
 */
int correspondContoursHistogram(list *slice1, list* slice2, double **dm, 
				int m, int n);

/**
 * connect contours based on kmeans clustering of the log minimum distance
 * to each other.
 */
int correspondContoursKmeans(list *slice1, list* slice2, double **dm, 
			     int m, int n);

/**
 * determine the optimal connection threshold for contours via kmeans
 */
double getGlobalConnectionThresholdKmeans(list *dms, list *ms, list *ns);

/**
 * determine the optimal connection threshold for contours via histogram
 */
double getGlobalConnectionThresholdHistogram(list *dms, list *ms, list *ns);

/**
 * gets distance to the nearest vertex on a contour from a point
 */
double getDistanceToNearestVertex(contour *c, vertex *v);

/**
 * gets the minimum distance to the contour from a vertex
 */
double getDistanceToContour(contour *c, vertex *v);

/** 
 * get the distance from a point to the nearest point on a line segment
 */
double getDistanceToSegment(vertex *p1, vertex *p2, vertex *v);

/**
 * get the distance between two contours in feature space
 */
double getContoursDist(contour *cont1, contour *cont2, double *mindist);

/**
 * clear all correspondence relations
 */
int clearAllCorrespondence(list *slices);

/**
 * free a created distance matrix
 */
int freeDistanceMatrix(double **dm, int m, int n);

/**
 * build data array from multiple distance matrices
 */
void buildDataArray(list *dms, list* ms, list *ns, double **data, int *numDP);

/**
 * mask a data array, discarding data outside of the given percentiles
 */
void maskPercentiles(double *data, int *numDatapoints,
		     double minPercentile, double maxPercentile);

/**
 * build a histogram of multiple matrices entries 
 */
histogram *buildHistogramMultiDM(double *data, int numDP, int numBins);

/**
 * build a histogram of a matrix's entries 
 */
histogram *buildHistogramDM(double **A, int m, int n, int numBins);

/**
 * frees a histogram 
 */
int deleteHistogram(histogram *hist);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/correspondence.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
