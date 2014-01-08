/*****************************************************************************
 * branch.h is the header file for an implementation of a solution to the
 * branching problem for libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#ifndef BRANCH_H
#define BRANCH_H

#include"libsrTypes.h"
#include"libsrUtil.h"
#include"list.h"
#include"gsl/gsl_permutation.h"
#include"gsl/gsl_combination.h"

#include"branch.extern"

/**
 * take two sets of contours for adjacent slices that already have
 * correspondences identified. For each contour one slice one that is
 * connected to multiple contours on slice two, break it up.
 *
 * 1) identify possible branch locations by finding ends of contours that
 * aren't aligned with ends of corresponding contours.
 *
 * 2) create a new contour pair from those contours.
 * 3) repeat until all correspodences are one-to-one
 */
list *branchSlices(list *slice1, list *slice2);

/**
 * split contours via maximum likelihood estimation of boundaries
 */
void splitContoursMLE(contour *cont, list *adjacent,
                      contour ***closestContours);

/**
 * create contour pairs and add them to a list from a set of closest contours
 */
void addContourPairsFromClosest(contour *curCont, contour **closestContours,
                                list *contourPairs);

/**
 * take in a list of contour pairs that are connected to the split
 * contour and split the target contour into many smaller ones, each
 * of which corresponds to one of the contours in the list of contour
 * pairs passed.  returns a vector of contour assignments for each
 * vertex in the targetcontour
 *
 * NOTE! POTENTIAL ERROR: CLOSED CONTOURS SHOULD BE HANDLED DIFFERENTLY!
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
                              int *contourAssignments);

/**
 * take in a list of contour pairs that all have the same c2Origin and
 * split the target contour into many smaller ones, each of which
 * corresponds to one of the c1 contours in the list of contour pairs
 * passed
 */
void splitContoursManyToOne(list *contourPairs);

/**
 * take in a list of contour pairs that all have the same c1Origin and
 * split the target contour into many smaller ones, each of which
 * corresponds to one of the c2 contours in the list of contour pairs
 * passed
 */
void splitContoursOneToMany(list *contourPairs);

/**
 * use an assignment vector to create actual contour pairs
 * corresponding to a branch where all original c2 contours were the
 * same. so break up the c2 origin and create new contours
 */
void splitC2Origin(list *contourPairs, int *contourAssignments);

/**
 * use a contour assignment vector to create actual contour pairs
 * corresponding to a branch where all original c2 contours were the
 * same. so break up the c1 origin and create new contours
 */
void splitC1Origin(list *contourPairs, int *contourAssignments);

/**
 * get the closest vertex in a to a contour to a given vertex
 */
vertex *getClosestVertexInContour(contour *cont, vertex *v);

/**
 * get the closest vertex in a to a contour to a given vertex, return index
 */
int getClosestVertexInContourInd(contour *cont, vertex *v);

/**
 * get the closest contour to a vertex
 */
contour *getClosestContour(vertex *v, list *contours);

/**
 * get the minimum distance from a vertex to a set of vertices
 */
float getMinDistance(vertex *v, list *vertices);

/**
 * evaluate a particular order and location of constant conoutr assignment
 */
int getNumMisclassifiedVertices(list *contours,
                                int numVertices,
                                contour **closestContours,
                                gsl_permutation *order,
                                list *possibleTransitions,
                                gsl_combination *transitions);

/**
 * searches two branched contours forvertices that were adjacent on the
 * original contour
 */
int matchBranchedContours(contour *c1, contour *c2,
                          int *c10, int *c11, int *c20, int *c21);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/branch.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
