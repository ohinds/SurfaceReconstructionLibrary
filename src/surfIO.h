/*****************************************************************************
 * surfIO.h is the header file for surface input/output functions for
 * libsr 
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 * 
 *
 *****************************************************************************/

#ifndef SURF_IO_H
#define SURF_IO_H

#include<time.h>
#include"libsrTypes.h"
#include"libsrUtil.h"
#include"list.h"
#include"ioUtil.h"
#include"surfUtil.h"


/**
 * reads a surface as a given format
 */
surface *readSurface(char *filename, enum SURFACE_FORMAT format);

/**
 * writes a surface as a given format
 */
int writeSurface(surface *surf, char *filename, enum SURFACE_FORMAT format);

/**
 * reads a surface from an off file 
 * returns a surface structure or null if reading flibvps
 */
surface *readOFF(char* filename);

/**
 * writes a surface in the off file formate
 * returns SUCCESS or FAILURE
 */
int writeOFF(surface *surf, char *filename);

/**
 * reads a surface from an obj file 
 * returns a surface structure or null if reading flibvps
 */
surface *readOBJ(char* filename);

/**
 * writes a surface in the obj file formate
 * returns SUCCESS or FAILURE
 */
int writeOBJ(surface *surf, char *filename);

/**
 * reads a surface in a freesurfer triangle surface binary format
 * returns a surface structure or null if reading fails
 */
surface *readMGHsurf(char* filename);

/**
 * write a surface in a freesurfer triangle surface binary format
 */
int writeMGHsurf(surface *surf, char *filename);

/** my label and curvature formats **/

/**
 * read a label file
 */
int readLabelFile(surface *surf, char *filename);

/**
 * write a label file
 */
int writeLabelFile(surface *surf, char *filename);

/**
 * read a curvature file
 */
int readCurvatureFile(surface *surf, char *filename);

/**
 * write a curvature file
 */
int writeCurvatureFile(surface *surf, char *filename);


/** freesurfer label and curvature formats **/

/**
 * read a freesurfer label file
 */
int readMGHLabelFile(surface *surf, char *filename);

/**
 * write a freesurfer label file
 */
int writeMGHLabelFile(surface *surf, char *filename);

/**
 * read a freesurfer curvature file
 */
int readMGHCurvatureFile(surface *surf, char *filename);

/**
 * write a freesurfer curvature file
 */
int writeMGHCurvatureFile(surface *surf, char *filename);


/** slice contour io */

/**
 * read a set of slice contours from a file
 */
list *readSliceContours(char *filename);

/**
 * write a set of slice contours to a file
 */
int writeSliceContours(list *slices, char *filename);

/**
 * delete a list of slice contours
 */
void deleteSliceContours(list *slices);

/**
 * write a set of slice contour labels to a file
 */
int writeSliceContourLabels(list *slices, char *filename);

/**
 * read a set of slice contour labels from a file
 */
int readSliceContourLabels(list *slices, char *filename);

/**
 * reads a nuages file into a set of slices
 */
list *readNuagesSlices(char *filename);

/**
 * looks for a valid surface format name in a string,
 * returns both the format 
 */
int getSurfFormatFromString(char *str, enum SURFACE_FORMAT *format);

#endif

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/surfIO.h,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
