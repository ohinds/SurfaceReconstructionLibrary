/******************************************************************************
 * functions for testing inside/outside and some simple intersections.
 * Taken largely from code by Joseph O'Rourke, with contributions by Min Xu.
 * The original code is available at
 * http://orion.math.iastate.edu/burkardt/c_src/orourke/
 * slight modifications by Oliver Hinds <oph@bu.edu> but i didnt make
 * anything better, believe me.
 *
 * see inhedron.c for more information
 *****************************************************************************/

#ifndef INHEDRON_H
#define INHEDRON_H
#include <stdio.h>
#include <math.h>
#include"libsrTypes.h"
#define EXIT_FAILURE 1
#define X 0
#define Y 1
#define Z 2
#define MAX_INT   2147483647
#define DIM 3                  /* Dimension of points */
typedef int tPointi[DIM];   /* Type integer point */
typedef double tPointd[DIM];   /* Type double point */
extern tPointd Vertices[];

/*---------------------------------------------------------------------
  Function prototypes.
  ---------------------------------------------------------------------*/
char    InPolyhedron( int F, tPointd q, tPointd bmin, tPointd bmax, int radius );
char    SegPlaneInt( tPointi Triangle, tPointd q, tPointd r, tPointd p, int *m )
    ;
int     PlaneCoeff( tPointi T, tPointd N, double *d );
void    Assigndi( tPointd p, tPointi a );
int     ReadVertices( void );
int     ReadFaces( void );
void    NormalVec( tPointd q, tPointd b, tPointd c, tPointd N );
double  Dot( tPointd q, tPointd d );
void    SubVec( tPointd q, tPointd b, tPointd c );
char    InTri3D( tPointi T, int m, tPointd p );
char    InTri2D( tPointd Tp[3], tPointd pp );
int     AreaSign( tPointd q, tPointd b, tPointd c );
char    SegTriInt( tPointi Triangle, tPointd q, tPointd r, tPointd p );
char    InPlane( tPointi Triangle, int m, tPointd q, tPointd r, tPointd p);
int     VolumeSign( tPointd a, tPointd b, tPointd c, tPointd d );
char    SegTriCross( tPointi Triangle, tPointd q, tPointd r );
int     ComputeBox( int F, tPointd bmin, tPointd bmax );
void    RandomRay( tPointd ray, int radius );
void    AddVec( tPointd q, tPointd ray );
int     InBox( tPointd q, tPointd bmin, tPointd bmax );
char    BoxTest ( int n, tPointd a, tPointd b );
void    PrintPoint( tPointd q );
/*-------------------------------------------------------------------*/

#endif
