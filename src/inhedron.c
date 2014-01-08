/*
  This code is described in "Computational Geometry in C" (Second Edition),
  Chapter 7.  It is not written to be comprehensible without the
  explanation in that book.

  Compile:    gcc -o inhedron inhedron.c -lm
  Run (e.g.): inhedron < i.8

  Written by Joseph O'Rourke, with contributions by Min Xu.
  Last modified: April 1998
  Questions to orourke@cs.smith.edu.
  --------------------------------------------------------------------
  This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
  redistributed in its entirety provided that this copyright notice is
  not removed.
  --------------------------------------------------------------------
  *
  * functions for testing inside/outside and some simple intersections.
  * Taken largely from code by Joseph O'Rourke, with contributions by Min Xu.
  * The original code is available at
  * http://orion.math.iastate.edu/burkardt/c_src/orourke/
  * slight modifications by Oliver Hinds <oph@bu.edu> but i didnt make
  * anything better, believe me.
  */
#include"inhedron.h"
#define PMAX 10000             /* Max # of pts */
tPointd Vertices[PMAX];        /* All the points */
tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
int check = 0;
tPointi Box[PMAX][2];          /* Box around each face */


#include"libsrTypes.h"

int ComputeBox( int F, tPointd bmin, tPointd bmax )
{
  int i, j;
  double radius;

  for( i = 0; i < F; i++ )
    for( j = 0; j < DIM; j++ ) {
      if( Vertices[i][j] < bmin[j] )
        bmin[j] = Vertices[i][j];
      if( Vertices[i][j] > bmax[j] )
        bmax[j] = Vertices[i][j];
    }

  radius = sqrt( pow( (double)(bmax[X] - bmin[X]), 2.0 ) +
                 pow( (double)(bmax[Y] - bmin[Y]), 2.0 ) +
                 pow( (double)(bmax[Z] - bmin[Z]), 2.0 ) );

  return rint( radius +1 ) + 1;
}


void AddVec( tPointd q, tPointd ray )
{
  int i;

  for( i = 0; i < DIM; i++ )
    ray[i] = q[i] + ray[i];
}

int InBox( tPointd q, tPointd bmin, tPointd bmax )
{
  if( ( bmin[X] <= q[X] ) && ( q[X] <= bmax[X] ) &&
      ( bmin[Y] <= q[Y] ) && ( q[Y] <= bmax[Y] ) &&
      ( bmin[Z] <= q[Z] ) && ( q[Z] <= bmax[Z] ) )
    return TRUE;
  return FALSE;
}


/*---------------------------------------------------------------------
  'p': The segment lies wholly within the plane.
  'q': The q endpoint is on the plane (but not 'p').
  'r': The r endpoint is on the plane (but not 'p').
  '0': The segment lies strictly to one side or the other of the plane.
  '1': The segement intersects the plane, and 'p' does not hold.
  ---------------------------------------------------------------------*/
char	SegPlaneInt( tPointi T, tPointd q, tPointd r, tPointd p, int *m)
{
  tPointd N; double d;
  tPointd rq;
  double num, denom, t;
  int i;

  *m = PlaneCoeff( T, N, &d );
  /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],d);*/
  num = d - Dot( q, N );
  SubVec( r, q, rq );
  denom = Dot( rq, N );
  /*printf("SegPlaneInt: num=%lf, denom=%lf\n", num, denom );*/

  if ( denom == 0.0 ) {  /* Segment is parallel to plane. */
    if ( num == 0.0 )   /* q is on plane. */
      return 'p';
    else
      return '0';
  }
  else
    t = num / denom;
  /*printf("SegPlaneInt: t=%lf \n", t );*/

  for( i = 0; i < DIM; i++ )
    p[i] = q[i] + t * ( r[i] - q[i] );

  if ( (0.0 < t) && (t < 1.0) )
    return '1';
  else if ( num == 0.0 )   /* t == 0 */
    return 'q';
  else if ( num == denom ) /* t == 1 */
    return 'r';
  else return '0';
}
/*---------------------------------------------------------------------
  Computes N & D and returns index m of largest component.
  ---------------------------------------------------------------------*/
int	PlaneCoeff( tPointi T, tPointd N, double *d )
{
  int i;
  double t;              /* Temp storage */
  double biggest = 0.0;  /* Largest component of normal vector. */
  int m = 0;             /* Index of largest component. */

  NormalVec( Vertices[T[0]], Vertices[T[1]], Vertices[T[2]], N );
  /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
  *d = Dot( Vertices[T[0]], N );

  /* Find the largest component of N. */
  for ( i = 0; i < DIM; i++ ) {
    t = fabs( N[i] );
    if ( t > biggest ) {
      biggest = t;
      m = i;
    }
  }
  return m;
}
/*---------------------------------------------------------------------
  Reads in the number and coordinates of the vertices of a polyhedron
  from stdin, and returns n, the number of vertices.
  ---------------------------------------------------------------------*/
int     ReadVertices( void )
{
  int   i, n;

  do {
    scanf( "%d", &n );
    if ( n <= PMAX )
      break;
    printf("Error in read_vertex:  too many points; max is %d\n", PMAX);
  }
  while ( 1 );

  for ( i = 0; i < n; i++ ) {
    scanf( "%lf %lf %lf", &Vertices[i][X], &Vertices[i][Y], &Vertices[i][Z] );
  }
  printf("n = %3d vertices read\n",n);
  putchar('\n');

  return n;
}

/*---------------------------------------------------------------------
  a - b ==> c.
  ---------------------------------------------------------------------*/
void    SubVec( tPointd a, tPointd b, tPointd c )
{
  int i;

  for( i = 0; i < DIM; i++ )
    c[i] = a[i] - b[i];
}

/*---------------------------------------------------------------------
  Returns the dot product of the two input vectors.
  ---------------------------------------------------------------------*/
double	Dot( tPointd a, tPointd b )
{
  int i;
  double sum = 0.0;

  for( i = 0; i < DIM; i++ )
    sum += a[i] * b[i];

  return  sum;
}

/*---------------------------------------------------------------------
  Compute the cross product of (b-a)x(c-a) and place into N.
  ---------------------------------------------------------------------*/
void	NormalVec( tPointd a, tPointd b, tPointd c, tPointd N )
{
  N[X] = ( c[Z] - a[Z] ) * ( b[Y] - a[Y] ) -
      ( b[Z] - a[Z] ) * ( c[Y] - a[Y] );
  N[Y] = ( b[Z] - a[Z] ) * ( c[X] - a[X] ) -
      ( b[X] - a[X] ) * ( c[Z] - a[Z] );
  N[Z] = ( b[X] - a[X] ) * ( c[Y] - a[Y] ) -
      ( b[Y] - a[Y] ) * ( c[X] - a[X] );
}

/* Reads in the number of faces of the polyhedron and their indices from stdin,
   and returns the number n. */
int ReadFaces( void )
{
  int	i,j,k, n;
  int   w; /* temp storage for coordinate. */

  do {
    scanf( "%d", &n );
    if ( n <= PMAX )
      break;
    printf("Error in read_vertex:  too many points; max is %d\n", PMAX);
  }
  while ( 1 );

  for ( i = 0; i < n; i++ ) {
    scanf( "%d %d %d", &Faces[i][0], &Faces[i][1], &Faces[i][2] );
    /* Compute bounding box. */
    /* Initialize to first vertex. */
    for ( j=0; j < 3; j++ ) {
      Box[i][0][j] = Vertices[ Faces[i][0] ][j];
      Box[i][1][j] = Vertices[ Faces[i][0] ][j];
    }
    /* Check k=1,2 vertices of face. */
    for ( k=1; k < 3; k++ )
      for ( j=0; j < 3; j++ ) {
        w = Vertices[ Faces[i][k] ][j];
        if ( w < Box[i][0][j] ) Box[i][0][j] = w;
        if ( w > Box[i][1][j] ) Box[i][1][j] = w;
      }
    /* printf("Bounding box: (%d,%d,%d);(%d,%d,%d)\n",
       Box[i][0][0],
       Box[i][0][1],
       Box[i][0][2],
       Box[i][1][0],
       Box[i][1][1],
       Box[i][1][2] );
    */
  }
  printf("n = %3d faces read\n",n);
  putchar('\n');

  return n;
}

/* Assumption: p lies in the plane containing T.
   Returns a char:
   'V': the query point p coincides with a Vertex of triangle T.
   'E': the query point p is in the relative interior of an Edge of triangle T.
   'F': the query point p is in the relative interior of a Face of triangle T.
   '0': the query point p does not intersect (misses) triangle T.
*/

char    InTri3D( tPointi T, int m, tPointd p )
{
  int i;           /* Index for X,Y,Z           */
  int j;           /* Index for X,Y             */
  int k;           /* Index for triangle vertex */
  tPointd pp;      /* projected p */
  tPointd Tp[3];   /* projected T: three new vertices */

  /* Project out coordinate m in both p and the triangular face */
  j = 0;
  for ( i = 0; i < DIM; i++ ) {
    if ( i != m ) {    /* skip largest coordinate */
      pp[j] = p[i];
      for ( k = 0; k < 3; k++ )
        Tp[k][j] = Vertices[T[k]][i];
      j++;
    }
  }
  return( InTri2D( Tp, pp ) );
}

char    InTri2D( tPointd Tp[3], tPointd pp )
{
  int area0, area1, area2;

  /* compute three AreaSign() values for pp w.r.t. each edge of the face in 2D */
  area0 = AreaSign( pp, Tp[0], Tp[1] );
  area1 = AreaSign( pp, Tp[1], Tp[2] );
  area2 = AreaSign( pp, Tp[2], Tp[0] );

  if ( (( area0 == 0 ) && ( area1 > 0 ) && ( area2 > 0 )) ||
       (( area1 == 0 ) && ( area0 > 0 ) && ( area2 > 0 )) ||
       (( area2 == 0 ) && ( area0 > 0 ) && ( area1 > 0 )) )
    return 'E';

  if ( (( area0 == 0 ) && ( area1 < 0 ) && ( area2 < 0 )) ||
       (( area1 == 0 ) && ( area0 < 0 ) && ( area2 < 0 )) ||
       (( area2 == 0 ) && ( area0 < 0 ) && ( area1 < 0 )) )
    return 'E';

  if ( (( area0 >  0 ) && ( area1 > 0 ) && ( area2 > 0 )) ||
       (( area0 <  0 ) && ( area1 < 0 ) && ( area2 < 0 )) )
    return 'F';

  if ( ( area0 == 0 ) && ( area1 == 0 ) && ( area2 == 0 ) )
    return '0';

  if ( (( area0 == 0 ) && ( area1 == 0 )) ||
       (( area0 == 0 ) && ( area2 == 0 )) ||
       (( area1 == 0 ) && ( area2 == 0 )) )
    return 'V';

  else
    return '0';
}

int     AreaSign( tPointd a, tPointd b, tPointd c )
{
  double area2;

  area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
      ( c[0] - a[0] ) * (double)( b[1] - a[1] );

  /* The area should be an integer. */
  if      ( area2 >  SR_TOL ) return  1;
  else if ( area2 < -SR_TOL ) return -1;
  else                     return  0;
}

char    SegTriInt( tPointi T, tPointd q, tPointd r, tPointd p )
{
  int code = '?';
  int m = -1;

  code = SegPlaneInt( T, q, r, p, &m );


  if      ( code == '0')
    return '0';
  else if ( code == 'q')
    return InTri3D( T, m, q );
  else if ( code == 'r')
    return InTri3D( T, m, r );
  else if ( code == 'p' )
    return InPlane( T, m, q, r, p );
  else if ( code == '1' )
    return SegTriCross( T, q, r );
  else /* Error */
    return code;
}

// return whether the segment intersects a tri edge
// not complete! the c1 == '0' and c2 == '0' case is not handled!
char	InPlane( tPointi T, int m, tPointd q, tPointd r, tPointd p)
{
  tPointd tmp;
  int c1,c2;
  c1 = InTri3D(T,m,q);
  c2 = InTri3D(T,m,r);

  // yes
  if(c1 == 'E' && c2 == 'E') return 'f';
  if(c1 == 'E' && c2 == 'F') return 'f';
  if(c1 == 'F' && c2 == 'E') return 'f';
  if(c1 == 'F' && c2 == 'F') return 'f';
  if(c1 == '0' && c2 == 'F') return 'f';
  if(c1 == 'F' && c2 == '0') return 'f';

  // no
  if(c1 == 'V' && c2 == 'V') return '0';

  // maybe
  if(c1 == '0' && (c2 == 'V' || c2 == 'E')) {
    tmp[0] = r[0] + 0.001*(p[0]-r[0]);
    tmp[1] = r[1] + 0.001*(p[1]-r[1]);
    tmp[2] = r[2] + 0.001*(p[2]-r[2]);

    if(InTri3D(T,m,tmp) == 'F') return 'f';
    else return '0';
  }

  if((c1 == 'V' || c1 == 'E') && c2 == '0') {
    tmp[0] = p[0] + 0.001*(r[0]-p[0]);
    tmp[1] = p[1] + 0.001*(r[1]-p[1]);
    tmp[2] = p[2] + 0.001*(r[2]-p[2]);

    if(InTri3D(T,m,tmp) == 'F') return 'f';
    else return '0';
  }

  return '0';
}

/*---------------------------------------------------------------------
  The signed volumes of three tetrahedra are computed, determined
  by the segment qr, and each edge of the triangle.
  Returns a char:
  'v': the open segment includes a vertex of T.
  'e': the open segment includes a point in the relative interior of an edge
  of T.
  'f': the open segment includes a point in the relative interior of a face
  of T.
  '0': the open segment does not intersect triangle T.
  ---------------------------------------------------------------------*/

char SegTriCross( tPointi T, tPointd q, tPointd r )
{
  int vol0, vol1, vol2;

  vol0 = VolumeSign( q, Vertices[ T[0] ], Vertices[ T[1] ], r );
  vol1 = VolumeSign( q, Vertices[ T[1] ], Vertices[ T[2] ], r );
  vol2 = VolumeSign( q, Vertices[ T[2] ], Vertices[ T[0] ], r );


  /* Same sign: segment intersects interior of triangle. */
  if ( ( ( vol0 > 0 ) && ( vol1 > 0 ) && ( vol2 > 0 ) ) ||
       ( ( vol0 < 0 ) && ( vol1 < 0 ) && ( vol2 < 0 ) ) )
    return 'f';

  /* Opposite sign: no intersection between segment and triangle */
  if ( ( ( vol0 > 0 ) || ( vol1 > 0 ) || ( vol2 > 0 ) ) &&
       ( ( vol0 < 0 ) || ( vol1 < 0 ) || ( vol2 < 0 ) ) )
    return '0';

  else if ( ( vol0 == 0 ) && ( vol1 == 0 ) && ( vol2 == 0 ) )
    return '0';

  /* Two zeros: segment intersects vertex. */
  else if ( ( ( vol0 == 0 ) && ( vol1 == 0 ) ) ||
            ( ( vol0 == 0 ) && ( vol2 == 0 ) ) ||
            ( ( vol1 == 0 ) && ( vol2 == 0 ) ) )
    return 'v';

  /* One zero: segment intersects edge. */
  else if ( ( vol0 == 0 ) || ( vol1 == 0 ) || ( vol2 == 0 ) )
    return 'e';

  else
    fprintf( stderr, "Error 2 in SegTriCross\n" );

  return '0';
}

int     VolumeSign( tPointd a, tPointd b, tPointd c, tPointd d )
{
  double vol, tol = 0.00000001;
  double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
  double bxdx, bydy, bzdz, cxdx, cydy, czdz;

  ax = a[X];
  ay = a[Y];
  az = a[Z];
  bx = b[X];
  by = b[Y];
  bz = b[Z];
  cx = c[X];
  cy = c[Y];
  cz = c[Z];
  dx = d[X];
  dy = d[Y];
  dz = d[Z];

  bxdx=bx-dx;
  bydy=by-dy;
  bzdz=bz-dz;
  cxdx=cx-dx;
  cydy=cy-dy;
  czdz=cz-dz;
  vol =   (az-dz) * (bxdx*cydy - bydy*cxdx)
      + (ay-dy) * (bzdz*cxdx - bxdx*czdz)
      + (ax-dx) * (bydy*czdz - bzdz*cydy);


  /* The volume should be an integer. */
  if(fabs(vol) < tol) return 0;
  else if      ( vol > 0 )   return  1;
  else                    return  -1;

}

/*
  This function returns a char:
  '0': the segment [ab] does not intersect (completely misses) the
  bounding box surrounding the n-th triangle T.  It lies
  strictly to one side of one of the six supporting planes.
  '?': status unknown: the segment may or may not intersect T.
*/
char BoxTest ( int n, tPointd a, tPointd b )
{
  int i; /* Coordinate index */
  int w;

  for ( i=0; i < DIM; i++ ) {
    w = Box[ n ][0][i]; /* min: lower left */
    if ( (a[i] < w) && (b[i] < w) ) return '0';
    w = Box[ n ][1][i]; /* max: upper right */
    if ( (a[i] > w) && (b[i] > w) ) return '0';
  }
  return '?';
}
