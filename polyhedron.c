/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  

Compile:    gcc -o inhedron inhedron.c -lm (or simply: make)
Run (e.g.): inhedron < i.8

Written by Hsien-Yi Liu, refer to Joseph O'Rourke, with contributions by Min Xu.

--------------------------------------------------------------------
This code is Copyright 1998 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h> 
#define EXIT_FAILURE 1
#define X 0
#define Y 1
#define Z 2
#define MAX_INT   2147483647 
typedef enum { FALSE, TRUE } bool;

#define DIM 3                  /* Dimension of points */
typedef int    tPointi[DIM];   /* Type integer point */
typedef double tPointd[DIM];   /* Type double point */
#define PMAX 10000             /* Max # of pts */
tPointd Vertices[PMAX];        /* All the points */
tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
tPointd com_Vertices[PMAX];    
tPointi com_Faces[PMAX];         
int check = 0;
tPointi Box[PMAX][2];          /* Box around each face */
int n_facets, n_vertices;      /* Original polyhedron*/
int com_facets, com_vertices;  /* Original polyhedron*/

/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/
char  InPolyhedron( int F, tPointd q, tPointd bmin, tPointd bmax, int radius );
char  SegPlaneInt( tPointi Triangle, tPointd q, tPointd r, tPointd p, int *m );
int   PlaneCoeff( tPointi T, tPointd N, double *D );
void  Assigndi( tPointd p, tPointi a );
int   ReadVertices( void );
int   ReadFaces( void );
void  NormalVec( tPointd q, tPointd b, tPointd c, tPointd N );
double Dot( tPointd q, tPointd d );
void  SubVec( tPointd q, tPointd b, tPointd c );
char  InTri3D( tPointi T, int m, tPointd p );
char  InTri2D( tPointd Tp[3], tPointd pp );
int   AreaSign( tPointd q, tPointd b, tPointd c );
char  SegTriInt( tPointi Triangle, tPointd q, tPointd r, tPointd p );
char  InPlane( tPointi Triangle, int m, tPointd q, tPointd r, tPointd p);
int   VolumeSign( tPointd a, tPointd b, tPointd c, tPointd d );
char  SegTriCross( tPointi Triangle, tPointd q, tPointd r );
int   ComputeBox( int F, tPointd bmin, tPointd bmax );
void  RandomRay( tPointd ray, int radius );
void 	AddVec( tPointd q, tPointd ray );
int  	InBox( tPointd q, tPointd bmin, tPointd bmax );
char 	BoxTest ( int n, tPointd a, tPointd b );
void 	PrintPoint( tPointi q );
int	irint( double x);
void  read_ori(void);
void  read_com(void);
void  init_bounding(int n);
//double get_time_elapsed(struct timeval &t1);
/*-------------------------------------------------------------------*/
int main(){
    int n, F, i;
    tPointd q, bmin, bmax;
    int radius;

    //srandom( (int) time( (long *) 0 ) ); 
    read_ori();
    read_com();
    n = n_vertices;
    F = n_facets;
    init_bounding(n_facets);

    /* Initialize the bounding box */
    for ( i = 0; i < DIM; i++ ){
        bmin[i] = bmax[i] = Vertices[0][i];
        printf("bmin=%lf\n", Vertices[0][i]);
    }
    // bmin --> doublemax
    // bmax --> doublemin
    radius = ComputeBox( n, bmin, bmax );
    printf("radius=%d\n", radius);
    int counter = com_vertices - 1;
    time_t begin = time(NULL);
    for(int i = 0; i < n_vertices; i++){
       printf( "In = %fn %f %f \n", Vertices[i][X],Vertices[i][Y], Vertices[i][Z]);
    }
    /*
    while( counter >= 0 ) {
        q[X] = com_Vertices[counter][X];
        q[Y] = com_Vertices[counter][Y];
        q[Z] = com_Vertices[counter][Z];
        printf( "\n %d -------->q = %lf %lf %lf\n", counter, q[X], q[Y], q[Z] );
        printf( "In = %c\n", InPolyhedron( F, q, bmin, bmax, radius ) );
        counter--;
    }
    time_t end = time(NULL); 
    printf("Time elpased is %ld seconds", (end - begin));
    */
     while( scanf( "%lf %lf %lf", &q[X], &q[Y], &q[Z] ) != EOF ) {
         printf( "\n----------->q = %lf %lf %lf\n", 
            q[X], q[Y], q[Z] );
         printf( "In = %c\n", InPolyhedron( F, q, bmin, bmax, radius ) );
  }
}

/*
  This function returns a char:
    'V': the query point a coincides with a Vertex of polyhedron P.
    'E': the query point a is in the relative interior of an Edge of polyhedron P.
    'F': the query point a is in the relative interior of a Face of polyhedron P.
    'i': the query point a is strictly interior to polyhedron P.
    'o': the query point a is strictly exterior to( or outside of) polyhedron P.
*/
char InPolyhedron( int F, tPointd q, tPointd bmin, tPointd bmax, int radius ){
   tPointd r;  /* Ray endpoint. */
   tPointd p;  /* Intersection point; not used. */
   int f, k = 0, crossings = 0;
   char code = '?';
 
   /* If query point is outside bounding box, finished. */
   if ( !InBox( q, bmin, bmax ) )
      return 'o';
  
   LOOP:
   while( k++ < F ) {
      crossings = 0;
  
      RandomRay( r, radius ); 
      AddVec( q, r ); 
      printf("Ray endpoint: (%lf,%lf,%lf)\n", r[0],r[1],r[2] );
  
      for ( f = 0; f < F; f++ ) {  /* Begin check each face */
         if ( BoxTest( f, q, r ) == '0' ) {
              code = '0';
              //printf("BoxTest = 0!\n");
         }
         else code = SegTriInt( Faces[f], q, r, p );
         //printf( "Face = %d: BoxTest/SegTriInt returns %c\n\n", f, code );

         /* If ray is degenerate, then goto outer while to generate another. */
         if ( code == 'p' || code == 'v' || code == 'e' ) {
            printf("Degenerate ray\n");
            goto LOOP;
         }
   
         /* If ray hits face at interior point, increment crossings. */
         else if ( code == 'f' ) {
            crossings++;
            printf( "crossings = %d\n", crossings );
         }

         /* If query endpoint q sits on a V/E/F, return that code. */
         else if ( code == 'V' || code == 'E' || code == 'F' )
            return( code );

         /* If ray misses triangle, do nothing. */
         else if ( code == '0' )
            ;

         else 
            fprintf( stderr, "Error, exit(EXIT_FAILURE)\n" ), exit(1);

      } /* End check each face */

      /* No degeneracies encountered: ray is generic, so finished. */
      break;

   } /* End while loop */
 
   printf( "Crossings = %d\n", crossings );
   /* q strictly interior to polyhedron iff an odd number of crossings. */
   if( ( crossings % 2 ) == 1 )
      return   'i';
   else return 'o';
}

int ComputeBox( int F, tPointd bmin, tPointd bmax ){
    int i, j, k;
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
    printf("radius = %lf\n", radius);

    return irint( radius +1 ) + 1;
}

/* Return a random ray endpoint */
void RandomRay( tPointd ray, int radius ){
    double x, y, z, w, t;

    /* Generate a random point on a sphere of radius 1. */
    /* the sphere is sliced at z, and a random point at angle t
       generated on the circle of intersection. */
    z = 2.0 * (double) rand() / MAX_INT - 1.0;
    t = 2.0 * M_PI * (double) rand() / MAX_INT;
    w = sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
  
    ray[X] = radius * x;
    ray[Y] = radius * y;
    ray[Z] = radius * z;
  
     /*printf( "RandomRay returns %6d %6d %6d\n", ray[X], ray[Y], ray[Z] );*/
}

void AddVec( tPointd q, tPointd ray ){
    int i;
  
    for( i = 0; i < DIM; i++ )
        ray[i] = q[i] + ray[i];
}

int InBox( tPointd q, tPointd bmin, tPointd bmax )
{
    int i;

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
    tPointd N; double D;
    tPointd rq;
    double num, denom, t;
    int i;

    *m = PlaneCoeff( T, N, &D );
    /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],D);*/
    num = D - Dot( q, N );
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
    for( i = 0; i < DIM; i++ ){
        //printf("PPP -1: t=%lf \n", p[i] );
        p[i] = q[i] + t * ( r[i] - q[i] );
    }

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
int	PlaneCoeff( tPointi T, tPointd N, double *D )
{
    int i;
    double t;              /* Temp storage */
    double biggest = 0.0;  /* Largest component of normal vector. */
    int m = 0;             /* Index of largest component. */

    NormalVec( Vertices[T[0]], Vertices[T[1]], Vertices[T[2]], N );
    /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
    *D = Dot( Vertices[T[0]], N );

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
a - b ==> c.
---------------------------------------------------------------------*/
void SubVec( tPointd a, tPointd b, tPointd c )
{
    int i;

    for( i = 0; i < DIM; i++ )
       c[i] = a[i] - b[i];
}

/*---------------------------------------------------------------------
Returns the dot product of the two input vectors.
---------------------------------------------------------------------*/
double  Dot( tPointd a, tPointd b )
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
/* Assumption: p lies in the plane containing T.
    Returns a char:
     'V': the query point p coincides with a Vertex of triangle T.
     'E': the query point p is in the relative interior of an Edge of triangle T.
     'F': the query point p is in the relative interior of a Face of triangle T.
     '0': the query point p does not intersect (misses) triangle T.
*/

char 	InTri3D( tPointi T, int m, tPointd p )
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

char 	InTri2D( tPointd Tp[3], tPointd pp )
{
    int area0, area1, area2;

    /* compute three AreaSign() values for pp w.r.t. each edge of the face in 2D */
    area0 = AreaSign( pp, Tp[0], Tp[1] );
    area1 = AreaSign( pp, Tp[1], Tp[2] );
    area2 = AreaSign( pp, Tp[2], Tp[0] );
    printf("area0=%d  area1=%d  area2=%d\n",area0,area1,area2);

    if ( ( area0 == 0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
         ( area1 == 0 ) && ( area0 > 0 ) && ( area2 > 0 ) ||
         ( area2 == 0 ) && ( area0 > 0 ) && ( area1 > 0 ) ) 
        return 'E';

    if ( ( area0 == 0 ) && ( area1 < 0 ) && ( area2 < 0 ) ||
         ( area1 == 0 ) && ( area0 < 0 ) && ( area2 < 0 ) ||
         ( area2 == 0 ) && ( area0 < 0 ) && ( area1 < 0 ) )
        return 'E';                 
   
    if ( ( area0 >  0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
         ( area0 <  0 ) && ( area1 < 0 ) && ( area2 < 0 ) )
        return 'F';

    if ( ( area0 == 0 ) && ( area1 == 0 ) && ( area2 == 0 ) )
        fprintf( stderr, "Error in InTriD\n" ), exit(EXIT_FAILURE);

    if ( ( area0 == 0 ) && ( area1 == 0 ) ||
         ( area0 == 0 ) && ( area2 == 0 ) ||
         ( area1 == 0 ) && ( area2 == 0 ) )
         return 'V';

    else  
        return '0';  
}

int  AreaSign( tPointd a, tPointd b, tPointd c )  
{
    double area2;

    area2 = ( b[0] - a[0] ) * ( c[1] - a[1] ) -
            ( c[0] - a[0] ) * ( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}                            

char  SegTriInt( tPointi T, tPointd q, tPointd r, tPointd p )
{
    int code = '?';
    int m = -1;

    code = SegPlaneInt( T, q, r, p, &m );
    printf("SegPlaneInt code=%c, m=%d; p=(%lf,%lf,%lf)\n", code,m,p[X],p[Y],p[Z]
);

    if  ( code == '0')
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

char	InPlane( tPointi T, int m, tPointd q, tPointd r, tPointd p)
{
    /* NOT IMPLEMENTED */
    return 'p';
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
 
   printf( "SegTriCross:  vol0 = %d; vol1 = %d; vol2 = %d\n", 
       vol0, vol1, vol2 ); 
     
   /* Same sign: segment intersects interior of triangle. */
   if ( ( ( vol0 > 0 ) && ( vol1 > 0 ) && ( vol2 > 0 ) ) || 
        ( ( vol0 < 0 ) && ( vol1 < 0 ) && ( vol2 < 0 ) ) )
       return 'f';
   
   /* Opposite sign: no intersection between segment and triangle */
   if ( ( ( vol0 > 0 ) || ( vol1 > 0 ) || ( vol2 > 0 ) ) &&
        ( ( vol0 < 0 ) || ( vol1 < 0 ) || ( vol2 < 0 ) ) )
       return '0';

   else if ( ( vol0 == 0 ) && ( vol1 == 0 ) && ( vol2 == 0 ) )
       fprintf( stderr, "Error 1 in SegTriCross\n" ), exit(EXIT_FAILURE);
   
   /* Two zeros: segment intersects vertex. */
   else if ( ( ( vol0 == 0 ) && ( vol1 == 0 ) ) || 
             ( ( vol0 == 0 ) && ( vol2 == 0 ) ) || 
             ( ( vol1 == 0 ) && ( vol2 == 0 ) ) )
       return 'v';

   /* One zero: segment intersects edge. */
   else if ( ( vol0 == 0 ) || ( vol1 == 0 ) || ( vol2 == 0 ) )
       return 'e';
   
   else
       fprintf( stderr, "Error 2 in SegTriCross\n" ), exit(EXIT_FAILURE);
}

int  VolumeSign( tPointd a, tPointd b, tPointd c, tPointd d )
{ 
    double vol;
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

    bxdx = bx - dx;
    bydy = by - dy;
    bzdz = bz - dz;
    cxdx = cx - dx;
    cydy = cy - dy;
    czdz = cz - dz;
    vol =   (az-dz) * (bxdx*cydy - bydy*cxdx)
         + (ay-dy) * (bzdz*cxdx - bxdx*czdz)
         + (ax-dx) * (bydy*czdz - bzdz*cydy);


    /* The volume should be an integer. */
    if      ( vol > 0.5 )   return  1;
    else if ( vol < -0.5 )  return -1;
    else                    return  0;
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
        if ( ((int)a[i] < w ) && ((int)b[i] < w) ) return '0';
            w = Box[ n ][1][i]; /* max: upper right */
        if ( ((int)a[i] > w) && ((int)b[i] > w) ) return '0';
    }
    return '?';
}


/* irint not available in some libraries, so... */

int irint( double x )
{
    return (int) rint( x );
}
void read_ori(void)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    float a,b,c;
    fp = fopen("big.off", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
        count++;
        char *token = strtok(line, " "); 
        int token_count = 0;
        while (token != NULL) {  
            // init facets and vertices
            if(count <= 2){
                printf("setting of file  %s\n", token);
                if(token_count == 0){
                    n_vertices = atoi(token);
                }else if(token_count == 1){
                    n_facets = atoi(token);
                }
                token_count++;

            }else if(count > 3 || count <  n_vertices+ 4){
                if(token_count == 0){
                    Vertices[count - 4][X] = atof(token);
                }else if(token_count == 1){
                    Vertices[count - 4][Y] = atof(token);
                }else{
                    Vertices[count - 4][Z] = atof(token);
                }
                token_count++;
            } else{
                if(token_count == 1){
                    Faces[count - 144][X] = atoi(token);
                }else if(token_count == 2){
                    Faces[count - 144][Y] = atoi(token);
                }else if(token_count == 3){
                    Faces[count - 144][Z] = atoi(token);
                }
                token_count++;
            } 
            token = strtok(NULL, " "); 
        } 
    }
    if (line)
        free(line);
}
void read_com(void)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    float a,b,c;
    fp = fopen("small.off", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
        count++;
        char *token = strtok(line, " "); 
        int token_count = 0;
        while (token != NULL) {  
            // init facets and vertices
            if(count <= 2){
                printf("setting of file  %s\n", token);
                if(token_count == 0){
                    com_vertices = atoi(token);
                }else if(token_count == 1){
                    com_facets = atoi(token);
                }
                token_count++;
            }else if(count > 3 && count < n_vertices + 4){
                if(token_count == 0){
                    com_Vertices[count - 4][X] = atof(token);
                }else if(token_count == 1){
                    com_Vertices[count - 4][Y] = atof(token);
                }else{
                    com_Vertices[count - 4][Z] = atof(token);
                }
                token_count++;
            } else{
                if(token_count == 1){
                    com_Faces[count - 144][X] = atoi(token);
                }else if(token_count == 2){
                    com_Faces[count - 144][Y] = atoi(token);
                }else if(token_count == 3){
                    com_Faces[count - 144][Z] = atoi(token);
                }
                token_count++;
            } 
            token = strtok(NULL, " "); 
        } 
    }
    if (line)
        free(line);
}
void init_bounding(int n){
   int i,j,k,w;
   for ( i = 0; i < n; i++ ) {
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
    /*
    printf("Bounding box: (%d,%d,%d);(%d,%d,%d)\n",
       Box[i][0][0],
       Box[i][0][1],
       Box[i][0][2],
       Box[i][1][0],
       Box[i][1][1],
       Box[i][1][2] );
    */
  }
}