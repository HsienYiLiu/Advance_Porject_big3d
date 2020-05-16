#include <curand_kernel.h>
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
//typedef enum { FALSE, TRUE } bool;

#define DIM 3                  /* Dimension of points */
typedef int    tPointi[DIM];   /* Type integer point */
typedef double tPointd[DIM];   /* Type double point */
#define PMAX 1000000             /* Max # of pts */
typedef enum boolean{ FALSE, TRUE } boolean;
tPointd Vertices[PMAX];        /* All the points */
tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
tPointd com_Vertices[PMAX];
tPointi com_Faces[PMAX];
int check = 0;
tPointi Box[PMAX][2];          /* Box around each face */
int n_facets, n_vertices;      /* Original polyhedron*/
int com_facets, com_vertices;  /* Original polyhedron*/

void read_ori(void);
void read_com(void);
int ComputeBox( int F, tPointd bmin, tPointd bmax );
int irint( double x );
char BoxTest ( int n, tPointd a, tPointd b );
__device__ int InBox( tPointd q, tPointd bmin, tPointd bmax );
__device__ void RandomRay( tPointd ray, int radius );
__device__ void AddVec( tPointd q, tPointd ray );
__global__ void InPolyhedron( int F,  tPointd * bmin, tPointd * bmax,int radius, tPointd * c_com_V, int * out);
//read_ori();
int main(){
    int n, F, i;
    tPointd q, bmin, bmax;
    int radius;
    read_ori();
    read_com();
    n = n_vertices;
    F = n_facets;
    // Allocate memory
    for ( i = 0; i < DIM; i++ ){
        bmin[i] = bmax[i] = Vertices[0][i];
    }
    radius = ComputeBox( n, bmin, bmax );
    int counter = com_vertices - 1;
    //printf("counter %d\n,", counter);
    // setting for judge
    tPointd *d_bmin, *d_bmax, *c_com_V;
    int *out,*result;
    
    //char out[counter];
    //char result[counter];
    result = (int *)malloc(sizeof(int)*counter);
   
    cudaMalloc(&c_com_V,sizeof(tPointd)*counter);
    cudaMalloc(&d_bmax,sizeof(tPointd)*3);
    cudaMalloc(&d_bmin,sizeof(tPointd)*3);
    cudaMalloc(&out,sizeof(int)*counter);

    cudaMemcpy(c_com_V, com_Vertices, sizeof(tPointd)*counter, cudaMemcpyHostToDevice);
    cudaMemcpy(d_bmin, bmin, sizeof(tPointd)*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(d_bmax, bmax, sizeof(tPointd)*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(out, result, sizeof(int)*counter, cudaMemcpyHostToDevice);
    InPolyhedron<<<counter, 1>>>(n,d_bmin, d_bmax,radius,c_com_V, out);
    printf("111111\n");
    cudaMemcpy(result,out, sizeof(int)*counter, cudaMemcpyDeviceToHost);
    printf("2222222\n");
    printf("resultFFFF %d\n",result[1]);

    cudaFree(d_bmax);
    cudaFree(d_bmin);
    cudaFree(c_com_V);
    cudaFree(out);
    free(result);   

    //InPolyhedron( F, q, bmin, bmax, radius )
    /*
    while( counter >= 0 ) {
        q[X] = com_Vertices[counter][X];
        q[Y] = com_Vertices[counter][Y];
        q[Z] = com_Vertices[counter][Z];
        printf( "\n %d -------->q = %lf %lf %lf\n", counter, q[X], q[Y], q[Z] );
        //printf( "In = %c\n", InPolyhedron( F, q, bmin, bmax, radius ) );
        counter--;
    }*/
    return 0;
}
__global__ void check_each(int F,tPointd * ori_F,tPointd * r, tPointd *q){

}

__global__ void InPolyhedron( int F,  tPointd * bmin, tPointd * bmax,int radius, tPointd * c_com_V, int * out )
{
   //volatile bool *found = FALSE;
   volatile __shared__ bool FoundIt;
   // initialize shared status
    FoundIt = false;
   //__syncthreads();
   tPointd r;  /* Rayendpoint. */
   tPointd p,q;  /* Intersection point; not used. */
   int f, k = 0, crossings = 0;
   int code = -1;
   int i = blockIdx.x;
   printf("i %d\n",i);
   q[0] = c_com_V[i][0];
   q[1] = c_com_V[i][1];
   q[2] = c_com_V[i][2];
   /* If query point is outside bounding box, finished. */
   if ( !InBox( q, *bmin, *bmax ) ){
      out[i] = 3;
      FoundIt = true;
      printf("wpwowow %d\n", out[i]);
   }
      //return 'o';
   
   LOOP:
   while( k++ < F && FoundIt == false) {
      crossings = 0;
  
      RandomRay( r, radius ); 
      AddVec( q, r ); // add the ray with the point to create end point
      printf("Ray endpoint: (%lf,%lf,%lf)\n", r[0],r[1],r[2] );
      /*
      for ( f = 0; f < F; f++ ) {  // Begin check each face 
         if ( BoxTest( f, q, r ) == '0' ) {
              out[i] = '0';
              //code = '0';
              printf("BoxTest = 0!\n");
         }
         //else code = SegTriInt( Faces[f], q, r, p );
         printf( "Face = %d: BoxTest/SegTriInt returns %c\n\n", f, code );

         //If ray is degenerate, then goto outer while to generate another.
         if ( code == 'p' || code == 'v' || code == 'e' ) {
            printf("Degenerate ray\n");
            goto LOOP;
         }
   
         //If ray hits face at interior point, increment crossings.
         else if ( code == 'f' ) {
            crossings++;
            printf( "crossings = %d\n", crossings );
         }

         //If query endpoint q sits on a V/E/F, return that code.
         else if ( code == 'V' || code == 'E' || code == 'F' )
            //return code;
            out[i] = code;

         //If ray misses triangle, do nothing. 
         else if ( code == '0' )
            ;

         else 
            fprintf( stderr, "Error, exit(EXIT_FAILURE)\n" ), exit(1);

      } */
      break;

   }  
   printf( "Crossings = %d\n", crossings );
   /* q strictly interior to polyhedron iff an odd number of crossings. */
   if( ( crossings % 2 ) == 1 )
      //return   'i';
      out[i] = 1;
   //else return 'o';
   else out[i] = 9;
   free(r);
   printf("result -->  %d\n", out[i]);
}
__device__ int InBox( tPointd q, tPointd bmin, tPointd bmax )
{
  int i;

  if( ( bmin[X] <= q[X] ) && ( q[X] <= bmax[X] ) &&
      ( bmin[Y] <= q[Y] ) && ( q[Y] <= bmax[Y] ) &&
      ( bmin[Z] <= q[Z] ) && ( q[Z] <= bmax[Z] ) )
    return TRUE;
  return FALSE;
}
/* Return a random ray endpoint */
__device__ void RandomRay( tPointd ray, int radius )
{
  double x, y, z, w, t;
  int tId = threadIdx.x + (blockIdx.x * blockDim.x);
  curandState state;
  curand_init((unsigned long long)clock() + tId, 0, 0, &state);

  double rand1 = curand_uniform_double(&state);
  double rand2 = curand_uniform_double(&state);
  /* Generate a random point on a sphere of radius 1. */
  /* the sphere is sliced at z, and a random point at angle t
     generated on the circle of intersection. */
  z = 2.0 * (double) 0 / MAX_INT - 1.0;
  t = 2.0 * M_PI * (double) rand2 / MAX_INT;
  printf("check %lf\n",rand1);
  w = sqrt( 1 - z*z );
  x = w * cos( t );
  y = w * sin( t );
  
  ray[X] = radius * x;
  ray[Y] = radius * y;
  ray[Z] = radius * z;
  
  /*printf( "RandomRay returns %6d %6d %6d\n", ray[X], ray[Y], ray[Z] );*/
}
__device__ void AddVec( tPointd q, tPointd ray )
{
  int i;
  
  for( i = 0; i < DIM; i++ )
    ray[i] = q[i] + ray[i];
}
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
__global__ void cal(tPointd *bmin, tPointd *bmax,tPointd *V,int F){

    int i = blockIdx.x; // will give you X block Index at that particular thread
    int j = blockIdx.y; // will give you Y block Index at that particular thread. 
    if(i < F){
        //j = j%3;
        for(j = 0; j < 3; j++){
            if( V[i][j] < *bmin[j] )
                *bmin[j] = V[i][j];
            if( V[i][j] > *bmax[j] ){
                *bmax[j] = V[i][j];
                //printf("V %lf\n",V[i][j]);
            }
            printf("Check i = %d, j = %d, F = %d\n",i,j,F);
        }
    }
    printf("bmax %lf, bmin %lf \n",*bmax[Y],*bmin[Y]);
}
int ComputeBox( int n, tPointd bmin, tPointd bmax ){
  int i, j;
  double radius;
  tPointd *d_bmin, *d_bmax, *d_a, *max, *min;
  max = (tPointd *)malloc(sizeof(tPointd)*DIM); // Allocate array1 on host
  min = (tPointd *)malloc(sizeof(tPointd)*DIM); // Allocate array2 on host 

  cudaMalloc(&d_a,sizeof(tPointd)*n);
  cudaMalloc(&d_bmax,sizeof(tPointd)*3);
  cudaMalloc(&d_bmin,sizeof(tPointd)*3);

  cudaMemcpy(d_a, Vertices, sizeof(tPointd)*n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bmin, bmin, sizeof(tPointd)*DIM, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bmax, bmax, sizeof(tPointd)*DIM, cudaMemcpyHostToDevice);

  //dim3 blockSize(256);
  //dim3 gridSize((n + blockSize.x) / blockSize.x);
  cal<<<n, 1>>>(d_bmin, d_bmax, d_a, n);
  cudaMemcpy(max,d_bmax, sizeof(tPointd)*DIM, cudaMemcpyDeviceToHost);
  cudaMemcpy(min,d_bmin, sizeof(tPointd)*DIM, cudaMemcpyDeviceToHost);
  printf("------------------------\n");
  printf("bmax %lf bmin %lf \n",*max[X],*min[X]);
  printf("bmax %lf, bmin %lf \n",*max[Y],*min[Y]);
  printf("bmax %lf, bmin %lf \n",*max[Z],*min[Z]);
  radius = sqrt( pow( (double)(*max[X] - *min[X]), 2.0 ) +
                 pow( (double)(*max[Y] - *min[Y]), 2.0 ) +
                 pow( (double)(*max[Z] - *min[Z]), 2.0 ) );
  printf("radius = %lf\n", radius);
  cudaFree(d_bmax);
  cudaFree(d_bmin);
  cudaFree(d_a);
  free(max);
  free(min);

  return irint( radius +1 ) + 1;
}

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
    int i = 0;
    int j,k,n,w;

    if (fp == NULL)
        exit(EXIT_FAILURE);
    while ((read = getline(&line, &len, fp)) != -1) {
        count++;
        char *token = strtok(line, " ");
        int token_count = 0;
        while (token != NULL ) {
            // init facets and vertices
            if(count <= 2){
                if(token_count == 0){
                    n_vertices = atoi(token);
                }else if(token_count == 1){
                    n_facets = atoi(token);
                }
                token_count++;
            }else if(count > 3 && count <  n_vertices + 4){
                if(token_count == 0){
                    Vertices[count - 4][X] = atof(token);
                }else if(token_count == 1){
                    Vertices[count - 4][Y] = atof(token);
                }else{
                    Vertices[count - 4][Z] = atof(token);
                }
                token_count++;
            } else if(count >= n_vertices + 4){
                i = count - n_vertices - 4;

                if(token_count == 1){
                    Faces[i][X] = atoi(token);
                }else if(token_count == 2){
                    Faces[i][Y] = atoi(token);
                    //printf("->>>>  %d\n",Faces[count - 144][X]);
                }else if(token_count == 3){
                    Faces[i][Z] = atoi(token);
                    for ( j=0; j < 3; j++ ) {
                        Box[i][0][j] = Vertices[ Faces[i][0] ][j];
                        Box[i][1][j] = Vertices[ Faces[i][0] ][j];
                  }

               for ( k=1; k < 3; k++ )
               for ( j=0; j < 3; j++ ) {
                  w = Vertices[ Faces[i][k] ][j];
                  //printf("->>>>  %d\n",Faces[i][k]);
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
                token_count++;
            }
            token = strtok(NULL," ");
            //free(token);
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
    fp = fopen("t.off", "r");
    int i ;
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
            }else if(count > 3 && count <  n_vertices + 4){
               if(token_count == 0){
                    com_Vertices[count - 4][X] = atof(token);
                }else if(token_count == 1){
                    com_Vertices[count - 4][Y] = atof(token);
                }else{
                    com_Vertices[count - 4][Z] = atof(token);
                }
                token_count++;
            }else if(count >= n_vertices + 4){
                i = count - n_vertices - 4;
                if(token_count == 1){
                    com_Faces[i][X] = atoi(token);
                }else if(token_count == 2){
                    com_Faces[i][Y] = atoi(token);
                }else if(token_count == 3){
                    com_Faces[i][Z] = atoi(token);
                }
                token_count++;
            }
            token = strtok(NULL, " ");
        }
    }
    if (line)
        free(line);
}
