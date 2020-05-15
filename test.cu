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
    printf("counter %d\n,", counter);
    while( counter >= 0 ) {
        q[X] = com_Vertices[counter][X];
        q[Y] = com_Vertices[counter][Y];
        q[Z] = com_Vertices[counter][Z];
        printf( "\n %d -------->q = %lf %lf %lf\n", counter, q[X], q[Y], q[Z] );
        //printf( "In = %c\n", InPolyhedron( F, q, bmin, bmax, radius ) );
        counter--;
    }
    return 0;
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
  /*
  for( i = 0; i < n; i++ )
    for( j = 0; j < DIM; j++ ) {
      if( Vertices[i][j] < bmin[j] )
              bmin[j] = Vertices[i][j];
      if( Vertices[i][j] > bmax[j] ) 
              bmax[j] = Vertices[i][j];
    }
  */
  dim3 blockSize(256);
  dim3 gridSize((n + blockSize.x) / blockSize.x);
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
  //free(bmin);
  //free(bmax);
  printf("test");
  /*
  radius = sqrt( pow( (double)(bmax[X] - bmin[X]), 2.0 ) +
                 pow( (double)(bmax[Y] - bmin[Y]), 2.0 ) +
                 pow( (double)(bmax[Z] - bmin[Z]), 2.0 ) );
  printf("radius = %lf\n", radius);*/
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
               
               printf("Bounding box: (%d,%d,%d);(%d,%d,%d)\n",
                  Box[i][0][0],
                  Box[i][0][1],
                  Box[i][0][2],
                  Box[i][1][0],
                  Box[i][1][1],
                  Box[i][1][2] );
               
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
