#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#define X 0
#define Y 1
#define Z 2
#define DIM 3                  /* Dimension of points */
typedef int    tPointi[DIM];   /* Type integer point */
typedef double tPointd[DIM];   /* Type double point */
#define PMAX 100000             /* Max # of pts */
tPointd Vertices[PMAX];        /* All the points */
tPointi Faces[PMAX];           /* Each triangle face is 3 indices */
int main(void)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    float a,b,c;
    fp = fopen("0.off", "r");
    int n_facets, n_vertices;
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
            }else if(count > 3 && count < 144){
                if(token_count == 0){
                    Vertices[count - 4][X] = atof(token);
                }else if(token_count == 1){
                    Vertices[count - 4][Y] = atof(token);
                }else{
                    Vertices[count - 4][Z] = atof(token);
                }
                token_count++;
                //printf("vertice --> %s\n", token);
            } else{
                if(token_count == 1){
                    Faces[count - 144][X] = atoi(token);
                }else if(token_count == 2){
                    Faces[count - 144][Y] = atoi(token);
                }else if(token_count == 3){
                    Faces[count - 144][Z] = atoi(token);
                }
                //printf("faces--> %s\n", token);
                token_count++;
            } 
            
            //printf("XDDD --> %s\n", token); 
            
            token = strtok(NULL, " "); 
        } 
        
        //printf("Retrieved line of length %zu:\n", read);
        //printf("%f \n", a);
    }
    for(int i = 0; i < n_vertices; i++){
        printf("No %d, facets --> %f, %f, %f \n", i,Vertices[i][X],Vertices[i][Y],Vertices[i][Z]);
    }
    printf("n_vertices --> %d\n",n_vertices);
    printf("n_f--> %d\n",n_facets);
    fclose(fp);
    if (line)
        free(line);
    exit(EXIT_SUCCESS);
}