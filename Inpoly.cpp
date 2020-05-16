char BoxTest(int n, tPointi a, tPointi b){
    int i;
    int w;
    for(i = 0; i < DIM; i++){
        w = Box[n][0][i];
        if( (a[i] < w) && (b[i] < w)) return '0';
        w = Box[n][1][i];
        if( (a[i] > w) && (b[i] > w)) return '0';
    }
    return '?';
}
char InPolyhedron(int F, tPointi q, tPointi bmin, tPointi bmax, int radius){
    tPointi r;
    tPointi p;
    int f, k = 0, crossing = 0;
    char code = '?';
    /* If query point is outside bounding box, finished*/
    if( !Inbox(q, bmin, bmax)){
        return 'o';
    }
    LOOP:
    while(k++ < F){
        crossing = 0;
        RandomRay(r, radius);
        AddVec(q,r,r);
        for(f = 0; f < F; f++){
            if(BoxTest(f, q, r) == 'o')
                code = 'o';
            else code = SegTriInt(Faces[f], q, r, p);
            if( code == 'p' || code == 'v' || code == 'e'){
                printf("Degenerate Ray\n");
                goto LOOP;
            }else if(code == 'f'){
                crossing++;
                printf("crossing number = %d\n", crossing);
            }else if(code == 'V' || code == 'E' || code == 'F' ){
                return(code);
            }else if(code =='O');
            else{
                fprintf(stderr, "ERROR, exit(EXIT_FAILURE)\n"), exit(1);
            }
        }
        break;
    }
    if(crossing % 2 == 1)
        return 'i';
    else return 'o';
}