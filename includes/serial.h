#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#define _GNU_SOURCE

struct group{	
	int groupSize;
	int * groupAtoms;		
};



int doSerial(int k, struct group** g1 , struct group** g2, char * dcd_file){	
	int natoms=0;
    void *raw_data = open_dcd_read(dcd_file, "dcd", &natoms);
    if (!raw_data) {
        printf("File named %s not found\n",dcd_file);
        return 1;
    }

	dcdhandle *dcd = (dcdhandle *) raw_data;
    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);
    int read_failed = 0;
	int subA[k];
	int subB[k];
	float dists[k];
	printf("%f\n",FLT_MAX );
	for (int i = 0; i < k; ++i)
	{
		dists[i]=FLT_MAX;
	}
	float curSlow=FLT_MAX;
    double start = omp_get_wtime();
	for (int i = 0; i < dcd->nsets; ++i) {
        int rc = read_next_timestep(raw_data, natoms, &timestep);
        if (rc) {
            read_failed = 1;
            break;
        }
        int a= sizeof(timestep.coords)/sizeof(int);
        //printf("%d %d %d\n",a,sizeof(int),sizeof(timestep.coords));
        /* At this point 
           dcd contains
           dcd->x = Array of X coordinates of all atoms for timestep i
           dcd->y = Array of Y coordinates of all atoms for timestep i
           dcd->z = Array of Z coordinates of all atoms for timestep i
           
           timestep contains
           timestep.coords = Array of packed XYZ coordinates of all atoms for timestep i
                             [X1, Y1, Z1, X2, Y2, Z2, ..., Xn, Yn, Zn] where n = natoms
          
           Both are overwritten next loop 
        */
        //printf("Timestep %d\n", i);
        //printf("i: x    y      z\n");
        int n1 = (*g1)->groupSize;// > 10 ? 10 : natoms;
        int n2 = (*g2)->groupSize;
        int *p1=(*g1)->groupAtoms;
        int *p2=(*g2)->groupAtoms;
        float ax, ay, az,bx,by,bz, *current,*current2,curDis;
        


        for (int j = 0; j < n1; ++j) {
            current = timestep.coords + 3 * p1[j];
            ax=current[0];
            ay=current[1];
            az=current[2];
            for (int f = 0; f < n2; ++f)
            {
            	current2 = timestep.coords + 3 * p2[f];
				bx=current2[0];
            	by=current2[1];
            	bz=current2[2];
            	curDis=sqrt((ax-bx)*(ax-bx) +(ay-by)*(ay-by)+(az-bz)*(az-bz));
            	//printf("%f\n",curDis);
            	if (curDis<20.97)
            	{
            		printf("%d %d %d %3.10f \n",i,p1[j],p2[f],curDis);
            	}
            }
            //printf("%d\n", p1[j]);
            //printf("%d: %3.2f %3.2f %3.2f\n", j, current[0], current[1], current[2]);
        }
    	if (i>10)
    	{
    		break;
    	}
        //printf("\n");
    }
    double end=omp_get_wtime();
    printf("Time: \t %f \n", end-start);
    free(timestep.coords);    
    close_file_read(raw_data);
    if (read_failed) {
        return 2;
    }

    printf("done in series \n");
    return 0;


	return 0;
}