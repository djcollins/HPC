
#include <stdio.h>
#include <stdlib.h>
#include "includes/dcdplugin.c"
#include <math.h>
#include <time.h>
#include <omp.h> //NB
#define _GNU_SOURCE
#include "serial.h"
//#include "group.h"

static int getRange(char *indec){
	//printf("indeces are %s\n",indeces);
	int size=0;
	char l[256];
	strcpy(l,indec);
	char *e;
	int ind;
	e=strchr(l,'-');
	ind=(int) (e-l);
	char left[100];
	char right[100];
	memcpy(left,&l[0],ind);
	memcpy(right,&l[ind]+1,strlen(l));
	//int low=atoi(strtok(NULL,"-"));
	//int hi=atoi(strtok(NULL,","));
    //printf("found a range from %s to %s \n",left,right);
	return (atoi(right)-atoi(left)+1);//(upper-lower);

}
int * getAts(char *ind,int size){
	int *atoms=malloc(sizeof(int)*size);
	int counter=0;
	char line[256];
	ind[strcspn(ind, "\r\n")] = 0;
	strcpy(line, ind);
	char *token=strtok(line,",");
    while(token!=NULL){
    	if (strchr(token,'-')!=NULL)
    	{
    		//printf("dash foune %s\n",token);
    		int size=0;
			char l[256];
			strcpy(l,token);
			char *e;
			int ind;
			e=strchr(l,'-');
			ind=(int) (e-l);
			char left[100];
			char right[100];
			memcpy(left,&l[0],ind);
			memcpy(right,&l[ind]+1,strlen(l));
			for(int s=atoi(left);s<atoi(right)+1;s++){
				atoms[counter]=s;
    			counter=counter+1;	
			}
    	}
    	else{
    		atoms[counter]=atoi(token);
    		counter=counter+1;
    	}
    	token=strtok(NULL,",");
    	//printf("line is %s token is %s \n",line,token);
    }
	return atoms;
}

struct group * getList(char *indeces){
	int groupSize=0; //we predetermine the size of the array needed to hold the atom indexes
	char line[256];
	indeces[strcspn(indeces, "\r\n")] = 0;
	strcpy(line, indeces);
	char *token=strtok(line,",");
    while(token!=NULL){


    	if (strchr(token,'-')!=NULL)
    	{

    		//printf("dash foune %s\n",token);
    		groupSize=groupSize+getRange(token);
    	}
    	else{
    		groupSize=groupSize+1;
    	}
    	token=strtok(NULL,",");
    	//printf("line is %s token is %s \n",line,token);
    }
    int * atoms=malloc(sizeof(int)*groupSize);
    atoms= getAts(indeces,groupSize);
    struct group* g=(struct group*)malloc(sizeof(struct group));
    g->groupSize=groupSize;
    g->groupAtoms=atoms;
    //free(atoms);
    return g;
}

int main(int argc, char **argv){
	//char *fileName="input.txt";
	char *tst="input.txt";
	//reading from file
	int getopt(int argc, char * const argv[], const char *optstring);
  	extern char *optarg;
  	extern int optind, opterr, optopt;
	int ich;
	char *fileName;
  	while ((ich = getopt (argc, argv, "ai:c")) != EOF) {
    	switch (ich) {
      		case 'a': /* Flags/Code when -a is specified */
        		break;
      		case 'i': 
      				fileName=optarg;
      				printf("found file name %s\n",fileName);
      				/* Flags/Code when -b is specified */
                	/* The argument passed in with b is specified */
					/* by optarg */
        		break;
      		case 'c': /* Flags/Code when -c is specified */
        		break;
      		default: /* Code when there are no parameters */
        		printf("default file name\n");
        		fileName="input";
        		break;
    	}
  	}
 	char dcd_file[256];
	char tmp[256];
	char a[256];
	char b[256];
	int k;
  	fileName=tst;
  	//printf("files name is %s\n",fileName);
	if(fileName!=NULL){
   		FILE *file = fopen ( fileName, "r" );
   		//printf("got here\n");
   		if ( file != NULL ){

			if(fgets(dcd_file, 100, file) && fgets(tmp, 100, file) && fgets(a, 100, file) && fgets(b, 100, file)){
				dcd_file[strcspn(dcd_file,"\r\n")]=0;
				tmp[strcspn(tmp,"\r\n")]=0;
				a[strcspn(a,"\r\n")]=0;
				b[strcspn(b,"\r\n")]=0;
				k=atoi(tmp);
				//printf("dcd is %s %s %s %s\n",dcd_file,tmp,a,b);	
				fclose ( file );
			}
			else{
				printf("input error in text file");
				return 0;	
			}	
   		}
   		else{
      		printf("file named %s was not found\n",fileName);
      		return 0;
   		}
		
	}
	else{
		printf("file named %s is empty?\n",fileName);
		return 0;
	}
    // int natoms=0;
    // void *raw_data = open_dcd_read(dcd_file, "dcd", &natoms);
    // if (!raw_data) {
    //     printf("File named %s not found\n",dcd_file);
    //     return 1;
    // }


    //now we call a function to return group structures,
    //parsing the indeces from the textfile
    struct group* g1 =getList(a);
    struct group *g2=getList(b);

    //group 1 must be bigger, do check and swap if needed

    if(g1->groupSize<g2->groupSize){
    	struct group temp=*g1;
    	*g1=*g2;
    	*g2=temp;
    }
    int wow=doSerial(k,&g1,&g2,dcd_file);
 //    dcdhandle *dcd = (dcdhandle *) raw_data;
 //    molfile_timestep_t timestep;
 //    timestep.coords = (float *) malloc(3 * sizeof(float) * natoms);
 //    int read_failed = 0;
 //    double start = omp_get_wtime();
	// for (int i = 0; i < dcd->nsets; ++i) {
 //        int rc = read_next_timestep(raw_data, natoms, &timestep);
 //        if (rc) {
 //            read_failed = 1;
 //            break;
 //        }
 //        int a= sizeof(timestep.coords)/sizeof(int);
 //        printf("%d %d %d\n",a,sizeof(int),sizeof(timestep.coords));
 //        if (i<5)
 //        {
 //        	printf("%d %d %d \n",dcd->x[313],dcd->y[313],dcd->z[313]);
 //        	float ax=dcd->x[313];
 //        	float ay=dcd->y[313];
 //        	float az=dcd->z[313];
 //        	printf("%d %d %d \n",dcd->x[168052],dcd->y[168052],dcd->z[168052]);
 //        	float bx=dcd->x[168052];
 //        	float by=dcd->y[168052];
 //        	float bz=dcd->z[168052];
 //        	float dis=sqrt((ax-bx)*(ax-bx) +(ay-by)*(ay-by)+(az-bz)*(az-bz));
 //        	printf("step %i distance is %f \n",i,dis);
 //        }
 //         At this point 
 //           dcd contains
 //           dcd->x = Array of X coordinates of all atoms for timestep i
 //           dcd->y = Array of Y coordinates of all atoms for timestep i
 //           dcd->z = Array of Z coordinates of all atoms for timestep i
           
 //           timestep contains
 //           timestep.coords = Array of packed XYZ coordinates of all atoms for timestep i
 //                             [X1, Y1, Z1, X2, Y2, Z2, ..., Xn, Yn, Zn] where n = natoms
          
 //           Both are overwritten next loop 
        
 //        //printf("Timestep %d\n", i);
 //        //printf("i: x    y      z\n");
 //        int n = natoms > 10 ? 10 : natoms;
 //        float tx, ty, tz, *current;
 //        #pragma omp parallel private(tx,ty,tz,current) shared(timestep)
 //        {
 //        for (int j = 0; j < n+100000; ++j) {
 //            current = timestep.coords + 3 * j;
 //            tx=current[0];
 //            ty=current[1];
 //            tz=current[2];
 //            //printf("%d: %3.2f %3.2f %3.2f\n", j, current[0], current[1], current[2]);
 //        }
 //    }
 //        //printf("\n");
 //        if (i >= 10) break;
 //    }
 //    double end=omp_get_wtime();
 //    printf("Time: \t %f \n", end-start);
 //    free(timestep.coords);    
 //    close_file_read(raw_data);
 //    if (read_failed) {
 //        return 2;
 //    }

 //    printf("done\n");
 //    return 0;
    }
    //use this to read in certain indexes
//int read_dcdstep(fio_fd fd, int N, float *X, float *Y, float *Z, 
  //                       float *unitcell, int num_fixed,
    //                     int first, int *indexes, float *fixedcoords, 
      //                   int reverseEndian, int charmm) {