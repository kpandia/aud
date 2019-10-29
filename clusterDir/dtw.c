#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>

int nModels,*n,d,*indexArray;
struct mat {
	int n;
	float **val;
	};struct mat *inMat;

float computeDTW(int x, int y, int d) {
	int n1,n2,i,j,k;
	float **distMat,**dtwMat,total,tempf,top,mid,bot,cheapest;
	n1=n[x];n2=n[y];
	dtwMat=(float **)malloc(n1*sizeof(float *));
	distMat=(float **)malloc(n1*sizeof(float *));
	for(i=0;i<n1;i++) {
		distMat[i]=(float *)malloc(n2*sizeof(float));
		dtwMat[i]=(float *)malloc(n2*sizeof(float));
		}
	for(i=0;i<n1;i++) {
		for(j=0;j<n2;j++) {
			total=0;
			for(k=0;k<d;k++) {
				total+=pow(((inMat[y].val[j][k])-(inMat[x].val[i][k])),2);
				}
			distMat[i][j]=sqrt(total)/d;
			}
		}


	dtwMat[0][0]=distMat[0][0];
	for(i=1;i<n1;i++)
		dtwMat[i][0]=dtwMat[i-1][0]+distMat[i][0];
	for(i=1;i<n2;i++)
		dtwMat[0][i]=dtwMat[0][i-1]+distMat[0][i];
	for(i=1;i<n1;i++) {
		for(j=1;j<n2;j++) {
			top=dtwMat[i][j-1];                                                                     
			mid=dtwMat[i-1][j-1];                                                                   
            bot=dtwMat[i-1][j];                                                                     
            if(top<mid && top<bot) {                                                                  
                cheapest=top;
                }                                                                                     
            else if(mid<top && mid<bot) {                                                             
                cheapest=mid;
                }                                                                                     
            else {                                                                                    
                cheapest=bot;
                }                                                                                     
            dtwMat[i][j]=cheapest+distMat[i][j];                                                     
			}
		}
	total=dtwMat[i-1][j-1]/n1;
	for(i=0;i<n1;i++){
		free(distMat[i]);
		free(dtwMat[i]);
		}
	return(total);
	}

void main(int argc, char *argv[]) {
	if(argc != 5) {
		printf("Usage: %s list #models model_index outFile\n",argv[0]);
		exit(0);
		}
	int i,j,k,rc,X;
	float **outMat;
	FILE *listFP, *tempFP, *outFile;
	listFP=fopen(argv[1],"r");
	nModels=atoi(argv[2]);
	X=atoi(argv[3]);
	outFile=fopen(argv[4],"w");
	struct inMat;


	outMat=(float **)malloc(nModels*sizeof(float *));
	indexArray=(int *)malloc(nModels*sizeof(int ));
	inMat=(struct mat *)malloc(nModels*sizeof(struct mat ));
	n=(int *)malloc(nModels*sizeof(int **));
	for(i=0;i<nModels;i++) {
		outMat[i]=(float *)malloc(nModels*sizeof(float));
		}
	printf("Loading files ...\n");
	for(i=0;i<nModels;i++) {
		char fileName[100];
		fscanf(listFP,"%s\n",fileName);
		tempFP=fopen(fileName,"r");
		fscanf(tempFP,"%d %d\n",&d,&n[i]);
		inMat[i].n=n[i];
		inMat[i].val=(float **)malloc(n[i]*sizeof(float *));
		for(j=0;j<n[i];j++) {
			inMat[i].val[j]=(float *)malloc(d*sizeof(float));
			for(k=0;k<d;k++) {
				fscanf(tempFP,"%f",&inMat[i].val[j][k]);
				}
			}
		fclose(tempFP);	
		}
	printf("Computing DTW distance\n");
	for(i=0;i<nModels;i++) {
		fprintf (outFile,"%f ",computeDTW(X,i,d));
		}
	fprintf(outFile,"\n");
	printf("Completed !!\n");
		
	}
