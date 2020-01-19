/*
    This file contains the functions for clustering using VQ

    Copyright (C) 1998-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>

    This file is part of SpeakerID-IITM.

    SpeakerID-IITM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SpeakerID-IITM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SpeakerID-IITM.  If not, see <http://www.gnu.org/licenses/>. 
*/


#include "fe/FrontEndDefs.h"
#include "fe/FrontEndTypes.h"
//#include "fe/InitAsdf.h"
#include "fe/DspLibrary.h"
#include "stdlib.h"
#include "math.h"
/****************************************************************************
 *   Function             : initVQ - initialises the codebook with a set
 *                        : of vectors chosen randomly from the input codebook
 *   Input args           : asdf (front-end structure), 
 *                        : vfv - vector of featureVectors
 *                        : numClusters : codebook size
 *   Output args          : clusters - vector of codebook vectors     	  
 *****************************************************************************/



void InitVQ (VECTOR_OF_F_VECTORS *vfv,int numVectors, 
              VECTOR_OF_F_VECTORS *clusterMeans, 
	     VECTOR_OF_F_VECTORS *clusterVars, int numClusters, int seed) {
  int                          index;
  int                          i, j;
  int                          random;
  float                        rmax;
  int                         *indexArray;
  int                          flag;
  srand(seed);
  fflush(stdout);
  indexArray = (int *) calloc (numClusters, sizeof(int));
  for (i = 0; i < numClusters; i++) {
    flag = 1;
    while (flag) {
      random = rand();
      rmax = RAND_MAX;
      index = (int) ((float) (random)/rmax*numVectors);
      index = index;
      flag = (int) FindIndex(indexArray, i-1, index);
      if (!flag) {
	flag = (int) FindMatch(vfv, numVectors, indexArray, i-1, index);
	//	index = (index+5)%numVectors;
	//flag = (int) FindMatch(vfv, numVectors, indexArray, i-1, index);
      }
    }
    for (j = 0; j <vfv[0]->numElements; j++) {
      clusterMeans[i]->array[j] = vfv[index]->array[j];
      clusterVars[i]->array[j] = 1.0;
    }
    indexArray[i] = index;
    
  }
  free(indexArray);
}



/****************************************************************************
 *   Function             : ComputeDiscriminant - computes euclidean distance
 *                        : distance between two vectors
 *   Input args           : clustVector, fvect : input vectors 
 *   Outputs              : ComputeDiscrimnant - distance      	  
 *****************************************************************************/

float ComputeDiscriminant(F_VECTOR *clusterMean, F_VECTOR *clusterVar, F_VECTOR *fvect, int varianceNormalize) {

  int                     i;
  float                   sum = 0;

  for (i = 0; i < fvect->numElements; i++) {
    if (varianceNormalize) 
      sum = sum + (clusterMean->array[i] - fvect->array[i])*
      (clusterMean->array[i] - fvect->array[i])/
      (clusterVar->array[i]);
    else
      sum = sum + (clusterMean->array[i] - fvect->array[i])*
	(clusterMean->array[i] - fvect->array[i]);
  } 
  return(sum);
}

/****************************************************************************
 *   Function             : DecideWhichCluster - determines index of codebook
 *                        : to which a given vector belongs
 *   Input args           : fvect : input vector, clusters - codebook 
 *   Outputs              : DecideWhichCluster - codebook index      	  
 *****************************************************************************/

int DecideWhichCluster(F_VECTOR *fvect, VECTOR_OF_F_VECTORS *clusterMeans,
		       VECTOR_OF_F_VECTORS *clusterVars, 
			 int numClusters, int varianceNormalize) {
  int                     i;
  float                   tempDesc;
  int                     index;
  float                   Discriminant;

  Discriminant  = ComputeDiscriminant(clusterMeans[0], clusterVars[0], 
				      fvect, varianceNormalize);
  index = 0;
  for (i = 1; i < numClusters; i++) {
    tempDesc = ComputeDiscriminant(clusterMeans[i], clusterVars[i], 
				   fvect, varianceNormalize);
    if (tempDesc < Discriminant) {
      index = i;
      Discriminant = tempDesc;
    }
  }
  return(index); 
}

/****************************************************************************
 *   Function             : ComputeVQ - compute codebook
 *                        : for the given set of vectors
 *   Input args           : vfv - input vectors,
 *                        : numVectors - number of input vectors
 *                        : numClusters - codebook size 
 *                        : numIterations - number of iterations
 *   Outputs		  : clusterMeans - array of cluster means
 *			  : clusterVars - array of cluster vars
 *			  : clusterElemCnt - number of elements in
 *			    each cluster
 *****************************************************************************/



void ComputeVQ(VECTOR_OF_F_VECTORS *vfv, int numVectors, 
		VECTOR_OF_F_VECTORS *clusterMeans, 
		VECTOR_OF_F_VECTORS *clusterVars, 
		float *clusterElemCnt, int numClusters, 
	       int varianceNormalize, float ditherMean, 
               int iterations, float seed) {
 
  int                            i,j,k;
  static VECTOR_OF_F_VECTORS     *tempMeanClusters, *tempVarClusters;
  int                            clusterNumber;
  int                            featLength;
  int                            flag;

  featLength = vfv[0]->numElements;
  tempMeanClusters = (VECTOR_OF_F_VECTORS *) calloc (numClusters, sizeof(VECTOR_OF_F_VECTORS));
  tempVarClusters = (VECTOR_OF_F_VECTORS *) calloc (numClusters, sizeof(VECTOR_OF_F_VECTORS));
  for (i = 0; i < numClusters; i++) {
    tempMeanClusters[i] = (F_VECTOR *) AllocFVector(featLength);
    tempVarClusters[i] = (F_VECTOR *) AllocFVector(featLength);
  }
//  printf("temp clusters allocated\n");
  fflush(stdout);
  if (ditherMean == 0) ditherMean = 1.0; 
  InitVQ (vfv, numVectors, clusterMeans, clusterVars, numClusters, seed);

  for ( k = 0; k < iterations; k++) {
    for (i = 0; i < numClusters; i++) {
      for (j = 0; j < featLength; j++) {
	tempMeanClusters[i]->array[j] = 0.0;
	tempVarClusters[i]->array[j] = 0.0;
      }
      clusterElemCnt[i] = 0;
    }
//    printf(" VQ iteration number = %d \n",k); 
    fflush(stdout);
    for (i = 0; i < numVectors; i++) {
      clusterNumber = DecideWhichCluster(vfv[i], 
					 clusterMeans, clusterVars, 
					 numClusters, varianceNormalize);
      clusterElemCnt[clusterNumber] = clusterElemCnt[clusterNumber]+1;
//for(j=0;j<numClusters;j++)
//	printf("%d ",clusterElemCnt[j]);
//printf("\n");
      for (j = 0; j < featLength; j++){
	tempMeanClusters[clusterNumber]->array[j] = 
	  tempMeanClusters[clusterNumber]->array[j] + vfv[i]->array[j];
	//  tempVarClusters[clusterNumber]->array[j] = 
	//  tempVarClusters[clusterNumber]->array[j] + 
	// (vfv[i]->array[j] - clusterMeans[clusterNumber]->array[j]) * 
	//  (vfv[i]->array[j] - clusterMeans[clusterNumber]->array[j]);
      }
   }  

    for (i = 0; i < numClusters; i++) {
       for (j = 0; j < featLength; j++){
	 if (clusterElemCnt[i] != 0)
	   tempMeanClusters[i]->array[j] = 
	     tempMeanClusters[i]->array[j]/clusterElemCnt[i];
	 else
	   tempMeanClusters[i]->array[j] = clusterMeans[i]->array[j];
      }
     }
     
     for (i = 0; i < numVectors; i++) {
       clusterNumber = DecideWhichCluster(vfv[i], clusterMeans, 
					  clusterVars, numClusters, varianceNormalize);
       for (j = 0; j < featLength; j++){
	 tempVarClusters[clusterNumber]->array[j] = 
	   tempVarClusters[clusterNumber]->array[j] + 
	   (vfv[i]->array[j] - tempMeanClusters[clusterNumber]->array[j]) * 
	   (vfv[i]->array[j] - tempMeanClusters[clusterNumber]->array[j]);
       }
     }
//printf("\ncount\n");
//for(i=0;i<numClusters;i++)
//printf("%f ",clusterElemCnt[i]); 
     for (i = 0; i < numClusters; i++){
       flag = ZeroFVector(tempMeanClusters[i]);
       flag = (flag || ZeroFVector(tempVarClusters[i]));
       for (j = 0; j < featLength; j++) {	     
	 if ((clusterElemCnt[i] > 1) && !flag ) {
	   clusterMeans[i]->array[j] = 
	     tempMeanClusters[i]->array[j];
	   clusterVars[i]->array[j] = 
	     tempVarClusters[i]->array[j]/clusterElemCnt[i];
	 } else {
	   clusterElemCnt[i] = 1;
	   clusterVars[i]->array[j] = 0.0001;
	   printf("Warning: ElemCnt=0, Flooring variance\n");
	   printf("mean %d %d = %e var %d %d = %e\n",i,j, 
		  clusterMeans[i]->array[j], i,j, clusterVars[i]->array[j]);
	   fflush(stdout);
	 }
       }
     }
  }
}













/*-------------------------------------------------------------------------
 * $Log: VQ_Modified.c,v $
 * Revision 1.1  2002/12/20 07:38:29  hema
 * Initial revision
 *
 *
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of VQ_Modified.c
 -------------------------------------------------------------------------*/
