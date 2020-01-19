/*
    This file contains a collection of procedures for Gaussian
    Mixture Modelling

    Copyright (C) 2002-2016 Speech and Music Technology Lab,
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


#include "FrontEndDefs.h"
#include "FrontEndTypes.h"
#include "InitAsdf.h"
#include "DspLibrary.h"
#include "stdlib.h"
#include "math.h"
#include "VQ.h"
#include "GMM.h"
/****************************************************************************
 *   Function             : InitGMM - initialises the GMMs with a set
 *                        : of mean vectors and variance vectors
 *   Input args           : seed - seed for randomising
 *                        : vfv - vector of featureVectors
 *                        : numMixtures : number of Vectors
 *   Output args          : mixtureMeans - vector of GMM mean vectors
 *                        : mixtureVars  - vector of GMM variance vectors 
 *****************************************************************************/

void InitGMM (VECTOR_OF_F_VECTORS *vfv,int numVectors, 
              VECTOR_OF_F_VECTORS *mixtureMeans,               
	      VECTOR_OF_F_VECTORS *mixtureVars, 
	     int numMixtures, int seed) {

  int                          index;
  int                          i, j;
  int                          random;
  float                        rmax;


  srand(seed);
  printf("seed = %d\n",seed);
  fflush(stdout);
  for (i = 0; i < numMixtures; i++) {
      random = rand();
      rmax = RAND_MAX;
      index = (int) ((float) (random/rmax*numVectors));
     for (j = 0; j <vfv[0]->numElements; j++)
      mixtureMeans[i]->array[j] = vfv[index]->array[j];
    for (j = 0; j <vfv[0]->numElements; j++)
      mixtureVars[i]->array[j] = 1.0;
  }
}

/****************************************************************************
 *   Function             : ComputeProbability - computes euclidean distance
 *                        : distance between two vectors
 *   Input args           : mixtureMean, mixtureVar, priorProb, 
 *                        : probScaleFactor, fvect : input vector 
 *   Outputs              : ComputeDiscrimnant - distance      	  
 *****************************************************************************/

float ComputeProbability(F_VECTOR *mixtureMean, 
			  F_VECTOR *mixtureVar, float priorProb, 
			 F_VECTOR *fvect, float probScaleFactor) {

  int                     i;

  float                   sumProb = 0;
  float                   scale = 0.0;
  float                   floorValue = 0.0;

  // if (prevMean != mixtureMean) {
    for (i = 0; i < fvect->numElements; i++)
      {
	if (mixtureVar->array[i] != 0.0)
	  scale = scale + log(mixtureVar->array[i]);
      } /*  for (..i < fvect->numElements..)  */
    scale = 0.5*(fvect->numElements*log(2.0*PI) + scale);
    scale = log(priorProb) - scale;
    //    prevMean = mixtureMean;
  /*    printf("scale = %f \n", scale);
	scanf("%*c"); */
    //}

sumProb = scale;
for (i = 0; i < fvect->numElements; i++) {
  floorValue =  floorValue +
    (mixtureMean->array[i] - fvect->array[i])*
    (mixtureMean->array[i] - fvect->array[i])
    /(2*mixtureVar->array[i]);
}
  sumProb = sumProb - floorValue; //fvect->numElements;
  return(sumProb + log(probScaleFactor));
}  

/****************************************************************************
 *   Function             : ComputeMixtureContribution - determines index of Mixture
 *                        : to which a given vector belongs
 *   Input args           : fvect : input vector to be classified
 *                                  mixtureMeans,  mixtureVars,numMixtures,
 *                                  mixtureWeights, numVectors, probScaleFactor
 *   Outputs              : the contribution for fvect to each mixture      	  
 *****************************************************************************/

float *ComputeMixtureContribution(F_VECTOR *fvect, 
                      VECTOR_OF_F_VECTORS *mixtureMeans, 
                      VECTOR_OF_F_VECTORS *mixtureVars, 
                      int numMixtures, float *mixtureWeights, 
		      float probScaleFactor, float *mixtureContribution) {
  int                 i;
  float               evidence;

  for (i = 0; i < numMixtures; i++)
    if (mixtureWeights[i] == 0)
      mixtureWeights[i] = 1.0E-50;

  mixtureContribution[0] = ComputeProbability (mixtureMeans[0], mixtureVars[0], 
					       mixtureWeights[0], fvect, probScaleFactor);
  evidence = mixtureContribution[0];
  for (i = 1; i < numMixtures; i++) {
      mixtureContribution[i] = ComputeProbability(mixtureMeans[i], 
						  mixtureVars[i], mixtureWeights[i], 
						  fvect, probScaleFactor);
      evidence = LogAdd(evidence, mixtureContribution[i]); 
  }
  for (i = 0; i < numMixtures; i++) {
    mixtureContribution[i] = mixtureContribution[i] - evidence;
    if (mixtureContribution [i] < -100) mixtureContribution[i] = -100;
  }
//for (i = 0; i < numMixtures; i++)
//printf("%f ",mixtureContribution[i]);
  return(mixtureContribution);
}

/****************************************************************************
 *   Function             : ComputeGMM - compute GMMs
 *                        : for the given set of vectors
 *   Input args           : vfv - input vectors,
 *                        : numVectors - number of input vectors
 *                        : numMixtures - codebook size 
 *                        : VQIter - number of VQ iterations
 *                        : GMMIter - number of GMM iterations
 *   Outputs		  : mixtureMeans - array of mixture means
 *			  : mixtureVars - array of mixture vars
 *			  : mixtureWeights - number of elements in
 *			    each mixture
 *****************************************************************************/



void ComputeGMM(VECTOR_OF_F_VECTORS *vfv, int numVectors, 
		VECTOR_OF_F_VECTORS *mixtureMeans, 
		VECTOR_OF_F_VECTORS *mixtureVars, 
		float *mixtureWeights, int numMixtures, 
		int VQIter, int GMMIter, float probScaleFactor,
                int ditherMean, int varianceNormalize, float varianceFloor, int seed) {
probScaleFactor=1;
//pass by value fails here -- CHECK
//printf("nmix = %d no. of vect = %d scale factor = %f\n", numMixtures,numVectors,probScaleFactor);
  int                            i,j,k;
  static VECTOR_OF_F_VECTORS     *tempMeans, *tempVars;
  static float                   *tempMixtureWeights;
  float                          *mixtureContribution;
  float                          expValue;
  int                            mixtureNumber;
  int                            featLength;
  int                            flag;
  //int                            total;
  //int                            minIndex, maxIndex;
  //int                            minMixtureSize, maxMixtureSize;
  

  featLength = vfv[0]->numElements;
  tempMeans = (VECTOR_OF_F_VECTORS *) calloc (numMixtures, 
					      sizeof(VECTOR_OF_F_VECTORS));
  tempVars = (VECTOR_OF_F_VECTORS *) calloc (numMixtures, 
					      sizeof(VECTOR_OF_F_VECTORS));
  mixtureContribution = (float *) AllocFloatArray(mixtureContribution, numMixtures);
  for (i = 0; i < numMixtures; i++) { 
    tempMeans[i] = (F_VECTOR *) AllocFVector(featLength);
    tempVars[i] = (F_VECTOR *) AllocFVector(featLength);
  }
  tempMixtureWeights = (float *) AllocFloatArray(tempMixtureWeights, 
						 numMixtures);
//  printf("temp mixtures allocated\n");
//  fflush(stdout);
  //InitGMM (vfv, numVectors, mixtureMeans, mixtureVars, numMixtures, seed);
  ComputeVQ(vfv, numVectors, mixtureMeans, mixtureVars,
      mixtureWeights, numMixtures, varianceNormalize, 
      ditherMean, VQIter, seed);
  for (i = 0; i < numMixtures; i++) 
    mixtureWeights[i] = mixtureWeights[i]/numVectors;

  for ( k = 0; k < GMMIter; k++) {
    printf(" GMM iteration number = %d\n", k);
    fflush(stdout);
    for (i = 0; i < numMixtures; i++) {
      for (j = 0; j < vfv[0]->numElements; j++) {
	tempMeans[i]->array[j] = 0;
	tempVars[i]->array[j] = 0;
      }
      tempMixtureWeights[i] = 0;
    }

    /* Compute the mean of each mixture */
    for (i = 0; i < numVectors; i++) {
      mixtureContribution = (float *) ComputeMixtureContribution(vfv[i], 
								 mixtureMeans, 
								 mixtureVars, numMixtures, 
								 mixtureWeights,
								 probScaleFactor, mixtureContribution);
//printf("%f ",mixtureContribution);
      for (mixtureNumber = 0; mixtureNumber < numMixtures; mixtureNumber++){
	 expValue = expf(mixtureContribution[mixtureNumber]);
//printf("%f %f ",mixtureContribution[mixtureNumber],expValue);
        if (isnan(expValue)) expValue = 0;
	  tempMixtureWeights[mixtureNumber] = tempMixtureWeights[mixtureNumber]+expValue;
	for (j = 0; j < featLength; j++){
	  tempMeans[mixtureNumber]->array[j] = 
	    tempMeans[mixtureNumber]->array[j] + 
	    vfv[i]->array[j]*expValue;
	}
      }
    }  
    
    for (i = 0; i < numMixtures; i++) 
      for (j = 0; j < vfv[0]->numElements; j++) {
      tempMeans[i]->array[j] = tempMeans[i]->array[j]/tempMixtureWeights[i];
    }

    /* Compute the variance and weights of each mixture */

      
    for (i = 0; i < numVectors; i++) {
      mixtureContribution = (float *) ComputeMixtureContribution(vfv[i], 
								 mixtureMeans, 
								 mixtureVars, numMixtures, 
								 mixtureWeights,
								 probScaleFactor, 
								 mixtureContribution);
     for (mixtureNumber = 0; mixtureNumber < numMixtures; mixtureNumber++) {
       expValue = expf(mixtureContribution[mixtureNumber]);
       if (isnan(expValue)) expValue = 0;
 	for (j = 0; j < featLength; j++){
	  tempVars[mixtureNumber]->array[j] = 
	    tempVars[mixtureNumber]->array[j] + 
	  (vfv[i]->array[j] - tempMeans[mixtureNumber]->array[j]) * 
	  (vfv[i]->array[j] - tempMeans[mixtureNumber]->array[j])
	    *expValue;
	}
      }
    }
    for (i = 0; i < numMixtures; i++) {
      for (j = 0; j < vfv[0]->numElements; j++)
	tempVars[i]->array[j] = tempVars[i]->array[j]/tempMixtureWeights[i];
      tempMixtureWeights[i] = tempMixtureWeights[i]/numVectors;
      

//      	for (i = 0; i < numMixtures; i++)
//      printf("mixWeights %d = %f\n", i, tempMixtureWeights[i]);
//      scanf("%*c");
      }

    /* Floor variance, adjust means if necessary */

    for (i = 0; i < numMixtures; i++){
      flag = ZeroFVector(tempMeans[i]);
      flag = (flag || ZeroFVector(tempVars[i]));
      if (!flag){
	for (j = 0; j < featLength; j++) {
	  mixtureMeans[i]->array[j] = tempMeans[i]->array[j];
	  mixtureVars[i]->array[j] = tempVars[i]->array[j];
	  if (mixtureVars[i]->array[j] < varianceFloor) 
	    mixtureVars[i]->array[j] = varianceFloor;
	}
	mixtureWeights[i] = tempMixtureWeights[i];
      } else {
	//	     mixtureMeans[i]->array[j] = tempMeans[i]->array[j];
	for (j = 0; j < featLength; j++)
	  mixtureVars[i]->array[j] = varianceFloor;
	mixtureWeights[i] = 1.0/(float)numVectors;
	printf("Flooring Mixture %d %d mean = %e var = %e mixWeights = %e\n",i,j, 
	       mixtureMeans[i]->array[j],
	       mixtureVars[i]->array[j], mixtureWeights[i]);
	fflush(stdout);
	  
      }
    }
  }

      for (i = 0; i < numMixtures; i++)
    for (j = 0; j < featLength; j++) {
	printf("mean %d %d = %f var %d %d = %f\n",i,j, 
	       mixtureMeans[i]->array[j], i, j, mixtureVars[i]->array[j]);
	fflush(stdout);
	}

}
