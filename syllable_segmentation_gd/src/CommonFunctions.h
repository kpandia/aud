/*
    This file is a header for CommonFunctions.c

    Copyright (C) 2009-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Srikanth Madikeri, Hema A Murthy <hema@cse.iitm.ac.in>

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

#ifndef __COMMON_FUNK_
#define __COMMON_FUNK_

#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "sp/sphere.h"
#include "constants.h"
#include "FrontEndDefs.h"
#include "FrontEndTypes.h"
#include "DspLibrary.h"
#include "InitAsdf.h"
#include "SphereInterface.h"
#include "BatchProcessWaveform.h"
#include "VQ.h"
#include "GMM.h"
#include "math.h"
#include "QuickSort.h"
#include "string.h"
#include<getopt.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
#define MAX_FILE_NAME_LENGTH 512
#define VFV VECTOR_OF_F_VECTORS
#define FVec F_VECTOR
#define finc(i,j,k) for(i=0;i<j;i+=k)
typedef unsigned int uint;

int compare_float (const float*, const float *);
float ComputeThreshold(ASDF *asdf, float threshold);

VECTOR_OF_F_VECTORS *ComputeFeatureVectors(ASDF *asdf, FILE *speakerFile, 
					   char *featureName, unsigned long *numVectors,
					   float thresholdScale); 


VECTOR_OF_F_VECTORS* 
ExtractFeatureVectors (ASDF *asdf, char *wavname,
   			           char *featureName,
							  unsigned long  *numVectors, 
							  float thresholdScale) ;


float 
MixtureProbability(F_VECTOR *fvect, 
					    VECTOR_OF_F_VECTORS *speakerMean,
		             VECTOR_OF_F_VECTORS *speakerVar, 
						 F_VECTOR *speakerWts, 
		             int numClusters, 
						 float probScaleFactor);

float*  
ComputeLikelihood (VECTOR_OF_F_VECTORS *vfv, 
                   unsigned long  numVectors, 
						 int numSpeakers, 
						 VECTOR_OF_F_VECTORS **speakerMeans, 
						 VECTOR_OF_F_VECTORS **speakerVars, 
						 VECTOR_OF_F_VECTORS *speakerWts, 
						 int *numClusters, float *Distortion, 
						 float probScaleFactor) ;

						 
F_VECTOR*
DivideFVectorElements (F_VECTOR *fvec1, F_VECTOR *fvec2, F_VECTOR *fvec3);

F_VECTOR*
ComputeMeanVfv (VECTOR_OF_F_VECTORS *vfv, unsigned long numVectors);

F_VECTOR*
ComputeVarVfv (VECTOR_OF_F_VECTORS *vfv, unsigned long numVectors);

F_VECTOR*
SubtractFVectorElements (F_VECTOR *fvec1, F_VECTOR *fvec2, F_VECTOR *fvec3);



VECTOR_OF_F_VECTORS*
CenterMeanToZero (VECTOR_OF_F_VECTORS *vfv, unsigned long numVectors, F_VECTOR *meanCache);


VECTOR_OF_F_VECTORS*
MakeZeroMeanUnitVar (VECTOR_OF_F_VECTORS *vfv, unsigned long numVectors);


VECTOR_OF_F_VECTORS*
ExtractFeatureVectorsVAD (ASDF *asdf, char *wavname,
                         char *featureName, 
								 unsigned long *numVectors,
								 VECTOR_OF_F_VECTORS *ergMeans,
								 VECTOR_OF_F_VECTORS *ergVars,
								 float *ergWts);
								 
VFV*
ExtractFeatureVectorsFromASR (ASDF * asdf, char *wavname,
			  char *featureName,
			  unsigned long *numVectors,
			  char *vadFileName);
VECTOR_OF_F_VECTORS*
ComputeDeltaVectors (VECTOR_OF_F_VECTORS *vfv,
						  const unsigned long numVectors,
						  const unsigned int deltaDifference);

						  
VECTOR_OF_F_VECTORS*
ComputeAcclVectors (VECTOR_OF_F_VECTORS *vfv,
								  const unsigned int numVectors,
								  const unsigned int deltaDifference,
								  const unsigned int deltaDeltaDifference);

F_VECTOR*
ConcantenateFVector (F_VECTOR *fvect1, F_VECTOR *fvect2);

VECTOR_OF_F_VECTORS*
ConcatVfv (VECTOR_OF_F_VECTORS *vfv1,
		    VECTOR_OF_F_VECTORS *vfv2,
			 unsigned int numVectors);

int
CopyFVector (F_VECTOR *src, F_VECTOR *dest);

F_VECTOR*
CloneFVector (F_VECTOR *src);

VECTOR_OF_F_VECTORS*
ReadVfvFromFile (char *fileName, unsigned int *numVectors);

VECTOR_OF_F_VECTORS*
JoinVfv (VECTOR_OF_F_VECTORS *vfv1, unsigned int nv1,
	VECTOR_OF_F_VECTORS *vfv2, unsigned int nv2,
	unsigned int *numVectors,
	unsigned int cleanUp);

VECTOR_OF_F_VECTORS*
BatchReadVfvFromFile (const char *listFileName, unsigned long *numVectors);

VECTOR_OF_F_VECTORS *
BatchExtractFeatureVectorsVAD (ASDF * asdf, char *batchFileName,
			  char *featureName,
			  unsigned long *numVectors,
			  VECTOR_OF_F_VECTORS * ergMeans,
			  VECTOR_OF_F_VECTORS * ergVars, float *ergWts);

int
ReadTriGaussianFile (char *filename, 
                    VECTOR_OF_F_VECTORS **ergMeans,
                    VECTOR_OF_F_VECTORS **ergVars,
                    float **ergWts);
int
ReadGMMFile (char *filename,
            VECTOR_OF_F_VECTORS **means,
            VECTOR_OF_F_VECTORS **vars,
            float **wts,
	          int nMix,
            int dim);

VECTOR_OF_F_VECTORS*
BatchExtractFeatureVectors (ASDF * asdf, char *fileName,
		       char *featureName,
		       unsigned long *numVectors, float ts);

int
SaveModelToFile (VECTOR_OF_F_VECTORS *ubmMeans,
                VECTOR_OF_F_VECTORS *ubmVars,
                float *ubmWeights,
                int numMix, char *fileName);

VFV* 
Gaussianization (VFV* v, unsigned long nv, uint wnd);

float*
flinspace (float s, float e,uint len);

float erfinv(float);

int
WriteVfvToFile (VFV *vfv, uint nvec, char *filename);

VECTOR_OF_F_VECTORS*
ReadVfvFromBinFile (char *fileName, unsigned int *numVectors);

VECTOR_OF_F_VECTORS*
BatchReadVfvFromBinFile (const char *listFileName, unsigned long *numVectors);
unsigned int*
GetVADDecisionsEnergy (ASDF *asdf, 
                       char *wavname,
                       unsigned long *numVectors, 
                       float thresholdScale);
unsigned int*
GetVADDecisionsTGMM (ASDF *asdf, char *wavname, 
                     unsigned long *numVectors,
                     char *tgmmfilename
                     );
unsigned int*
GetVADDecisionsTrans(ASDF *asdf, 
                     char *transfile,
                     unsigned long *numVectors
                    );

unsigned int*
GetVADDecisions (ASDF *asdf, char *wavname, unsigned int vadtype, 
                 float thresholdScale, char *vadfile, 
                 unsigned long *numVectors);
int
FreeVfv (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors);
int
set_select_channel(char x);
int**
ReadTopCValuesFromFile(char *filename,int cValue,
                       unsigned int numVectors);
VFV*
NormalizeVariance(VFV *vfv, unsigned long numVectors, F_VECTOR *var);
#endif
