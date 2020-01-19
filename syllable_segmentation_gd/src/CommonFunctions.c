/*
    This file contains the functions required to compute feature vectors

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


#include "CommonFunctions.h"


float
ComputeThreshold (ASDF * asdf, float threshold)
{
  int i;
  float ave = 0;
  if (threshold <= 0.0) return 0.0;
  FVec* fv;
  for (i = 0; i < asdf->numFrames; i++) {
    fv = (F_VECTOR *) GsfRead (asdf, i, "frameEnergy");
    ave += fv->array[0];
    FreeFVector (fv);
  }
  ave = ave / asdf->numFrames;
  return (ave * threshold);
}


/****************************************************
 * Function:	ComputeFeatureVectors
 * Details:		Given a file that contains a list
 * 			   of wave files, this func zero-centers
 * 			   and variance normalizes 
 * 			   each wave file individually and then
 * 			   returns the concatenated vfvs  			   
 ***************************************************/


VECTOR_OF_F_VECTORS *
ComputeFeatureVectors (ASDF * asdf, FILE * speakerFile,
		       char *featureName, unsigned long *numVectors,
		       float thresholdScale)
{
  F_VECTOR *fvect, *evec;
  VFV *vfv, *nvfv;
  char line[500];
  unsigned long totalFrames;
  int i, j, k, s, l, gflag,gwnd;
  int frameNo, prevFrameNo;
  char wavname[200];
  float energy, threshold;
  char select_channel_option;

  totalFrames = 0;
  rewind (speakerFile);
  while (fgets (line, 500, speakerFile) != NULL)
    {
      sscanf (line, "%s %c", wavname, &select_channel_option);
      printf ("file read :%s\n", wavname);
      fflush (stdout);
      set_select_channel (select_channel_option);
      GsfOpen (asdf, wavname);
      totalFrames = totalFrames + asdf->numFrames;
      printf
	("[DEBUG]wavname = %s totFrames =  %lu numSamp =  %lu numFrames = %u\n",
	 wavname, totalFrames, asdf->numSamples, asdf->numFrames);
      fflush (stdout);
      GsfClose (asdf);
    }
  printf ("total no frames = %lu\n", totalFrames);
  fflush (stdout);
  vfv =
    (VFV*) malloc (totalFrames * sizeof (VFV));
  if (vfv == NULL)
    {
      printf ("unable to allocate space for vfv \n");
      exit (-1);
    } 
  rewind (speakerFile);
  frameNo = 0;
  prevFrameNo = -1;
  gflag = GetIAttribute (asdf,"stGauss");
  if (gflag) 
    gwnd = GetIAttribute (asdf, "stGaussWnd");
  
  while (EOF != fscanf (speakerFile, "%s %c", wavname, &select_channel_option))
    {
      GsfOpen (asdf, wavname);
      printf ("numFrames = %d\n", asdf->numFrames);
      threshold = (float) ComputeThreshold (asdf, thresholdScale);
      for (i = 0; i < asdf->numFrames; i++)
      {
        evec = (F_VECTOR *) GsfRead (asdf, i, "frameEnergy");
        energy = evec->array[0];
        if (energy >= threshold)
          {
            fvect = (F_VECTOR *) GsfRead (asdf, i, featureName);
            if (fvect == NULL)
              {
                printf ("problems fvect\n");
                fflush (stdout);
                exit (-1);
              }
            vfv[frameNo] = fvect;
            frameNo++;
          }
        FreeFVector (evec);
      }
      if (GetIAttribute (asdf, "zeroMean") == 0)
        {
          printf ("zeroMean is set to 0\n");
          s = prevFrameNo + 1;
          l = frameNo - s;
          MakeZeroMeanUnitVar (&vfv[s], l);          
          prevFrameNo = frameNo - 1;
        }
      GsfClose (asdf);
    }
  *numVectors = frameNo;
  printf ("total no of frame processed = %ld\n", *numVectors);
  return (vfv);
}

VFV* ExtractFeatureVectors (ASDF * asdf, char *wavname,
		       char *featureName,
		       unsigned long *numVectors, float thresholdScale) {

  F_VECTOR *fvect, *evec = NULL, *tempfvec = NULL;
  VECTOR_OF_F_VECTORS *vfv, *nvfv, *tempvfv;
  int i,gwnd, dovad;
  uint gflag = 0;
  unsigned long frameNum;
  float threshold, energy;

  GsfOpen (asdf, wavname);
  dovad = 0;
  gflag = (int) GetIAttribute (asdf,"stGauss");
  if (gflag) {
    gwnd = GetIAttribute (asdf,"stGaussWnd");

//    printf ("-----------------\n");
//    printf ("%d %d\n", gflag, gwnd);
//    printf ("-----------------\n");

  }

  if (!dovad) {
      vfv = (VFV *) calloc (asdf->numFrames,sizeof (VFV));
      threshold = (float) ComputeThreshold (asdf, thresholdScale);
      frameNum = 0;
//      printf ("DEBUG: thresholdInExtFeat = %f\n", threshold);
      evec = AllocFVector(1);
      finc(i,asdf->numFrames,1) {
          evec = FrameComputeEnergy(asdf,i,evec);
        energy = evec->array[0];
        if (energy >= threshold) {
          fvect = (F_VECTOR *) GsfRead (asdf, i, featureName);
            if (fvect == NULL) {
              printf ("problems fvect\n");
              fflush (stdout);
              exit (-1);
            }
            vfv[frameNum++] = fvect;
        }      
      }
      FreeFVector (evec);
  }
  else {
      vfv = (VFV *) calloc (asdf->numVoicedFrames, sizeof (VFV));
      frameNum = 0;
      finc (i, asdf->numFrames, 1)  
        if ((tempfvec = GsfRead (asdf, i, featureName)) != NULL)
            vfv[frameNum++] = tempfvec ;
  }
  //GsfClose (asdf);  
  //MakeZeroMeanUnitVar (vfv, frameNum);
  if (gflag) {
//    printf ("Calling Gaussianization\n");
    nvfv = Gaussianization (vfv,frameNum,gwnd);
    vfv = nvfv;
  }
// TODO: free vfv  
  *numVectors = frameNum;
  return (vfv);
}



float
MixtureProbability (F_VECTOR * fvect,
		    VECTOR_OF_F_VECTORS * speakerMean,
		    VECTOR_OF_F_VECTORS * speakerVar,
		    F_VECTOR * speakerWts,
		    int numClusters, float probScaleFactor)
{

  int i;
  float maxProb = 0, probValue;
  int featLength;

  featLength = fvect->numElements;
  maxProb = ComputeProbability (speakerMean[0],
				speakerVar[0],
				speakerWts->array[0], fvect, probScaleFactor);
  for (i = 1; i < numClusters; i++)
    {
      probValue = ComputeProbability (speakerMean[i],
				      speakerVar[i],
				      speakerWts->array[i],
				      fvect, probScaleFactor);
      maxProb = LogAdd (maxProb, probValue);
    }

  return (maxProb);

}


float *
ComputeLikelihood (VECTOR_OF_F_VECTORS * vfv,
		   unsigned long numVectors,
		   int numSpeakers,
		   VECTOR_OF_F_VECTORS ** speakerMeans,
		   VECTOR_OF_F_VECTORS ** speakerVars,
		   VECTOR_OF_F_VECTORS * speakerWts,
		   int *numClusters, float *Distortion, float probScaleFactor)
{
  int i, j;
  float mixProbValue;
  float evidence;

  for (i = 0; i < numSpeakers; i++)
    {
      Distortion[i] = 0;
      for (j = 0; j < numVectors; j++)
	{
	  mixProbValue = MixtureProbability (vfv[j], speakerMeans[i],
					     speakerVars[i],
					     speakerWts[i],
					     numClusters[i], probScaleFactor);
	  Distortion[i] = Distortion[i] + mixProbValue;
	}
      printf ("Distortion %d = %f\n", i, Distortion[i]);
      fflush (stdout);
      Distortion[i] = Distortion[i] / numVectors;
    }
  evidence = Distortion[0];
  for (i = 1; i < numSpeakers; i++)
    evidence = LogAdd (evidence, Distortion[i]);
  for (i = 0; i < numSpeakers; i++)
    Distortion[i] = Distortion[i] - evidence;

  return (Distortion);
}


F_VECTOR *
ComputeMeanVfv (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors)
{
  unsigned int i, j, d;
  F_VECTOR *meanFv;

  if (numVectors == 0)
    return NULL;


  d = vfv[0]->numElements;
  meanFv = (F_VECTOR *) AllocFVector (d);
  for (j = 0; j < d; j++)
    meanFv->array[j] = 0.0;

  for (i = 0; i < numVectors; i++)
    {
      for (j = 0; j < d; j++)
	{
	  meanFv->array[j] = meanFv->array[j] + vfv[i]->array[j];
	}
    }

  for (j = 0; j < d; j++)
    meanFv->array[j] = meanFv->array[j] / numVectors;
  return meanFv;
}

F_VECTOR *
ComputeVarVfv (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors)
{
  unsigned int i, j, d;
  F_VECTOR *meanFv, *varFv;

  if (numVectors == 0)
    return NULL;

  d = vfv[0]->numElements;
  meanFv = ComputeMeanVfv (vfv, numVectors);
  varFv = (F_VECTOR *) AllocFVector (d);
  for (j = 0; j < d; j++)
    varFv->array[j] = 0.0;

  for (i = 0; i < numVectors; i++)
    {
      for (j = 0; j < d; j++)
	{
	  varFv->array[j] =
	    varFv->array[j] +
	    ((vfv[i]->array[j] - meanFv->array[j]) * (vfv[i]->array[j] -
						      meanFv->array[j]));
	}
    }

  for (j = 0; j < d; j++)
    varFv->array[j] = sqrt (varFv->array[j]) / numVectors;
  free(meanFv->array);
  free(meanFv);
  return varFv;
}

F_VECTOR *
SubtractFVectorElements (F_VECTOR * fvec1, F_VECTOR * fvec2, F_VECTOR * fvec3)
{
  int i, d;
  d = fvec1->numElements;

  for (i = 0; i < d; i++)
    fvec3->array[i] = fvec1->array[i] - fvec2->array[i];

  return fvec3;
}

F_VECTOR *
DivideFVectorElements (F_VECTOR * fvec1, F_VECTOR * fvec2, F_VECTOR * fvec3)
{
  int i, d;
  float floor_value = 1E-5;

  d = fvec1->numElements;

  for (i = 0; i < d; i++)
    {
      if (fvec2->array[i] == 0.0)
          fvec2->array[i] = floor_value ;
      fvec3->array[i] = fvec1->array[i] / fvec2->array[i];
    }
  return fvec3;
}

VECTOR_OF_F_VECTORS *
CenterMeanToZero (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors,
		  F_VECTOR * meanVectorCache)
{
  unsigned int i, j;
  F_VECTOR *meanVector;

  if (numVectors == 0 || vfv == NULL)
    return NULL;

  if (meanVectorCache == NULL)
    meanVector = ComputeMeanVfv (vfv, numVectors);
  else
    meanVector = meanVectorCache;

  for (i = 0; i < numVectors; i++)
    SubtractFVectorElements (vfv[i], meanVector, vfv[i]);

  return vfv;
}

VECTOR_OF_F_VECTORS *
MakeZeroMeanUnitVar (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors)
{
  unsigned int i, j;
  F_VECTOR *meanVector, *varVector;

  meanVector = ComputeMeanVfv (vfv, numVectors);
  varVector = ComputeVarVfv (vfv, numVectors);

  CenterMeanToZero (vfv, numVectors, meanVector);
  for (i = 0; i < numVectors; i++)         
    DivideFVectorElements (vfv[i], varVector, vfv[i]);
  FreeFVector (meanVector);
  FreeFVector (varVector);
  return vfv;
}

VFV*
NormalizeVariance(VFV *vfv, unsigned long numVectors, F_VECTOR *var) 
{
    int i;
    if (var == NULL) {
         var = (F_VECTOR *) ComputeVarVfv(vfv, numVectors);
    }
    for (i = 0; i < numVectors; i++)
      DivideFVectorElements (vfv[i], var, vfv[i]);
    return vfv;
}


unsigned int*
GetVADDecisionsEnergy (ASDF *asdf, 
                       char *wavname,
                       unsigned long *numVectors, 
                       float thresholdScale) {
    F_VECTOR *energyVector;
    int zm;
    unsigned long numSpeechVectors;
    int i,j;
    unsigned int *idx;
    float probScaleFactor, thresholdValue;      
    
    GsfOpen (asdf, wavname);
    numSpeechVectors = asdf->numFrames;
    printf ("numSpeechVectors = %lu", numSpeechVectors);
    zm = asdf->zeroMean;
    asdf->zeroMean = 0;
    probScaleFactor = (float) GetFAttribute (asdf, "probScaleFactor");
    thresholdValue = ComputeThreshold (asdf, thresholdScale);
    i = -1;
    idx = (unsigned int *) calloc (numSpeechVectors, sizeof(unsigned int));
    energyVector = (F_VECTOR*) AllocFVector (1);
    if (energyVector == NULL) {
        printf ("Unable to allocate FVector\n");
        exit(0);
    }
    finc(j,numSpeechVectors,1) {
        energyVector = FrameComputeEnergy (asdf, j,energyVector);    
        if (energyVector->array[0] >= thresholdValue) idx[++i] = j;
    }
    FreeFVector (energyVector);
    GsfClose (asdf);
    asdf->zeroMean = zm;
    numSpeechVectors = i+1;
    printf ("numSpeechVectors = %lu\n", numSpeechVectors);
    *numVectors = numSpeechVectors;
    idx = realloc (idx,numSpeechVectors*sizeof(unsigned int));
    return idx;
} 

unsigned int*
GetVADDecisionsTGMM (ASDF *asdf, char *wavname, 
                     unsigned long *numVectors,
                     char *tgmmfilename
                     ) {
    VFV *energyVector;
    VFV *ergMeans, *ergVars;
    float *ergWts, threshold;
    int  zm, i,j;
    unsigned long numSpeechVectors;
    int silMixNo = 0, speechMixNo = 2 , noiseMixNo = 1;
    int *idx;
    float temp_pr[3], probScaleFactor;

    ReadTriGaussianFile(tgmmfilename, &ergMeans, &ergVars, &ergWts);
    GsfOpen (asdf, wavname);
    numSpeechVectors = asdf->numFrames;
    zm = asdf->zeroMean;
    asdf->zeroMean = 0;
    probScaleFactor = (float) GetFAttribute (asdf, "probScaleFactor");
    energyVector = (VFV*) calloc (numSpeechVectors, sizeof(VFV));
    finc(j,numSpeechVectors,1) {
      energyVector[j] = (F_VECTOR*) AllocFVector (1);
      if (energyVector[j] == NULL) {
        printf ("Unable to allocate FVector\n");
        exit(0);
      }
      energyVector[j] = FrameComputeLogEnergy (asdf, j,energyVector[j]);    
    }
    GsfClose (asdf);
    MakeZeroMeanUnitVar (energyVector,numSpeechVectors);
        
    idx = (unsigned int*) calloc (numSpeechVectors,sizeof(unsigned int));
    numSpeechVectors = 0;
    threshold = ergMeans[2]->array[0] - sqrt(ergVars[2]->array[0]);
    finc(i,asdf->numFrames,1) {
        if(energyVector[i]->array[0] > threshold)
               idx[numSpeechVectors++] = i;
    }
    idx = realloc (idx,numSpeechVectors * sizeof(unsigned int));
    *numVectors = numSpeechVectors;
    asdf->zeroMean = zm;
    return idx;
}

unsigned int*
GetVADDecisionsTrans(ASDF *asdf, 
                     char *transfile,
                     unsigned long *numVectors
                    ) {
    FILE *vadFile = NULL;
    unsigned int lframeno, rframeno, *idx = NULL;       
    int samplingRate = asdf->samplingRate, 
        frameRate = asdf->frameAdvanceSamples;
    float ratio = ((float) samplingRate)/frameRate;
    float lval, rval;
    unsigned long numSpeechVectors;
    int i,j;

    vadFile = fopen (transfile, "r");
    if (vadFile == NULL) return NULL;
    numSpeechVectors = 0;
    i = -1;
    idx = (unsigned int*) malloc (sizeof(unsigned int)*asdf->numFrames);
    while (!feof (vadFile)) {
        fscanf (vadFile, "%f %f\n", &lval,&rval);
        lframeno =  ceil (lval * ratio);
        rframeno =  ceil (rval * ratio);
        for (j=lframeno;j<=rframeno;j++) idx[++i] = j;
    }
    numSpeechVectors = i+1;
    *numVectors = numSpeechVectors;
    idx = realloc (idx,numSpeechVectors*sizeof(unsigned int));
    fclose (vadFile);
    return idx;
}

unsigned int*
GetVADDecisions (ASDF *asdf, char *wavname, unsigned int vadtype, 
                 float thresholdScale, char *vadfile, 
                 unsigned long *numVectors) {
    if (!vadtype || vadtype > 3) return NULL;
    else if (vadtype == 1) return GetVADDecisionsEnergy(asdf,
                                                        wavname,
                                                        numVectors,
                                                        thresholdScale);
    else if (vadtype == 2) return GetVADDecisionsTGMM(asdf,
                                                      wavname,
                                                      numVectors,
                                                      vadfile
                                                     );
    else return GetVADDecisionsTrans (asdf,vadfile,numVectors);
                                                      
}

VECTOR_OF_F_VECTORS *
ExtractFeatureVectorsVAD (ASDF * asdf, char *wavname,
			  char *featureName,
			  unsigned long *numVectors,
			  VECTOR_OF_F_VECTORS * ergMeans,
			  VECTOR_OF_F_VECTORS * ergVars, float *ergWts)
{

  VECTOR_OF_F_VECTORS *energyVector, *vfv, *newVfv, *nvfv;
  int numSpeechVectors, speechMixtureNumber, silenceMixtureNumber,
    noiseMixtureNumber;
  int i, j, gflag, gwnd, zm;
  float tempProbabilities[3], probScaleFactor;
  float thresholdScale = 0.0;

  probScaleFactor = (float) GetFAttribute (asdf, "probScaleFactor");
  
  GsfOpen (asdf, wavname);
  numSpeechVectors = asdf->numFrames;
  zm = asdf->zeroMean;
  // zero mean is ignored assuming the input utterance is very
  // long. this means, if filter bank energies were used in 
  // front-end-dsp, they would have been normalized w.r.t the whole
  // utterance. this would prove detrimental. For mfcc, unsetting zero
  // mean did not change much of the input
  //
  asdf->zeroMean = 0;
  energyVector = (VFV*) calloc (numSpeechVectors, sizeof(VFV));
  finc(j,numSpeechVectors,1) {
    energyVector[j] = (F_VECTOR*) AllocFVector (1);
    if (energyVector[j] == NULL) {
      printf ("Unable to allocate FVector\n");
      exit(0);
    }
    energyVector[j] = FrameComputeLogEnergy (asdf, j,energyVector[j]);    
  }
  
  MakeZeroMeanUnitVar (energyVector, numSpeechVectors);
  asdf->zeroMean = zm;
  gflag = (int) GetIAttribute(asdf,"stGauss");
  //PutIAttribute(asdf,"stGauss",0);
  vfv = (VECTOR_OF_F_VECTORS *)
    ExtractFeatureVectors (asdf, wavname,
			   featureName, numVectors, thresholdScale);
  GsfClose(asdf);
  numSpeechVectors = 0;
  // Pre-allocating space. Includes extra space. Should be freed
  // by caller
  newVfv =
    (VECTOR_OF_F_VECTORS *) calloc (*numVectors,
				    sizeof (VECTOR_OF_F_VECTORS));

  silenceMixtureNumber = 0;
  noiseMixtureNumber = 1;
  speechMixtureNumber = 2;
  finc(i,*numVectors,1) {
      if (energyVector[i]->array[0] > ergMeans[noiseMixtureNumber]->array[0])
              newVfv[numSpeechVectors++] = vfv[i];

    }
  printf ("numSpeechVectors is %d while numVectors is %lu\n", numSpeechVectors,
	  *numVectors);
  FreeVfv(energyVector,*numVectors);
  *numVectors = numSpeechVectors;

  if (gflag) {
    gwnd = GetIAttribute (asdf,"stGaussWnd");
    nvfv = Gaussianization (newVfv,numSpeechVectors,gwnd);
    FreeVfv (newVfv, numSpeechVectors);    
    newVfv = nvfv;
  }
  return newVfv;
}

VFV* ExtractFeatureVectorsFromASR (ASDF * asdf, char *wavname,
                                   char *featureName,
                                   unsigned long *numVectors,
                                   char *vadFileName) {
    VFV *vfv, *newVfv, *nvfv;
    FILE *vadFile = NULL;
    int i, j, gflag, gwnd, zm, nboundaries, idx, k;
    int numSpeechVectors, samplingRate, frameRate;
    float probScaleFactor, *endPoints = NULL, tempf;
    float thresholdScale = 0.0;

    probScaleFactor = (float) GetFAttribute (asdf, "probScaleFactor");
    samplingRate = (int) asdf->samplingRate;
    frameRate = (int) asdf->frameAdvanceSamples;
    
    gflag = GetIAttribute(asdf, "stGauss");
    asdf->stGauss = 0;
    vfv = (VFV *)
      ExtractFeatureVectors (asdf, wavname,
           featureName, numVectors, thresholdScale);
    
    vadFile = fopen (vadFileName,"r");
    if (vadFile == NULL) return NULL;
    nboundaries = 0;
    while (!feof (vadFile)) {
      fscanf (vadFile, "%f %f\n", &tempf,&tempf);
      nboundaries++;
    }
    rewind (vadFile);
    endPoints = (float *) calloc (nboundaries * 2, sizeof(float));
    numSpeechVectors = 0;
    for (i = 0; i < nboundaries; i++) {
        j = i*2;
        fscanf (vadFile, "%f %f\n", &endPoints[j],&endPoints[j+1]);      
        endPoints[j] = (float) ceil (endPoints[j] * samplingRate / frameRate);
        endPoints[j+1] = (float) ceil (endPoints[j+1] * samplingRate / frameRate);
        numSpeechVectors += (endPoints[j+1] - endPoints[j] + 1);
    }
    fclose (vadFile);
    // Pre-allocating space. Includes extra space. Should be freed
    // by caller
    newVfv =
      (VECTOR_OF_F_VECTORS *) calloc (numSpeechVectors,
              sizeof (VECTOR_OF_F_VECTORS));
    printf ("numSpeechVectors = %d\n", numSpeechVectors);
    idx = 0;
    for (i = 0; i < nboundaries; i++) {
        j = i*2;
        for (k = endPoints[j]; k <= endPoints[j+1]; k++) 
            newVfv[idx++] = vfv[(int)k];
    }
    // at this point idx == numSpeechVectors. verified!
    *numVectors = numSpeechVectors;

    asdf->stGauss = gflag;
    if (gflag) {
        gwnd = GetIAttribute (asdf,"stGaussWnd");
        nvfv = Gaussianization (newVfv,(unsigned long) numSpeechVectors,gwnd);
        FreeVfv (newVfv, numSpeechVectors);    
        newVfv = nvfv;
    }
    return newVfv;
}

VFV* ComputeDeltaVectors (VECTOR_OF_F_VECTORS * vfv,
                          const unsigned long numVectors,
                          const unsigned int deltaDifference)
{
    int i, j, frameIndex, featLength;
    VECTOR_OF_F_VECTORS *deltaVfv;
    F_VECTOR *fvect, *temp, *prev, *next;
    float normalizingConst = 0.0;

    deltaVfv = (VFV *) calloc (numVectors, sizeof (VFV));
    if (deltaVfv == NULL)
    {
        printf ("Unable to allocate memory\n");
        return NULL;
    }
    for (i=1; i <= deltaDifference; i++)
        normalizingConst = normalizingConst + (i * i * 2);

    featLength = vfv[0]->numElements;
    for (j = 0; j < numVectors; j++)
    {
        prev = (F_VECTOR *) AllocFVector (featLength);
        next = (F_VECTOR *) AllocFVector (featLength);
        temp = (F_VECTOR *) AllocFVector (featLength);
        fvect = (F_VECTOR *) AllocFVector (featLength);

        for (i = 0; i < featLength; i++)
            fvect->array[i] = 0.0;
        for (i = 1; i <= deltaDifference; i++) 
        {
            if (j - i >= 0)
                CopyFVector (vfv[j-i], prev);
        //	    prev = vfv[j - i];
            else
                InitFVector (prev);
            if (j + i < numVectors)
                CopyFVector (vfv[j + i], next);
        //	    next = vfv[j + i];
            else
                InitFVector (next);
            LinearVectorDifference (prev, next, temp);
            LinearVectorScalarMultiply ((float) (i), temp, temp);
            LinearVectorAddition (temp, fvect, fvect);
        }

        LinearVectorScalarDivide (normalizingConst, fvect, fvect);
        deltaVfv[j] = (F_VECTOR *) fvect;
        FreeFVector (prev);
        FreeFVector (next);
    }
    return deltaVfv;
}

VECTOR_OF_F_VECTORS *
ComputeAcclVectors (VECTOR_OF_F_VECTORS * vfv,
		    const unsigned int numVectors,
		    const unsigned int deltaDifference,
		    const unsigned int deltaDeltaDifference)
{
  VECTOR_OF_F_VECTORS *deltaVfv, *acclVfv;
  int i;
  deltaVfv = ComputeDeltaVectors (vfv, numVectors, deltaDifference);

  acclVfv = ComputeDeltaVectors (deltaVfv, numVectors, deltaDeltaDifference);

  for (i = 0; i < numVectors; i++)
    {
      free (deltaVfv[i]->array);
      free (deltaVfv[i]);
    }
  free (deltaVfv);
}

F_VECTOR *
ConcantenateFVector (F_VECTOR * fvect1, F_VECTOR * fvect2)
{
  F_VECTOR *newFVector;
  int totalElements;
  int numElements2 = fvect2->numElements;
  int i,j;

  totalElements = fvect1->numElements + fvect2->numElements;
  newFVector = (F_VECTOR *) AllocFVector (totalElements);
  if (newFVector == NULL)
    return NULL;

  for (i = 0; i < fvect1->numElements; i++)
    newFVector->array[i] = fvect1->array[i];

  j = 0;
  for (i = fvect1->numElements; j < fvect2->numElements;)
    newFVector->array[i++] = fvect2->array[j++];
  
  return newFVector;
}

VECTOR_OF_F_VECTORS *
ConcatVfv (VECTOR_OF_F_VECTORS * vfv1,
	   VECTOR_OF_F_VECTORS * vfv2, unsigned int numVectors)
{
  VECTOR_OF_F_VECTORS *newVfv;
  unsigned int numElements;
  int i;

  if (vfv1 == NULL || vfv2 == NULL)
	return NULL;
  newVfv = (VECTOR_OF_F_VECTORS *)
    calloc (numVectors, sizeof (VECTOR_OF_F_VECTORS));

  if (newVfv == NULL)
    printf ("Unable to allocate memory in ConcatVfv\n");
  for (i = 0; i < numVectors; i++)
    newVfv[i] = (F_VECTOR *) ConcantenateFVector (vfv1[i], vfv2[i]);

  return newVfv;
}

int
FreeFVector (F_VECTOR * fvect)
{
  if (fvect == NULL)
	return 1;
  if (fvect->array != NULL)
      free (fvect->array);
  free (fvect);
  fvect = NULL;
  return 0;
}

int
FreeVfv (VECTOR_OF_F_VECTORS * vfv, unsigned long numVectors)
{
  int i;
  if (vfv == NULL) return 1;
  for (i = 0; i < numVectors; i++)
    {
      if (vfv[i] != NULL) {
          free (vfv[i]->array);
          free (vfv[i]);
          vfv[i] = NULL;
      }
    }
  free (vfv);
  vfv = NULL;
  return 1;
}

int
CopyFVector (F_VECTOR *src, F_VECTOR *dest)
  {
	int i;
	if (src == NULL || dest == NULL)
		return 1;    
        if (src->array == NULL || dest->array == NULL)
		return 2;
	
	if (src->numElements != dest->numElements)
		return 3;

	for (i = 0; i < src->numElements; i++)
		dest->array[i] = src->array[i];
	return 0;	
  }

F_VECTOR*
CloneFVector (F_VECTOR *src)
  { 
 	F_VECTOR *fvect;
 	if (src == NULL || src->array == NULL)
		return NULL;
	fvect = (F_VECTOR *) AllocFVector (src->numElements);
	if (CopyFVector (src, fvect))
		return NULL;
	return fvect;
  }

VECTOR_OF_F_VECTORS*
ReadVfvFromFile (char *fileName, unsigned int *numVectors)
  {
	int dim, length, i, j;
    int num_elem_read = 0;
	VECTOR_OF_F_VECTORS *vfv = NULL;
	FILE *dataFile;

	dataFile = fopen (fileName, "r");
	if (dataFile == NULL)
		return NULL;
	num_elem_read = fscanf (dataFile, "%d %d", &dim, &length);
	
	if (num_elem_read < 2 || length == 0)
		return NULL;
	*numVectors = (unsigned int) length;

	vfv = (VECTOR_OF_F_VECTORS *)
		calloc (*numVectors, sizeof(VECTOR_OF_F_VECTORS));

	if (vfv == NULL)
		return NULL;

	i = 0;
	while (i < *numVectors)
	  {
		length = *numVectors;
		vfv[i] = (F_VECTOR*) AllocFVector (dim);
		j = 0;
		while (j < dim)		
			fscanf (dataFile, " %f", &vfv[i]->array[j++]);
		i++;
	  }
	if (i != *numVectors)
		printf ("In ReadVfvFromFile: reading stopped abruptly\n");
	fclose (dataFile);
	return vfv;
  }

VECTOR_OF_F_VECTORS*
ReadVfvFromBinFile (char *fileName, unsigned int *numVectors) {
	int dim, length, i, j,k;
	VECTOR_OF_F_VECTORS *vfv = NULL;
	FILE *dataFile;
  unsigned char *uctemp;

	dataFile = fopen (fileName, "rb");
	if (dataFile == NULL)
		return NULL;
    uctemp = (unsigned char *) &dim;
    finc(k,sizeof(int),1)
        uctemp[k] = (unsigned char) fgetc (dataFile);
    
    uctemp = (unsigned char *) &length;
    finc(k,sizeof(int),1)
        uctemp[k] = (unsigned char) fgetc (dataFile);
    
	
	if (length == 0)
		return NULL;
	*numVectors = (unsigned int) length;

	vfv = (VECTOR_OF_F_VECTORS *)
		calloc (*numVectors, sizeof(VECTOR_OF_F_VECTORS));

	if (vfv == NULL)
		return NULL;

	i = 0;
	while (i < length)
	  {
		vfv[i] = (F_VECTOR*) AllocFVector (dim);
		j = 0;
		while (j < dim) {		
            uctemp = (unsigned char *) &vfv[i]->array[j++];
            finc(k,sizeof(float),1) 
                uctemp[k] = (unsigned char) fgetc (dataFile);
        }
		i++;
	  }
	if (i != *numVectors)
		printf ("In ReadVfvFromFile: reading stopped abruptly\n");
	fclose (dataFile);
	return vfv;
}
/* 
 * Write a stream to a file. returns 0 on success, 1 on failure
 */

int
WriteVfvToFile (VFV *vfv, uint nvec, char *filename) {
    FILE *f = NULL;
    int i,j,d;
    
    f = fopen (filename, "w");
    if (f == NULL) return 0;
    
    if (nvec == 0) return 0;
    d = vfv[0]->numElements;
    fprintf (f, "%d %d\n", d, nvec);

    finc(i,nvec,1) {
        finc(j,d,1) 
            fprintf (f, " %e", vfv[i]->array[j]);
        fprintf(f, "\n");
    }
    fclose (f);
    return 0;
}

int
WriteVfvToBinFile (VFV *vfv, uint nvec, char *filename) 
{
    int             i, j, k, 
                    dim;
    FILE            *opf = NULL;  
    size_t          sizeofint = sizeof (int),
                    sizeoffloat = sizeof (float);

    if (nvec == 0) 
        return 1;

    opf = fopen (filename, "wb");
    if (opf == NULL) {
        return 2;
    }

    // read dim value
    dim = vfv[0]->numElements;
    finc(k,sizeofint,1)  
       putc (((unsigned char *) &dim)[k], opf);
    
    finc(k,sizeofint,1)  
       putc (((unsigned char *) &nvec)[k], opf);
    
    for (j = 0; j < nvec; j++) {
        for (i = 0; i < dim; i++) {        
            finc(k,sizeoffloat,1)                         // write each byte to file
                putc (*(((unsigned char *) &vfv[i]->array[j]) + k), opf);
        }
    }
    fclose (opf);
    return 0;
}


VECTOR_OF_F_VECTORS*
JoinVfv (VECTOR_OF_F_VECTORS *vfv1, unsigned int nv1,
	VECTOR_OF_F_VECTORS *vfv2, unsigned int nv2,
	unsigned int *numVectors,
	unsigned int cleanUp)
  {
	unsigned int nv;
	VECTOR_OF_F_VECTORS *vfv = NULL;
	int i,j;

	if (vfv1 == NULL || vfv2 == NULL)
		return NULL;

	nv = nv1 + nv2;
	vfv = (VECTOR_OF_F_VECTORS*)
		 	calloc (nv,
			sizeof (VECTOR_OF_F_VECTORS));
	  
	
	for (i = 0; i < nv1; i++)
		vfv[i] = vfv1[i];		

	j = 0;
	for (i = nv1; i < nv; i++)
		vfv[i] = vfv2[j++];

	*numVectors = nv;	
	if (!cleanUp)
		return vfv;

	else if (cleanUp == 1) free (vfv1);
	else if (cleanUp == 2) free (vfv2);
	else if (cleanUp == 3)
	  {
		free (vfv1);
		free(vfv2);
	  }
	return vfv;
  }

int
ReadTriGaussianFile (char *filename, 
                    VECTOR_OF_F_VECTORS **ergMeans,
                    VECTOR_OF_F_VECTORS **ergVars,
                    float **ergWts)
  {
					return ReadGMMFile (filename, ergMeans,
													ergVars, ergWts,
													3, 1);
  }

int
ReadGMMFile (char *filename,
            VECTOR_OF_F_VECTORS **rmeans,
            VECTOR_OF_F_VECTORS **rvars,
            float **rwts,
	    	int nMix,
            int dim)
{
	  FILE *gmmFile;
		VECTOR_OF_F_VECTORS *means;
		VECTOR_OF_F_VECTORS *vars;
		float *wts;
	  int i,j,k;
	 if ( dim < 0 || nMix < 0)
		 return 1;	
  means =
    (VECTOR_OF_F_VECTORS *) calloc (nMix, sizeof (VECTOR_OF_F_VECTORS));
  vars =
    (VECTOR_OF_F_VECTORS *) calloc (nMix, sizeof (VECTOR_OF_F_VECTORS));

  for (i = 0; i < nMix; i++)
    {
      means[i] = (F_VECTOR *) AllocFVector (dim);
      vars[i] = (F_VECTOR *) AllocFVector (dim);
    }

  wts  = (float *) calloc (nMix, sizeof(float));

  gmmFile = fopen (filename, "r");
  for (i = 0; i < nMix; i++)
    {
      fscanf (gmmFile, "%f\n", &wts[i]);
      for (j = 0; j < dim; j++)
	{
	  fscanf (gmmFile, " %f %f", &means[i]->array[j],
		  &vars[i]->array[j]);
	}
      fscanf (gmmFile, "\n");
    }
  fclose (gmmFile);
	*rmeans = means;
	*rvars = vars;
	*rwts = wts;
  return 0;
}

VECTOR_OF_F_VECTORS *
BatchExtractFeatureVectorsVAD (ASDF * asdf, char *batchFileName,
			  char *featureName,
			  unsigned long *numVectors,
			  VECTOR_OF_F_VECTORS * ergMeans,
			  VECTOR_OF_F_VECTORS * ergVars, float *ergWts)
  {
     FILE *batchFile;
     VECTOR_OF_F_VECTORS *allVfv, *vfv;
     unsigned int tempNV, tempNV2, nv;
     char wavname[256];

     batchFile = fopen (batchFileName, "r");
     if (batchFile == NULL)
        {
          printf ("BatchExtractFeatureVectorsVAD: Unable to open ");
          printf ("%s ...\n", batchFileName);
        }
     tempNV = 0;
     nv = 0;
     allVfv = NULL;
 
     while (!feof (batchFile))
       {
          fscanf (batchFile, "%s", wavname);
	  vfv = ExtractFeatureVectorsVAD (asdf, wavname, featureName,
	                             (unsigned long *) &tempNV,
	                             ergMeans, ergVars, ergWts);
          if (allVfv == NULL)            
	    {
	      allVfv = vfv;
	      nv = tempNV;
	      continue;
	    }
	 allVfv = JoinVfv (allVfv, nv, vfv, tempNV,&tempNV2, 3);
	 nv += tempNV2;
       }         
       *numVectors = nv;
       fclose (batchFile);
       return allVfv;
  }

VECTOR_OF_F_VECTORS*
BatchReadVfvFromBinFile (const char *listFileName, unsigned long *numVectors)
  {
    FILE *listFile, *vfvfile;
    char *wavname, *fn;
    VECTOR_OF_F_VECTORS *allVfv, *vfv, *vfv2;
    unsigned int tempNV2, nv, i, nvacc;
    int tempNV, temp;
    unsigned char *uctemp;

    listFile = fopen (listFileName, "r");
    if (listFile == NULL)
      {
        printf ("Unable to open list file...\n");
	return NULL;
      }
    
     tempNV = 0;
     nv = 0;
     allVfv = NULL;
     // open each file and find out size
    while (!feof (listFile)) {
        fn = (char *) calloc (256, sizeof(char));
        fscanf (listFile, "%s\n",fn);
        vfvfile = fopen (fn, "rb");
        uctemp = (unsigned char *) &temp;
        for (i = 0; i < sizeof (int); i++) 
            uctemp[i] = (unsigned char) fgetc (vfvfile);
        uctemp = (unsigned char *) &tempNV;
        for (i = 0; i < sizeof (int); i++) 
            uctemp[i] = (unsigned char) fgetc (vfvfile);
        nv += tempNV;
        fclose (vfvfile);
        free (fn);
    }
    allVfv = (VECTOR_OF_F_VECTORS*) calloc (nv, sizeof(VECTOR_OF_F_VECTORS));
    tempNV = 0;
    nvacc = 0;
    rewind (listFile);
    while (!feof (listFile))
      { 
        fn = (char *) calloc (256, sizeof(char));
        fscanf (listFile, "%s",fn);        
        vfv = ReadVfvFromBinFile (fn, &tempNV);
        if (vfv == NULL) continue;
        for (i = 0; i < tempNV; i++) {
          allVfv[nvacc++] = vfv[i];
        }
        free (fn);
        free (vfv);
      }
      fclose (listFile);
      *numVectors = nv;
      return allVfv;
  }

VECTOR_OF_F_VECTORS*
BatchReadVfvFromFile (const char *listFileName, unsigned long *numVectors)
  {
    FILE *listFile, *vfvfile;
    char *wavname, *fn;
    VECTOR_OF_F_VECTORS *allVfv, *vfv, *vfv2;
    unsigned int tempNV, tempNV2, nv, i, temp, nvacc;

    listFile = fopen (listFileName, "r");
    if (listFile == NULL)
      {
        printf ("Unable to open list file...\n");
	return NULL;
      }
    
     tempNV = 0;
     nv = 0;
     allVfv = NULL;
     // open each file and find out size
    while (!feof (listFile)) {
        fn = (char *) calloc (256, sizeof(char));
        fscanf (listFile, "%s\n",fn);
        vfvfile = fopen (fn, "r");
        fscanf (vfvfile, "%u %u\n", &temp, &tempNV);            
        nv += tempNV;
        fclose (vfvfile);
        free (fn);
    }
    allVfv = (VECTOR_OF_F_VECTORS*) calloc (nv, sizeof(VECTOR_OF_F_VECTORS));
    tempNV = 0;
    nvacc = 0;
    rewind (listFile);
    while (!feof (listFile))
      { 
        fn = (char *) calloc (256, sizeof(char));
        fscanf (listFile, "%s",fn);        
        vfv = ReadVfvFromFile (fn, &tempNV);
        if (vfv == NULL) continue;
        for (i = 0; i < tempNV; i++) {
          allVfv[nvacc++] = vfv[i];
        }
        free (fn);
        free (vfv);
      }
      fclose (listFile);
      *numVectors = nv;
      return allVfv;
  }

VECTOR_OF_F_VECTORS*
BatchExtractFeatureVectors (ASDF * asdf, char *fileName,
		       char *featureName,
		       unsigned long *numVectors, float ts)
  {
    
    FILE *listFile;
    char *wavname, fn[256];
    VECTOR_OF_F_VECTORS *allVfv, *vfv;
    unsigned int tempNV, tempNV2, nv;

    listFile = fopen (fileName, "r");
    if (listFile == NULL)
      {
        printf ("Unable to open list file...\n");
	return NULL;
      }
     tempNV = 0;
     nv = 0;
     allVfv = NULL;
    while (!feof (listFile))
      { 
        fscanf (listFile, "%s",fn);
	vfv = ExtractFeatureVectors (asdf, featureName, fn, (unsigned long *) &tempNV,ts) ;
          if (allVfv == NULL)            
	    {
	      allVfv = vfv;
	      nv = tempNV;
	      continue;
	    }
	 allVfv = JoinVfv (allVfv, nv, vfv, tempNV,&tempNV2, 3);
	 nv += tempNV2;
      }

      fclose (listFile);
      *numVectors = nv;
       return allVfv;
}

int
SaveModelToFile (VECTOR_OF_F_VECTORS *ubmMeans,
                VECTOR_OF_F_VECTORS *ubmVars,
                float *ubmWeights,
                int numMix, char *fileName)
  {
    FILE *outFile;
    int featLength ,i,j,k;
    outFile = fopen (fileName, "w");
    if (outFile == NULL)
      {
         printf ("Unable to open new file\n");
         return 1;
      } 
    
    if (numMix < 1 || ubmMeans == NULL ||
       ubmVars == NULL || ubmWeights == NULL)
      {
        printf ("SaveModelToFile: Invalid arguments\n");
        return 2;
      }
    featLength = ubmMeans[0]->numElements; 
    for (k = 0; k < numMix; k++)
      {
		fprintf (outFile, "%e\n", ubmWeights[k]);
		for (j = 0; j < featLength; j++)		  
			fprintf (outFile, " %e %e", 
                                ubmMeans[k]->array[j], 
                                ubmVars[k]->array[j]);
		  
		fprintf (outFile,"\n");
      }
      fclose (outFile);
      return 0;
  }

/****************************************************************
 * Function: erfinv
 * Desc    : returns t for which erf(t) = int_0_x 2/sqrt(pi) *
 *                                        exp(-t^2) dt
 * Ip arg  : x val
 * Op      : t val
****************************************************************/

float erfinv (float x)
{
  float a = 0.147;
  float pi = 3.14159265359;
  float t1 = (2.0 / (pi * a)); 
  float t2 = (log(1- x*x)) ;
  float sign = x > 0.0 ? 1.0 : -1.0;
  t1 = t1 + t2/2.0;
  return sign * sqrt(-t1 + sqrt ((t1*t1) - (t2/a)));
}

/***************************************************************
 * Function: flinspace
 * Desc    : creates an float array of length 'len' from 's'
 *           to 'e' like an arithmetic progression
 * Ip args : s - start
 *           e - end
 *           len - length
 * Op      : floating pt array starting from s, ending at e with
 *           a length of len
 ***************************************************************/           

float*
flinspace (float s, float e,uint len)
{
  float step = (float) (e-s)/(float)(len-1);
  int i;
  float *x;
  x = (float *) calloc(len, sizeof(float));
  finc(i,len,1) {
    x[i] = s+(i*step);
  }
  return x;
}

int
compare_float (const float *a, const float *b)
{
  float x = *a ;
  float y = *b;
  if (x>y)
    return 1;
  else if (x==y)
    return 0;
  else return -1;
}

VFV*
Gaussianization (VFV* v, unsigned long nv, uint wnd)
{
  int i,j,k,s,e,idx,clen;
  int sidx, eidx;
  int ngt, nless, neq;
  uint l,r,d;
  VFV* g;
  F_VECTOR *fv;
  float *fa;
  const float sqrt2 = sqrt(2);
  static float *x = NULL;
  if (x == NULL) {
    x = flinspace (0.0,1.0,wnd+3);
   finc(i,wnd+3,1)
      x[i] = erfinv(x[i]*2.0 - 1.0) * sqrt2;
  }

  l = wnd / 2;
  r = wnd - l;
        l=0;
        r=wnd-1;
  d = v[0]->numElements;
  g = (VFV*) calloc (nv, sizeof(VFV));
  finc(j,nv,1)
    g[j] = (F_VECTOR*) AllocFVector (d);

  fa = (float*) calloc (wnd,sizeof(float));
  if (fa == NULL)
    {
      printf ("Unable to allocate floating pt array\n");
      exit(1);
    }
  finc(i,d,1) {
    finc(j,nv,1) {
      s = j+1-l; 
      if (s < 0) s = 0;
      e = j+r-1;
      if (e >=nv) e = nv-1;
if(nv < 500) {
s=0;
e=nv-1;
}
      clen = e-s+1;
      idx = 0;
      ngt = 0;
      nless = 0;
      neq = 0;
      for (k = s; k <= e; k++) {
          if (v[k]->array[i] > v[j]->array[i]) ngt++;
          else if (v[k]->array[i] < v[j]->array[i]) nless++;
          else neq++;
      }
      sidx = nless;
      eidx = sidx + neq - 1;
      if (eidx > sidx) sidx = round((eidx + sidx)/2);
      if (e-s < wnd) sidx = round(sidx*(l+r+1)/(clen));
      sidx++;
      g[j]->array[i] = x[sidx];
    }
  }
  free (fa);

  return g;
}

int
set_select_channel(char x)
{
  if (x == 'a'  || x == '0'  || x == 'A') 
        select_channel=0;
  else if(x == 'b' || x == '1' || x == 'B') 
        select_channel=1;
  else {
        fprintf(stderr,"incorrect option for select_channel\n");
        return 1;
  }
  return 0;
}

int**
ReadTopCValuesFromFile(char *filename,int cValue,
                       unsigned int numVectors)
{
    int **mixids = NULL;
    int i,j;
    FILE *mixfile = NULL;

    mixfile = fopen (filename, "r");
    if (mixfile == NULL) return NULL;
    mixids = (int **) malloc (sizeof(int*) * numVectors);
    if (mixids == NULL) return NULL;
    for (i = 0; i < numVectors; i++) {
        mixids[i] = (int *) malloc (sizeof(int)*cValue);
        if (mixids[i] == NULL) return NULL;
        for (j = 0; j < cValue; j++) fscanf(mixfile, "%d", &mixids[i][j]);
        for (j = 0; j < cValue; j++) mixids[i][j] = mixids[i][j] + 1;
        fscanf(mixfile,"\n");
    }
    fclose (mixfile);
    return mixids;
}

void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}
