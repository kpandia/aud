/*                                                                         
    This file identifies the syllable boundaries using group-delay
    processing of energy contour
                          
    Copyright (C) 1998-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
                
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>
   
    This file is part of KWS.                         
   
    KWS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or   
    (at your option) any later version.
                
    KWS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.                             

    You should have received a copy of the GNU General Public License          
    along with KWS.  If not, see <http://www.gnu.org/licenses/>. 
*/


#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "math.h"
#include "malloc.h"
#include "fe/FrontEndDefs.h"
#include "fe/FrontEndTypes.h"
#include "fe/BatchProcessWaveform.h"
#include "fe/InitAsdf.h"
#include "fe/DspLibrary.h"

/*-------------------------------------------------------------------------
 *  Usage -- Command line arguments for executing the program
 *    Args:	None
 *    Returns:	None
 *    Bugs:	
 * -------------------------------------------------------------------------*/
void Usage() {
  printf("WordsWithSilenceRemoval ctrlFile waveFile boundaryFile thresSilence (ms) thresVoiced (ms) thresEnergy1 thresEnergy2\n");
}	/*  End of Usage		End of Usage   */


/*-------------------------------------------------------------------------
 *  main -- main function
 *    Args:	8 arguments
 *    Returns:	none
 *    Bugs:	none
 * -------------------------------------------------------------------------*/

int main(int argc, char *argv[]) {

  FILE             *ctrlFile=NULL,*boundaryFile=NULL;

  char             *ctrlFilename=NULL, 
                   *boundaryFileName=NULL,
                   *waveFileName=NULL;


  float            *signal=NULL,*sigEgy=NULL, *sigEgyInv=NULL, *sigCopy=NULL,
                   *sigNonZero=NULL,*sigSymmetric=NULL,
                   *ax=NULL,*ay=NULL,
                   *sigCep=NULL,
                   *amag=NULL,*phase=NULL,
                   *derv=NULL,*dervSmth=NULL,
                   *negDerv=NULL,*negDervSmth=NULL;
  int              *peakArr=NULL, *valleyArr=NULL;

  int              nfft, 
                   mfft,
                   frameAdvanceSamples, frameSize, 
                   startVFrameNumber, endVFrameNumber;
  int              i, j,  k, numUVoicedFrames, numVoicedFrames, nfby2;
  int              startUVFrameNumber, endUVFrameNumber, segBoundary;
  int              winlen, medianOrder, lpOrder, percentFrames;
  long             numFrames, numSamples, samplingRate;
  float            alfa, rMin, rMax, scaleFactor,
    rMaxDerv, rMaxNegDerv, energyValue;
  ASDF             *asdf=NULL;
  float            thresEnergy, thresZero, thresSpecFlatness, 
                   durationOfSilence, durationOfVoiced, rMaxEgy,
                   thresVoiced, thresSilence, thresEnergy1, thresEnergy2,
		   rmax, average;
  short            *vU=NULL, *vUSmth=NULL;
  F_VECTOR         *waveform=NULL, *energy=NULL; 

  if (argc != 8) {
    Usage();
    _exit(-1);
  }

  ctrlFilename = argv[1];
  waveFileName = argv[2];
  boundaryFileName = argv[3];
  sscanf(argv[4], "%f", &thresSilence);
  sscanf(argv[5], "%f", &thresVoiced);
  sscanf(argv[6], "%f", &thresEnergy1);
  sscanf(argv[7], "%f", &thresEnergy2);

  ctrlFile = fopen(ctrlFilename,"r");
  boundaryFile = fopen(boundaryFileName,"w");
  
  asdf = (ASDF *) malloc(sizeof(ASDF));
  InitializeStandardFrontEnd(asdf,ctrlFile);
  GsfOpen(asdf, waveFileName);
  
  alfa = GetFAttribute(asdf,"gamma");
  frameAdvanceSamples = GetIAttribute(asdf, "frameAdvanceSamples");
  medianOrder = (int) GetIAttribute(asdf, "medianOrder");
  frameSize = (int) GetIAttribute(asdf, "windowSize");
  thresEnergy = (float) GetFAttribute(asdf, "thresEnergy");
  thresZero = (float) GetFAttribute(asdf, "thresZero");
  thresSpecFlatness = (float) GetFAttribute(asdf, "thresSpecFlatness");
  numFrames = (int) GetIAttribute(asdf, "numFrames");
  numSamples = (int)GetIAttribute(asdf, "numSamples");
  lpOrder = (int)GetIAttribute(asdf, "lpOrder");
  percentFrames = (int)GetIAttribute(asdf,"percentFrames");
  samplingRate = (int) GetIAttribute(asdf, "samplingRate");
  signal = (float *) AllocFloatArray(signal, numFrames+1); 
  sigNonZero = (float *) AllocFloatArray(sigNonZero, numFrames+1); 
  sigEgy = (float *) AllocFloatArray(sigEgy, numFrames+1); 
  sigEgyInv = (float *) AllocFloatArray(sigEgyInv, numFrames+1); 
  waveform = (F_VECTOR *)AllocFVector(frameSize);
  energy = (F_VECTOR *)AllocFVector(1);
  if ((thresSilence > (float) numSamples/samplingRate*1000.0)|| 
      (thresVoiced > (float) numSamples/samplingRate*1000.0)) 
    {printf("Thresholds too large \n");
      _exit(-1);
    }

  rmax = (float) (asdf->waveform[ImaxShort0(asdf->waveform, numSamples)]);
  for (i = 0; i < numFrames; i++) {
    waveform = (F_VECTOR *) FrameComputeWaveform(asdf, i, waveform);
    energyValue = ((F_VECTOR *) (FrameComputeEnergy(asdf, i, energy)))->array[0];
      signal[i] = energyValue;
  }

  j = 0;
  for (i = 0; i < numFrames; i++) {
    if (signal[i] != 0.0) {
      sigNonZero[j] = signal[i];
      j++;
    }
  }
  rMin  = sigNonZero[Imin0(signal, j)];
  for (i = 0; i < numFrames; i++) {
    if (signal[i] == 0.0) 
      signal[i] = rMin;
  }
  rMaxEgy = signal[Imax0(signal, numFrames)];
  for (i = 0; i < numFrames; i++) {
    signal[i] = signal[i]/rMaxEgy;
    sigEgy[i] = signal[i];
    //if (signal[i] != 0.0) 
      //signal[i] = expf(alfa*log(signal[i]));
    if ((signal[i] != 0.0) && (signal[i] != 1.0)) {
      sigEgyInv[i] = 1.0/signal[i];
      signal[i] = expf(alfa*log(1.0/signal[i]));
    }
  }


  /* The following code is not much use -- except to
     ensure that the FFT size is set to the largest possible */

  /* begin FFT size fix */
  mfft = ceil(log((float)numFrames*2)/log(2.0)); 
//  mfft=14;
  nfft = (int) pow(2,mfft); 
  nfby2 = nfft/2;
  Cstore(nfft);
  /* end FFT size fix */

 /* Allocating various arrays required for plotting/group delay
   computation */

  /* Allocate arrays to the largest possible FFT size */

  peakArr = (int *) AllocIntArray(peakArr, nfby2+1);
  valleyArr = (int *) AllocIntArray(valleyArr, nfby2+1);
  sigCopy = (float*) AllocFloatArray(sigCopy, nfft+1);
  sigSymmetric = (float*) AllocFloatArray(sigSymmetric, nfft+1);
  sigCep = (float *) AllocFloatArray(sigCep,nfft+1);
  ay = (float *) AllocFloatArray(ay,nfft+1);
  sigCopy[1] = 100.0;
  Rfft(sigCopy, sigCep, ay, mfft, nfft, 1); 
  ax = (float *) AllocFloatArray(ax,nfft+1);
  derv  = (float *) AllocFloatArray(derv,nfby2+1);
  negDerv = (float *) AllocFloatArray(negDerv,nfby2+1);
  dervSmth = (float *) AllocFloatArray(dervSmth,nfby2+1);
  negDervSmth = (float *) AllocFloatArray(negDervSmth,nfby2+1);
  amag = (float *) AllocFloatArray(amag,nfft+1);
  phase = (float *) AllocFloatArray(phase,nfft+1);
  scaleFactor = GetFAttribute(asdf,"winScaleFactor");
  PutIAttribute(asdf,"fftSize",nfft);
  PutIAttribute(asdf,"fftOrder",mfft);
FILE *spectraFile=fopen("temp.spec","w");

  //printf("nfft = %d mfft =  %d numFrames %ld\n",nfft, mfft, numFrames);

PutIAttribute(asdf,"windowSize",30);

  for (i = 0; i < numFrames; i++) {
    waveform = (F_VECTOR *) FrameComputeWaveform(asdf, i, waveform);
  }


  vU = (short *) calloc(numFrames+1, sizeof(short));
  vUSmth = (short *) calloc(numFrames+1, sizeof(short));
  VoicedUnvoiced(asdf->waveform, numSamples, vU, frameAdvanceSamples,
                 frameSize, lpOrder, thresZero, thresEnergy, thresSpecFlatness, percentFrames);
  vU[numFrames] = vU[numFrames-1];
  vUSmth[numFrames] = vUSmth[numFrames-1];
  MedianSmoothShortArray(vU, numFrames, 3, vUSmth);
  i = 0;
  /* If the duration of silence is smaller than an expected
     duration -- change the silences to voiced */
  startVFrameNumber= 0;
  endVFrameNumber = -1;
  startUVFrameNumber = 0;
  endUVFrameNumber = -1;
  while (i < numFrames) {
    if (vUSmth[i] == 0){
      while ((i <  numFrames) && (vUSmth[i] == 0))
	i++;
      startUVFrameNumber = endVFrameNumber+1;
      endUVFrameNumber = i-1;
      durationOfSilence = (endUVFrameNumber - startUVFrameNumber+1)*frameAdvanceSamples/
	(float)samplingRate*1000;
      if (durationOfSilence < thresSilence){
	for (k = startUVFrameNumber; k <= endUVFrameNumber; k++)
	  vUSmth[k] = 1;
	endUVFrameNumber = startVFrameNumber-1;
	i = startVFrameNumber;
      }
    }else {
      while ((i <  numFrames) && (vUSmth[i] == 1))
	i++;
      startVFrameNumber = endUVFrameNumber+1;
      endVFrameNumber = i-1;
      numVoicedFrames = endVFrameNumber - startVFrameNumber + 1;
      durationOfVoiced = numVoicedFrames*frameAdvanceSamples/(float)samplingRate*1000;
      if (durationOfVoiced < thresVoiced) {
	for (k = startVFrameNumber; k <= endVFrameNumber; k++)
	  vUSmth[k] = 0;
	endVFrameNumber = startUVFrameNumber-1;
	i = startUVFrameNumber;
      }
    }
  }
  /* Handle final Voiced/Unvoiced segment separately 
  if (endVFrameNumber > endUVFrameNumber) {
    if (durationOfVoiced < thresVoiced) {
      for (k = startVFrameNumber; k <= endVFrameNumber; k++)
	vU[k] = 0;
    }
  } else {
    if (durationOfSilence < thresSilence) {
      for (k = startUVFrameNumber; k <= endUVFrameNumber; k++)
	vU[k] = 1;
    }
    }*/
  i = 0;
  startVFrameNumber= 0;
  endVFrameNumber = -1;
  startUVFrameNumber = 0;
  endUVFrameNumber = -1; 
  segBoundary = 0;
  while (i < numFrames) {
    if (vUSmth[i] == 0) {
      while ((i < numFrames) && (vUSmth[i] == 0))
	i++;
      endUVFrameNumber = i-1;
      startVFrameNumber = i;
      numUVoicedFrames = endUVFrameNumber - startUVFrameNumber+ 1;
      printf("startUVFrameNumber = %d endUVFrameNumber = %d \n", startUVFrameNumber, endUVFrameNumber);
      printf ("numUVFrames = %d \n", numUVoicedFrames);

for (k =1; k < numUVoicedFrames; k++) {
          if (signal[k+startUVFrameNumber-1] != 0.0)
            fprintf(spectraFile,"%f %f %f %f %f %f %f %f\n",
                  (k+startUVFrameNumber-1)*(float)frameAdvanceSamples/(float) samplingRate,
                    signal[k+startUVFrameNumber-1],sigEgy[k+startUVFrameNumber-1], 0.0,0.0,0.0,0.0,0.0);
          else
            fprintf(spectraFile,"%f %f %f %f %f %f %f %f\n",
                  (k+startUVFrameNumber-1)*(float)frameAdvanceSamples/(float) samplingRate,
                    0.0, 0.0, 0.0,0.0,0.0,0.0,0.0);
}
	fprintf(boundaryFile, "%d %d %d\n", startUVFrameNumber,  endUVFrameNumber-1, -1);

    } else {
      while ((i <  numFrames) && (vUSmth[i] == 1))
	i++;
      endVFrameNumber = i-1;
      startUVFrameNumber = i;
      numVoicedFrames = endVFrameNumber - startVFrameNumber+ 1;
     printf("startVFrameNumber = %d endVFrameNumber = %d\n", startVFrameNumber, endVFrameNumber);
	printf("numVoiceFrames = %d \n", numVoicedFrames);
      if (numVoicedFrames > 0) {
	mfft = ceil(log((float)numVoicedFrames*2)/log(2.0)); 
	nfft = (int) pow(2,mfft); 
	nfby2 = nfft/2;
	Cstore(nfft);
	printf("nfft = %d mfft = %d \n", nfft, mfft);
	 
	sigCopy[1] = signal[startVFrameNumber];
	for (k = 2; k <= numVoicedFrames; k++) {
	  sigCopy[k] = signal[startVFrameNumber+k-1];
	  sigCopy[nfft-k+2] = sigCopy[k];
	}
	rMax = sigCopy[Imax(sigCopy, numVoicedFrames)];
	//for (k = 1; k <= numVoicedFrames; k++)
	for (k = 1; k <= nfft; k++)
	  sigCopy[k] = sigCopy[k]/rMax;
	rMin = sigCopy[Imin(sigCopy,numVoicedFrames)];
	for (k = numVoicedFrames + 1; k <= nfft-numVoicedFrames + 1; k++)
	  sigCopy[k] = rMin;
	for (k = 1; k <= nfft; k++)
	  sigSymmetric[k] = sigCopy[k];
	winlen = (int) ceil((float)numVoicedFrames/scaleFactor);
	Rfft(sigCopy,sigCep,ay,mfft,nfft,1);
	
	for (k = 1; k <= winlen; k++){
	  sigCopy[k] = sigCep[k]*HanW(k,winlen);
	}
	for (k = winlen+1; k <= nfft; k++)
	  sigCopy[k] = 0;
	
	Rfft(sigCopy,ax,ay,mfft,nfft,-1);
	SpectrumReal(nfft,ax,ay,amag,phase);
	
	for (k = 1; k < nfby2-1; k++){
	  derv[k] = phase[k] - phase[k+1];  
	  negDerv[k] = phase[k+1] - phase[k];
	}
	derv[nfby2] = derv[nfby2-1];
	negDerv[nfby2] = negDerv[nfby2-1];
//	RemoveAverage(derv, nfby2, &average);
//	RemoveAverage(negDerv, nfby2, &average);
	dervSmth = (float *) MedianSmoothArray(derv, nfby2, 
					       medianOrder, dervSmth);
	negDervSmth = (float *) MedianSmoothArray(negDerv, nfby2, 
						medianOrder, negDervSmth);
	rMaxDerv = derv[Imax(derv,nfby2)];
	rMaxNegDerv = negDerv[Imax(negDerv, nfby2)];



for (k =1; k < numVoicedFrames; k++) {
          if (signal[k+startVFrameNumber-1] != 0.0)
            fprintf(spectraFile,"%f %f %f %f %f %f %f %f\n",
                  (k+startVFrameNumber-1)*(float)frameAdvanceSamples/(float) samplingRate,
                    signal[k+startVFrameNumber-1],sigEgy[k+startVFrameNumber-1],sigCopy[k], phase[k], amag[k],
                  dervSmth[k]/fabs(rMaxDerv), negDervSmth[k]/fabs(rMaxNegDerv));
          else
            fprintf(spectraFile,"%f %f %f %f %f %f %f %f\n",
                  (k+startVFrameNumber-1)*(float)frameAdvanceSamples/(float) samplingRate,
                    0.0, 0.0, 0.0, phase[k],amag[k],
                  dervSmth[k]/fabs(rMaxDerv), negDervSmth[k]/fabs(rMaxNegDerv));
}


	FindPeaks(derv,peakArr, nfby2);
	FindPeaks(negDerv,valleyArr,nfby2);
	segBoundary = 0;
	for (j = 1; j <= numVoicedFrames; j++){
	  if(peakArr[j] == 1) { 
	     if(dervSmth[j]/fabs(rMaxDerv) > thresEnergy1) {
		if(segBoundary+startVFrameNumber !=0 )
			if (j+startVFrameNumber-2 != segBoundary+startVFrameNumber -1 )
		       		fprintf(boundaryFile, "%d %d %d\n", segBoundary+startVFrameNumber -1,  j+startVFrameNumber-2, 1);
			else {}
		else
			if(j+startVFrameNumber-2 != 0)
			       fprintf(boundaryFile, "%d %d %d\n", 0,  j+startVFrameNumber-2, 1);
			else {}
	       segBoundary = j;
	     } else if(dervSmth[j]/fabs(rMaxDerv) > thresEnergy2) {
		if(segBoundary+startVFrameNumber !=0)
			if (j+startVFrameNumber-2 != segBoundary+startVFrameNumber -1 )
			      fprintf(boundaryFile, "%d %d %d\n", segBoundary+startVFrameNumber -1,j+startVFrameNumber-2, 2);
			else {}
		else
			if(j+startVFrameNumber-2 !=0 )
			      fprintf(boundaryFile, "%d %d %d\n", 0,j+startVFrameNumber-2, 2);
			else {}
	       segBoundary = j;
	     } else {
		if(segBoundary+startVFrameNumber !=0 )
			if (j+startVFrameNumber-2 !=segBoundary+startVFrameNumber -1 )
		               fprintf(boundaryFile, "%d %d %d\n", segBoundary+startVFrameNumber -1, j+startVFrameNumber-2, 3);
			else {}
		else
			if(j+startVFrameNumber-2 !=0 )
		               fprintf(boundaryFile, "%d %d %d\n", 0, j+startVFrameNumber-2, 3);
			else {}
               segBoundary = j;
             }
	  }
	}
	/* Fix end of voiced segment */
	if (segBoundary < numVoicedFrames) {
		if(segBoundary + startVFrameNumber != 0 )
			if(segBoundary + startVFrameNumber-1 !=  startVFrameNumber+numVoicedFrames-2)
				  fprintf(boundaryFile, "%d %d %d\n",segBoundary + startVFrameNumber-1, startVFrameNumber+numVoicedFrames-2,1);
			else {}
		else
			if(startVFrameNumber+numVoicedFrames != 0)
				fprintf(boundaryFile,"%d %d %d\n",0,startVFrameNumber+numVoicedFrames-2,1);
			else {}
	}
      }
    }
  }


  free(ax);
  free(ay);
  free(amag);
  free(phase);
  free(derv);
  free(negDerv);
  free(sigCopy);
  free(sigNonZero);
  free(signal);
  free(peakArr);
  free(sigCep);
  free(dervSmth);
  free(negDervSmth);
  free(waveform->array);
  free(waveform);
  free (energy->array);
  free(energy);			    
  //  fprintf(boundaryFile,"\n");
  fclose(ctrlFile);
  fclose(boundaryFile);
  free(valleyArr);
  return(0);
  }		 
	/*  End of main		End of main   */

















/*-------------------------------------------------------------------------
 * $Log: WordsWithSilenceRemoval.c,v $
 * Revision 1.1  2013/02/08 10:23:34  hema
 * Initial revision
 *
 *
 *
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of WordsWithSilenceRemoval.c
 -------------------------------------------------------------------------*/
