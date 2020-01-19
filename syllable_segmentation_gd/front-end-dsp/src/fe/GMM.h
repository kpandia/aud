/*
    This is a header file for GMM.c

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

#ifndef GMM_H
#define GMM_H
void InitGMM (VECTOR_OF_F_VECTORS *vfv,int numVectors, 
              VECTOR_OF_F_VECTORS *mixtureMeans,               
	      VECTOR_OF_F_VECTORS *mixtureVars, 
	     int numMixtures, int seed);
float ComputeProbability(F_VECTOR *meanVector, F_VECTOR *varVector, 
			 float priorProb, F_VECTOR *fvect, 
			 float probScaleFactor);
float *ComputeMixtureContribution(F_VECTOR *fvect, VECTOR_OF_F_VECTORS *clusterMeans, 
		       VECTOR_OF_F_VECTORS *clusterVars, int numClusters,
		       float *clusterWeight,
				  float probScaleFactor, float *mixtureContribution);
void ComputeGMM(VECTOR_OF_F_VECTORS *vfv, int numVectors, 
		VECTOR_OF_F_VECTORS *clusterMeans, 
		VECTOR_OF_F_VECTORS *clusterVars, 
		float *clusterWeight, int numClusters, 
		int VQIter, int GMMIter, float probScaleFactor, 
                int ditherMean, int varianceNormalize, float varianceFloor, int seed);
#endif


/*-------------------------------------------------------------------------
 * $Log$
 *
 * Local Variables:
 * time-stamp-active: t
 * time-stamp-line-limit: 20
 * time-stamp-start: "Last modified:[ 	]+"
 * time-stamp-format: "%3a %02d-%3b-%:y %02H:%02M:%02S by %u"
 * time-stamp-end: "$"
 * End:
 *                        End of GMM.h
 -------------------------------------------------------------------------*/
