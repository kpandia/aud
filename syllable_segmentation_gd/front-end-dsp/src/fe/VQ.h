/*
    This is a header file for VQ.c

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

#ifndef VQ_H
#define VQ_H
#include "FrontEndDefs.h"
#include "FrontEndTypes.h"
void InitVQ (VECTOR_OF_F_VECTORS *vfv, 
	     VECTOR_OF_F_VECTORS *clusterMeans, 
	     VECTOR_OF_F_VECTORS *clusterVars,int num_clusters, float seed); 
float ComputeDiscriminant(F_VECTOR *clusterMean, 
			  F_VECTOR *clusterVar, F_VECTOR *fvect, int varianceNormalize);
int DecideWhichCluster(F_VECTOR *fvect, VECTOR_OF_F_VECTORS *clusterMeans,
		       VECTOR_OF_F_VECTORS *clusterVars, int num_clusters,
		       int varianceNormalize);
void ComputeVQ(VECTOR_OF_F_VECTORS *vfv, int num_vectors, VECTOR_OF_F_VECTORS *cluster_means, VECTOR_OF_F_VECTORS *cluster_vars, 
float *cluster_elem_cnt, int num_clusters, int varianceNormalize, 
float ditherMean, int iterations, float seed);
#endif




