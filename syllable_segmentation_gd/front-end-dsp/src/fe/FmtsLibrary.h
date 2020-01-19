/*
    This is a header file for FmtsLibrary.c

    Copyright (C) 1996-2016 Speech and Music Technology Lab,
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

#ifndef FMTS_LIB_H
#define FMTS_LIB_H
float 	*FmtsMinGd(float *signal,int npts,int nfft,
		   int mfft, int winLen,float alfa, 
		   float *freq,int *num);
float *FmtsLpMag(float *signal,int npts, int frameShift, int nfft,
	      int mfft, int LPOrder,float *freq,int *num);

float *FmtsLpPhase(float *signal,int npts, int frameShift, int nfft,
	      int mfft, int LPOrder,float *freq,int *num);
float *FmtsCepstrum(float *signal,int npts,int nfft, 
		    int mfft, int winLen,float *freq,int *num);
float *FmtsModGd(float *signal,int npts,int nfft,
		 int mfft, int winLen, int  gdSmthWinSize, float gamma,
		 float gdPosScale, float gdNegScale,
		 float *freq,int *num);
float *FmtsModGdLP(float *signal,int npts,int nfft,
		 int mfft, int winLen, int  lpOrder, int medianOrder, float gamma,
		 float gdPosScale, float gdNegScale,
		 float *freq,int *num);
#endif



