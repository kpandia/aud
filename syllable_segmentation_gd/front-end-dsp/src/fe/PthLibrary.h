/*
    This is a header file for PthLibrary.c

    Copyright (C) 1997-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>, Chaitanya

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


#ifndef PTH_LIBRARY_H
#define PTH_LIBRARY_H
float PitchMinGd(float *Signal,int Npts,int Nfft,int Mfft, int PthLow,int PthHgh,float Alfa);
float 	PitchCepstrum(float *Signal,int Npts,int Nfft,int Mfft, int PthLow, int PthHgh);
float PitchModifiedGd(float *Signal,int Npts,int Nfft,int Mfft, int PthLow, int PthHgh, int WinLen, float winScaleFactor, int gdSmthWinSize, int medOrder, float gamma, float gdPosScale, float gdNegScale);
float PitchModifiedGdLP(float *Signal,int Npts,int Nfft,int Mfft, int PthLow, int PthHgh, int WinLen, float winScaleFactor, int lpOrder, 
			int medianOrder, float gamma, float gdPosScale, float gdNegScale);
float 	PitchLP(float *Signal,int Npts, int frameShift, int PthLow, int PthHgh);
#endif
