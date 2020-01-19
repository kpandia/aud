/*
    This is a header file containing some definitions

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


#ifndef CONSTANTS_H
#define CONSTANTS_H
#define SAMPLINGFREQ  11025
#define MAXWAVESIZE   SAMPLINGFREQ*18  /* 18 secs buffer length */
#define MAXFRAMESIZE  512 
#define MAXLPORDER    25
#define PI            3.1415926535898
#define PI2           PI*2
#define MAXFRAMES     MAXWAVESIZE/64   /* Assuming a lowest shift of 64 */

/* Approximate formant ranges for formant extraction algo */
/* The exact values used in a program depends on the FFT resolution */
#define F1_LOW   300.0
#define F1_HIGH  900.0
#define F2_LOW   900.0
#define F2_HIGH  2500.0
#define F3_LOW   2500.0
#define F3_HIGH  3500.0
#endif
