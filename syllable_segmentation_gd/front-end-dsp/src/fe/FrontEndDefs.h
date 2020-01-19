/*
    This is a header file containing some definitions

    Copyright (C) 2001-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>

    This file is part of VAD-IITM.

    VAD-IITM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    VAD-IITM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with VAD-IITM.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef FRONT_END_DEFS
#define FRONT_END_DEFS
#define PI            3.1415926535898
/* definitions for complex number arithmetic */

/* complex number type */

 typedef struct cmplx { float re,im ; }complex;

/* temporary variable defenitions for complex arithmetic */

 float  rp_a,im_a,rp_b,im_b;

/* add complex no's a and b and store the result in c */

# define cadd(c,a,b) rp_a = a.re; im_a = a.im; \
                     rp_b = b.re; im_b = b.im; \
                     c.re = rp_a + rp_b;       \
                     c.im = im_a + im_b

/* conjugate f complex number a stored in c */

# define conjg(c,a) rp_a = a.re; im_a = a.im; \
                   c.re = rp_a; \
                   c.im = -im_a
 
/* subtract b from a and store the result in c */ 

# define csub(c,a,b) rp_a = a.re; im_a = a.im; \
                     rp_b = b.re; im_b = b.im; \
                     c.re = rp_a - rp_b;       \
                     c.im = im_a - im_b

/* multiply a and b and store in c */

# define cmul(c,a,b) rp_a = a.re; im_a = a.im;     \
                     rp_b = b.re; im_b = b.im;     \
                     c.re = rp_a*rp_b - im_a*im_b; \
                     c.im = rp_a*im_b + im_a*rp_b

/* divide a by b and store the result in c */

# define cdiv(c,a,b) rp_a = a.re; im_a = a.im; \
                     rp_b = b.re; im_b = b.im; \
                     c.re = ( rp_a*rp_b + im_a*im_b ) \
                           /( rp_b*rp_b + im_b*im_b );\
                     c.im = ( im_a*rp_b - rp_a*im_b ) \
                           /( rp_b*rp_b + im_b*im_b )

# define cabs(b) ((float)sqrt((double)(b.re*b.re+b.im*b.im)))

# define cabs2(b) (float) (b.re*b.re+b.im*b.im)
#define LOG_ZERO  (-1.0E20)
#define LOG_ONE  (0.0)
#define LOG_SMALL (-0.5E10)
#define MINLOGARG (-708.3)
#endif














/*-------------------------------------------------------------------------
 * $Log: front-end-defs.h,v $
 * Revision 1.1  2001/10/25 12:28:05  hema
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
 *                        End of front-end-defs.h
 -------------------------------------------------------------------------*/
