/*
    This file does the Quicksort of N floating point numbers

    Copyright (C) 2001-2016 Speech and Music Technology Lab,
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
#include "string.h"

/*-------------------------------------------------------------------------
 *  Swap -- Swaps two float elements
 *    Args:	int, int
 *    Returns:	int
 *    Throws:	
 *    See:	
 *    Bugs:	
 -------------------------------------------------------------------------*/
void Swap (float *a, float *b) {
  float temp; 
  temp = *a;
  *a = *b;
  *b = temp;
  
}	/*  End of Swap		End of Swap   */

/*-------------------------------------------------------------------------
 *  QuickSort -- Recursive sort using CAR Hoare's Algorithm
 *    Args:	floa*, int, int
 *    Returns:	Nothing
 *    Throws:	
 *    See:	
 *    Bugs:	
 -------------------------------------------------------------------------*/

void QuickSort(float *v, int left, int right) {
  int   i,j,k;
  float x, w;
  i = left;
  j = right;
  //  printf("i = %d right = %d left = %d j = %d\n", i, right, left, j
  //);
  x = v[(left+right)/2];
  do {
    while (v[i] < x)
      i++;
    while (x < v[j]) 
      j--;
    //   printf("swap ? i = %d %f j = %d %f\n", i, v[i], j, v[j]); 
     if (i <= j) {
       // printf("before swapping i = %d %f j = %d %f\n", i, v[i], j, v[j]); 
       Swap(&v[i],&v[j]);
       // printf("after swapping i = %d %f j = %d %f\n", i, v[i], j, v[j]); 
       i++;
       j--;
       //printf("new vals of i and j i = %d j = %d\n", i, j);
       /*       for (k = 1; k <= 8; k++)
	 printf("%f ", v[k]);
	 printf("\n");*/
     }
 }   while (i < j);
  /* for(k = left; k<= right; k++)
     printf("%d ", v[k]); */
//printf("\n");

  if (left < j) 
    QuickSort(v,left,j);
  if (i < right)
    QuickSort(v,i,right);
  

}	/*  End of QuickSort		End of QuickSort   */


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
 *                        End of QuickSort.c
 -------------------------------------------------------------------------*/







