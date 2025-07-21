/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/
/*
		The ideas and code in this program are based on
		
		SIAM journal on scientific and statistical computing.
		Society for Industrial and Applied Mathematics.
		SELF-SORTING IN-PLACE FAST FOURIER-TRANSFORMS		
		1991-07-01  Clive Temperton.
		
		I've transcribed the original FORTRAN code to simple C++ with a few obvious optimizations.
		I've retained the original variable names and loop numbers for easy back reference.
*/

#include <iostream>
#include <math.h>
#include <algorithm>
#include "simpleRadix2IOIPFFTindex.h"
#include "simple32NIOIPFFT.h"
/*
      
*/
void initBuffer(u32 N,  FFTType * buffer )
{
	for (u32 i = 0; i < N; i++) buffer[i] = (FFTType) 0.0 ;
}

#define FORMAT "4.3f"

void testsimple32N(u32 p, u32 f)
{
	simple32NIOIPFFTtype  *fft = new simple32NIOIPFFTtype();
  u32 *indexMap;
	FFTType *re;
	FFTType *im;

	fft->SetLength(p);
	u32 N = fft->Status();
	if( N ){
		 indexMap = new u32[N];
		 re = new FFTType[N];
		 im = new FFTType[N];
  	 initBuffer(N, re); 
	 	 initBuffer(N, im); 
     for(u32 i = 0; i < N; i++) indexMap[i] = i;
		 
     unsigned int ix = f;
  	 im[ix	]  =  0.5;   
	   im[N-ix]  = -0.5;   
     fft->InverseFFT(re, im, indexMap);
		 
 		 for(u32 i = 0 ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("\n");

	}
	
}

int main(int argc, char **argv)
{
	  testsimple32N(9, 48);
	  return 0;
}
