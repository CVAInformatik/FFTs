/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/

#include "simpleRadix2IOIPFFT.h"

#define FORMAT "4.3f"

void initBuffer(u32 N,  FFTType * buffer )
{
	for (u32 i = 0; i < N; i++) buffer[i] = (FFTType) 0.0 ;
}

void putValue(u32 N, FFTType *s, u32 x, u32 y, FFTType val) 
{
	if(( x < N ) && (y < N))
		 s[(y*N) + x] = val;
	else
		std::cout << "invalid arguments N, x, y " << N << " " << x  << " " << y << std::endl; 
}

FFTType getValue(u32 N, FFTType *s, u32 x, u32 y)
{
 	if(( x < N ) && (y < N))
		   return s[(y*N) + x] ;
	else
		std::cout << "invalid arguments N, x, y " << N << " " << x  << " " << y << std::endl; 

  return ( FFTType) 0.0;
}

void test(u32 P )
{
	simpleRadix2IOIPFFTtype* fft = new simpleRadix2IOIPFFTtype();
	
	fft->SetLength(P);
  u32 N = fft->Status();
 	printf("N  %3d\n", N);
  
  if( N ){
		FFTType *sre1 = new FFTType[N * N];
		FFTType *sim1 = new FFTType[N * N];
		FFTType *sre2 = new FFTType[N * N];
		FFTType *sim2 = new FFTType[N * N];
	
    initBuffer(N*N, sre1); initBuffer(N*N, sim1);
    initBuffer(N*N, sre2); initBuffer(N*N, sim2);
    
    putValue(N,sre1, 0, 0, 1.0);
    putValue(N,sre1, 1, 0, 1.0);
    putValue(N,sre1, 0, 1, 1.0);

    printf("\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre1 + (i*N), sim1+(i*N),1 );

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre1+i, sim1+i, N );

    printf("\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    // note: We reverse the order of the axes, when we invert !

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1+i, sim1+i, N );

    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1 + (i*N), sim1+(i*N),1 );


    printf("\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");
  
  	delete [] sre1;  
  	delete [] sim1;  
  	delete [] sre2;  
  	delete [] sim2;  
  }
 
}


int main(int argc, char **argv)
{
	test(3);
	return 0;
}