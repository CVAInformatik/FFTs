/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/

#include "simpleRadix2IOIPFFT.h"

#define FORMAT "8.3f"

void test(u32 P )
{
	simpleRadix2IOIPFFTtype* fft = new simpleRadix2IOIPFFTtype();
	
	fft->SetLength(P);
  u32 N = fft->Status();
 	printf("N  %3d\n", N);
  
  if( N ){
		FFTType *sre = new FFTType[N];
		FFTType *sim = new FFTType[N];
	
  	printf("\n\n");
  
#define f (10.0)
#define f2 (19.0)
	  for(u32 i = 0; i < N ; i++) { 
  		FFTType t = (2.0*PI*f)/N;
  		FFTType t2 = (2.0*PI*f2)/N;
  		sre[i] = cos(i*t)+sin(i*t2);
  		sim[i] = 0;
    	printf("(sre, sim) i = %3d ( %2.3f, %2.3f)  : %2.3f\n", i, sre[i], sim[i], sqrt(( sre[i]*sre[i]) + (sim[i] *sim[i])));
  	}
  	printf("\n\n");

		fft->ForwardFFT(sre, sim);
	  printf("\n\n");
		for( u32 j = 0; j < N; j++) 
	      printf(" %3d : ( %" FORMAT ", %" FORMAT ")   %" FORMAT "\n", j, sre[j], sim[j], sqrt(( sre[j]*sre[j]) + (sim[j] *sim[j])));
  	printf("\n\n");

  	fft->InverseFFT(sre,sim);

		for( u32 j = 0; j < N; j++) 
	      printf(" %3d : ( %" FORMAT ", %" FORMAT ") : %" FORMAT "\n", j, sre[j], sim[j], sqrt(( sre[j]*sre[j]) + (sim[j] *sim[j])));

   	delete [] sre;
   	delete [] sim;
  }
}


int main(int argc, char **argv)
{
	test(6);
	return 0;
}