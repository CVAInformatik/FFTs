


#include <iostream>
#include <math.h>
#include <algorithm>
#include "FFTBase.h"
#include "simpleRadix2IOIPFFTindex.h"
#include "simple32NIOIPFFT.h"
#include "FFTClass.h"
/*
      
*/
void initBuffer(u32 N,  FFTType * buffer )
{
	for (u32 i = 0; i < N; i++) buffer[i] = (FFTType) 0.0 ;
}



#define FORMAT "4.3f"

void testVirtual(u32 p, u32 f)
{
	FFTBaseClass  *fft = new IndexedSimpleRadix2IOIPFFTtype();
	//FFTBaseClass  *fft = new simple32NIOIPFFTtype();
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
  	 im[ix	]  =  0.5*N;   
	   im[N-ix]  = -0.5*N;   
     fft->InverseFFT(re, im, indexMap);
     
#define PCOUNT  100 	 
 		 for(u32 i = 0 ; i < PCOUNT; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("................\n");
 		 for(u32 i = N-PCOUNT ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("\n");

     fft->ForwardFFT(re, im, indexMap);
 		 for(u32 i = 0 ; i < PCOUNT; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("................\n");
 		 for(u32 i = N-PCOUNT ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("\n");

	}
	
}

void testFFTClass(u32 l, u32 f)
{
	FFTClass  *fft = new FFTClass();

	FFTType *re;
	FFTType *im;

	fft->Length(l);
	u32 N = fft->Status();
	if( N ){
		 re = new FFTType[N];
		 im = new FFTType[N];
  	 initBuffer(N, re); 
	 	 initBuffer(N, im); 

		 
     unsigned int ix = f;
  	 im[ix	]  =  0.5*N;   
	   im[N-ix]  = -0.5*N;   
     fft->InverseFFT(re, im);
     
#define PCOUNT  100 	 
 		 for(u32 i = 0 ; i < PCOUNT; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("................\n");
 		 for(u32 i = N-PCOUNT ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("\n");

     fft->ForwardFFT(re, im);
 		 for(u32 i = 0 ; i < PCOUNT; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("................\n");
 		 for(u32 i = N-PCOUNT ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , re[i], im[i] );
     printf("\n");

	}
	
}

void testFindLength()
{
		FFTClass  *fft = new FFTClass();
    for(int i = 0 ; i < 257 ; i++){
    	 printf(" min %4d  length  %4d\n", i, fft->FindLength(i) );
    }

	
}


int main(int argc, char **argv)
{
	  //testVirtual(10, 48);
	  //testFFTClass(3072, 48);
	  //testFFTClass(1536, 48);
	  //testFFTClass(1024, 48);
	  testFindLength();
	  return 0;
}
