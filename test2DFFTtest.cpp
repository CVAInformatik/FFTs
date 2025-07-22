/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/

//#include "simpleRadix2IOIPFFT.h"
#include "simpleRadix2IOIPFFTindex.h"

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
	IndexedSimpleRadix2IOIPFFTtype* fft = new IndexedSimpleRadix2IOIPFFTtype();
	
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
    
    putValue(N,sre1, 4, 2, 1.0);
    putValue(N,sre1, 4+1, 2, 1.0);
    putValue(N,sre1, 4, 2+1, 1.0);

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

     /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1 + (i*N), sim1+(i*N),1 );

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1+i, sim1+i, N );


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

/*
     trying to make 2D cross correlation 
*/
void test2(u32 P )
{
	IndexedSimpleRadix2IOIPFFTtype* fft = new IndexedSimpleRadix2IOIPFFTtype();
	
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

#define OX 5
#define OY 5
    
    putValue(N,sre1, OX, OY, 2.0);
    putValue(N,sre1, OX+1, OY, 1.0);
    putValue(N,sre1, OX,OY+1, 1.0);

    printf("SR1: \n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");


#define OXd 3
#define OYd 6

    putValue(N,sre2, OX+OXd,   OY+OYd,   2.0);
    putValue(N,sre2, OX+1+OXd, OY+OYd,   1.0);
    putValue(N,sre2, OX+OXd,   OY+1+OYd, 1.0);

    printf(" SR2: \n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre2[y*N+x],sim2[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre1 + (i*N), sim1+(i*N),1 );

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre1+i, sim1+i, N );

    
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre2 + (i*N), sim2+(i*N),1 );

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->ForwardFFT(sre2+i, sim2+i, N );


    printf("FFT af S1 :\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    printf("FFT af S2 :\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre2[y*N+x],sim2[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    // now we multiply pointwise the elements of S1 
    // with the complex conjugates of S2
    // it is the 'reference' signal we conjugate !
    //  (a-ib)*(c+id) = (ac+bd)+i(-bc +ad)
    // we do it in a linear sweep
    for(u32 i = 0; i < N*N; i++){
        FFTType tre = sre1[i]*sre2[i]+sim1[i]*sim2[i];
        FFTType tim = -(sim1[i]*sre2[i])+sre1[i]*sim2[i];
        sre1[i] = tre;
        sim1[i] = tim;
		}
		
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1 + (i*N), sim1+(i*N),1 );

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++)
       fft->InverseFFT(sre1+i, sim1+i, N );


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

void test3(){
	IndexedSimpleRadix2IOIPFFTtype* fft = new IndexedSimpleRadix2IOIPFFTtype();
	
	fft->SetLength(8);  // 256 punkter
  u32 N = fft->Status();
 	printf("N  %3d\n", N);
  
  if( N ){
		FFTType *sre1 = new FFTType[N];
		FFTType *sim1 = new FFTType[N];
		FFTType *sre2 = new FFTType[N];
		FFTType *sim2 = new FFTType[N];
		
		initBuffer(N, sre1);
		initBuffer(N, sim1);
		initBuffer(N, sre2);
		initBuffer(N, sim2);

		sre1[40] = 1.0;
		sre1[41] = 2.0;
		sre1[42] = 1.0;

		sre2[21] = 1.0;
		sre2[22] = 2.0;
		sre2[23] = 1.0;

		fft->ForwardFFT(sre1, sim1);
		fft->ForwardFFT(sre2, sim2);
		
 		// now we multiply pointwise the elements of S2 
    // with the complex conjugates of S1
    //  (a-ib)*(c+id) = (ac+bd)+i(-bc + ad)
    for(u32 i = 0; i < N; i++){
        FFTType tre = sre1[i]*sre2[i]+sim1[i]*sim2[i];
        FFTType tim = -(sim1[i]*sre2[i])+(sre1[i]*sim2[i]);
        sre1[i] = tre; // we use s1 to get the result
        sim1[i] = tim;
		}

		fft->InverseFFT(sre1, sim1);
		for(u32 i = 0 ; i < N; i++)
    		printf("%-4d: (% " FORMAT ",% " FORMAT") \n",i , sre1[i],sim1[i] );
		
	}
}

void testRotation(uint length)
{
	IndexedSimpleRadix2IOIPFFTtype* fft = new IndexedSimpleRadix2IOIPFFTtype();

	fft->SetLength(length);  // 2^length punkter
	u32 N = fft->Status();
 	printf("N  %3d\n", N);
	if( N ){
		FFTType *sre1 = new FFTType[N];
		FFTType *sim1 = new FFTType[N];
		u32     *indexMap = new u32[N];
		for(int r = 1 ; r < N; r+=2){
			std::cout << "begin rotation " << r << std::endl;
			fft->rotate(r, sre1, sim1, indexMap);
			std::cout << "done rotation " << r << std::endl;
		}
	}
}

#undef FORMAT
#define FORMAT "5.2f"
void test2index(u32 P )
{
	IndexedSimpleRadix2IOIPFFTtype* fft = new IndexedSimpleRadix2IOIPFFTtype();
	
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
    
    u32 *indexMap =  new u32[N];
    
#undef OX
#undef OY
#define OX 5
#define OY 5
    
    putValue(N,sre1, 0, 0, 1.0);
    putValue(N,sre1, 1, 0, 1.0);
    putValue(N,sre1, 0, 1, 1.0);
    putValue(N,sre1, 1, 1, 1.0);

    putValue(N,sre1, 5, 0, 1.0);
    putValue(N,sre1, 6, 0, 1.0);
    putValue(N,sre1, 7, 0, 1.0);
    putValue(N,sre1, 5, 1, 1.0);
    putValue(N,sre1, 6, 1, 1.0);
    putValue(N,sre1, 7, 1, 1.0);
    putValue(N,sre1, 5, 2, 1.0);
    putValue(N,sre1, 6, 2, 1.0);
    putValue(N,sre1, 7, 2, 1.0);


    printf("SR1: \n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(%" FORMAT ",%" FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");

#undef   OXd
#undef   OYd
#define OXd 1
#define OYd 1

    putValue(N,sre2, OX+OXd,   OY+OYd,   1.0);
    putValue(N,sre2, OX+1+OXd, OY+OYd,   1.0);
    putValue(N,sre2, OX+OXd,   OY+1+OYd, 1.0);
    putValue(N,sre2, OX+1+OXd, OY+1+OYd, 1.0);
    printf(" SR2: \n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre2[y*N+x],sim2[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = (i*N) + ix;
       fft->ForwardFFT(sre1, sim1,indexMap );
     }

    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = i + (ix*N);
       fft->ForwardFFT(sre1, sim1, indexMap );
    }
    
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = (i*N) + ix;
       fft->ForwardFFT(sre2, sim2,indexMap );
		}
    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = i + (N*ix);
       fft->ForwardFFT(sre2, sim2,indexMap );
		}

    printf("FFT af S1 :\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre1[y*N+x],sim1[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    printf("FFT af S2 :\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(% " FORMAT ",% " FORMAT") ", sre2[y*N+x],sim2[y*N+x] );
      printf("\n");
  	}
    printf("\n");

    // now we multiply pointwise the elements of S1 
    // with the complex conjugates of S2
    // it is the 'reference' signal we conjugate !
    //  (a-ib)*(c+id) = (ac+bd)+i(-bc +ad)
    // we do it in a linear sweep
    // we use sre1/sim1 for the result
    for(u32 i = 0; i < N*N; i++){
        FFTType tre = sre1[i]*sre2[i]+sim1[i]*sim2[i];
        FFTType tim = -(sim1[i]*sre2[i])+sre1[i]*sim2[i];
        sre1[i] = tre;
        sim1[i] = tim;
		}
		
    /* along the x-axis  N ffts stride 1 */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = (i*N) + ix;
       fft->InverseFFT(sre1 , sim1,indexMap);
		}
    /* along the y-axis  N ffts stride N */
    for(int i = 0 ; i < N; i++){
    	 for(int ix = 0; ix < N; ix++)
    	     indexMap[ix] = i + (ix*N) ;
       fft->InverseFFT(sre1, sim1, indexMap );
		}

    printf("\n");
    for(int y = 0 ; y < N; y++){
    	for(int x = 0 ; x < N; x++)
    		printf("(%" FORMAT ")", (sre1[y*N+x]*sre1[y*N+x])+(sim1[y*N+x]*sim1[y*N+x] ));
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
	//test(3);
	//test2(4);
  //test3();
  //testRotation(3);
  //testRotation(5);
  test2index(5);
	return 0;
}