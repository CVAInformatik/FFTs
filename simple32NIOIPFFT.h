#pragma once
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
		
*/

#include <iostream>
#include <math.h>
#include <algorithm>
#include "FFTBase.h"
#include "simpleRadix2IOIPFFTindex.h"

/*
      this is a class which implements in-order-in-place 
      FFTs of length  3*(2^P), P = 3,4,5,6.......
      (24,48,96,192,384,768,1536, 3072.....)
      
      it is a combination of a length 3 FFT and 
      IndexedSimpleRadix2IOIPFFTtype using Chinese Remainder Mapping
      and explicit Rotation (which is a method in IndexedSimpleRadix2IOIPFFTtype).
      
      if you just want a power of two ( 8,16..)
      you should use IndexedSimpleRadix2IOIPFFTtype class
      instead,
      
*/

class simple32NIOIPFFTtype: public FFTBaseClass {
	

public: 
	
	simple32NIOIPFFTtype() {
				length = 0;
				indexMap0 = 0;
				indexMap1 = 0;
				indexMap2 = 0;
				activeMap = 0;
				strideMap = 0;
				fft2  = new IndexedSimpleRadix2IOIPFFTtype();
	};
	~simple32NIOIPFFTtype(){
		  delete fft2 ;
		  if(indexMap0) delete [] indexMap0 ;
		  if(indexMap1) delete [] indexMap1 ;
		  if(indexMap2) delete [] indexMap2 ;
		  if(activeMap) delete [] activeMap ;
		  if(strideMap) delete [] activeMap ;
	};
	
	// 
	//  SetLength(P)  ->  length = (2**P)*3
	//  Status() will return the actual length
	void SetLength(u32 P)
	{
		 fft2->SetLength(P);
		 length = fft2->Status();
 		 if (length > 2) { 
 		 	   // we need to find the inverses for the CRT mapping
 		 	   u32  P;
 		 	   u32  Q;

    		 for( u32 r = 1; r <length; r++){
		 		 	   for(P  = 1; P < 3  ; P++) 		 	   
  	           if ( (P*length) == ((r * 3	)+1)) ROT3 = P  ;
        	           	
 		 	   };
 		 	   
    		 for( u32 s = 1; s < 3; s++ ){
		 		 	   for(Q  = 1; Q < length ; Q++) 		 	     	
  	           if ( (Q*3) == ((s*length)+1)) ROT2 = Q  ;
 		 	   };

   		   if(indexMap0) delete [] indexMap0 ;
		  	 if(indexMap1) delete [] indexMap1 ;
		  	 if(indexMap2) delete [] indexMap2 ;
		  	 if(activeMap) delete [] activeMap ;       		 	
		  	 if(strideMap) delete [] strideMap ;       		 	
         indexMap0 = new u32[length];
         indexMap1 = new u32[length];
         indexMap2 = new u32[length];
         activeMap = new u32[length];
         strideMap = new u32[3*length];
         
         /* fill in the mapping tables, using CRT mapping */
         indexMap0[0] = 0;
         if (ROT3 == 1){
            indexMap1[0] = (length)%( length *3);
            indexMap2[0] = (2*length)%( length *3);  ;
         }
         else  {
            indexMap2[0] = (length)  %( length *3);
            indexMap1[0] = (2*length)%( length *3);  
         };

         for(int i = 1; i < length; i++)
         {
            indexMap0[i] = ( indexMap2[i-1] +1) % (length*3);
            indexMap1[i] = ( indexMap0[i-1] +1) % (length*3);
            indexMap2[i] = ( indexMap1[i-1] +1) % (length*3);         
         }
     }   	 
	};
	
	u32 Status() { return 3*length;}
	
	 	void ForwardFFT(FFTType *re, FFTType*im, unsigned int stride )
	 	{
	 		  if( strideMap == 0){
	 		  	strideMap = new u32[3*length];	 		  	
	 			  for(int ix = 0; ix < (3*length); ix++) strideMap[ix] = stride *ix;
	 			}
	 			ForwardFFT( re, im, strideMap);
	 	}
	 		  
	
 	void ForwardFFT(FFTType *re, FFTType*im, unsigned int *indexMap ){
 	 	for( u32 i = 0; i < length ; i++) activeMap[i] = indexMap[indexMap0[i]];
 	 	fft2->ForwardFFT(re,im,activeMap ); 
 	 	fft2->rotate( ROT2, re, im,activeMap);
 	 	for( u32 i = 0; i < length ; i++) activeMap[i] = indexMap[indexMap1[i]];
 	 	fft2->ForwardFFT(re,im,activeMap ); 
 	 	fft2->rotate( ROT2, re, im,activeMap);
 	 	for( u32 i = 0; i < length ; i++) activeMap[i] = indexMap[indexMap2[i]];
 	 	fft2->ForwardFFT(re,im,activeMap ); 
 	 	fft2->rotate( ROT2, re, im,activeMap);
   	unsigned int row[3];	 
   	for(int i = 0; i < length ; i++){
   	 	 row[0]= indexMap[indexMap0[i]];
   	   row[1]= indexMap[indexMap1[i]];
   	   row[2]= indexMap[indexMap2[i]];
   	   FFT3(ROT3, re,im, row);
   	}
 	}

	void UnscaledInverseFFT(FFTType *re, FFTType*im, unsigned int stride){
	  if( strideMap == 0){
	 		  	strideMap = new u32[3*length];	 		  	
	 			  for(int ix = 0; ix < (3*length); ix++) strideMap[ix] = stride *ix;
		}
		UnscaledInverseFFT( re, im, strideMap);
	}
 	
	void UnscaledInverseFFT(FFTType *re, FFTType*im, unsigned int *indexMap){
		ForwardFFT(im, re, indexMap);
	};

	void InverseFFT(FFTType *re, FFTType*im, unsigned int stride )
	{
		  if( strideMap == 0){
	 		  	strideMap = new u32[3*length];	 		  	
	 			  for(int ix = 0; ix < (3*length); ix++) strideMap[ix] = stride *ix;
 			}
 			InverseFFT( re, im, strideMap);
	}
	
	void InverseFFT(FFTType *re, FFTType*im, unsigned int *indexMap ){
	  	ForwardFFT(im, re, indexMap);
		  if( length > 0 ) { 
	    	for(u32 i = 0 ; i < 3*length; i++){
	    			re[indexMap[i]] = re[indexMap[i]]/(3*length);
		  	 		im[indexMap[i]] = im[indexMap[i]]/(3*length);
			  }
			}

	};

  private:
  
  	u32 *indexMap0 ;
  	u32 *indexMap1 ;
  	u32 *indexMap2 ;
  	u32 *activeMap;
  	u32 *strideMap ;
  	u32 ROT2;
  	u32 ROT3;
  	u32 length ;  // the length of the radix 2 FFT *3
  	
  	IndexedSimpleRadix2IOIPFFTtype* fft2;
  	
		/* a simple length 3 FFT after temperton, 
			(complex arithmetic) 		
          t1 = z1 + z2  t2 = z0 - 1/2*t1 t3 = sin60 * (z1-z2)
          x0 = z0 + t1    x1 = t2 +it3   x2 = t2 - it3        
		*/
		void FFT3(u32 rot, FFTType *re, FFTType*im, unsigned int *indexMap )
		{
        FFTType t1re = re[indexMap[1]] +re[indexMap[2]];
        FFTType t1im = im[indexMap[1]] +im[indexMap[2]];
        FFTType t2re = re[indexMap[0]] - (0.5 * t1re) ;
        FFTType t2im = im[indexMap[0]] - (0.5 * t1im) ;
        FFTType t3re = std::sin( PI/3 )*(re[indexMap[1]] - re[indexMap[2]]);
        FFTType t3im = std::sin( PI/3 )*(im[indexMap[1]] - im[indexMap[2]]);

				
				re[indexMap[0]] = re[indexMap[0]]+ t1re;
				im[indexMap[0]] = im[indexMap[0]]+ t1im;
        if(rot == 1){
						re[indexMap[1]] = t2re - t3im;	
						im[indexMap[1]] = t2im + t3re;
						re[indexMap[2]] = t2re + t3im;
						im[indexMap[2]] = t2im - t3re;
				}
				else {
						re[indexMap[2]] = t2re - t3im;	
						im[indexMap[2]] = t2im + t3re;
						re[indexMap[1]] = t2re + t3im;
						im[indexMap[1]] = t2im - t3re;					
				}

		}
};
