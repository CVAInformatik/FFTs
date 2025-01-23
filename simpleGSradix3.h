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
#include <vector>

typedef uint32_t u32 ;
typedef double FFTType ;
const FFTType PI = 3.141592653589793238463;

class simpleGSRadix3FFTtype {

public: 
	simpleGSRadix3FFTtype(){ Pw = 0; length = 1;};
	~simpleGSRadix3FFTtype(){	
		    ScrambleTable.clear();
		 	  delete [] wre;
		 	  delete [] wim;
	};
	
	// 
	//  SetLength(P)  ->  length = 3**P, 
	//  Status() will return the actual length
	void SetLength(u32 P)
	{
		 if (length > 1){
		 	  length = 1;
		 	  Pw = 0;
		 	  ScrambleTable.clear();
		 	  delete [] wre;
		 	  delete [] wim;
		 }
		 if( P > 0 ) {
		 		length = 1; 
		 		Pw = P; 
		 		for(int i = 0; i < P ; i++) length = length + length + length;  
		 		wre = new FFTType[length];
		 		wim = new FFTType[length];
 		    for(u32 i = 0; i < length; i++) {
					wre[i] = std::cos(( 2*PI *i )/length);
					wim[i] = std::sin(( 2*PI *i )/length);
    		}
		 }
	};
	
	u32 Status() { return length;};
	
	void ForwardFFT(FFTType *re, FFTType*im ){ 		 if( length > 1 ) GSradix3FFT(length,re, im, wre, wim  );		}

	void UnscaledInverseFFT(FFTType *re, FFTType*im ){ if( length > 1 ) GSradix3FFT(length,im, re, wre, wim  );	}

	void InverseFFT(FFTType *re, FFTType*im ){ 
		     if( length > 1 ) GSradix3FFT(length,im, re, wre, wim  );	
		     for(u32 i = 0 ; i < length; i++){
		     	 re[i] = re[i]/length;
		     	 im[i] = im[i]/length;
		     }
  }
	
#define RADIX 3

	void Unscramble(FFTType *re, FFTType*im ){ Unscramble(RADIX , length,  re, im ); }

	private: 
		u32 length ;
		u32 Pw     ;
		FFTType *wre ;
		FFTType *wim ;
    std::vector< u32> ScrambleTable ;

    void simpleRadix3(FFTType *re, FFTType* im, u32 start, u32 stride, u32 tw2, u32 tw3)  
    {    	
		//
		// DFT length 3
		//
		//Cplx t1 = data[ind[1]] +data[ind[2]];
		//printf("faster DFT3  Temperton PFA \n");
			FFTType t1re = re[start+ stride]+ re[start + stride+ stride];
			FFTType t1im = im[start+ stride]+ im[start + stride+ stride];
		//Cplx y0 = data[ind[0]] + t1;
			FFTType y0re = re[start] + t1re;
			FFTType y0im = im[start] + t1im;
		//Cplx y1 = data[ind[0]]- 0.5 * t1;
			FFTType y1re = re[start]- (0.5 *t1re);
			FFTType y1im = im[start]- (0.5 *t1im);
		//Cplx y2 = C1*( data[ind[1]] - data[ind[2]]);
			FFTType y2re = std::sin((2*PI)/3 )*( re[start+stride]- re[start+stride+stride]);
			FFTType y2im = std::sin((2*PI)/3 )*( im[start+stride]- im[start+stride+stride]);
		//data[ind[0]] = y0;
			re[start] = y0re;
			im[start] = y0im;
		//data[ind[1]] = y1 + I*y2;
			re[start+stride] = y1re - y2im;
			im[start+stride] = y1im + y2re;		
		//data[ind[2]] = y1 - I*y2;
			re[start+stride+stride] = y1re + y2im;
			im[start+stride+stride] = y1im - y2re;		    	
  	
    	// then the twiddlefactors
    	if(tw2 != 0){
    			FFTType t = (re[start+stride] * wre[tw2]) - (im[start+stride] * wim[tw2]);
		    	im[start+stride] = (im[start+stride] * wre[tw2]) + (re[start+stride] * wim[tw2]);
		    	re[start+stride] = t;
    	}
    	if(tw3 != 0){
    			FFTType t = (re[start+stride+stride] * wre[tw3]) - (im[start+stride+stride] * wim[tw3]);
		    	im[start+stride+stride] = (im[start+stride+stride] * wre[tw3]) + (re[start+stride + stride] * wim[tw3]);
		    	re[start+stride+stride] = t;
    	}
    } 
     	
		void GSradix3FFT(u32 N,FFTType *re, FFTType*im, FFTType *wre, FFTType* wim  )
    {

       u32 k = 1;
    	 for( u32 p = 0; p < Pw ; p++){
    	 		u32 dj = 1;

    	 		for( u32 _p = 1; _p < (Pw-p+1); _p++) dj = dj*RADIX ;
    	 	
    	 		for( u32 base = 0; base < N; base += dj)
    	 		{
    	 		  u32 stride = dj/3;
	 	 	      for( u32 m = 0; m < stride; m++) 	simpleRadix3(		re, im , base+m, stride, k*m ,2*k*m);  	 		
    	 		}    	 	
	     	  k = RADIX*k;
    	 	
    	 }
	 	}	
	 	
	 		
	/*  
	    a general radix P (P is a prime) unscrambler for 
	    the Decimation in Frequency FFT GS algorithm
	    
      the first time after a change of length, it will be a little slower.  
  */
	void Unscramble(u32 radix, u32 length, FFTType *re, FFTType* im ){

		if(length >1){
		  
		  if(ScrambleTable.size() == 0 ){
		  	
		  		for( u32 i = 0 ; i < length; i++)
						if(ScrambleTable.end() == find(ScrambleTable.begin(), ScrambleTable.end(), i))
						{
							u32 jn = i ;
							u32 index =1;
							u32 jr = 0;
							while (index < length ) {
								jr = jr *radix ;
								jr += (jn % radix );
								jn = jn /radix;
	  	  				index = index + index + index;
		    			}				    			
			      	if( i != jr) { ScrambleTable.push_back(i); ScrambleTable.push_back(jr);};	  
    				}		
		  }

    	for(u32 t = 0; t < ScrambleTable.size(); t = t+2)
    	{
    		FFTType temp  = re[ScrambleTable[t]]; 
    		re[ScrambleTable[t]] = re[ScrambleTable[t+1]]; re[ScrambleTable[t+1]] = temp;
	  		temp  = im[ScrambleTable[t]]; im[ScrambleTable[t]] = im[ScrambleTable[t+1]]; 
	  		im[ScrambleTable[t+1]] = temp;
    	}
		}
	
	}
};

