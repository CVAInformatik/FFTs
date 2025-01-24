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

#define RADIX 5


class simpleGSRadix5FFTtype {

public: 
	simpleGSRadix5FFTtype(){ Pw = 0; length = 1;};
	~simpleGSRadix5FFTtype(){	
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
		 		for(int i = 0; i < P ; i++) length = 5*length ;  
		 		wre = new FFTType[length];
		 		wim = new FFTType[length];
 		    for(u32 i = 0; i < length; i++) {
					wre[i] = std::cos(( 2*PI *i )/length);
					wim[i] = std::sin(( 2*PI *i )/length);
    		}
		 }
	};
	
	u32 Status() { return length;};
	
	void ForwardFFT(FFTType *re, FFTType*im ){ 		 if( length > 1 ) GSradix5FFT(length,re, im, wre, wim  );		}

	void UnscaledInverseFFT(FFTType *re, FFTType*im ){ if( length > 1 ) GSradix5FFT(length,im, re, wre, wim  );	}

	void InverseFFT(FFTType *re, FFTType*im ){ 
		     if( length > 1 ) GSradix5FFT(length,im, re, wre, wim  );	
		     for(u32 i = 0 ; i < length; i++){
		     	 re[i] = re[i]/length;
		     	 im[i] = im[i]/length;
		     }
  }
	

	void Unscramble(FFTType *re, FFTType*im ){ Unscramble(RADIX , length,  re, im ); }

	private: 
		u32 length ;
		u32 Pw     ;
		FFTType *wre ;
		FFTType *wim ;
    std::vector< u32> ScrambleTable ;

    void simpleRadix5(FFTType *re, FFTType* im, u32 start, u32 stride, u32 tw2, u32 tw3, u32 tw4, u32 tw5)  
    { 
    /*	based on 	Clive Temperton Journal of Computational Physics 75, 190-198 (1988) +
				Corrigendum Journ. of Comp.Physics 79, 486(1988) 
				"A New Set of Minimum-Add small-n Rotated DFT Modules"   	
		*/
		//
		// DFT length 5
		//printf("start %3d  stride %3d  tw2 %3d tw3 %3d tw4 %3d tw5 %3d\n",start, stride, tw2, tw3, tw4, tw5 );
		//
		FFTType Phi = 2.0 * M_PI/5.0;
		FFTType	C1 = 0.5*(std::cos(Phi) - std::cos(2.0 * Phi));
		FFTType C2 = std::sin( Phi );
		FFTType C3 = std::sin(2 * Phi);

		//Cplx t1 = data[ind[1]] + data[ind[4]];
		FFTType t1re = re[start+stride] +re[start+(4*stride)];
		FFTType t1im = im[start+stride] +im[start+(4*stride)];
		//Cplx t2 = data[ind[2]] + data[ind[3]];
		FFTType t2re = re[start+(2*stride)] +re[start+(3*stride)];
		FFTType t2im = im[start+(2*stride)] +im[start+(3*stride)];
		//Cplx t3 = data[ind[1]] - data[ind[4]];
		FFTType t3re = re[start+stride] - re[start+(4*stride)];
		FFTType t3im = im[start+stride] -im[start+(4*stride)];
		//Cplx t4 = data[ind[2]] - data[ind[3]];
		FFTType t4re = re[start+(2*stride)] - re[start+(3*stride)];
		FFTType t4im = im[start+(2*stride)] - im[start+(3*stride)];
		//Cplx t5 = t1+t2;
		FFTType t5re = t1re + t2re;
		FFTType t5im = t1im + t2im;
		//Cplx t6 = C1*(t1-t2);
		FFTType t6re = C1 * (t1re - t2re);
		FFTType t6im = C1 * (t1im - t2im);		
		//Cplx t7 = data[ind[0]] - 0.25*t5;
		FFTType t7re = re[start] - 0.25*t5re ;
		FFTType t7im = im[start] - 0.25*t5im ;
		//Cplx y0 = data[ind[0]] + t5;
		FFTType y0re = re[start] + t5re ;
		FFTType y0im = im[start] + t5im ;
		//Cplx y1 = t7+t6;
		FFTType y1re = t7re + t6re ;
		FFTType y1im = t7im + t6im ;
		//Cplx y2 = t7-t6;
		FFTType y2re = t7re - t6re ;
		FFTType y2im = t7im - t6im ;
		//Cplx y3 = C3*t3 - C2*t4;
		FFTType y3re = (C3 * t3re) - (C2 * t4re) ;
		FFTType y3im = (C3 * t3im) - ( C2 * t4im) ;		
		//Cplx y4 = C2*t3 + C3*t4;
		FFTType y4re = (C2 * t3re) + (C3 * t4re) ;
		FFTType y4im = (C2 * t3im)  + (C3 * t4im) ;		
		re[start] =  y0re;
		im[start] =  y0im;
		re[start+stride] =  y1re - y4im;
		im[start+stride] =  y1im + y4re;
		re[start+(2*stride)] =  y2re -  y3im;
		im[start+(2*stride)] =  y2im +  y3re;
		re[start+(3*stride)] =	y2re +  y3im;
		im[start+(3*stride)] =	y2im - y3re;
		re[start+(4*stride)] =  y1re + y4im;
		im[start+(4*stride)] =  y1im - y4re;
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
    	if(tw4 != 0){
    			FFTType t = (re[start+(3*stride)] * wre[tw4]) - (im[start+(3*stride)] * wim[tw4]);
		    	im[start+(3*stride)] = (im[start+(3*stride)] * wre[tw4]) + (re[start+(3*stride)] * wim[tw4]);
		    	re[start+(3*stride)] = t;
    	}
    	if(tw5 != 0){
    			FFTType t = (re[start+(4*stride)] * wre[tw5]) - (im[start+(4*stride)] * wim[tw5]);
		    	im[start+(4*stride)] = (im[start+(4*stride)] * wre[tw5]) + (re[start+(4*stride)] * wim[tw5]);
		    	re[start+(4*stride)] = t;
    	}

    } 
     	
		void GSradix5FFT(u32 N,FFTType *re, FFTType*im, FFTType *wre, FFTType* wim  )
    {

       u32 k = 1;
    	 for( u32 p = 0; p < Pw ; p++){
    	 		u32 dj = 1;

    	 		for( u32 _p = 1; _p < (Pw-p+1); _p++) dj = dj*RADIX ;
    	 	
    	 		for( u32 base = 0; base < N; base += dj)
    	 		{
    	 		  u32 stride = dj/RADIX;
	 	 	      for( u32 m = 0; m < stride; m++) 	simpleRadix5(		re, im , base+m, stride, k*m ,2*k*m, 3*k*m,4*k*m);  	 		
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
	  	  				index = radix*index ;
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
	  		//printf("swapping %3d:    %3d and %3d\n",t/2, ScrambleTable[t], ScrambleTable[t+1] );
    	}
		}
	
	}
};

