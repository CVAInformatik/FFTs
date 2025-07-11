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
//#include <vector>

typedef uint32_t u32 ;
typedef double FFTType ;
const FFTType PI = 3.141592653589793238463;

class simpleRadix2IOIPFFTtype {

public: 
	simpleRadix2IOIPFFTtype(){ length = 0;};
	~simpleRadix2IOIPFFTtype(){	
		 	  delete [] wre;
		 	  delete [] wim;
	};
	
	// 
	//  SetLength(P)  ->  length = 2**P, 
	//  Status() will return the actual length
	void SetLength(u32 P)
	{
		 if (length > 0){
		 	  length = 0;
		 	  delete [] wre;
		 	  delete [] wim;
		 }
		 if( P > 0 ) {
		 		length = 1; 
		 		for(int i = 0; i < P ; i++) length = length + length;  
		 		wre = new FFTType[length/2];
		 		wim = new FFTType[length/2];
 		    for(u32 i = 0; i < length/2; i++) {
					wre[i] = std::cos(( 2*PI *i )/length);
					wim[i] = std::sin(( 2*PI *i )/length);
    		}
		 }
	}
	
	u32 Status() { return length;}
	
	void ForwardFFT(FFTType *re, FFTType*im, unsigned int stride = 1 ){ if( length > 0 ) radix2IOIPFFT(length,re, im, wre, wim, stride  );	}

	void UnscaledInverseFFT(FFTType *re, FFTType*im, unsigned int stride = 1 ){ if( length > 0 ) radix2IOIPFFT(length,im, re, wre, wim, stride  );	}

	void InverseFFT(FFTType *re, FFTType*im, unsigned int stride = 1 ){ 
		     if( length > 0 ) radix2IOIPFFT(length,im, re, wre, wim, stride  );	
		     for(u32 i = 0 ; i < length; i++){
		     	 re[i*stride] = re[i*stride]/length;
		     	 im[i*stride] = im[i*stride]/length;
		     }
  }
	
	private: 
		u32 length ;
		FFTType *wre ;
		FFTType *wim ;

		void radix2IOIPFFT(u32 N,FFTType *re, FFTType*im, FFTType *wre, FFTType* wim, unsigned int stride = 1 )
    {
 	    u32 NH = N/2;
      u32 M = 0 ;
      u32 mtemp = 1 ;
			while (mtemp < N) {mtemp = mtemp+ mtemp; M++;}
  		
  
			for( u32 L = 1; L < ((M+3)/2) ; L++ ) //30  Temperton's FORTRAN Label
			{
		 		u32 LA = 1;
		 		for( u32 i = 0; i < (L-1); i++) LA = LA+LA;
		 		u32 IA = 0;
		 		u32 IB = NH/LA;
		 		u32 KK = 0;
		 		for( u32 K = 0; K < IB; K++)//20
		 		{
				 	for(u32 I = K ; I < N; I += N/LA)//10
				 	{
				 		 if( KK == 0){ // W = (1,0)
			 	  			FFTType Zre = 0.0, Zim = 0.0;
				 	  		Zre = ( re[(IA+I)*stride] - re[(IB+I)*stride]) ;
				 	  		Zim = ( im[(IA+I)*stride] - im[(IB+I)*stride]);
	  		 	  		re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
				 	  		im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
				 	  		re[(IB+I)*stride] = Zre;
				 	  		im[(IB+I)*stride] = Zim;		 			  	
			 			  }
							else if( KK == (N/4) ){ //W = (0,1)
			 	  			FFTType Zre = 0.0, Zim = 0.0;
				 	  		Zre = -( im[(IA+I)*stride] - im[(IB+I)*stride]);
				 	  		Zim =  re[(IA+I)*stride] - re[(IB+I)*stride];
	  		 	  		re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
				 	  		im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
				 	  		re[(IB+I)*stride] = Zre;
				 	  		im[(IB+I)*stride] = Zim;
							}
							else {
		 	  				FFTType Wimag =   wim[KK];
		 	  				FFTType Wreal =   wre[KK];
		 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  		Zre = (Wreal*( re[(IA+I)*stride] - re[(IB+I)*stride])) - (Wimag* ( im[(IA+I)*stride] - im[(IB+I)*stride]));
				 	  		Zim = (Wimag*( re[(IA+I)*stride] - re[(IB+I)*stride])) + (Wreal* ( im[(IA+I)*stride] - im[(IB+I)*stride]));
	  		 	  		re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
				 	  		im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
				 	  		re[(IB+I)*stride] = Zre;
				 	  		im[(IB+I)*stride] = Zim;
				 			}
 				  }// 10
			  	KK = KK + LA;
		 		}//20
			}//30
	   for( u32 L = (M+3)/2; L <= M ; L++ ) // 70
		 	{
				u32 LA = 1;
	  		for( u32 i = 0; i < (L-1); i++) LA = LA+LA;
				u32 IA = 0;
				u32 NH = N/2;
				u32 IB = NH/LA;
				u32 IC = LA;
				u32 ID = IB + LA;
				u32 KK = 0;
				for(u32 K = 0; K < IB; K++ ) //60
				{
		 	  	for(u32 J = K; J < LA; J = J+(N/LA) ) //50
		 	  	{
		 	  		for(u32 I = J; I < N; I = I + (2*LA)) //40
		 	  		{
		 	  			if( KK == 0){ // W = (1,0)
		 	  				FFTType Zre =  re[(IA+I)*stride] - re[(IB+I)*stride];
		 	  				FFTType Zim =  im[(IA+I)*stride] - im[(IB+I)*stride];
		 	  				re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
		 	  				im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
		 	  				re[(IB+I)*stride] = re[(IC+I)*stride] + re[(ID+I)*stride];
		 	  				im[(IB+I)*stride] = im[(IC+I)*stride] + im[(ID+I)*stride];
		 	  				FFTType tre, tim ;
		 	  				tre =  re[(IC+I)*stride] - re[(ID+I)*stride];
		 	  				tim =  im[(IC+I)*stride] - im[(ID+I)*stride];
		 	  				re[(ID+I)*stride] = tre;
		 	  				im[(ID+I)*stride] = tim;
          			re[(IC+I)*stride] = Zre ;
          			im[(IC+I)*stride] = Zim ;
  			  		}
							else if( KK== (N/4)){ // W = (0,1)
		 	  				FFTType Zre = - ( im[(IA+I)*stride] - im[(IB+I)*stride]);
		 	  				FFTType Zim =  re[(IA+I)*stride] - re[(IB+I)*stride];
		 	  				re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
		 	  				im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
		 	  				re[(IB+I)*stride] = re[(IC+I)*stride] + re[(ID+I)*stride];
		 	  				im[(IB+I)*stride] = im[(IC+I)*stride] + im[(ID+I)*stride];
		 	  				FFTType tre, tim ;
		 	  				tre =  - ( im[(IC+I)*stride] - im[(ID+I)*stride]);
		 	  				tim =   re[(IC+I)*stride] - re[(ID+I)*stride];
		 	  				re[(ID+I)*stride] = tre;
		 	  				im[(ID+I)*stride] = tim;
          			re[(IC+I)*stride] = Zre ;
          			im[(IC+I)*stride] = Zim ;
							}
							else {
					 	  	FFTType Wimag =   wim[KK];
		 				  	FFTType Wreal =   wre[KK];
		 	  				FFTType Zre = (Wreal*( re[(IA+I)*stride] - re[(IB+I)*stride])) - (Wimag* ( im[(IA+I)*stride] - im[(IB+I)*stride]));
		 	  				FFTType Zim = (Wimag*( re[(IA+I)*stride] - re[(IB+I)*stride])) + (Wreal* ( im[(IA+I)*stride] - im[(IB+I)*stride]));
		 	  				re[(IA+I)*stride] = re[(IA+I)*stride] + re[(IB+I)*stride];
		 	  				im[(IA+I)*stride] = im[(IA+I)*stride] + im[(IB+I)*stride];
		 	  				re[(IB+I)*stride] = re[(IC+I)*stride] + re[(ID+I)*stride];
		 	  				im[(IB+I)*stride] = im[(IC+I)*stride] + im[(ID+I)*stride];
		 	  				FFTType tre, tim ;
		 	  				tre = (Wreal*( re[(IC+I)*stride] - re[(ID+I)*stride])) - (Wimag* ( im[(IC+I)*stride] - im[(ID+I)*stride]));
		 	  				tim =  (Wimag*( re[(IC+I)*stride] - re[(ID+I)*stride])) + (Wreal* ( im[(IC+I)*stride] - im[(ID+I)*stride]));
		 	  				re[(ID+I)*stride] = tre;
		 	  				im[(ID+I)*stride] = tim;
          			re[(IC+I)*stride] = Zre ;
          			im[(IC+I)*stride] = Zim ;
          		}
		 	  		}//40
		 	  	}// 50			
		 	  	KK = KK + LA;
				}//60
			}//70
		}	
};

