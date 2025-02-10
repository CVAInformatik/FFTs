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
		
		this version, where Tempertons two for loops have been collected into one with a test for which part
		we are in, happened during porting the FFT to Ocaml, 
		
		The stride parameter, with default value 1, can be used in various ways
		
		Stride = 1  :  Just a complex signal stored as two independent vectors
		Stride = 2  :  A complex signal, real and imaginary stored interleaved, as you would get with a struct
		strides = M and N  :  A complex 2D MxN signal, stored as two linearized independent vectors, e.g.
		                N = 16 M = 32 
		                0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17.....
		                interpreted as the picture 
		                 0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
		                16,17,.....                         ,31, 
		                32,33,.....                         ,47, 
		                .........
		                496,497............................ ,511

                    in this case you would make M= 32 length 16= N FFTs with stride  1 and
                                                N =16 length 32 =M FFTs with stride 16 = N
                                                
*/

#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>

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
	
	void ForwardFFT(FFTType *re, FFTType*im, u32 stride = 1 ){ 
		   if( length > 0 ) 
		   	  if( stride == 1) radix2IOIPFFT(length,re, im, wre, wim  );
		   	  else radix2IOIPFFT(stride, length,re, im, wre, wim  );		
		   	}

	void UnscaledInverseFFT(FFTType *re, FFTType*im, u32 stride = 1 ){ 
		      if( length > 0 ) 
		      	if( stride == 1) radix2IOIPFFT(length,im, re, wre, wim  );	
		      	else radix2IOIPFFT(stride, length,im, re, wre, wim  );	
		     }

	void InverseFFT(FFTType *re, FFTType*im, u32 stride = 1 ){ 
		     if( length > 0 ) 
		     	  if ( stride == 1)
		     	  	 radix2IOIPFFT(length,im, re, wre, wim  );
		     	  else 	
		     	  	 radix2IOIPFFT(stride, length,im, re, wre, wim  );
		     	  if ( stride == 1)
					     for(u32 i = 0 ; i < length; i++){
		  			   	 	re[i] = re[i]/length;
		     	 				im[i] = im[i]/length;
		     			 }
		     		 else 
		     		 	 for(u32 i = 0 ; i < length; i++){
		  			   	 	re[stride * i] = re[stride * i]/length;
		     	 				im[stride * i] = im[stride * i]/length;
		     			 }
  }
	
	private: 
		u32 length ;
		FFTType *wre ;
		FFTType *wim ;

		void radix2IOIPFFT(u32 N,FFTType *re, FFTType*im, FFTType *wre, FFTType* wim  )
    {
 	    u32 NH = N/2;
      u32 M = 0 ;
      u32 mtemp = 1 ;
			while (mtemp < N) {mtemp = mtemp+ mtemp; M++;}

	 		u32 LA = 1;
     	for( u32 L_ = 1 ; L_ <= M ; L_++ ) //30+70  Temperton's FORTRAN Label
     	{     	
		 		if( L_ - 1 ) LA = LA+LA;
		 		u32 IA = 0;
		 		u32 IB = NH/LA;
		 		u32 KK = 0;
		 		
     		if( (L_) <  ((M+3)/2))
				{
		 			for( u32 K = 0; K < IB; K++)//20
		 			{
					 	for(u32 I = K ; I < N; I += N/LA)//10
				 		{
				 			 	if( KK == 0){ // W = (1,0)
			 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = ( re[IA+I] - re[IB+I]) ;
				 	  			Zim = ( im[IA+I] - im[IB+I]);
	  		 	  			re[IA+I] = re[IA+I] + re[IB+I];
				 	  			im[IA+I] = im[IA+I] + im[IB+I];
				 	  			re[IB+I] = Zre;
				 	  			im[IB+I] = Zim;		 			  	
			 			   	}
							 	else if( KK == (N/4) ){ //W = (0,1)
			 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = -( im[IA+I] - im[IB+I]);
				 	  			Zim =  re[IA+I] - re[IB+I];
	  		 	  			re[IA+I] = re[IA+I] + re[IB+I];
				 	  			im[IA+I] = im[IA+I] + im[IB+I];
				 	  			re[IB+I] = Zre;
				 	  			im[IB+I] = Zim;
								}
								else {
		 	  					FFTType Wimag =   wim[KK];
		 	  					FFTType Wreal =   wre[KK];
		 	  					FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = (Wreal*( re[IA+I] - re[IB+I])) - (Wimag* ( im[IA+I] - im[IB+I]));
				 	  			Zim = (Wimag*( re[IA+I] - re[IB+I])) + (Wreal* ( im[IA+I] - im[IB+I]));
	  		 	  			re[IA+I] = re[IA+I] + re[IB+I];
				 	  			im[IA+I] = im[IA+I] + im[IB+I];
				 	  			re[IB+I] = Zre;
				 	  			im[IB+I] = Zim;
				 				}
 				  	}// 10
			  		KK = KK + LA;
		 			}//20
				}//30
	   		else // 70
		 		{
					u32 IC = LA;
					u32 ID = IB + LA;
					for(u32 K = 0; K < IB; K++ ) //60
					{
			 	  	for(u32 J = K; J < LA; J = J+(N/LA) ) //50
			 	  	{
			 	  		for(u32 I = J; I < N; I = I + (2*LA)) //40
			 	  		{
			 	  			if( KK == 0){ // W = (1,0)
			 	  				FFTType Zre =  re[IA+I] - re[IB+I];
		 		  				FFTType Zim =  im[IA+I] - im[IB+I];
		 	  					re[IA+I] = re[IA+I] + re[IB+I];
		 	  					im[IA+I] = im[IA+I] + im[IB+I];
		 	  					re[IB+I] = re[IC+I] + re[ID+I];
		 	  					im[IB+I] = im[IC+I] + im[ID+I];
		 	  					FFTType tre, tim ;
		 	  					tre =  re[IC+I] - re[ID+I];
		 	  					tim =  im[IC+I] - im[ID+I];
		 	  					re[ID+I] = tre;
		 	  					im[ID+I] = tim;
          				re[IC+I] = Zre ;
          				im[IC+I] = Zim ;
  			  			}
								else if( KK== (N/4)){ // W = (0,1)
		 	  					FFTType Zre = - ( im[IA+I] - im[IB+I]);
		 	  					FFTType Zim =  re[IA+I] - re[IB+I];
		 	  					re[IA+I] = re[IA+I] + re[IB+I];
		 	  					im[IA+I] = im[IA+I] + im[IB+I];
		 	  					re[IB+I] = re[IC+I] + re[ID+I];
		 	  					im[IB+I] = im[IC+I] + im[ID+I];
		 	  					FFTType tre, tim ;
		 	  					tre =  - ( im[IC+I] - im[ID+I]);
		 	  					tim =   re[IC+I] - re[ID+I];
		 	  					re[ID+I] = tre;
		 	  					im[ID+I] = tim;
          				re[IC+I] = Zre ;
          				im[IC+I] = Zim ;
								}
								else {
					 	  		FFTType Wimag =   wim[KK];
		 				  		FFTType Wreal =   wre[KK];
		 	  					FFTType Zre = (Wreal*( re[IA+I] - re[IB+I])) - (Wimag* ( im[IA+I] - im[IB+I]));
		 	  					FFTType Zim = (Wimag*( re[IA+I] - re[IB+I])) + (Wreal* ( im[IA+I] - im[IB+I]));
		 	  					re[IA+I] = re[IA+I] + re[IB+I];
		 	  					im[IA+I] = im[IA+I] + im[IB+I];
		 	  					re[IB+I] = re[IC+I] + re[ID+I];
		 	  					im[IB+I] = im[IC+I] + im[ID+I];
		 	  					FFTType tre, tim ;
		 	  					tre = (Wreal*( re[IC+I] - re[ID+I])) - (Wimag* ( im[IC+I] - im[ID+I]));
		 	  					tim =  (Wimag*( re[IC+I] - re[ID+I])) + (Wreal* ( im[IC+I] - im[ID+I]));
		 	  					re[ID+I] = tre;
		 	  					im[ID+I] = tim;
          				re[IC+I] = Zre ;
          				im[IC+I] = Zim ;
          			}
		 	  			}//40
		 	  		}// 50			
		 	  		KK = KK + LA;
					}//60
				}//70
			}	
  	}
  	
  	
		void radix2IOIPFFT(u32 stride, u32 N, FFTType *re, FFTType*im, FFTType *wre, FFTType* wim  )
    {
 	    u32 NH = N/2;
      u32 M = 0 ;
      u32 mtemp = 1 ;
			while (mtemp < N) {mtemp = mtemp+ mtemp; M++;}

	 		u32 LA = 1;
     	for( u32 L_ = 1 ; L_ <= M ; L_++ ) //30+70  Temperton's FORTRAN Label
     	{     	
		 		if( L_ - 1 ) LA = LA+LA;
		 		u32 IA = 0;
		 		u32 IB = NH/LA;
		 		u32 KK = 0;
		 		
     		if( (L_) <  ((M+3)/2))
				{
		 			for( u32 K = 0; K < IB; K++)//20
		 			{
					 	for(u32 I = K ; I < N; I += N/LA)//10
				 		{
				 			 	if( KK == 0){ // W = (1,0)
			 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = ( re[stride*(IA+I)] - re[stride*(IB+I)]) ;
				 	  			Zim = ( im[stride*(IA+I)] - im[stride*(IB+I)]);
	  		 	  			re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
				 	  			im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
				 	  			re[stride*(IB+I)] = Zre;
				 	  			im[stride*(IB+I)] = Zim;		 			  	
			 			   	}
							 	else if( KK == (N/4) ){ //W = (0,1)
			 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = -( im[stride*(IA+I)] - im[stride*(IB+I)]);
				 	  			Zim =  re[stride*(IA+I)] - re[stride*(IB+I)];
	  		 	  			re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
				 	  			im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
				 	  			re[stride*(IB+I)] = Zre;
				 	  			im[stride*(IB+I)] = Zim;
								}
								else {
		 	  					FFTType Wimag =   wim[KK];
		 	  					FFTType Wreal =   wre[KK];
		 	  					FFTType Zre = 0.0, Zim = 0.0;
				 	  			Zre = (Wreal*( re[stride*(IA+I)] - re[stride*(IB+I)])) - (Wimag* ( im[stride*(IA+I)] - im[stride*(IB+I)]));
				 	  			Zim = (Wimag*( re[stride*(IA+I)] - re[stride*(IB+I)])) + (Wreal* ( im[stride*(IA+I)] - im[stride*(IB+I)]));
	  		 	  			re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
				 	  			im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
				 	  			re[stride*(IB+I)] = Zre;
				 	  			im[stride*(IB+I)] = Zim;
				 				}
 				  	}// 10
			  		KK = KK + LA;
		 			}//20
				}//30
	   		else // 70
		 		{
					u32 IC = LA;
					u32 ID = IB + LA;
					for(u32 K = 0; K < IB; K++ ) //60
					{
			 	  	for(u32 J = K; J < LA; J = J+(N/LA) ) //50
			 	  	{
			 	  		for(u32 I = J; I < N; I = I + (2*LA)) //40
			 	  		{
			 	  			if( KK == 0){ // W = (1,0)
			 	  				FFTType Zre =  re[stride*(IA+I)] - re[stride*(IB+I)];
		 		  				FFTType Zim =  im[stride*(IA+I)] - im[stride*(IB+I)];
		 	  					re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
		 	  					im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
		 	  					re[stride*(IB+I)] = re[stride*(IC+I)] + re[stride*(ID+I)];
		 	  					im[stride*(IB+I)] = im[stride*(IC+I)] + im[stride*(ID+I)];
		 	  					FFTType tre, tim ;
		 	  					tre =  re[stride*(IC+I)] - re[stride*(ID+I)];
		 	  					tim =  im[stride*(IC+I)] - im[stride*(ID+I)];
		 	  					re[stride*(ID+I)] = tre;
		 	  					im[stride*(ID+I)] = tim;
          				re[stride*(IC+I)] = Zre ;
          				im[stride*(IC+I)] = Zim ;
  			  			}
								else if( KK== (N/4)){ // W = (0,1)
		 	  					FFTType Zre = - ( im[stride*(IA+I)] - im[stride*(IB+I)]);
		 	  					FFTType Zim =  re[stride*(IA+I)] - re[stride*(IB+I)];
		 	  					re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
		 	  					im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
		 	  					re[stride*(IB+I)] = re[stride*(IC+I)] + re[stride*(ID+I)];
		 	  					im[stride*(IB+I)] = im[stride*(IC+I)] + im[stride*(ID+I)];
		 	  					FFTType tre, tim ;
		 	  					tre =  - ( im[stride*(IC+I)] - im[stride*(ID+I)]);
		 	  					tim =   re[stride*(IC+I)] - re[stride*(ID+I)];
		 	  					re[stride*(ID+I)] = tre;
		 	  					im[stride*(ID+I)] = tim;
          				re[stride*(IC+I)] = Zre ;
          				im[stride*(IC+I)] = Zim ;
								}
								else {
					 	  		FFTType Wimag =   wim[KK];
		 				  		FFTType Wreal =   wre[KK];
		 	  					FFTType Zre = (Wreal*( re[stride*(IA+I)] - re[stride*(IB+I)])) - (Wimag* ( im[stride*(IA+I)] - im[stride*(IB+I)]));
		 	  					FFTType Zim = (Wimag*( re[stride*(IA+I)] - re[stride*(IB+I)])) + (Wreal* ( im[stride*(IA+I)] - im[stride*(IB+I)]));
		 	  					re[stride*(IA+I)] = re[stride*(IA+I)] + re[stride*(IB+I)];
		 	  					im[stride*(IA+I)] = im[stride*(IA+I)] + im[stride*(IB+I)];
		 	  					re[stride*(IB+I)] = re[stride*(IC+I)] + re[stride*(ID+I)];
		 	  					im[stride*(IB+I)] = im[stride*(IC+I)] + im[stride*(ID+I)];
		 	  					FFTType tre, tim ;
		 	  					tre = (Wreal*( re[stride*(IC+I)] - re[stride*(ID+I)])) - (Wimag* ( im[stride*(IC+I)] -  im[stride*(ID+I)]));
		 	  					tim =  (Wimag*( re[stride*(IC+I)] - re[stride*(ID+I)])) + (Wreal* ( im[stride*(IC+I)] - im[stride*(ID+I)]));
		 	  					re[stride*(ID+I)] = tre;
		 	  					im[stride*(ID+I)] = tim;
          				re[stride*(IC+I)] = Zre ;
          				im[stride*(IC+I)] = Zim ;
          			}
		 	  			}//40
		 	  		}// 50			
		 	  		KK = KK + LA;
					}//60
				}//70
			}	
  	}

};

