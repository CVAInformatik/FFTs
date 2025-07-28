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
		
		I've transcribed the original FORTRAN code to simple C++ with a few obvious optimizations.
		I've retained the original variable names and loop numbers for easy back reference.
*/

#include <iostream>
#include <math.h>
#include <algorithm>
//#include <vector>

#include "FFTBase.h"

class IndexedSimpleRadix2IOIPFFTtype : public FFTBaseClass {

public: 
	IndexedSimpleRadix2IOIPFFTtype(){ 
		    length = 0;
		 	  wre = 0;
		 	  wim = 0;
		 	  indexMap = 0;
		 	  rotationMap = 0;		    
		    };
	~IndexedSimpleRadix2IOIPFFTtype(){	
		 	  if( wre) delete [] wre;
		 	  if( wim) delete [] wim;
		 	  if(indexMap) delete [] indexMap;
		 	  if( rotationMap) delete [] rotationMap;
	};
	
	// 
	//  SetLength(P)  ->  length = 2**P, 
	//  Status() will return the actual length
	void SetLength(u32 P)
	{
  	if(wre) delete [] wre;
		if(wim) delete [] wim;
		if(indexMap) delete [] indexMap;
		if(rotationMap) delete [] rotationMap;
			
		if (length > 0){
		 	  length = 0;
		}
		if( P > 2 ) {
		 		length = 1; 
		 		for(int i = 0; i < P ; i++) length = length + length;  

		 		wre = new FFTType[length/2];
		 		wim = new FFTType[length/2];
 		    for(u32 i = 0; i < length/2; i++) {
					wre[i] = std::cos((  2*PI *i )/length);
					wim[i] = std::sin((  2*PI *i )/length);
    		}
    		indexMap = new unsigned int[length];
    		for(u32 i = 0; i < length; i++) indexMap[i] = i;
    		rotationMap = new unsigned int[length];
		}
	}
	
	u32 Status() { return length;}
	
	void ForwardFFT(FFTType *re, FFTType*im, unsigned int stride = 1 )
	{ 
		if( length > 0 ) {
			for (int i = 0; i < length; i++) indexMap[i] = i*stride  ;
			radix2IOIPFFT(length,re, im, wre, wim, indexMap  );	
		}
	}
	
	void ForwardFFT(FFTType *re, FFTType*im, unsigned int *indexMap ){
		 if( length > 0 ) radix2IOIPFFT(length,re, im, wre, wim, indexMap );
	}

	void UnscaledInverseFFT(FFTType *re, FFTType*im, unsigned int stride = 1 )
	{ 
		if( length > 0 ) {
			for (int i = 0; i < length; i++) indexMap[i] = i*stride  ;
			radix2IOIPFFT(length,im, re, wre, wim, indexMap  );	
		}
	}
	
	void UnscaledInverseFFT(FFTType *re, FFTType*im, unsigned int *indexMap){
		 if( length > 0 ) radix2IOIPFFT(length,im, re, wre, wim, indexMap  );
	}

	void InverseFFT(FFTType *re, FFTType*im, unsigned int stride = 1 )
	{ 
		for (int i = 0; i < length; i++) indexMap[i] = i*stride  ;
  	InverseFFT(re, im,indexMap  );	
  }
  
	void InverseFFT(FFTType *re, FFTType*im, unsigned int *indexMap )
	{ 
     if( length > 0 ){ radix2IOIPFFT(length,im, re, wre, wim, indexMap  );	
		     for(u32 i = 0 ; i < length; i++){
		     	 re[indexMap[i]] = re[indexMap[i]]/length;
		     	 im[indexMap[i]] = im[indexMap[i]]/length;
		     }
		}
  }

 	void rotate(u32 rotation, FFTType *re, FFTType*im, unsigned int stride = 1 )
	{ 
		for (int i = 0; i < length; i++) indexMap[i] = i*stride  ;
	 	rotate(rotation, re, im, indexMap )	;	
	}

 	void rotate(u32 rotation, FFTType *re, FFTType*im, unsigned int *indexMap )
	{ 
		if(rotation == 1) return; //noop
		if(0 == (rotation % 2)) {
			std::cerr << "rotation must be odd " << rotation << std::endl;
		}
		else {
			u32 rot = rotation % length;
			// build map
			for( u32 ix = 0 ; ix < length; ix++) 
			      rotationMap[ix] = (rot * ix) % length;
			 	
			// use (and destroy) map      
			for( u32 ix = 0 ; ix < length; ix++) 
			{ 
				
				if( rotationMap[ix] == 0 ) continue;
				else if (rotationMap[ix] == ix) {
					  rotationMap[ix] = 0;
					  continue;
					}
				else if (rotationMap[rotationMap[ix]] == ix) { 
						// a simple swap 
					FFTType tre ;
					FFTType tim ;
					tre =  re[indexMap[rotationMap[ix]]];
					tim =  im[indexMap[rotationMap[ix]]];
          re[indexMap[rotationMap[ix]]] = re[indexMap[ix]] ;     
          im[indexMap[rotationMap[ix]]] = im[indexMap[ix]] ;     
          re[indexMap[ix]] = tre;
          im[indexMap[ix]] = tim;
					rotationMap[rotationMap[ix]] = 0;
					rotationMap[ix] = 0;
				}
				else {
					FFTType tre ;
					FFTType tim ;
					u32 startIndex = ix;					
					u32 destinationIndex = startIndex;
					u32 sourceIndex = rotationMap[destinationIndex];
					tre =  re[indexMap[startIndex]];
					tim =  im[indexMap[startIndex]];
					do {
					  re[indexMap[destinationIndex]] = re[indexMap[sourceIndex]];
					  im[indexMap[destinationIndex]] = im[indexMap[sourceIndex]];
						rotationMap[destinationIndex] = 0;
						destinationIndex = sourceIndex ;
						sourceIndex = rotationMap[sourceIndex];
					} while (destinationIndex && (sourceIndex != startIndex ));
					re[indexMap[destinationIndex]] = tre; 
					im[indexMap[destinationIndex]] = tim;
					rotationMap[destinationIndex] = 0;
				}
			}
		}
	}

	
	private: 
		u32 length ;

		FFTType *wre ;
		FFTType *wim ;
		unsigned int *indexMap ;
		unsigned int *rotationMap ;

		void radix2IOIPFFT(u32 N,FFTType *re, FFTType*im, FFTType *wre, FFTType* wim, unsigned int *indexMap)
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
				 	  		Zre = ( re[indexMap[(IA+I)] ] - re[indexMap[(IB+I)]]) ;
				 	  		Zim = ( im[indexMap[(IA+I)] ] - im[indexMap[(IB+I)]]) ;
	  		 	  		re[indexMap[IA+I]] = re[indexMap[(IA+I)]] + re[indexMap[(IB+I)]];
				 	  		im[indexMap[IA+I]]  = im[indexMap[(IA+I)]] + im[indexMap[(IB+I)]];
				 	  	  re[indexMap[IB+I]] = Zre;
				 	  		im[indexMap[IB+I]] = Zim;		 			  	
			 			  }
							else if( KK == (N/4) ){ //W = (0,1)
			 	  			FFTType Zre = 0.0, Zim = 0.0;
				 	  		Zre = -( im[indexMap[IA+I]] - im[indexMap[IB+I]]);
				 	  		Zim =  re[indexMap[IA+I]] - re[indexMap[IB+I]];
	  		 	  	re[indexMap[(IA+I)]] = re[indexMap[(IA+I)]] + re[indexMap[IB+I]];
				 	  		im[indexMap[(IA+I)]] = im[indexMap[IA+I]] + im[indexMap[IB+I]];
				 	  		re[indexMap[(IB+I)]] = Zre;
				 	  		im[indexMap[(IB+I)] ] = Zim;
							}
							else {
		 	  				FFTType Wimag =   wim[KK];
		 	  				FFTType Wreal =   wre[KK];
		 	  				FFTType Zre = 0.0, Zim = 0.0;
				 	  		Zre = (Wreal*( re[indexMap[IA+I]] - re[indexMap[IB+I]])) - (Wimag* ( im[indexMap[IA+I]] - im[indexMap[IB+I]]));
				 	  		Zim = (Wimag*( re[indexMap[IA+I]] - re[indexMap[IB+I]])) + (Wreal* ( im[indexMap[IA+I]] - im[indexMap[IB+I]]));
	  		 	  		re[indexMap[(IA+I)] ] = re[indexMap[(IA+I)] ] + re[indexMap[(IB+I) ]];
				 	  		im[indexMap[(IA+I)] ] = im[indexMap[(IA+I)] ] + im[indexMap[(IB+I) ]];
				 	  		re[indexMap[(IB+I)] ] = Zre;
				 	  		im[indexMap[(IB+I)] ] = Zim;
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
		 	  				FFTType Zre =  re[indexMap[IA+I]] - re[indexMap[IB+I]];
		 	  				FFTType Zim =  im[indexMap[IA+I]] - im[indexMap[IB+I]];
		 	  				re[indexMap[IA+I]] = re[indexMap[IA+I]] + re[indexMap[IB+I]];
		 	  				im[indexMap[IA+I]] = im[indexMap[IA+I]] + im[indexMap[IB+I]];
		 	  				re[indexMap[IB+I]] = re[indexMap[IC+I]] + re[indexMap[ID+I]];
		 	  				im[indexMap[IB+I]] = im[indexMap[IC+I]] + im[indexMap[ID+I]];
		 	  				FFTType tre, tim ;
		 	  				tre =  re[indexMap[IC+I]] - re[indexMap[ID+I]];
		 	  				tim =  im[indexMap[IC+I]] - im[indexMap[ID+I]];
		 	  				re[indexMap[ID+I]] = tre;
		 	  				im[indexMap[ID+I]] = tim;
          			re[indexMap[IC+I]] = Zre ;
          			im[indexMap[IC+I]] = Zim ;
  			  		}
							else if( KK== (N/4)){ // W = (0,1)
		 	  				FFTType Zre = - ( im[indexMap[IA+I]] - im[indexMap[IB+I]]);
		 	  				FFTType Zim =     re[indexMap[IA+I]] - re[indexMap[IB+I]];
		 	  				re[indexMap[(IA+I)]] = re[indexMap[IA+I]] + re[indexMap[IB+I]];
		 	  				im[indexMap[(IA+I)]] = im[indexMap[IA+I]] + im[indexMap[IB+I]];
		 	  				re[indexMap[(IB+I)]] = re[indexMap[IC+I]] + re[indexMap[ID+I]];
		 	  				im[indexMap[(IB+I)]] = im[indexMap[IC+I]] + im[indexMap[ID+I]];
		 	  				FFTType tre, tim ;
		 	  				tre =  - ( im[indexMap[IC+I]] - im[indexMap[ID+I]]);
		 	  				tim =   re[indexMap[IC+I]] - re[indexMap[ID+I]];
		 	  				re[indexMap[ID+I]] = tre;
		 	  				im[indexMap[ID+I]] = tim;
          			re[indexMap[IC+I]] = Zre ;
          			im[indexMap[IC+I]] = Zim ;
							}
							else {
					 	  	FFTType Wimag =   wim[KK];
		 				  	FFTType Wreal =   wre[KK];
		 	  				FFTType Zre = (Wreal*( re[indexMap[IA+I]] - re[indexMap[IB+I]])) - (Wimag* ( im[indexMap[IA+I]] - im[indexMap[IB+I]]));
		 	  				FFTType Zim = (Wimag*( re[indexMap[IA+I]] - re[indexMap[IB+I]])) + (Wreal* ( im[indexMap[IA+I]] - im[indexMap[IB+I]]));
		 	  				re[indexMap[IA+I]] = re[indexMap[IA+I]] + re[indexMap[IB+I]];
		 	  				im[indexMap[IA+I]] = im[indexMap[IA+I]] + im[indexMap[IB+I]];
		 	  				re[indexMap[IB+I]] = re[indexMap[IC+I]] + re[indexMap[ID+I]];
		 	  				im[indexMap[IB+I]] = im[indexMap[IC+I]] + im[indexMap[ID+I]];
		 	  				FFTType tre, tim ;
		 	  				tre = (Wreal*( re[indexMap[IC+I]] - re[indexMap[ID+I]])) - (Wimag* ( im[indexMap[IC+I]] - im[indexMap[ID+I]]));
		 	  				tim = (Wimag*( re[indexMap[IC+I]] - re[indexMap[ID+I]])) + (Wreal* ( im[indexMap[IC+I]] - im[indexMap[ID+I]]));
		 	  				re[indexMap[ID+I]] = tre;
		 	  				im[indexMap[ID+I]] = tim;
          			re[indexMap[IC+I]] = Zre ;
          			im[indexMap[IC+I]] = Zim ;
          		}
		 	  		}//40
		 	  	}// 50			
		 	  	KK = KK + LA;
				}//60
			}//70
		}	

};

