#pragma once
/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/
#include "FFTBase.h"
#include "simpleRadix2IOIPFFTindex.h"
#include "simple32NIOIPFFT.h"


class FFTClass
{
	public:   
		   FFTClass(){ FFTBaseClass  *fft = new IndexedSimpleRadix2IOIPFFTtype();};
		   ~FFTClass() {delete fft ;};
	     u32 Status() {return fft->Status();};
	     void ForwardFFT(FFTType *re, FFTType *im, unsigned int stride = 1 )
	     { fft->ForwardFFT(re, im, stride ); };
	     void ForwardFFT(FFTType *re, FFTType *im, unsigned int *indexMap )
	     {  fft->ForwardFFT(re, im, indexMap ); };
	     void UnscaledInverseFFT(FFTType *re, FFTType *im, unsigned int stride = 1)
	     { fft->UnscaledInverseFFT(re, im, stride);};
	     void UnscaledInverseFFT(FFTType *re, FFTType *im, unsigned int *indexMap) 
	     { fft->UnscaledInverseFFT(re, im, indexMap);};
	     void InverseFFT(FFTType *re, FFTType *im, unsigned int stride = 1)
	     { fft->InverseFFT(re, im, stride); };
	     void InverseFFT(FFTType *re, FFTType *im, unsigned int *indexMap)
	     { fft->InverseFFT(re, im, indexMap);};

        // to set the length of the FFT
        //  NB this is the actual length not some power
        //  valid values  8,16,24,32,48,64,96,128,192,256,....
        u32 Length( u32 _len){
        	 u32 P = 0;
        	 u32 len = _len;
        	 delete fft;
        	 while(len && (len& 1) == 0){  P++, len = len>> 1; }
        	 if (len == 3) {
        	 	  fft = new simple32NIOIPFFTtype();
        	 	  fft-> SetLength(P); 
        	 } else if ( len == 1){
        	 	  fft = new IndexedSimpleRadix2IOIPFFTtype();
        	 	  fft-> SetLength(P);
        	 }
        	 else fft = new IndexedSimpleRadix2IOIPFFTtype();
        	 return fft->Status();
        };
	   
	private:
		FFTBaseClass  *fft ; // will always point at something 
};
