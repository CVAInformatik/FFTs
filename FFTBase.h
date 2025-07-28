#pragma once
/*
Copyright  © 2025 Claus Vind-Andreasen

This program is free software; you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 - 1307 USA
*/

typedef uint32_t u32 ;
typedef double FFTType ;
const FFTType PI = 3.141592653589793238463;


class FFTBaseClass {
public:
	  virtual void SetLength(u32 p)	= 0;
	  virtual u32 Status() = 0;
	  virtual void ForwardFFT(FFTType *re, FFTType *im, unsigned int stride ) = 0;
	  virtual void ForwardFFT(FFTType *re, FFTType *im, unsigned int *indexMap ) = 0;
	  virtual void UnscaledInverseFFT(FFTType *re, FFTType *im, unsigned int stride) = 0;
	  virtual void UnscaledInverseFFT(FFTType *re, FFTType *im, unsigned int *indexMap) = 0;
	  virtual void InverseFFT(FFTType *re, FFTType *im, unsigned int stride) = 0;
	  virtual void InverseFFT(FFTType *re, FFTType *im, unsigned int *indexMap) = 0;
	  
};
