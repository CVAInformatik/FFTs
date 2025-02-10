

(* #use "fft.mli";;*)
#use "fft.ml";;


let rec addSineAux _N freq phase amplitude i v =
        let arg = phase +.(( freq *. 2.0 *. (float_of_int i) )  /. (float_of_int _N ))  *. Float.pi
        in  v +.  (amplitude *. sin(arg))
          
let rec  addSine (arr:float array) freq phase amplitude  =      
      let f =  addSineAux  (Array.length arr) freq phase amplitude 
			in Array.mapi_inplace ( fun i a -> (f i  a) ) arr 
        
let _N = 16384 ;;

let (twre, twim) = FFT2.mkTwiddles _N  (_N/2)  ;;

let testSignalRe = Array.make _N 0.0;;

let testSignalIm = Array.make _N 0.0;;

addSine testSignalRe 2.0 (Float.pi /. 2.) 2. ;;

testSignalRe;;

addSine testSignalRe 5.0 (Float.pi /. 2.) 1. ;;

addSine testSignalRe 25.0 0. 1. ;;

FFT2.fft2 _N testSignalRe testSignalIm twre twim ;;

testSignalRe;;

testSignalIm;;

FFT2.inverse_fft2 _N testSignalRe testSignalIm twre twim ;;

testSignalRe;;


