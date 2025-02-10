
module FFT2 = struct
(*
		The ideas and code in this program are based on
		
		SIAM journal on scientific and statistical computing.
		Society for Industrial and Applied Mathematics.
		SELF-SORTING IN-PLACE FAST FOURIER-TRANSFORMS		
		1991-07-01  Clive Temperton.
		
 this code is basically Tempertons Fortan -> C++ -> non-functional ocaml (that is: using the IMPERATIVE  features in the Ocaml language) *)

(*  standard classic butterfly *)
let fft2_2  n c1 c2 tw (rea:float array) (ima:float array) (twre:float array) (twim: float array) = 
						try
						   if( tw = 0)  
						   then  
				 	  		while true do
					 	  		let    zre =  rea.(c1) -. rea.(c2) 
					 	  		in let zim =  ima.(c1) -. ima.(c2)
		  		 	  		in 
		  		 	  		rea.(c1) <-  rea.(c1) +. rea.(c2);
				 	  		 	ima.(c1) <-  ima.(c1) +. ima.(c2);
				 	  		 	rea.(c2) <- zre;
				 	  		 	ima.(c2) <- zim;
				 	  		 	raise Exit	 	
				 	  		done
							else if(  tw = (n/4) ) then 
							  while true do 
				 	  		let zre = -1. *. ( ima.(c1) -. ima.(c2))
				 	  		in let zim =  rea.(c1) -. rea.(c2)
	  		 	  		in rea.(c1) <- rea.(c1) +. rea.(c2);
				 	  		   ima.(c1) <- ima.(c1) +. ima.(c2);
				 	  		   rea.(c2) <- zre;
				 	  		   ima.(c2) <- zim;
				 	  	  raise Exit	 	
				 	  		done
							else 
							  while true do 
		 	  				let  wimag =   twim.(tw)
		 	  				in let wreal = twre.(tw)
				 	  		in let zre = (wreal *.( rea.(c1) -. rea.(c2))) -. (wimag *. ( ima.(c1) -. ima.(c2)))
				 	  		in let zim = (wimag *.( rea.(c1) -. rea.(c2))) +. (wreal *. ( ima.(c1) -. ima.(c2)))
	  		 	  		in rea.(c1) <- rea.(c1) +. rea.(c2);
				 	  		ima.(c1) <- ima.(c1) +. ima.(c2);
				 	  		rea.(c2) <-  zre ;
				 	  		ima.(c2) <-  zim ;
								raise Exit	 	
				 	  		done		 		
					with Exit -> ()


(*  tempertons coupled  butterfly (mixing to standard butterflies, used in the second half to get an In-Order-In-Place  FFT *)
(* this could probably be fixed to use the previous butter in some way, but for the moment I'll leave it alone *)
let fft2_4  n c1 c2 c3 c4 tw (rea:float array) (ima:float array) (twre:float array) (twim: float array) = 
						try
						   if( tw = 0) then 
				 	  		while true do
		 	  				let zre =  rea.(c1) -. rea.(c2) in 
		 	  				let zim =  ima.(c1) -. ima.(c2) in
		 	  				let tre =  rea.(c3) -. rea.(c4) in
		 	  				let tim =  ima.(c3) -. ima.(c4) in
		 	  				rea.(c1) <- rea.(c1) +. rea.(c2);
		 	  				ima.(c1) <- ima.(c1) +. ima.(c2);
		 	  				rea.(c2) <- rea.(c3) +. rea.(c4);
		 	  				ima.(c2) <- ima.(c3) +. ima.(c4);
		 	  				rea.(c4) <- tre;
		 	  				ima.(c4) <- tim;
          			rea.(c3) <- zre;
          			ima.(c3) <- zim;
			 	  		 	raise Exit	 	
				 	  		done
							else if(  tw = (n/4) ) then 
							  while true do 
		 	  				let zre = -1. *. ( ima.(c1) -. ima.(c2)) in 
		 	  				let zim =  rea.(c1) -. rea.(c2) in 
		 	  				let tre =  -1. *. ( ima.(c3) -. ima.(c4)) in
		 	  				let tim =   rea.(c3) -. rea.(c4) in
		 	  				rea.(c1) <- rea.(c1) +. rea.(c2);
		 	  				ima.(c1) <- ima.(c1) +. ima.(c2);
		 	  				rea.(c2) <- rea.(c3) +. rea.(c4);
		 	  				ima.(c2) <- ima.(c3) +. ima.(c4);
		 	  				rea.(c4) <- tre;
		 	  				ima.(c4) <- tim;
          			rea.(c3) <- zre ;
          			ima.(c3) <- zim ; 
				 	  	  raise Exit	 	
				 	  		done
	 						else 
							  while true do 
					 	  	let wimag = twim.(tw)  in
		 				  	let wreal = twre.(tw)  in
		 	  				let tre =  (wreal *. ( rea.(c3) -. rea.(c4))) -. (wimag *. ( ima.(c3) -. ima.(c4))) in
		 	  				let tim =  (wimag *. ( rea.(c3) -. rea.(c4))) +. (wreal *. ( ima.(c3) -. ima.(c4))) in
		 	  				let zre =  (wreal *. ( rea.(c1) -. rea.(c2))) -.  (wimag *. ( ima.(c1) -. ima.(c2))) in
		 	  				let zim =  (wimag *. ( rea.(c1) -. rea.(c2))) +.  (wreal *. ( ima.(c1) -. ima.(c2))) in
		 	  				rea.(c1) <- rea.(c1) +. rea.(c2);
		 	  				ima.(c1) <- ima.(c1) +. ima.(c2);
		 	  				rea.(c2) <- rea.(c3) +. rea.(c4);
		 	  				ima.(c2) <- ima.(c3) +. ima.(c4);
		 	  				rea.(c4) <- tre;
		 	  				ima.(c4) <- tim;
          			rea.(c3) <- zre;
          			ima.(c3) <- zim;
          			raise Exit	 	
				 	  		done;		 		
					with Exit -> ()
        


let fft2 _N rea ima twre twim =
    let _NH = ref  (_N/2) in
    let _M = ref 0 in
    let  mtemp = ref 1 in
		while !mtemp < _N do 
				mtemp := !mtemp+ !mtemp; 
				_M := !_M + 1 ;
		done;
		let _LA = ref 1 in
		for _L = 1 to !_M do
		   (* Printf.printf "_L  %i \n" _L ;*)
		   if _L = 1 then _LA  := 1 else _LA := !_LA + !_LA ;
		   let _IA =  ref 0 in
		   let _IB = ref ((!_NH)/(!_LA)) in
		   let _KK = ref 0 in
		   if _L < ((!_M+3)/2) 
		   then
		      for _K = 0 to (!_IB-1) do
		        let _I = ref _K in
		        while !_I < _N do
		        	fft2_2 _N (!_IA+(!_I)) (!_IB+(!_I)) (!_KK) rea ima twre twim ;	
		        	_I := !_I + ((_N)/(!_LA));
		        done;
		        _KK := !_KK+ !_LA;
		      done
		   else
		      let _IC = ref  !_LA in
		      let _ID = ref (!_IB + !_LA) in
		      for _K = 0 to (!_IB-1) do 
		        let _J = ref _K in 
		        while !_J  < !_LA do 
			        let _I = ref !_J in 
		        	while !_I < _N do 
		        	   fft2_4 _N (!_IA+ (!_I)) (!_IB+(!_I)) (!_IC+(!_I))  (!_ID+(!_I)) (!_KK) rea ima twre twim ;	
		        		_I := !_I + (2 * (!_LA));
		        	done;
		        _J := !_J+ (_N/(!_LA));
		        done;	 
		      	_KK := !_KK+ !_LA;
	      	done 
		done

let unscaled_inverse_fft2 _N rea ima twre twim =  fft2 _N ima rea twre twim

let inverse_fft2 _N rea ima twre twim =  
							for i = 0 to 0 do
                unscaled_inverse_fft2 _N rea ima twre twim ;
								Array.map_inplace (fun x-> x /. (float_of_int _N) ) rea ;
								Array.map_inplace (fun x-> x /. (float_of_int _N) ) ima ;
							done
							

        
let mkTwiddles _N n  =
         let twre = Array.make n 0.0  in
         let twim = Array.make n 0.0  in         
         for i = 0 to n - 1 do 
            let arg = ((float_of_int ( i * 2) ) /. (float_of_int _N ))  *. Float.pi in
         		twre.(i)<- cos(arg);
         		twim.(i)<- sin(arg);
         done ;
         (twre, twim) 
        			
end;;