lmax := 300; //initialising range in which we want to apply 1st elimination
pmax := 3000; //initialising range for powers p up to which we want to search

load "IntegralFrobeniusMatrix.m";

function WeierstrassDiscriminant(a1,a2,a3,a4,a6) //compute the discriminant of the curve
    b2 := a1^2 + 4*a2;
    b4 := 2*a4 + a1*a3;
    b6 := a3^2 + 4*a6;
    b8 := a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;
    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6;
end function;

function F(A, B, l) //return Frey curve over F_l in case c even
	K := GF(l); //define finite field F_l
	//compute coefficients in the long Weierstrass equation for given values A, B over F_l
	a1 := K!1;
	a2 := (K!3*K!(B-A)-K!2)/K!8;
	a3 := K!0;
	a4 := K!3*K!(A+B)^2/K!64;
	a6 := (K!9*K!(B-A)*K!(A+B)^2)/K!512;
	//check if the curve obtained is non-singular
	if WeierstrassDiscriminant(a1, a2, a3, a4, a6) ne K!0 then 	//if yes, return curve model
		return true, EllipticCurve([GF(l) | a1, a2, a3, a4, a6]);
	else								//if no, return 0
		return false, 0;
	end if;
end function;

function FF(A, B, l) //return Frey curve over F_l in case c odd
	K := GF(l); //define finite field F_l
	//compute coefficients in the long Weierstrass equation for given values A, B over F_l
	a1 := K!0;
	a2 := K!0;
	a3 := K!0;
	a4 := K!3*K!A*K!B;
	a6 := K!B^3-K!A^3;
	//check if the curve obtained is non-singular
	if WeierstrassDiscriminant(a1, a2, a3, a4, a6) ne K!0 then	//if yes, return curve model
		return true, EllipticCurve([GF(l) | a1, a2, a3, a4, a6]);
	else								//if no, return 0
		return false, 0;
	end if;
	end function;

//Function computing B_l(F) 
function Bl(lvl, l, f)
	b := ((l+1)^2 - Coefficient (f, l)^2); //Initialise the product with ((l+1)^2 - a_l(E_f)^2) 
	traces := {}; //Initialise the set S_l
	//we now compute traces of Frobenius at l for all possible a, b in the finite field F_l
	for A in [0..l-1] do
		for B in [0..l-1] do
			if lvl eq 126 then 			//case c even
				flag, F1:= F(A, B, l);
			else					//case c odd
				flag, F1:= FF(A, B, l);
			end if;
			if flag then				//for every non-singular curve obtained
				i := TraceOfFrobenius(F1);	//compute its trace of Frobenius at l
        			Include(~traces, i);		//add it to the list of traces
        			//function Include() omits duplicates, so only unique traces are kept
    			end if;
    		end for;
	end for;
	for t in traces do
		b*:= (t - Coefficient(f,l)); //compute the desired product as in Proposition 2.1.1
	end for;
	return b;
end function;

function reps(EC) //for each of the curves that survived the initial sieve, find representatives fulfilling condition of proposition 2.1.2
	ECrep := [* *];
	fail := 0; //flag indicating if a class failed to deliver the desired representative
	for pair in EC do
		found := #ECrep; //the number of representatives found so far
		r :=CremonaReference(pair[1]);
		while r[#r] in "0123456789" do //retrieve the isogeny class of the curve
			Prune(~r); 
		end while;
		temp := EllipticCurves(CremonaDatabase(), r); //list all representatives in the isogeny class
		for index in [1..#temp] do
			E := MinimalModel(temp[index]);	//for each representative find its minimal model
			Delta := Discriminant(E);
			if Valuation(Delta, 7) eq 2 then //check condition 2 of the Proposition 2.1.2
    				Frob := TracesOfFrobenius(E, Floor(pmax^2/4)); //tabelarise the traces of Frobenius needed for the local criterion //will possibly be removed
				Append(~ECrep, <E, pair[2], Frob>);
				break;//move on to the next class if a representative is found
			end if;
		end for;
		if(#ECrep eq found) then //if for one class no representative fulfilling condition 2 has been found, Theorem 2.1.3 can offer no further insight
			fail := 1;
			break;
		end if;
	end for;
	if(fail eq 1) then 
		ECrep := [];
	end if;
	return ECrep;  
end function;

// This function checks that p divides the order of Frobenius.
function CheckpdivFrob(E,l,p)
    FF:=GF(l);
    E:=ChangeRing(E,FF);
    G:=GL(2,p);
    F:=G!IntegralFrobenius(E);
    if Order(F) mod p eq 0 then
        return true;
    else
        return false;
    end if;
end function;

//This function checks that for all primes q dividing \Delta_l (q/p) is not -1.
function Checkcond4(Dl2, p)
    primeFactors := PrimeDivisors(Dl2);
    for q in primeFactors do
        if LegendreSymbol(q, p) eq -1 then
            return false;
        end if;
    end for;
    return true;
end function;

function local_elim(EC)
	ElimPrimes := [* *];
	for p in PrimesInInterval(17,pmax) do
		if p mod 4 ne 3 then //Assumption (1) of Theorem 2.1.3
			continue; //if not fulfilled, go to the next prime
		end if;
		p_good := 0; //counter for how many curves p fulfilled all conditions 1-6 of Theorem 2.1.3
		how_many := 0; //counter how many curves are being checked against the local criterion
		for pair in EC do
			E := pair[1];
			//check if the curve shouldn't have been eliminated due to Proposition 2.1.1 already
    			if pair[2] mod p ne 0 then 
    				print "The isogeny class of the curve ", CremonaReference(pair[1]), " has been eliminated during the first elimination with gcd of its Bl coefficients equal ", pair[2];
    				continue;
    			end if;
    			Frob := pair[3];//tabelarised values of the traces of Frobenius
    			how_many+:=1;
    			i := 1;//index for Frob
    			for l in PrimesInInterval(2,Floor(p^2/4)) do
		                a := Frob[i];
                		i+:=1;
                		if l eq 7 then 
                			continue; 
                		end if;
	                	if l mod p ne 1 and Abs(a) lt p - 2*Sqrt(l) then //conditions (6) and (5)
                    			Dl:=a^2-4*l;
                    			if Dl mod p eq 0 then
                    				is_sq, Dl2 := IsSquare(ExactQuotient((-Dl), p));
                 				if is_sq then 
	                    				if l ne p and Checkcond4(Dl2, p) and CheckpdivFrob(E,l,p) then //Conditions (2), (3), (4)
        	                			p_good:=p_good+1;
        	                			break; 
        	            				end if;
                    				end if;
                    			end if;
                    		end if;
                	end for;
    			
		end for;
		if p_good eq how_many then
		        Append(~ElimPrimes,p);
        	end if;
	end for;
	return ElimPrimes;
end function;

//initialising arrays for newforms and corresponding elliptic curves
NFs := Newforms(CuspForms(126)) cat Newforms(CuspForms(252)) cat Newforms(CuspForms(504));
EC := [* *];

//Main part of the program
for i in [1..#NFs] do
	f := NFs[i][1]; //retrieve newform
	lvl := Level(f); //and its level 
	gcd:=0; //initialise variable for storing gcd
	for l in PrimesInInterval(1,lmax) do
		if lvl mod l ne 0 then //Assumption of Prop 2.1.1: l \nmid N_E
			val := Integers()!Bl(lvl, l, f); //Compute B_l(F)
			gcd := GCD(gcd, val);
		end if;
	end for;
	E := EllipticCurve(f); //retrieve an elliptic curve associated with f 
	EC := Append(EC, <E, gcd>); //add to the list the elliptic curve together with gcd of its B_l(E)
end for;
survivors := [* *]; //initialising array for storing curves surviving initial sieve
for pair in EC do
	//Since the smallest p for which the criterion l \nmid p becomes relevant is 17, if all common divisors of Bls are smaller than 17, the non-zero Bls will remain non-zero for any p considered later; thus, the curve can be eliminated immediately. 
	if pair[2] eq 0 then
		Append(~survivors, pair);
	elif pair[2] eq 1 or Max(Set(PrimeDivisors(pair[2]))) lt 17 then
		print "The isogeny class of the curve ", CremonaReference(pair[1]), " has been eliminated during the first elimination with gcd of its Bl coefficients equal ", pair[2];
	else
		Append(~survivors, pair);	
	end if;
end for;
EC := reps(survivors); //update the list of the remaining curves
if IsEmpty(EC) then 
	print "The local criterion could not eliminate the remaining curves.";
else
	print "The remaining curves have been eliminated for the following primes: ", local_elim(EC);
end if;

