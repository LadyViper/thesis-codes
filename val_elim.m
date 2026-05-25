

EC := [ <126,"a">, <126,"b">, <252,"a">, <252,"b">, <504,"a">, <504,"b">, <504,"c">, <504,"d">, <504,"e">, <504,"f">, <504,"g"> ];
D := CremonaDatabase();
Ef :=[];

for pair in EC do
	N := pair[1];
	label := pair[2];
	for index in [1..NumberOfCurves(D, N, label)] do
		E := MinimalModel(EllipticCurve(D, N, label, index));
		Delta := Discriminant(E);
		if Valuation(Delta, 7) eq 2 then
			CremonaReference(E);
			Append(~Ef, E);
		end if;
	end for;
end for;

//This part has been taken from Diana and Nuno's code, modified for optimisation 

//This code verifies the claims in Section 10.

load "/home/viper/Desktop/Studium/M VI Semester/Master Thesis/IntegralFrobeniusMatrix.m";

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

ps:=[];
t0 := Cputime();
for p in PrimesInInterval(20,3000) do
    if p mod 4 eq 3 then //Assumption (1) of Theorem 5.3.
        okp:=0;
    for label in Ef do
    E:=EllipticCurve(label);
                for l in PrimesInInterval(5,Floor(p^2/16)) do //Assumption (5) of Theorem 5.3.
                    a:=TraceOfFrobenius(E,l);
                    Dl:=a^2-4*l;
                    if Dl mod p eq 0 then
                    is_sq, Dl2 := IsSquare(ExactQuotient((-Dl), p));
                    if is_sq then 
	                    if l ne p and Checkcond4(Dl2, p) and CheckpdivFrob(E,l,p) then //Assumptions (2), (3), (4) of Theorem 5.3.
        	                okp:=okp+1;
        	                break l;
        	            end if;
                    end if;
                    end if;
                end for;
    end for;
    	print "p works for okp curves", p, okp;    
        if okp eq #Ef then
        print "p is good for all curves",p;
        Append(~ps,p);
        end if;
    end if;
end for;
print "Total time:", Cputime(t0);


//Curves obtained:
//126a2
//126a6
//126b2
//126b6
//252a1
//252b1
//504a1
//504b2
//504c2
//504d1
//504e2
//504f2
//504g2

//594.900

