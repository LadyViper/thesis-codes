//NFs := Newforms(CuspForms(126));
//NFs;
//f:=NFs[1,1];
//f;
EC:=["126a1", "126b1", "252a1","252b1", "504a1", "504b1", "504c1","504d1","504e1","504f1", "504g1", "504h1" ]; 
//E := EllipticCurve("90c2");
//E;

function WeierstrassDiscriminant(a1,a2,a3,a4,a6)
    b2 := a1^2 + 4*a2;
    b4 := 2*a4 + a1*a3;
    b6 := a3^2 + 4*a6;
    b8 := a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;

    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6;
end function;

function F(A, B, l)
K := GF(l);
a1 := K!1;
a2 := (K!3*K!(B-A)+K!2)/K!8;
a3 := K!0;
a4 := K!3*K!(A+B)^2/K!64;
a6 := (K!9*K!(B-A)*K!(A+B)^2)/K!512;
if WeierstrassDiscriminant(a1, a2, a3, a4, a6) ne K!0 then
	return true, EllipticCurve([GF(l) | a1, a2, a3, a4, a6]);
else
	return false, 0;
end if;
end function;

function FF(A, B, l)
K := GF(l);
a1 := K!0;
a2 := K!0;
a3 := K!0;
a4 := K!3*K!A*K!B;
a6 := K!B^3-K!A^3;
if WeierstrassDiscriminant(a1, a2, a3, a4, a6) ne K!0 then
	return true, EllipticCurve([GF(l) | a1, a2, a3, a4, a6]);
else
	return false, 0;
end if;
end function;

function Bl(N, j, f)
b := ((j+1)^2 - Coefficient (f, j)^2);
K := GF(j);
traces := {};
for A in [0..j-1] do
	for B in [0..j-1] do
		if N eq 126 then
			flag, F1:= F(A, B, j);
			else
			flag, F1:= FF(A, B, j);
		end if;
		if flag then
		i := TraceOfFrobenius(F1);
        	Include(~traces, i);
    		end if;
    	end for;
end for;
for t in traces do
	b*:= (t - Coefficient(f,j));
end for;
return b;
end function;

g:=0;



for label in EC do
t0 := Cputime();
E:=EllipticCurve(label);
N:=Conductor(E);
f := ModularForm(E);
for k in [1..1000] do
	if N mod k ne 0 then
		if IsPrime(k) then
		val := Bl(N, k, f);
		//k; val;
		g:= GCD(g, val);
		end if;
	end if;
end for;
if g ne 0 then
	print "B is", Factorisation(g);
	print "ELIMINATE", label;
else print "For primes smaller or equal than 500 the curve", label, "cannot be eliminated";	
end if;
print "Total time:", Cputime(t0);
end for;


