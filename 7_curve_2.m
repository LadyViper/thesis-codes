NFs := Newforms(CuspForms(504));
NFs;
f:=NFs[8,1];
f;
E := EllipticCurve(f);
E;

function WeierstrassDiscriminant(a1,a2,a3,a4,a6)
    b2 := a1^2 + 4*a2;
    b4 := 2*a4 + a1*a3;
    b6 := a3^2 + 4*a6;
    b8 := a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2;

    return -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6;
end function;

function F(A, B, l)
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

function Bl(j)
b := ((j+1)^2 - Coefficient (f, j)^2);
K := GF(j);
traces := {};
for A in [0..j-1] do
	for B in [0..j-1] do
		flag, F1:= F(A, B, j);
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
for k in [1..1000] do
	if 126 mod k ne 0 then
	if IsPrime(k) then
	val := Bl(k);
	//k; val;
	if val ne 0 then print "Non-zero coefficient found for l = ", k; end if;
	g:= GCD(g, val);
	end if;
	end if;
end for;
if g ne 0 then Factorisation(g);
else print "Curve cannot be eliminated";
end if;
