function laguerre = av3_laguerre(r,s,x)
% AV3_LAGUERRE creates Laguerre Polynomials
%
% LAGUERRE = AV3_LAGUERRE(R, S, X) 
%
%   laguerre = av3_laguerre(r,s,x)
%
% SEE ALSO
%   H-atom

%design polynomial
k=0:r-s;
p=(-1).^(k+s) .* ((av3_factorial(r)).^2) ./ (av3_factorial(k).*av3_factorial(k+s).*av3_factorial(r-k-s));
%vgl. Schwabl, "Quantenmechanik" p. 128 
laguerre = polyval(p,x);



function fac=av3_factorial(n)
%   fac=av3_factorial(n)
%   wie factorial, aber fuer arrays
for ind=1:size(n,2)
    fac(ind) = prod(1:n(ind));
end;