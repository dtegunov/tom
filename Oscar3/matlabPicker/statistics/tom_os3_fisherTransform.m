function res=tom_os3_fisherTransform(ccc,q,n)
% Calculates Fisher Transform according to Pencek slides Image Assesment -2005 

    z = fisherTrans(ccc);
    q = 1.96;
    left = z - q/sqrt(n-3);
    right = z + q/sqrt(n-3);  


    %return interval borders scaled to r value
    res = [ccc (exp(2*left)-1)/(exp(2*left)+1) (exp(2*right)-1)/(exp(2*right)+1)];

%%
%formula for fisher transformation
function res = fisherTrans(r)

res = 0;
%avoid 0s in transformation
if(r == 1)
    return;
end;

res = 0.5*log((1+r)/(1-r));