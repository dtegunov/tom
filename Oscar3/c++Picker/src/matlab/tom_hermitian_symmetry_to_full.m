function v = tom_hermitian_symmetry_to_full(v2, sizez)


if (floor(sizez/2)+1 ~= size(v2,3))
    error('wrong sizez');
end;

v = nan(size(v2,1),size(v2,2), sizez);

v(:,:,1:size(v2,3)) = v2;


sizey = size(v2,2);
sizex = size(v2,1);

for (z= ((floor(sizez/2)+1):(sizez-1)))
    for (y= (0:(sizey-1)))
        for (x= (0:(sizex-1)))
            v(x+1,y+1,z+1) = conj(v2(mod(sizex-x,sizex)+1, mod(sizey-y,sizey)+1, mod(sizez-z,sizez)+1));
        end;
    end;
end;
    

