function T = tom_mark_cv_getTensor(P)


psize = size(P);
if (psize(1) ~= 3 || psize(2) ~= 4 || psize(3)<2 || psize(3)>4)
    error('The camera matrix must be a array of size 3x4xN (N=1~4)');
end;

switch (psize(3))
    case 2
        if (false)
            T = nan(3, 3);
            for (i=1:3)
                for (j=1:3)
                    T(j,i) = (-1)^(i+j)*det([P(setdiff([1:3], i), :, 1); P(setdiff([1:3], j), :, 2)]);
                end;
            end;
        else
            T = crossm(P(:,:,2)*null(P(:,:,1))) * P(:,:,2) * pinv(P(:,:,1)); 
        end;
        %T = T/norm(T)/sign(sum(T(:)));
    case 3
        T = nan(3, 3, 3);
        for (i=1:3)
            for (q=1:3)
                for (r=1:3)
                    T(q,r,i) = (-1)^(i+1) * det([P(setdiff([1:3], i), :, 1); P(q, :, 2); P(r, :, 3)]);
                end;
            end;
        end;
    case 4
        error('not yet implemented');
    otherwise
        error('Only supported for 2-, 3- or 4-view. The camera matrix must be a array of size 3x4xN (N=1~4)');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix of cross product cross(a,b) == crossm(a)*b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function am = crossm(a)

am = [[    0, -a(3),  a(2)]; ...
      [ a(3),     0, -a(1)]; ...
      [-a(2),  a(1),     0]];

