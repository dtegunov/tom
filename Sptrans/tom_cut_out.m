function out=tom_cut_out(in,pos,size_c,fill_flag)
%TOM_CUT_OUT cuts one array out of another one
%
%   out=tom_cut_out(in,pos,size_c,fill_flag)
%
%The array out will be cut out from the array in so that
%its upper left corner will become pixel (POS(1) POS(2)) of array in
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   pos                 ...
%   size_c              ...
%   fill_flag           ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%      out=tom_cut_out(in,[3 3],[5 5],'no-fill');                           
%
%       in=
%      -4     0    -1   -13   -14     7   -10         
%     -16     4     2     8     6    12    15        
%       2     2     1     1     1     1     1         
%       3    -1     1     1     1     1     1   
%     -11     8     1     1     1     1     1      
%      12    22     1     1     1     1     1        
%
%         out=
%         1  1  1  1  1
%         1  1  1  1  1
%         1  1  1  1  1 
%         1  1  1  1  1
%         1  1  1  1  1
%
% in=rand(128,128);
% out=tom_cut_out(in,'center',[64 64 64]);
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_PEAK, TOM_PASTE
%
%   created by FB hurz! 07/28/05
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


if nargin < 4
    fill_flag='no-fill';
end;


%kick out values < 1

if (~isnumeric(pos))
    pos=floor((size(in)-size_c)./2)+1;
end;

pos=((pos>0).*pos)+(pos<1);

num_of_dim=size(size(in),2);
if (size(in,1)==1)
    in_size=size(in,2);
    num_of_dim=1;
else
    in_size=size(in);
    num_of_dim=size(size(in),2);
end;
 
if ( (pos(1)+size_c(1))<=in_size(1))
    bound(1)=pos(1)+size_c(1)-1;
else
    bound(1)=in_size(1);
end;

if (num_of_dim > 1)
    if ((pos(2)+size_c(2))<=in_size(2) )
        bound(2)=pos(2)+size_c(2)-1;
    else
        bound(2)=in_size(2);
    end;
else
    pos(2)=1;
    bound(2)=1;
end;

if (num_of_dim > 2)
    if ((pos(3)+size_c(3))<=in_size(3) )
        bound(3)=pos(3)+size_c(3)-1;
    else
        bound(3)=in_size(3);
    end;
else
    pos(3)=1;
    bound(3)=1;
end;

pos=round(pos);
bound=round(bound);

% cut it
if  (num_of_dim>1)
    out=in(pos(1):bound(1),pos(2):bound(2),pos(3):bound(3));
else
    out=in(pos(1):bound(1));
end;


if (strcmp(fill_flag,'fill')== 1)

end;








