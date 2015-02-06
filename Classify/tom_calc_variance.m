function var=tom_calc_variance(input_data,classes,class)
%tom_calc_variance calculates the variance of a feature vector
%
%    var=tom_calc_variance(input_data,classes,class)
%
%   plots positions as a 2d grid
%   useful for ploting a som
%
%PARAMETERS
%
%  INPUT
%   input_data    data matrix
%   classes       vector containing the class number (optional)
%   class         class nr variance should be calculated (optional)
%  
%  OUTPUT
%   var           variance vector
%
%EXAMPLE
%    var=tom_calc_variance(input_data,classes,1)
%
%REFERENCES
%
%SEE ALSO
%   -
%
%   created by FB 11/17/06
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



if (nargin ==1);
    classes=ones(size(input_data,1),1);
    class=1;
end;

if ((sum(classes==class))<2)
    var=zeros(size(input_data,2),1);
    return;
end;

num=length(find(classes));
tmp=zeros(num,size(input_data,2));

zz=1;
for i=1:size(input_data,1)
    if (classes(i)==class)
        tmp(zz,:)=input_data(i,:);
        zz=zz+1;
    end;
end;

tmp=tmp(1:zz,:);

clear('input_data');

var=std(tmp).^2;



