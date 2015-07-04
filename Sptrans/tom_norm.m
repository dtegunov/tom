function c=tom_norm(a,scf,mask)
%TOM_NORM  Normalizes the image values between 0 and Scaling Factor (scf)
%
%   c=tom_norm(a,scf,mask)
%
%   C=TOM_NORM(A,SCF,mask) This function normalizes the values of the input image A between 0 and
%   Scaling Factor, SCF. The output is the image A with normalized values.
%
%   If SCF is 'phase' the output is normalized as: (A-mean(A))/mean(A)
%   
%   If SCF is '3std' or '2std' or '1std', the data is limited to 3,
%   respectively 2 or 1
%   standard deviations and the mean is subtracted.
%   
%   mean0+1std
%    
%   If a mask is given, the image values for normalization are only taken
%   from inside this mask
%
%
%  INPUT
%   a                   input image
%   scf                 scaling factor
%   mask                apply in mask
%  
%  OUTPUT
%   c                   output image
%
%EXAMPLE
%       c=tom_norm(a,10)
%
%           10  1  19  8   5
%           12  9  17  10  10
%       a = 3   3  2   0   8
%           15  8  0   7   8
%           0   8  6   15  8
%
%           5.263   0.526   10.00   4.210   2.631
%           6.315   4.736   8.947   5.263   5.263
%       c = 1.578   1.578   1.056   0       4.210   
%           7.894   4.210   0       3.684   4.210
%           0       4.219   3.157   7.894   4.210
%
%       c=tom_norm(a,'phase')
%       c=tom_norm(a,'2std')
%       c=tom_norm(a,'3std')
%       c=tom_norm(a,'mean0+1std')
%       
%       space
%
% Example for oscar usage:
%       for i=1:53;
%           in=tom_emread('particle_' num2str(i) '.em'external link); 
%           in.Value=tom_norm(in.Value,'3std');
%           in.Value=in.Value+100; 
%           in.Value=tom_norm(in.Value,'phase',mask); 
%           tom_emwrite('normed/particle_' num2str(i) '.em'external link,in);
%       end;
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_PEAK, TOM_LIMIT, TOM_FILTER
%
%   created by AL 08/17/02
%   updated by SN 10/13/05
%   updated by AK 16/02/06
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

if nargin == 1
    scf = 'mean0+1std';
end;

if ischar(scf) && strcmp(scf,'mean0+1std')
    a = (a.*2) + 100;
end 

%norm only inside of mask
if nargin == 3
    if size(a) ~= size(mask)
        error('image and mask must have the same dimensions');
    end
    ind = find(mask>0);
    %a=a.*mask;
    [mea, ma, mi, st] = tom_dev(a(ind),'noinfo');
else
    [mea, ma, mi, st] = tom_dev(a,'noinfo');
end

%phase norm
if ischar(scf) && strcmp(scf,'phase')
    c=(a-mea)./(mea+0.000000000001);
    
    %c=(a-mea)./(mea);
    
%histogram equalisation    
elseif ischar(scf) && strcmp(scf,'3std')
    c = tom_limit(a,mea-3*st,mea+3*st)-mea;

elseif ischar(scf) && strcmp(scf,'2std')
    c = tom_limit(a,mea-2*st,mea+2*st)-mea;
elseif ischar(scf) && strcmp(scf,'1std')
    c = tom_limit(a,mea-st,mea+st)-mea;
elseif ischar(scf) && strcmp(scf,'mean0+1std')
    if (st ~= 0)
        c = (a-mea)./st;
    else
        c = (a-mea);
        warning('std = 0 !!!');
    end;
    
elseif ischar(scf) && strcmp(scf,'oscar')
    [mea2, ma2, mi2, st2] = tom_dev(a,'noinfo');
    %3std
    c = tom_limit(a,mea2-3*st2,mea2+3*st2)-mea2;
    %+100
    c = c + 100;
    
    %phase
    if nargin == 3
        if size(a) ~= size(mask)
            error('image and mask must have the same dimensions');
        end
        ind = find(mask);
        [mea, ma, mi, st] = tom_dev(c(ind),'noinfo');
    else
        [mea, ma, mi, st] = tom_dev(c,'noinfo');
    end

    c=(c-mea)./(mea+0.000000000001);
elseif isnumeric(scf)
    if st~=0 && mea~=scf
        a=a-mi;
        %c=scf*a./(ma+0.00000000000001);
        c=scf.*a./max(max(max(a)));    
    elseif st==0 && mea~=scf
        c=a./scf;
    else
        c=a;
    end
    
else
    error('Only numbers and keyword ''phase'', ''2std'' , ''3std'' or ''mean0+1std'' allowed.');
end;
