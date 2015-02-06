function varargout=tom_cut(varargin)
%TOM_CUT Filtering in Fourier Space
%
%   varargout=tom_cut(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%   [J h Hd]=tom_cut(I,type,order,cutoff,PROPERTIES) This function applies
%   a high- or a lowpass filter to the image. The user can specify
%   parameters such as the order of the filter or the cut off frequency. J
%   contains the filtered image while h is the frequency response and Hd is
%   the bandpass respone. Order of the filter means its dimensions. Usually
%   choose the dimensions of the image that is going to be filtered.
%
%
%
%
%   
%   PROPERTIES is a cell array containing one string. The only acceptable
%   value for this string is 'gaussian'. If the user chooses to apply a
%   gaussina window then he should also justify the deviation of the
%   gaussian distribution
%
%EXAMPLE
%   A=IMREAD('RICE.TIF');
%
%   J=TOM_CUT(A,'lowpass',16,0.4); User applies a lowpass filter to the
%   image A with a 0.4 cut off frequency.
%   
%   [J H]=TOM_CUT(A,'highpass',32,0.3); User applies a higpass filter to
%   the image A with a 0.3 cut off frequency
%
%   [J H HD]=TOM_CUT(A,'highpass',16,0.2,'GAUSSIAN',0.5); User applies a
%   highpass filter using a gaussian window (deviation 0.5) with a 0.2 cut
%   off frequency
%
%REFERENCES
%
%SEE ALSO
%   TOM_PEAK, TOM_FILTER, TOM_COMP, TOM_FIRDEMO
%
%   created by AL 09/01/02
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

if nargin<4
    error('Not Enough Inout Arguments');
    %break;
end

IMG=varargin{1};
art=varargin{2};
order=varargin{3};
frcut=varargin{4};
if nargin == 4    
    [h Hd]=DesignFilter(art,order,frcut);
    J=filter2(h,IMG);
elseif nargin == 6
    smooth=varargin{5};
    dev=varargin{6};
    [h Hd]=DesignFilter(art,order,frcut);
    h=AddWindow(Hd,smooth,dev);
    J=filter2(h,IMG);
else
    error('Not proper number of Arguments');
    %break;
end
if nargout == 1
    varargout{1}=J;
elseif nargout == 2
    varargout{1}=J;
    varargout{2}=h;
elseif nargout == 3
    varargout{1}=J;
    varargout{2}=h;
    varargout{3}=Hd;
end
    
    
    
        
    
%%%
%%%  Sub-Function - DesignFilter
%%%

function [h,Hd]=DesignFilter(tp,or,ctf)
    
[f1 f2]=freqspace(or,'meshgrid');
Hd=zeros(or);
Hd(sqrt(f1.^2+f2.^2)<ctf)=1;
if isequal(tp,'lowpass')    
elseif isequal(tp,'highpass')
    Hd=1-Hd;
else
    error('Unknown Type Of Filter');
    %break;
end
h=fsamp2(Hd);
%figure;freqz2(h);
return;


%%%
%%%  Sub-Function - AddWindow
%%%  

function h=AddWindow(Hd,tp,dev);

[s1 s2]=size(Hd);
win=fspecial(tp,s1,dev);
win=win./max(win(:));
h=fwind2(Hd,win);
%figure;freqz2(h);
return;



