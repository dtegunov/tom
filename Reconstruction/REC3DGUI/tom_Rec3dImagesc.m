function tom_Rec3dImagesc(varargin)
%TOM_REC3DIMAGESC is a module of TOM_REC3DGUI.
%
%   tom_Rec3dImagesc(varargin)
%
%   It displays a matrix in EM style.
%
%PARAMETERS
%
%  INPUT
%   tom_imagesc(input,'fixed') displays the real size of the matrix input. 
%   1 pixel of the image correspond to 1 pixel on the display.
%
%   tom_imagesc(input,'fixed noinfo') displays the real size of the matrix
%   input. No information image are displayed.
%
%   tom_imagesc(input,'noinfo') displays matrix input as an image. No
%   information image are displayed.
%
%   tom_imagesc(input,'header') displays matrix input as an image. A
%   created header is added to the input structure.
%
%   tom_imagesc(input,'fixed','header') displays the real size of the
%   matrix input. A created header is added to ths input structure.
%   'fixed' can be replaced by 'fixed noinfo' or 'noinfo'.
%
%   tom_imagesc(input,'header',expr_header) displays matrix input as an
%   image. The header expression expr_header is added to the input structure.
%
%   tom_imagesc(input,'range',[min max]) displays matrix input as an image
%   in the range [min max].
%
%   tom_imagesc(input,'fixed','header',expr_header) displays the real size of 
%   the image. The header expression expr_header is added to the input structure.
%   'fixed' can be replaced by 'fixed noinfo' or 'noinfo'.
%
%   tom_imagesc(input,'fixed','range',[min max]) displays the real size of the
%   matrix input in the range [min max]. 'fixed' can be replaced by 'fixed noinfo'
%   or 'noinfo'. 
%
%   tom_imagesc(input,'header','range',[min max]) displays the matrix input
%   in the range [min max]. A created header is added to the input structure.
%
%   tom_imagesc(input,'fixed','header','range',[min max]) displays the real
%   size of the matrix input in the range [min max]. A created header is
%   added to the input structure. 'fixed' can be replaced by 'fixed noinfo' or
%   'noinfo'.
%   
%   tom_imagesc(input,'header',expr_header,'range',[min max]) displays matrix 
%   input as an image in the range [min max]. A created header is added to the
%   input structure. 
%
%   tom_imagesc(input,'fixed','header',expr_Header,'range',[min max]) displays 
%   the real size of the matrix input in the range [min max]. The header
%   expression expr_Header is added to the input structure. 'fixed' can be 
%   replaced by 'fixed noinfo' or 'noinfo'.
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   created by SN 12/08/02
%   updated by WDN 02/27/04
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


%if nargin==1
input=varargin{1};
param.fixinfo='';
param.range='';

if nargin==2
    input=varargin{1};
    st=varargin{2};
    switch st
        case  {'fixed','fixed noinfo','noinfo'}
            param.fixinfo=st;
        case 'header'
            input=tom_emheader(input);%input.Value=input.Value';            
            set(gcf,'Userdata',input);
    end
elseif nargin==3
    input=varargin{1};
    st=varargin{2};
    switch st
        case {'fixed','fixed noinfo','noinfo'}
            param.fixinfo=st;
            if strfind('header',varargin{3})
                input=tom_emheader(input);%input.Value=input.Value';            
                set(gcf,'Userdata',input);
            end
        case 'header'
            tmp=input;
            clear input;
            input.Value=tmp;input.Header=varargin{3};            
            
            set(gcf,'Userdata',input);
        case 'range'
            param.range=varargin{3};                        
    end
elseif nargin==4
    input=varargin{1};
    st=varargin{2};
    switch st
        case {'fixed','fixed noinfo','noinfo'}
            param.fixinfo=st;
            head_or_ran=varargin{3};
            if strfind('header',head_or_ran)
                input.Value=input;input.Header=varargin{4};
                set(gcf,'Userdata',input);
            elseif strfind('range',head_or_ran)
                param.range=varargin{4};
            end
        case 'header'
            input=tom_emheader(input);%input.Value=input.Value';            
            set(gcf,'Userdata',input);
            param.range=varargin{4};                        
    end
    
elseif nargin==5
    input=varargin{1};
    st=varargin{2};
    switch st
        case {'fixed','fixed noinfo','noinfo'}
            param.fixinfo=st;
            input=tom_emheader(input);%input.Value=input.Value'; 
            set(gcf,'Userdata',input);            
        case 'header'
            input.Value=input;input.Header=varargin{3};
            set(gcf,'Userdata',input);            
    end  
    param.range=varargin{5};
elseif nargin==6
    input=varargin{1};
    param.fixinfo=varargin{2};
    input.Value=input;input.Header=varargin{4};
    set(gcf,'Userdata',input);
    param.range=varargin{6};  
end


wasstruct=0;
if isstruct(input) 
    in=input.Value';
    wasstruct=1;  clear input.Value;
else 
    in=input'; clear input;
end;    
in_red=imresize(in,.25);

if isempty(param.range)
    [meanv max min std]=tom_devinternal(in);
    if (meanv-2*std)>=(meanv+3*std)
        imagesc(in);
    else
        imagesc(in,[meanv-(3*std) meanv+(3*std)]);
    end;
else
    imagesc(in,[param.range]);
end
axis image;axis ij;
colormap gray;
if ~isempty(param.fixinfo)
    switch param.fixinfo
        case 'fixed'
            set(gca,'Units','pixels');
            pp=get(gca,'Position');sf=size(in);            
            set(gca,'Position',[pp(1) pp(2) sf(1) sf(2)]);
            %t=title(['Info: (' num2str(size(in)) ') pixel']);
            if wasstruct==1
                if input.Header.Objectpixelsize==0
                    dose=1000000; %This is to prevent an error message
                else
                    dose=(meanv./input.Header.CountsPerElectron)./(input.Header.Objectpixelsize.^2);
                end
                t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3) ' dose:' num2str(dose,3) ]);
            else
                meanv=mean2(mean2(in));
                std=std2(std2(in));
                t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3)]);                 
            end
        case 'fixed noinfo'
            set(gca,'Units','pixels');
            pp=get(gca,'Position');sf=size(in);            
            set(gca,'Position',[pp(1) pp(2) sf(1) sf(2)]);
        case 'noinfo'
    end
else
    if wasstruct==1
        if input.Header.Objectpixelsize==0
            dose=1000000; %This is to prevent an error message
        else
            dose=(meanv./input.Header.CountsPerElectron)./(input.Header.Objectpixelsize.^2);
        end
        t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3) ' dose:' num2str(dose,3) ]);
    else
        meanv=mean2(mean2(in));
        std=std2(std2(in));
        t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3)]);
    end
    %t=title(['Info: (' num2str(size(in,1)) ' x ' num2str(size(in,2)) ') pixel']);
end
%elseif nargin==1
%    axis image; axis ij; %colormap gray;
%    if wasstruct==1
%        if input.Header.Objectpixelsize==0
%            dose=1000000; %This is to prevent an error message
%        else
%            dose=(meanv./input.Header.CountsPerElectron)./(input.Header.Objectpixelsize.^2);
%        end
%        t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3) ' dose:' num2str(dose,3) ]);
%    else
%        t=title(['Info: (' num2str(size(in)) ') pixel, mean:' num2str(meanv,3) ' std:' num2str(std,3)]); 
%    end;
%end
clear in;
function [a,b,c,d,e]=tom_devinternal(A);

[s1,s2,s3]=size(A);
a=sum(sum(sum(A)))/(s1*s2*s3);
b=max(max(max(A)));
c=min(min(min(A)));
d=std2(A);
e=d^2;
clear A;

