function tom_imagesc(varargin)
%TOM_IMAGESC displays a matrix in EM style 
%    The EM display style goes from left to right an top bottom. 
%    Matlab convention is in the transposed way. To speed up    
%    the procedure a highly binned image is calculated and used
%    for the scaling of the image.
%
%    tom_imagesc(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin    ...
%  
%  OUTPUT
%
%   Syntax: 
%   tom_imagesc(input) displays matrix input as an image.    
%
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
%   load earth; tom_imagesc(X,'header','range',[10 100]);
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMBROWSE, TOM_EMHEADER, TOM_CTFFIT, IMAGESC
%
%   created by SN 08/12/02
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



% check "children and menu, thanks!' SN
%set(h,'Tag',num2str(h),'Userdata',input);
%im=get(1,'Children'); s=size(im,1);
%if s==0
h=gcf;z=0;
f=findobj('Label','Process');
p=get(f,'Parent');
if isempty(p)
    menu0 = uimenu (gcf,'label', 'Process');    	    
    uimenu (menu0,'label', 'Calculate PS','Callback','A = getimage(gcf);A=A'';A=tom_ps(A);A=tom_emheader(A);tom_imagesc(A);');
%    uimenu (menu0,'label', 'Contrast','Callback','n=gca;m=findobj(gcf);input=get(m(1),''Userdata'');set(m(1),''Doublebuffer'',''on'');input=tom_emheader(input);tom_isc(n,input)');
%    uimenu (menu0,'label', 'Fit CTF manually','Callback','m=findobj(gcf);input=get(m(1),''Userdata'');input=tom_emheader(input);tom_ctffit(input);');
%    uimenu (menu0,'label', 'Fit CTF manually','Callback','m=findobj(gcf);input=get(m,''Userdata'');tom_ctffit(input,input.Header.Defocus);');
%    uimenu (menu0,'label', 'Save as EM file','Callback','A = getimage(gcf);A=A''; A=tom_emheader(A);tom_emwrite(A);');	    
%    end
    menu1=uimenu(gcf,'label','Save as');
%    uimenu (menu1,'label','EM file','Callback','A=get(gcf,''Userdata'');if ~isstruct(A); A=A'';A=tom_emheader(A);end;tom_emwrite(A);'); 
    uimenu (menu1,'label','EM file','Callback','A=getimage(gcf);A=A'';tom_emwrite(A);');
else
    if size(p,1)==1 % p is a number with 1 argument
        p={p};% p is now a cell
    end
    for j=1:size(p,1)
        if cell2mat(p(j))==h% p is a cell when there is at leat 2 arguments.
            z=z+1;
            break;
        end
    end
    if j==size(p,1) & z==0
        menu0 = uimenu (gcf,'label', 'Process');    	    
        uimenu (menu0,'label', 'Calculate PS','Callback','A = getimage(gcf);A=A'';A=tom_ps(A);A=tom_emheader(A);tom_imagesc(A);');
        menu1=uimenu(gcf,'label','Save as');
        %uimenu (menu1,'label','EM file','Callback','A = get(gcf,''Userdata'');if ~isstruct(A); A=A'';A=tom_emheader(A);disp(''in loop if'');disp(A);end;disp(''Not in loop if'');disp(A);tom_emwrite(A);disp(''Not in loop if'');disp(A)');      
        uimenu (menu1,'label','EM file','Callback','A=getimage(gcf);A=A'';tom_emwrite(A);');
    end
end

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

