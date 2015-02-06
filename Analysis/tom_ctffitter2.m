function [pbest, success] = tom_ctffitter2(filename,filter_w1,filter_w2,decay_width,mask_w1,mask_w2)
%TOM_CTFFITTER creates ...
%
%   [Dz, success] = tom_ctffitter(filename, startval, filterval, split, demomode,norings,lowcutoff)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%   startval            ...
%   filterval           ...
%   split               ...
%   demomode            ...
%   norings             ...
%   lowcutoff           ...
%  
%  OUTPUT
%   Dz   		...
%   success		...
%
%EXAMPLE
%   ... = tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

error(nargchk(1, 6, nargin, 'struct'))

if nargin < 6
    mask_w2 = 0.3;
end

if nargin < 5
    mask_w1 = 0.035;
end

if nargin < 4
    decay_width = 0.02;
end

if nargin < 3
    filter_w2 = 0.3;
end

if nargin < 2
    filter_w1 = 0.02;
end

demomode = 1;

%open the image
if ischar(filename)
    try
        im = tom_emreadc(filename);
    catch
        error('Could not open file.');
    end

    %calculate and integrate power spectrum
    im.Value = tom_calc_periodogram(single(im.Value),256);
    im.Value = tom_psd_enhance(im.Value,true,true,filter_w1,filter_w2,decay_width,mask_w1,mask_w2);
    im.Header.Size = [256 256 1];
    ps = tom_cart2polar(single(im.Value));
    ps = sum(ps,2)./(size(ps,2));
    ps = tom_norm(ps,1)-.5;
else
    im = filename;
    im.Value = tom_psd_enhance(im.Value,true,true,filter_w1,filter_w2,decay_width,mask_w1,mask_w2);
    im.Header.Size = [256 256 1];
    ps = tom_cart2polar(single(im.Value));
    ps = sum(ps,2)./(size(ps,2));
    ps = tom_norm(ps,1)-.5;
end

%generate figure or reuse existing figure
if demomode == 1
    figure(34);cla;
end


%display the experimental power spectrum
if demomode == 1
    plot(ps);hold on;title('Original powerspectrum');set(gca,'XLim',[1 128]); 
end

pixs = im.Header.Size(1);
pix_size = im.Header.Objectpixelsize.*1e-10;
x = 0:1/(pixs*pix_size):1/(2*pix_size);% von, Increment, Nyqvist;

global header;
header = im.Header;

Dz(1) = im.Header.Defocus;


[pbest,fval,exitflag,output]=easyfit(x,ps,Dz,@sinfun,[],[]);

disp(['Dz: ' num2str(pbest(1))]);
disp(output.message);

if demomode == 1
    fitted = sinfun(pbest,x);
    hold on;
    plot(fitted,'r--');
end



function FittedCurve = sinfun(Dz,q)

global header;
voltage = header.Voltage;
Cs = header.Cs/1000;
Dz = Dz.*1e-10;

voltagest=voltage.*(1+voltage./1022000); %for relativistic calc
lambda=sqrt(150.4./voltagest)*10^-10;

FittedCurve = sin( (pi./2).* (Cs.*lambda.^3.*q.^4 - 2.*Dz(1).*lambda.*q.^2) ).^2;
FittedCurve = (FittedCurve-0.5);
FittedCurve = FittedCurve(1:end-1)';
