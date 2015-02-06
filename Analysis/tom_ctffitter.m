function [Dz, success,envparams,noiseparams] = tom_ctffitter(filename, startval, filterval, split, binning, demomode,norings,lowcutoff)
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

error(nargchk(1, 8, nargin, 'struct'))

if nargin < 8
    lowcutoff = 10;
end

if nargin < 7
    norings = 5;
end

if nargin < 6
    binning = 0;
end;

if nargin < 5
    demomode = 1;
end

if nargin < 4
    split = 1;
end

if nargin < 3
    filterval = 1;
end

if nargin < 2
    startval = 0;
end

%open the image
if ischar(filename)
    try
        %im = tom_emreadc(filename,'binning',binning);
        im = tom_emreadc(filename);
    catch
        error('Could not open file.');
    end

    %lowcutoff = lowcutoff ./ 2^binning ./ (split./2);
    
    if filterval > 1
        im.Value = tom_filter(im.Value,filterval,'quadr','real');
    end

    %calculate and integrate power spectrum
    im.Value = tom_calc_periodogram(single(im.Value),256);
    im.Value = tom_xmipp_psd_enhance(im.Value,true,true,0.02,0.4,0.02,0.02,0.4);
    im.Header.Size = [256 256 1];
    ps = tom_cart2polar(im.Value);
    ps = sum(ps,2)./(size(ps,2));
    %ps = smooth(ps,round(length(ps).*.05))';
%    ps = tom_image2ctf(im.Value,1,0,1,0);

else
    ps = filename;
end
%close(gcf);

%cutoff frequencies
%highcutoff = round(length(ps)-length(ps)./10);
%lowcutoff = round(lowcutoff);
%ps_fit = ps(lowcutoff:highcutoff);

ps_fit=ps;
highcutoff=128;
lowcutoff=1;

%generate figure or reuse existing figure
if demomode == 1
    figure;
    %     if isempty(findobj('Tag','ctffitterfigure'))
    %         figure('Tag','ctffitterfigure');
    %     else
    %         clf(findobj('Tag','ctffitterfigure'));
    %     end
end


%display the experimental power spectrum
if demomode == 1
    subplot(3,1,1);plot(ps);hold on;title('Original powerspectrum');set(gca,'XLim',[1 128]);
    set(gcf,'Position',[281   128   969   936]);
end

%fit the noise function to the experimental power spectrum
lb = [];
ub = [];
Aeq = [1 sqrt(length(ps)) length(ps) length(ps).^2];
beq = ps(end);
nonlcon = [];
%x0 = [10 1 1e-6 1e-6]';
x0 = [0 0 0 0]';
snoise = double([lowcutoff:highcutoff])';
A = [ones(size(snoise,1),1) sqrt(snoise) snoise snoise.^2];
b = double(ps_fit);
options = optimset('LargeScale','off','Display','off','TolCon',1e-4,'TolFun',1e-4,'MaxFunEvals',300,'MaxIter',30);
[noiseparams,fun,eflag] = fmincon(@noiseobjfun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options,A,b);
noisefunction = noiseplotfun([1:length(ps)]',noiseparams);

if demomode == 1
    plot(noisefunction,'--r');
end

%calculation of noise free power spectrum
ps_noisefree = ps - noisefunction;

%cutoff frequencies
lowcutoff = round((im.Header.Size(1)./11)./split);
highcutoff = round((im.Header.Size(1)./2.5)./split);
ps_fit = ps_noisefree(lowcutoff:highcutoff);
%ps_fit = ps(lowcutoff:highcutoff);
%fit the envelope function to the experimental power spectrum
lb = [];
ub = [];
Aeq = [1 1 1 1;1 sqrt(length(ps_noisefree)) length(ps_noisefree) length(ps_noisefree).^2];
beq = [ps_noisefree(1);ps_noisefree(end)];
nonlcon = [];
x0 = [max(ps_fit) 0 0 0]';
senv = double([lowcutoff:highcutoff])';
A = [ones(size(senv,1),1) sqrt(senv) senv senv.^2];
b = double(ps_fit');
options = optimset('LargeScale','off','Display','off','TolCon',1e-4,'TolFun',1e-4,'MaxFunEvals',300,'MaxIter',30);
[envparams,fun,eflag] = fmincon(@envobjfun, x0, -A, -b, Aeq, beq, lb, ub, nonlcon, options,A,b);
envfunction = envplotfun([1:length(ps)]',envparams);

envfunction = envfunction + noisefunction;
avgfunction = (envfunction + noisefunction)./2;

if demomode == 1
    plot(envfunction,'--r');
    plot(avgfunction,'--g');
end



%plot function without noise and envelope
finalps = (ps - noisefunction)./(envfunction-noisefunction);

if demomode == 1
    subplot(3,1,2);plot([lowcutoff:highcutoff],tom_norm(finalps(lowcutoff:highcutoff),1));hold on;
    plot(zeros(length(finalps),1),'Color',[0 0 0]);title('Powerspectrum without envelope functions');set(gca,'XLim',[1 512]);
end

%get zeros of experimental ctf
[xx,expzeros] = crossing(finalps(lowcutoff:highcutoff));
try
    expzeros = expzeros(1:norings.*2)+lowcutoff-1;
catch
    expzeros = expzeros+lowcutoff-1;
end

%calculate first order difference from corrected ctf
finalps_diff1 = diff(finalps);
finalps_diff1 = smooth(finalps_diff1,round(5*length(ps)./512))';

%get minima/maxima of experimental ctf
[xx,expminmax] = crossing(finalps_diff1(lowcutoff:highcutoff));
try
    expminmax = expminmax(1:norings.*2)+lowcutoff-1;
catch
    expminmax = expminmax+lowcutoff-1;
end

%plot zeros and minima/maxima
if demomode == 1
    %    plot(expzeros,zeros(length(expzeros),1),'ro');
    plot(expminmax,finalps(round(expminmax)),'go');
end


%plot first order difference
if demomode == 1
    subplot(3,1,3);plot([lowcutoff:highcutoff],finalps_diff1(lowcutoff:highcutoff));hold on;title('First order derivation of powerspectrum without envelope functions');set(gca,'XLim',[1 512]);
    plot(zeros(length(finalps),1),'Color',[0 0 0]);
    plot(expminmax,zeros(length(expminmax),1),'go');
end

%fit defocus
%try
    %highcutoff = round(expminmax(norings.*2));
%catch
    %highcutoff = 200;
%end
lowcutoff=1;
highcutoff = 100;

%startval=0; % changed by SN, use default image.header defocus setting
[Dz, Dnom, resid, exitflag] = fitcurve(im.Header,finalps,highcutoff,lowcutoff,split,demomode,startval);

%correct the output values for the split
%expzeros = expzeros .* split;
%expminmax = expminmax .* split;

if exitflag > 0
    %if abs(im.Header.Defocus-Dz) < 2000 & sum(envfunction(200:400)-noisefunction(200:400)) > thresh
    success = 1;
    disp([im.Header.Filename '  Nominal Dz: ' num2str(Dnom) ', fitted Dz: ' num2str(Dz) '      OK']);

else
    success = 0;
    disp([im.Header.Filename '  Nominal Dz: ' num2str(Dnom) ', fitted Dz: ' num2str(Dz) '      FAILED']);

end
drawnow;


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise fit function                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = noiseobjfun(x,A,b)
bhat = ones(length(b),1);
l = length(bhat);
bhat(round(l/2):end)=0;
val = norm((A*x-b).*bhat,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envelope fit function                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = envobjfun(x,A,b)
bhat = ones(length(b),1);
l = length(bhat);
bhat(round(l/2):end)=0;
val = norm((A*x-b).*bhat,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise function                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = noiseplotfun(x,n)
for i=1:size(x,1)
    y(i) = (n(1)+n(2).*(double(i).^(1/2))+n(3).*double(i)+n(4).*double(i).^2);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envelope function                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = envplotfun(x,n)
for i=1:size(x,1)
    y(i) = (n(1)+n(2).*(double(i).^(1/2))+n(3).*double(i)+n(4).*double(i).^2);
    %y(i) = (n(1).*double(i).^2);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find zeros                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind,t0] = crossing(S,t,level,par)
% CROSSING find the crossings of a given level of a signal
%   ind = CROSSING(S) returns an index vector ind, the signal
%   S crosses zero at ind or at between ind and ind+1
%   [ind,t0] = CROSSING(S,t) additionally returns a time
%   vector t0 of the zero crossings of the signal S. The crossing
%   times are linearly interpolated between the given times t
%   [ind,t0] = CROSSING(S,t,level) returns the crossings of the
%   given level instead of the zero crossings
%   ind = CROSSING(S,[],level) as above but without time interpolation
%   [ind,t] = CROSSING(S,t,level,par) allows additional parameters
%       par.interpolation = {'none'|'linear'}
%

% Steffen Brueckner, 2002-09-25
% Copyright (c) Steffen Brueckner, 2002
% brueckner@sbrs.net

% check the number of input arguments
error(nargchk(1,4,nargin));

if nargin < 2 || isempty(t)
    % if not given: time vector == index vector
    t = 1:length(S);
elseif length(t) ~= length(S)
    error('t and S must be of identical length!');
end
if nargin < 3
    level = 0;
end
if nargin < 4
    par.interpolation = 'linear';
end

S   = S - level;
S1  = S(1:end-1) .* S(2:end);
ind = find(S1 <= 0);

% get times of crossings
t0 = double(t(ind));

if strcmp(par.interpolation,'linear')
    % linear interpolation of crossing
    for ii=1:length(t0)
        if abs(S(ind(ii))) > eps
            % interpolate only when data point is not already zero
            t0(ii) = t0(ii) - S(ind(ii)) .* (t(ind(ii)+1) - t(ind(ii))) ./ (S(ind(ii)+1) - S(ind(ii)));
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [estimates, d, residuals,exitflag] = fitcurve(header,expfun,highcutoff,lowcutoff,split,demomode, startval)

if startval == 0
    d = header.Defocus;
else
    d =startval.*10;
end

options = optimset('TolX',1e-10,'TolFun',1e-10,'Diagnostics','on','LargeScale','off','Display','on','MaxIter',200,'MaxFunEvals',10000,'FunValCheck','on');
[estimates,resnorm,residuals,exitflag] = lsqnonlin(@sinfun,d,-inf,+inf,options);

[sse,curve] = sinfun(estimates);
if demomode == 1
    subplot(3,1,2);plot(curve,'--g');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    function [sse, FittedCurve] = sinfun(Dz)
        voltage = header.Voltage;
        Cs = header.Cs./1000;
        Dz = Dz.*1e-10;
        pixs = header.Size(1)./split;
        pix_size = header.Objectpixelsize.*1e-10;
        voltagest=voltage.*(1+voltage./1022000); %for relativistic calc
        lambda=sqrt(150.4./voltagest)*10^-10;
        q=0:1/(pixs*pix_size):1/(2*pix_size);% von, Increment, Nyqvist
        FittedCurve = sin( (pi./2).* (Cs.*lambda.^3.*q.^4 - 2.*Dz.*lambda.*q.^2) ).^2;

        ErrorVector = (FittedCurve(lowcutoff:highcutoff) - expfun(lowcutoff:highcutoff));
        sse  = ErrorVector;
    end
end
