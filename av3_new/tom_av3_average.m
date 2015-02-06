function average=tom_av3_average(motl,method,threshold,iclass,waitbarflag)
%TOM_AV3_AVERAGE averages subtomograms
%
%   average=tom_av3_average(motl,method,threshold,iclass,waitbarflag)
%
%   tom_av3_average is designed for averaging subtomogram 
%   (PARTICLEFILENAME = 'filename'_#no.em) using the parameters of the
%   MOTL. If THRESHOLD is specified only particle with a CCC >
%   thresold*mean(ccc) are included into the average. If ICLASS is
%   specified only particles of this class will be included into the
%   average. 
%
%   ALIGN                ALIGN structure (see tom_AV3_TRANS_ROT_ALIG for format)
%   method              all angles and shifts can be applied:
%                       'direct' as is, or
%                       'inverse' angles Phi and Psi are flipped and change
%                       its signum. Translation also changes its signum.
%
%PARAMETERS
%
%  INPUT
%   motl                ...
%   method              ...
%   threshold           ...
%   iclass              ...
%   waitbarflag         ...
%  
%  OUTPUT
%   average             ...
%
%EXAMPLE
%   ... = tom_av3_average(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 03/31/05
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


if nargin < 5
    waitbarflag = 0;
end
if nargin<4
    iclass = 0;
end;
if nargin<3
    threshold = 0;
end;
if nargin < 2
    method = 'inverse';
end;

if waitbarflag == 1
    h = waitbar(0,'Creating average particle');
end

indx = find (tom_substr2stream(motl,'CCC') >= 0); 
meanv = mean(tom_substr2stream(motl,'CCC',indx));% find ultimate solution!!!!


try
    tsize = motl(1).Tomogram.Header.Size';
catch
    prompt = {'Enter templatesize x:','Enter templatesize y:','Enter templatesize z:'};
    dlg_title = 'Missing template size';
    num_lines = 1;
    def = {'','',''};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
    tsize = [str2num(answer{1}), str2num(answer{2}), str2num(answer{3})];
end

%avmatrix = zeros(tsize);

itomo_old = 0;
icount = 0;
%loop over all particles in motl-list

for indpart = 1:size(motl,2) 
    if (motl(indpart).CCC)>=threshold*meanv && motl(indpart).Class == iclass
        icount = icount +1;
        itomo = motl(indpart).Filename;
        xshift = motl(indpart).Shift.X;
        yshift =motl(indpart).Shift.Y;
        zshift = motl(indpart).Shift.Z;
        tshift = [xshift yshift zshift];
        phi= motl(indpart).Angle.Phi;
        psi=motl(indpart).Angle.Psi;
        the=motl(indpart).Angle.Theta;
        ifile = indpart; %????
%       name = [particlefilename '' num2str(ifile) '.em'];
        name=motl(indpart).Filename;
        particle = tom_emreadc(name); particle = particle.Value;
       % disp(name);
        %       particle = tom_limit(particle,-3,4,'z'); % throw away the gold
        if icount == 1
            wei = zeros(size(particle,1),size(particle,2),size(particle,3));
            average = wei;
        end;
%        if itomo_old ~= itomo %wedge stuff - exact weighting according to list
            %xx = find(wedgelist(1,:) == itomo);
            %minangle= wedgelist(2,xx);
            %maxangle= wedgelist(3,xx);

            maxangle = motl(indpart).Tomogram.AngleMin;
            minangle = motl(indpart).Tomogram.AngleMax;
            %FIXME
            if maxangle == 0
                maxangle = -65;
            end

            if minangle == 0
                minangle = 65;
            end

            
            wedge = tom_av3_wedge(particle,minangle,maxangle);
            itomo_old = itomo;
%        end;
        if isequal(method,'inverse')
            particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi -phi -the]));
        elseif isequal(method,'direct')
            %rotation matrix
            %particle = double(tom_rotate(tom_shift(particle,tshift),motl(indpart).Angle.Rotmatrix));
            particle = double(tom_rotate(tom_shift(particle,tshift),[phi psi the]));
        elseif isequal(method,'direct_dd')
            particle = double(tom_shift(tom_rotate(particle,[phi psi the]),tshift));
        else %sum
            particle = double(tom_shift(tom_rotate(particle,motl(indpart).Angle.Rotmatrix),tshift));
        end;
        %tom_emwrite(['particle_aligned/particle_' num2str(indpart) '.em'],particle);
%        particle = tom_bandpass(particle,0,50,5);
        
        particle = tom_norm(particle,1);
        %avmatrix(:,:,indpart) = particle(:,:,round(size(particle,3)./2));
        average = average + particle;
       % tom_emwrite(['particle_aligned_' num2str(indpart) '.em'],particle);
        if isequal(method,'inverse')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        elseif isequal(method,'direct')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[phi psi the])),0.5,1,'z'),0,0.5);
         elseif isequal(method,'direct_dd')
              tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[phi psi the])),0.5,1,'z'),0,0.5);
        else
            %rotation matrix
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,motl(indpart).Angle.Rotmatrix)),0.5,1,'z'),0,0.5);
            %tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[phi psi the])),0.5,1,'z'),0,0.5);
        end;
        wei = wei + tmpwei;
        if waitbarflag == 1
            waitbar(indpart./size(motl,2),h,[num2str(indpart), ' of ', num2str(size(motl,2)), ' files done.']);
        elseif waitbarflag == 0
            disp(['Particle no ' num2str(ifile) ' added to average'  ]);
        end
    end;%if - threshold
end;
lowp = floor(size(average,1)/2)-3;
wei = 1./wei;
rind = find(wei > 100000);
wei(rind) = 0;% take care for inf
average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));

if waitbarflag == 0
    disp(['Averaging finished - ' num2str(icount) ' particles averaged ... '  ]);
elseif waitbarflag == 1
    close(h);
end