function [Align_new,new_average]=tom_av3_apply_alignment_particles(Align,iteration_nr,particle_list,new_particle_path,new_particle_extension,aligned_particle_name,aligned_particle_extension)
%TOM_AV3_APPLY_ALIGNMENT_PARTICLES creates an average. Apply alignment ... 
%to e.g. unbinned particles
%
%   [Align_new,new_average]=tom_av3_apply_alignment_particles(Align,iteration_nr,particle_list,new_particle_path,new_particle_extension,aligned_particle_name,aligned_particle_extension)
%
%preselect particles by their index: serie1_list=0;ii=1;for i=1:size(serie1,2);if serie1(i)==1;serie1_list(ii)=i;ii=ii+1;end;end;
%
%check for Tempfilename and Filename ... TODO
%
%PARAMETERS
%
%  INPUT
%   Align                       ...
%   iteration_nr                ...
%   particle_list               ...
%   new_particle_path           ...
%   new_particle_extension      ...
%   aligned_particle_name       ...
%   aligned_particle_extension  ...
%  
%  OUTPUT
%   Align_new                   ...
%   new_average                 ...
%
%EXAMPLE
% [anew,average]=tom_av3_apply_alignment_particles(Align,4,[1:140],'./bin1/part26S_','.em','./aligned/part26S_','.em');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN mm/dd/yy
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


%part_alig=tom_reademheader(Align(1,1).Tempfilename);
%part_new=tom_reademheader([new_particle_path num2str(particle_list(1)) new_particle_extension]);
binning_factor = 1;
%binning_factor=part_new.Header.Size(1)./part_alig.Header.Size(1);


ii=1;
for i=1:size(particle_list,2)
    
    idx=particle_list(i);
    
    unterstrich=findstr(Align(1,idx).Tempfilename,'_');

    try
        punkt=findstr(Align(1,idx).Tempfilename,'.');
    catch
        punkt=size(Align(1,idx).Tempfilename);
    end;
    if isempty(punkt); punkt=size(Align(1,idx).Tempfilename)+1; end;
    unterstrich_index=unterstrich(size(unterstrich,2));
    punkt_index=punkt(size(punkt,2));
    
    fn = Align(1,idx).Tempfilename;
    
    filenr=fn(unterstrich_index+1:punkt_index-1);
    
    
    Align_new(1,ii)=Align(iteration_nr,idx);
    Align_new(1,ii).Filename=[new_particle_path num2str(filenr) new_particle_extension];

    for iteration=1:iteration_nr
        angles(iteration,1)=Align(iteration,idx).Angle.Phi;
        angles(iteration,2)=Align(iteration,idx).Angle.Psi;
        angles(iteration,3)=Align(iteration,idx).Angle.Theta;
        shifts(iteration,1)=Align(iteration,idx).Shift.X.*binning_factor;
        shifts(iteration,2)=Align(iteration,idx).Shift.Y.*binning_factor;
        shifts(iteration,3)=Align(iteration,idx).Shift.Z.*binning_factor;
    end;    
    
    [euler_out shift_out rott]=tom_sum_rotation(angles,shifts);
    
    Align_new(1,ii).Angle.Phi=euler_out(1);
    Align_new(1,ii).Angle.Psi=euler_out(2);
    Align_new(1,ii).Angle.Theta=euler_out(3);
    Align_new(1,ii).Angle.Rotmatrix=rott;
    
    Align_new(1,ii).Shift.X=shift_out(1);
    Align_new(1,ii).Shift.Y=shift_out(2);
    Align_new(1,ii).Shift.Z=shift_out(3);
    
    
    
    ii=ii+1;
    
end;

new_average=tom_av3_particles_and_average(Align_new,'sum',0,0,0,aligned_particle_name,aligned_particle_extension);

function average=tom_av3_particles_and_average(motl,method,threshold,iclass,waitbarflag,aligned_particle_name,aligned_particle_extension)

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

indx = find (tom_substr2stream(motl,'CCC') >= 0); meanv = mean(tom_substr2stream(motl,'CCC',indx));% find ultimate solution!!!!
% indx = find (motl(1,:) > threshold*meanv);
% motl = motl(:,indx);

% define average in advance
%name = [particlefilename '_' num2str(1) '.em'];
%particle = tom_emread(name);particle = particle.Value;
%wei = zeros(size(particle,1),size(particle,2),size(particle,3));
%average = wei;

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
    if (motl(indpart).CCC)>=threshold*meanv & motl(indpart).Class == iclass
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

            
            wedge = av3_wedge(particle,minangle,maxangle);
            itomo_old = itomo;
%        end;
        if isequal(method,'inverse')
            particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi -phi -the]));
        elseif isequal(method,'direct')
            %rotation matrix
            %particle = double(tom_rotate(tom_shift(particle,tshift),motl(indpart).Angle.Rotmatrix));
            particle = double(tom_rotate(tom_shift(particle,tshift),[phi psi the]));
        else %sum
            particle = double(tom_shift(tom_rotate(particle,motl(indpart).Angle.Rotmatrix),tshift));
        end;
%        tom_emwrite(['./particle_aligned/part_' num2str(indpart) '.em'],particle);
%        particle = tom_bandpass(particle,0,50,5);
        tom_emwritec2([aligned_particle_name num2str(indpart) aligned_particle_extension],particle);
        
        particle = tom_norm(particle,1);
        %avmatrix(:,:,indpart) = particle(:,:,round(size(particle,3)./2));
        average = average + particle;
        if isequal(method,'inverse')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        elseif isequal(method,'direct')
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