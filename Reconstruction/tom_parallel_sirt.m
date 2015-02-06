function Rec = tom_parallel_sirt(Rec3dProject, NumOfWorkers)
%TOM_PARALLEL_SIRT calculates a SIRT reconstruction.
%
%   Rec = tom_parallel_sirt(Rec3dProject)
%
%   In case of singleaxis tilting the wedge is oriented parallel to the y-axis.
%   Furthermore the computation is completely parallel.
%
%PARAMETERS
%
%  INPUT
%   REC3DPROJECT         structure with all reconstruction settings: 
%                                     The following settings are needed: The project must be aligned! Reconstruction size and 
%                                     binning needed, Iterations and Relaxation needed! All other parameters for example 
%                                     AlignTiltaxis, Weighting, Filter NOT needed! Per default: Tiltaxis aligned parallel to y-axis, 
%                                     weighting off, filter off.
%   NUMOFWORKERS      number of workers
%
%  OUTPUT
%   REC   3-d reconstruction
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   created by ME, AK 01/12/07
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

% -------------------------------------------------------------------------
% Rec3dProject -> ParallelSIRT
% -------------------------------------------------------------------------
% Head
ProjectName = Rec3dProject.ProjectName;
% Project
ProjectionsNamePre = Rec3dProject.PARAMETER.ProjectionsName;
ProjectionsPath = Rec3dProject.PARAMETER.ProjectionsPath;
NumOfProj = Rec3dProject.PROJECT.NumOfProj;
% Prepare projections
ProjectionsName = PrepareProjectionsName(NumOfProj, ProjectionsNamePre);
% Alignment
Tiltangles = Rec3dProject.PARAMETER.Tiltangles;
Tiltaxis = Rec3dProject.PARAMETER.Tiltaxis;
ProjDir = Rec3dProject.PARAMETER.ProjDir;
tx = Rec3dProject.PARAMETER.tx;
ty = Rec3dProject.PARAMETER.ty;
% Volume
SizeX = Rec3dProject.VOLUME.SizeX;
SizeY = Rec3dProject.VOLUME.SizeY;
SizeZ = Rec3dProject.VOLUME.SizeZ;
recSize = [SizeX SizeY SizeZ];
% Binning
PreBinning = Rec3dProject.VOLUME.PreBinning;
tx = tx./2.^PreBinning;
ty = ty./2.^PreBinning;
Imdim = Rec3dProject.PROJECT.Imdim;
Imdim = Imdim./(2.^PreBinning);
% Sirtparameter
Iterations = Rec3dProject.METHOD.Iterations;
Relaxation = Rec3dProject.METHOD.Relaxation;
% Relaxation
Relaxation = 1 - Relaxation;
% -------------------------------------------------------------------------

disp(['Start SIRT reconstruction: ' ProjectName]);
disp(['working directory: ' ProjectionsPath{1}]);
disp(['Number of workers: ' num2str(NumOfWorkers)]);
disp(['Iterations: ' num2str(Iterations)]);
disp(['Relaxation: ' num2str(Relaxation)]);

% Align wedge parallel to the y-axis
ProjDir(1:NumOfProj) = 0;



%generate aligned measured projections
%-------------------------------------------
disp('generate aligned measured projections');

% Loop over projections
parfor (k=1:NumOfProj)

    % Load projection     
    Proj = tom_emreadc([ProjectionsPath{k} '/' ProjectionsName{k}]);

    % Convert to single
    Proj.Value = single(Proj.Value);

    % PreBinning
    if (PreBinning > 0)
        Proj.Value = tom_bin(Proj.Value, PreBinning);
    end
        
    % Normalization
    ProjDev = mean(mean(Proj.Value));
    Proj.Value = (Proj.Value-ProjDev)./ProjDev;

    % Move projection
    Proj.Value = tom_move(Proj.Value, [-fix(tx(k)) -fix(ty(k))]);
    Proj.Value = fftshift(fft2(Proj.Value));
    Proj.Value = shift_in_fs(Proj.Value, [-(tx(k)-fix(tx(k))) -(ty(k)-fix(ty(k)))]);
    Proj.Value = real(ifft2(ifftshift(Proj.Value)));
    
    % Rotate projection
    Proj.Value = tom_rotate(Proj.Value, -Tiltaxis(k), 'linear', [Imdim./2 Imdim./2]);

    % Save SIRT measured projections
    tom_emwrite([ProjectionsPath{k} '/SIRT_mProj_' ProjectName '_' num2str(k) '.em'], Proj.Value);

end



%calculate pathlength projections
%----------------------------------------------
disp('calculate pathlength projections');

% Calculate SIRT_pProjections
parfor (k=1:NumOfProj)

    % Initialize PathlengthVolume
    PathlengthVolume = ones(recSize(1), recSize(2), recSize(3),'single');
    
    % Pathlength 'ones' / Align Tiltaxis
    PathlengthProj = tom_proj3d2c(double(PathlengthVolume), [ProjDir(k) Tiltangles(k)]);
    PathlengthProj = tom_rotate(PathlengthProj, [-ProjDir(k) 0 0], 'linear', [Imdim./2 Imdim./2]);
    PathlengthProj = tom_limit(PathlengthProj,1,1000);
    
    % Save SIRT_pProjections
    tom_emwrite([ProjectionsPath{k} '/SIRT_pProj_' ProjectName '_' num2str(k) '.em'], PathlengthProj);
       
end



%calculate start volume
%---------------------------------------------
disp('calculate start volume');

% Initialize Volume
tom_emwritec([ProjectionsPath{1} '/Rec.vol'], [recSize(1), recSize(2), recSize(3)],'new');
unix(['chmod 666 ' ProjectionsPath{1} '/Rec.vol']);

% Calculate slices
slicez = recSize(3)./NumOfWorkers;
z_offset = 0:slicez:(recSize(3)-slicez);

% Backprojection loop
parfor (i=1:NumOfWorkers)
    Rec_worker = zeros(recSize(1), recSize(2), slicez, 'single');    
    for k=1:NumOfProj
         % Load processed projection
         Proj = tom_emreadc([ProjectionsPath{k} '/SIRT_mProj_' ProjectName '_' num2str(k) '.em']);
         z_offset2=z_offset(i) - recSize(3)/2 + slicez/2;
         % Backprojection / Align Tiltaxis
         Proj.Value = single(Proj.Value);
         tom_backproj3d(Rec_worker, Proj.Value, ProjDir(k), Tiltangles(k), [0 0 z_offset2]);
    end
    
%     % Normalize start volume
%     tom_emwritec([ProjectionsPath{1} '/Rec.vol'],Rec_worker./NumOfProj,'subregion',[1 1 z_offset(i)+1],[recSize(1) recSize(2) slicez]);
    
    % Normalize start volume ENHANCED
    tom_emwritec([ProjectionsPath{1} '/Rec.vol'],Rec_worker,'subregion',[1 1 z_offset(i)+1],[recSize(1) recSize(2) slicez]);
    
end

% -------------------------------------------------------------------------
% SIRT Reconstruction
% -------------------------------------------------------------------------

for i=1:Iterations
    disp(['iteration' num2str(i)]); 
    
    %-----------------dPROJECTIONOPERATOR----------------
    disp('calculating difference projections'); 
    parfor (k=1:NumOfProj)
        
        cRec = tom_emreadc([ProjectionsPath{1} '/Rec.vol']);
        cRec = cRec.Value;
        cRec = cRec.*Relaxation;
        
        % Load pathlength projection
         pProj = tom_emreadc([ProjectionsPath{k} '/SIRT_pProj_' ProjectName '_' num2str(k) '.em']);
         
         % Load measured projections
         mProj = tom_emreadc([ProjectionsPath{k} '/SIRT_mProj_' ProjectName '_' num2str(k) '.em']);
         
         % Calculate correction projection cProj / Align Tiltaxis
         cProj = tom_proj3d2c(single(cRec), [ProjDir(k) Tiltangles(k)]);
         %cProj = p_forward(cRec,ProjDir(k), Tiltangles(k));
         cProj = tom_rotate(cProj, [-ProjDir(k) 0 0], 'linear', [Imdim./2 Imdim./2]);
         
         % Norm correction projection with pathlength projection
         cProj = cProj./pProj.Value;
         
         % Calculate difference projection dProj
         dProj = mProj.Value-cProj;
         
         % Save difference projection dProj
         dProj = tom_emheader(dProj);
         dProj.Header.Tiltangle = Tiltangles(k);
         dProj.Header.Tiltaxis = ProjDir(k);
         tom_emwrite([ProjectionsPath{k} '/SIRT_dProj_' ProjectName '_' num2str(k) '.em'],dProj);
         
    end
    %----------------------------------------------------
    
    %-----------------DIFFERENCEVOLUME------------------
    disp('calculating difference volume');
    % Initialize difference volume
    
    
    % Calculate difference volume
    
    % Initialize Volume
    tom_emwritec([ProjectionsPath{1} '/dRec.vol'], [recSize(1), recSize(2), recSize(3)],'new');
    unix(['chmod 666 ' ProjectionsPath{1} '/dRec.vol']);

    % Backprojection loop
    slicez = recSize(3)/NumOfWorkers;
    z_offset = 0:slicez:recSize(3)-slicez;

    parfor (ii=1:NumOfWorkers)
        Rec_worker = zeros(recSize(1), recSize(2), slicez, 'single');    
        for k=1:NumOfProj
            % Load processed projection
            Proj = tom_emreadc([ProjectionsPath{k} '/SIRT_dProj_' ProjectName '_' num2str(k) '.em']);
            z_offset2=z_offset(ii) - recSize(3)/2 + slicez/2;
            % Backprojection / Align Tiltaxis
            Proj.Value = single(Proj.Value);
            tom_backproj3d(Rec_worker, Proj.Value, ProjDir(k), Tiltangles(k), [0 0 z_offset2]);
        end
        
        % Normalize differencevolume
        tom_emwritec([ProjectionsPath{1} '/dRec.vol'],Rec_worker./NumOfProj,'subregion',[1 1 z_offset(ii)+1],[recSize(1) recSize(2) slicez]);
    end
    
    % Apply difference volume
    Rec = tom_emreadc([ProjectionsPath{1} '/Rec.vol']);
    dRec = tom_emreadc([ProjectionsPath{1} '/dRec.vol']);
    Rec.Value = Rec.Value .* Relaxation + dRec.Value;
    
    disp('saving reconstruction');
    tom_emwrite([ProjectionsPath{1} '/Rec.vol'],Rec.Value);
    clear dRec;
end









% -------------------------------------------------------------------------
% PrepareProjectionsName
% -------------------------------------------------------------------------
function [ProjectionsName] = PrepareProjectionsName(NumOfProj, ProjectionsNamePre)

for k=1:NumOfProj
    
    UnderlineIndex = findstr(ProjectionsNamePre{k}, '_');
    WhichUnderline = size(UnderlineIndex, 2);
    UnderlineIndex = UnderlineIndex(1,WhichUnderline);
    
    WordLength = size(ProjectionsNamePre{k}, 2);
    
    if strcmp(ProjectionsNamePre{k}(UnderlineIndex+1), '0') == 1 && ...
       strcmp(ProjectionsNamePre{k}(UnderlineIndex+2), '0') == 1
   
         for i=1:UnderlineIndex
              ProjectionsName{k}(i) = ProjectionsNamePre{k}(i);
         end
         
         for i=(UnderlineIndex+3):WordLength
              ProjectionsName{k}(i-2) = ProjectionsNamePre{k}(i);
         end
    
    elseif strcmp(ProjectionsNamePre{k}(UnderlineIndex+1), '0') == 1 && ...
       strcmp(ProjectionsNamePre{k}(UnderlineIndex+2), '0') ~= 1
         
         for i=1:UnderlineIndex
              ProjectionsName{k}(i) = ProjectionsNamePre{k}(i);
         end
         
         for i=(UnderlineIndex+2):WordLength
              ProjectionsName{k}(i-1) = ProjectionsNamePre{k}(i);
         end
    
    else
        
         ProjectionsName{k} = ProjectionsNamePre{k};
    
    end

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Shift in Fourier Space
% -------------------------------------------------------------------------
function im_out=shift_in_fs(im,v)
[dimx,dimy]=size(im);
[x,y]=ndgrid( -floor(dimx/2):-floor(dimx/2)+(dimx-1),...
-floor(dimy/2):-floor(dimy/2)+dimy-1);
v = v./[dimx dimy];
x = v(1)*x + v(2)*y;clear y;
im_out=(im.*exp(-2*pi*i*x));
% -------------------------------------------------------------------------
