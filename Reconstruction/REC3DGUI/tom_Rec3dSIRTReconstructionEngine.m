function [Rec3dReconstruction, Rec3dProject] = tom_Rec3dSIRTReconstructionEngine...
(Rec3dControl, Rec3dProject, ProjNameConv, FeedbackStyle)
%TOM_REC3DSIRTRECONSTRUCTIONENGINE is a module of TOM_REC3DGUI.
%
%   [Rec3dReconstruction,Rec3dProject] = tom_Rec3dSIRTReconstructionEngine...
%   (Rec3dControl,Rec3dProject,ProjNameConv,FeedbackStyle)
%
%   This programme operates a SIRT 3d-reconstruction.
%
%PARAMETERS
%
%  INPUT
%   Rec3dControl               control-structure of Rec3dGUI
%   Rec3dProject               main-structure of Rec3dGUI
%   ProjNameConv               '001' or '1'
%   FeedbackStyle              'guistyle' or 'funstyle'
%
%  OUTPUT
%   Rec3dReconstruction        3d-reconstruction
%   Rec3dProject               aktualized main-structure of Rec3dGUI
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DALIGNMENTENGINE   TOM_PROJ3D2C   TOM_BACKPROJ3D
%
%   01/01/07 ME
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
% Rec3dControl -> ReconstructionEngine
% -------------------------------------------------------------------------
ProtocolMode = Rec3dControl.ProtocolMode;
DelMode = Rec3dControl.DelMode;
DemoMode = Rec3dControl.DemoMode;
% -------------------------------------------------------------------------
% Rec3dProject -> ReconstructionEngine
% -------------------------------------------------------------------------
% HEAD
ProjectName = Rec3dProject.ProjectName;
% PROJECT
TiltingGeometry = Rec3dProject.PROJECT.TiltingGeometry;
NumOfProj = Rec3dProject.PROJECT.NumOfProj;
Imdim = Rec3dProject.PROJECT.Imdim;
% RECONSTRUCTION
ReconstructionName = Rec3dProject.NAMEPATHEXT.ReconstructionName;
ReconstructionPath = Rec3dProject.NAMEPATHEXT.ReconstructionPath;
ReconstructionExt = Rec3dProject.NAMEPATHEXT.ReconstructionExt;
TempFilesName = Rec3dProject.NAMEPATHEXT.TempFilesName;
TempFilesPath = Rec3dProject.NAMEPATHEXT.TempFilesPath;
TempFilesExt = Rec3dProject.NAMEPATHEXT.TempFilesExt;
% RECONSTRUCTION / Parameter
ProjectionsName = Rec3dProject.PARAMETER.ProjectionsName;
ProjectionsPath = Rec3dProject.PARAMETER.ProjectionsPath;
Tiltangles = Rec3dProject.PARAMETER.Tiltangles;
Tiltaxis = Rec3dProject.PARAMETER.Tiltaxis;
ProjDir = Rec3dProject.PARAMETER.ProjDir;
tx = Rec3dProject.PARAMETER.tx;
ty = Rec3dProject.PARAMETER.ty;
% RECONSTRUCTION / Volume
SizeX = Rec3dProject.VOLUME.SizeX;
SizeY = Rec3dProject.VOLUME.SizeY;
SizeZ = Rec3dProject.VOLUME.SizeZ;
PreBinning = Rec3dProject.VOLUME.PreBinning;
PostBinning = Rec3dProject.VOLUME.PostBinning;
% RECONSTRUCTION / Method
ReconstructionMethod = Rec3dProject.METHOD.ReconstructionMethod;
Normalization = Rec3dProject.METHOD.Normalization;
SmoothBorders = Rec3dProject.METHOD.SmoothBorders;
AlignTiltaxis = Rec3dProject.METHOD.AlignTiltaxis;
Handedness = Rec3dProject.METHOD.Handedness;
ApplyWeighting = Rec3dProject.METHOD.ApplyWeighting;
WeightingMethod = Rec3dProject.METHOD.WeightingMethod;
ObjectThickness = Rec3dProject.METHOD.ObjectThickness;
Taper = Rec3dProject.METHOD.Taper;
Iterations = Rec3dProject.METHOD.Iterations;
Relaxation = Rec3dProject.METHOD.Relaxation;
Pathlength = Rec3dProject.METHOD.Pathlength;
% RECONSTRUCTION / Filter
FILTER = Rec3dProject.FILTER;
% RECONSTRUCTION / Detail
DetailMode = Rec3dProject.DETAIL.DetailMode;
OverviewSizeZ = Rec3dProject.DETAIL.OverviewSizeZ;
OverviewPreBinning = Rec3dProject.DETAIL.OverviewPreBinning;
OverviewPostBinning = Rec3dProject.DETAIL.OverviewPostBinning;
NumberOfDetails = Rec3dProject.DETAIL.NumberOfDetails;
DetailCoordinates = Rec3dProject.DETAIL.DetailCoordinates;
av3_alignstruct = Rec3dProject.DETAIL.av3_alignstruct;
% -------------------------------------------------------------------------


% Calculate PreBinning
tx = tx./(2.^PreBinning);
ty = ty./(2.^PreBinning);
Imdim = Imdim./(2.^PreBinning); %Imdim is now binned projection dimension

% Relaxation
Relaxation = 1 - Relaxation;


% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check ReconstructionMethod
if strcmp(ReconstructionMethod, 'SIRT') ~= 1
    % ERROR
    msgbox('ReconstructionMethod must be SIRT!', 'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check ReconstructionName, ReconstructionPath, ReconstructionExt
if strcmp(ReconstructionName, '') == 1 && ...
   strcmp(ReconstructionPath, '') == 1 && ...
   strcmp(ReconstructionExt, '') == 1
    
    ReconstructionName = ProjectName;
    ReconstructionPath = pwd;
    ReconstructionExt = '.vol';
end

% Check ReconstructionPath
if strcmp(ReconstructionPath, '') == 1
    % ERROR
    msgbox('No Reconstruction Path!', 'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check PostBinning
if isequal(PostBinning, 0) ~= 1
    % ERROR
    msgbox('PostBinning must be zero for SIRT reconstruction!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];      
    return;
end

% Check Imdim
if isequal(SizeX, SizeY) ~= 1 || isequal(SizeX, Imdim) ~= 1 || isequal(SizeY, Imdim) ~=1
    % ERROR
    msgbox('SizeX and SizeY and Imdim must be equal for SIRT reconstruction!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];      
    return;
end

% Check Normalization
if strcmp(Normalization, 'off') == 1
    % ERROR
    msgbox('Normalization is set off. This is not possible for SIRT reconstruction!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check SmoothBorders
if SmoothBorders > Imdim
    % ERROR
    msgbox('Smooth Borders is larger than resulting Imdim.', 'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check AlignTiltaxis
if strcmp(TiltingGeometry, 'dualaxis') == 1 && strcmp(AlignTiltaxis, 'Y-axis') == 1
    % ERROR
    msgbox('Align Tiltaxis is set Y-axis. This is not possible for dualaxis projects!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check Handedness
if strcmp(TiltingGeometry, 'dualaxis') == 1 && isequal(Handedness, 180) == 1
    % ERROR
    msgbox('Handedness is set 180. This is not possible for dualaxis projects!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% % Check Weighting
% if strcmp(ApplyWeighting, 'on') == 1
%     % ERROR
%     msgbox('ApplyWeighting is set on. This is not possible for SIRT reconstruction!', ...
%                    'Do Reconstruction', 'error');
%     Rec3dReconstruction = [];
%     return;
% end

% Check DetailMode
if strcmp(DetailMode, 'on') == 1 || strcmp(DetailMode, 'av3') == 1
    % ERROR
    msgbox('DetailMode is set on or av3. This is not possible for SIRT reconstruction!', ...
                  'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end
% -------------------------------------------------------------------------


% FeedbackStyle
if strcmp(FeedbackStyle, 'funstyle') == 1
    ProtocolMode = 'off';
    DemoMode = 'off';
end

% ProtocolMode
if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(ProtocolMode, 'on') == 1
    disp(' ');
    disp('Do Reconstruction');
    disp(' ');
end


% -------------------------------------------------------------------------
% Calculate SIRT_mProjections
% -------------------------------------------------------------------------
% AlignTiltaxis / Modify ProjDir
if strcmp(AlignTiltaxis, 'Y-axis') == 1
   ProjDir(1:NumOfProj) = 0;
end

% Handedness / Modify Tiltaxis
if isequal(Handedness, 180) == 1
    Tiltaxis = Tiltaxis + 180;
end

% Prepare ProjectionsName
if strcmp(ProjNameConv, '001') == 1
    ProjectionsName = PrepareProjectionsName(NumOfProj, ProjectionsName);
end

% Calculate SIRT_mProjections
[NumOfmProjProcessed] = SIRTmProjProcessor...
(FeedbackStyle, ProtocolMode, DemoMode, ...
NumOfProj, Imdim, ...
ReconstructionName, ReconstructionPath, ReconstructionExt, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsName, ProjectionsPath, ...
Tiltangles, Tiltaxis, ProjDir, tx, ty, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, AlignTiltaxis, ...
Normalization, SmoothBorders, Taper, FILTER);

% Check SIRTmProjProcessor
if isequal(NumOfmProjProcessed, NumOfProj) ~= 1
    
    % DelMode
    if strcmp(DelMode, 'on') == 1
         
         for k=1:NumOfmProjProcessed
              workingdir = pwd;
              cd (ReconstructionPath);
              delete(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em']);
              cd (workingdir);
         end
              
    end
    
    Rec3dReconstruction = [];
    return;
    
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Calculate SIRT_pProjections
% -------------------------------------------------------------------------

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    RecTempWaitbar = waitbar(0, 'Calculate pathlength projections', ...
                                                 'CreateCancelBtn', 'closereq');
    set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
    set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          Calculate pathlength projections');
    end
end
% -----Feedback-----

% Initialize Pathlength Volume
if strcmp(Pathlength, 'ones') == 1
    PathlengthVolume = ones(SizeX,SizeY,SizeZ,'double');
end

% Calculate SIRT_pProjections
for k=1:NumOfProj
    
    % Pathlength 'ones' / Align Tiltaxis
    if strcmp(Pathlength, 'ones') == 1
         PathlengthProj = tom_proj3d2c(double(PathlengthVolume), [ProjDir(k) Tiltangles(k)]);
         PathlengthProj = tom_rotate(PathlengthProj, [-ProjDir(k) 0 0], 'linear', [Imdim./2 Imdim./2]);
         PathlengthProj = tom_limit(PathlengthProj,1,100000);
    end
    
    % Pathlength 'slab'
    if strcmp(Pathlength, 'slab') == 1
         PathlengthProj = (SizeZ.*ones(Imdim,Imdim,'single'))./cosd(Tiltangles(k));
    end
    
    % Add and aktualize header
    PathlengthProj = tom_emheader(PathlengthProj);
    PathlengthProj.Header.Tiltangle = Tiltangles(k); 
    PathlengthProj.Header.Size = [size(PathlengthProj.Value,1),size(PathlengthProj.Value,2),1];
    PathlengthProj.Header.Magic(4) = 5;
    PathlengthProj.Header.Tiltaxis = ProjDir(k);
    
    % Save SIRT_pProjections
    for jjj=1:10
         try
              workingdir = pwd;
              cd (ReconstructionPath);
              tom_emwrite(['SIRT_pProj_' ReconstructionName '_' num2str(k) '.em'], PathlengthProj);
              cd (workingdir);
              break;
         catch
              pause(0.5);
              disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
         end 
    end
              
    % Clear PathlengthProj
    clear PathlengthProj;
    
                                
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1
    % Check RecTempWaitbar
    RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         % Cancel button clicked
         if isempty(RecTempWaitbar)
              % DelMode
              if strcmp(DelMode, 'on') == 1
                   % Delete SIRT_mProj
                   for j=1:NumOfProj
                        workingdir = pwd;
                        cd (ReconstructionPath);
                        delete(['SIRT_mProj_' ReconstructionName '_' num2str(j) '.em']);
                        cd (workingdir);
                   end
                   % Delete SIRT_pProj
                   for i=1:k
                        workingdir = pwd;
                        cd (ReconstructionPath);
                        delete(['SIRT_pProj_' ReconstructionName '_' num2str(i) '.em']);
                        cd (workingdir);
                   end
              end
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp(['                    Processed pathlength projection: ' num2str(k)]);
                   disp('          Process cancelled!');
                   disp(' ');
              end
              Rec3dReconstruction = [];
              return;
         end
         % Aktualize RecTempWaitbar
         if ~isempty(RecTempWaitbar)
              waitbar(k./NumOfProj, RecTempWaitbar); drawnow;
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp(['                    Processed pathlength projection: ' num2str(k)]);
              end
         end
    end
    % -----Feedback-----
    
       
end


% Clear PathlengthVolume
clear PathlengthVolume;


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
    if ~isempty(RecTempWaitbar)
         delete(RecTempWaitbar); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp('          Pathlength projections calculated!');
         disp(' ');
    end
end
% -----Feedback-----


% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Calculate SIRT Start Volume
% -------------------------------------------------------------------------
[Rec3dReconstruction, NumOfProjBackprojected] = SIRTCalculateStartVolume...
(FeedbackStyle, ProtocolMode, DemoMode, ...
NumOfProj, Imdim, ...
ReconstructionName, ReconstructionPath, ReconstructionExt, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsName, ProjectionsPath, ...
Tiltangles, Tiltaxis, ProjDir, tx, ty, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, AlignTiltaxis, ...
Normalization, SmoothBorders, FILTER, SizeX, SizeY, SizeZ);

% Check SIRT Start Volume
if isequal(NumOfProjBackprojected, NumOfProj) ~= 1
         
    if strcmp(DelMode, 'on') == 1
    
        for k=1:NumOfProj
            workingdir = pwd;
            cd (ReconstructionPath);
            delete(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em']);
            delete(['SIRT_pProj_' ReconstructionName '_' num2str(k) '.em']);
            cd (workingdir);
        end
    
    end
    
    Rec3dReconstruction = [];
    return;

end

% Norm start volume with NumOfProj
Rec3dReconstruction = Rec3dReconstruction./NumOfProj;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% SIRT Reconstruction
% -------------------------------------------------------------------------
% Zaehler for convergence values
zaehler = 1;


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    RecTempWaitbar = waitbar(0, 'Iterative SIRT loop', 'CreateCancelBtn', 'closereq');
    set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
    set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          Calculate SIRT iterations');
         disp(' ');
    end
end
% -----Feedback-----


for i=1:Iterations
    
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(ProtocolMode, 'on') == 1
         % ProtocolMode
         disp(['                    Iteration: ' num2str(i)]);
         disp(' ');
    end
    % -----Feedback-----
    
    
    %---------------DUBLICATION--------------------------
    cRec = Rec3dReconstruction;
    %----------------------------------------------------
    
    %---------------RELAXATION---------------------------
    cRec = cRec.*Relaxation;
    Rec3dReconstruction = Rec3dReconstruction.*Relaxation;
    %-----------------------------------------------------
    
    %-----------------dPROJECTIONOPERATOR----------------
    for k=1:NumOfProj
        
         % Load pathlength projection
         for jjj=1:10
              try
                   workingdir = pwd;
                   cd (ReconstructionPath);
                   pProj = tom_emread(['SIRT_pProj_' ReconstructionName '_' num2str(k) '.em']);
                   cd (workingdir);
                   break;
              catch
                   pause(0.5);
                   disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
              end 
         end
         
         % Convert to values
         pProj = pProj.Value;
         
         % Load measured projections
         for jjj=1:10
              try
                   workingdir = pwd;
                   cd (ReconstructionPath);
                   mProj = tom_emread(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em']);
                   cd (workingdir);
                   break;
              catch
                   pause(0.5);
                   disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
              end 
         end
         
         % Convert to values
         mProj = mProj.Value;
         
         % Calculate correction projection cProj / Align Tiltaxis
         cProj = tom_proj3d2c(double(cRec), [ProjDir(k) Tiltangles(k)]);
         cProj = tom_rotate(cProj, [-ProjDir(k) 0 0], 'linear', [Imdim./2 Imdim./2]);
         
         % Norm correction projection with pathlength projection
         cProj = cProj./pProj;
         
         % Calculate difference projection dProj
         dProj = mProj-cProj;
         
         % Calculate Euclidian distance per projection
         EuclidianDistancePerProjectionPre(k) = (1./(Imdim.*Imdim)).*sum(sum(dProj.^2));
         
         % Calculate difference projection convergence values
         [a_dp,b_dp,c_dp,d_dp,e_dp] = tom_dev(dProj,'noinfo');
         f_dp = 1;%entropy(double(dProj));
         g_dp = (1./e_dp).*sum(sum((dProj).^2));
         h_dp = tom_ccc(mProj,cProj,'norm');
         
         dProjVal(1,zaehler) = a_dp; %Mean
         dProjVal(2,zaehler) = b_dp; %Max
         dProjVal(3,zaehler) = c_dp; %Min
         dProjVal(4,zaehler) = d_dp; %Std
         dProjVal(5,zaehler) = e_dp; %Var
         dProjVal(6,zaehler) = f_dp; %Entropy
         dProjVal(7,zaehler) = g_dp; %Chi
         dProjVal(8,zaehler) = h_dp; %CCC
         zaehler = zaehler + 1;
         
         % Save difference projection dProj
         for jjj=1:10
              try
                   workingdir = pwd;
                   cd (ReconstructionPath);
                   dProj = tom_emheader(dProj);
                   dProj.Header.Tiltangle = Tiltangles(k);
                   dProj.Header.Tiltaxis = ProjDir(k);
                   tom_emwrite(['SIRT_dProj_' ReconstructionName '_' num2str(k) '.em'],dProj);
                   cd (workingdir);
                   break;
              catch
                   pause(0.5);
                   disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
              end 
         end
         
         % -----Feedback-----
         % ProtocolMode
         if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(ProtocolMode, 'on') == 1
              disp(['                    Mean of difference projection ................... : ' num2str(a_dp)]);
              disp(['                    Euclidian distance of difference projection: ' ...
                                            num2str(EuclidianDistancePerProjectionPre(k))]);
              disp(' ');
         end
         % -----Feedback-----
         
         
    end
    %----------------------------------------------------
    
    %-----------------DIFFERENCEVOLUME------------------
    % Initialize difference volume
    DifferenceVolume = zeros(SizeX,SizeY,SizeZ,'single');
    
    % Calculate difference volume
    for k=1:NumOfProj
    
         % Load difference projections
         for jjj=1:10
              try
                   workingdir = pwd;
                   cd (ReconstructionPath);
                   dProj = tom_emread(['SIRT_dProj_' ReconstructionName '_' num2str(k) '.em']);
                   cd (workingdir);
                   break;
              catch
                   pause(0.5);
                   disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
              end 
         end
         
         % Convert to values
         dProj = dProj.Value;
         
         % Convert to single
         dProj = single(dProj);
         
         % Backprojection / Align Tiltaxis
         tom_backproj3d(DifferenceVolume, dProj, ProjDir(k), Tiltangles(k), [0 0 0]);
    
    end
    
    % Norm difference volume with NumOfProj
    DifferenceVolume = DifferenceVolume./NumOfProj;
    %----------------------------------------------------
    
    % Calculate Euclidian distance
    EuclidianDistancePerProjection{i} = EuclidianDistancePerProjectionPre;
    EuclidianDistancePerIteration(i) = (1./(NumOfProj)).*...
                                                           sum(EuclidianDistancePerProjectionPre);
    clear EuclidianDistancePerProjectionPre;
                                                
    % Calculate difference volume convergence values
    [a_diff,b_diff,c_diff,d_diff,e_diff]=tom_dev(DifferenceVolume,'noinfo');
    f_diff = 1;%entropy(double(DifferenceVolume));
    DifVolVal(1,i) = a_diff; %Mean
    DifVolVal(2,i) = 0; %Max
    DifVolVal(3,i) = 0; %Min
    DifVolVal(4,i) = d_diff; %Std
    DifVolVal(5,i) = e_diff; %Var
    DifVolVal(6,i) = f_diff; %Entropy
    
    % Apply difference volume
    Rec3dReconstruction = Rec3dReconstruction + DifferenceVolume;
    
    % Clear difference volume
    clear DifferenceVolume;
    
    % Calculate reconstruction volume convergence values
    [a_rec,b_rec,c_rec,d_rec,e_rec]=tom_dev(Rec3dReconstruction,'noinfo');
    f_rec = 1;%entropy(double(Rec3dReconstruction));
    RecVolVal(1,i) = a_rec; %Mean
    RecVolVal(2,i) = 0; %Max
    RecVolVal(3,i) = 0; %Min
    RecVolVal(4,i) = d_rec; %Std
    RecVolVal(5,i) = e_rec; %Var
    RecVolVal(6,i) = f_rec; %Entropy
    
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1
         % ProtocolMode
         if strcmp(ProtocolMode, 'on') == 1
              disp(' ');
              disp(['                    Mean of difference volume ..... : ' num2str(a_diff)]);
              disp(['                    Mean of reconstruction volume: ' num2str(a_rec)]);
              disp(['                    Euclidian distance per iteration: ' ...
                                            num2str(EuclidianDistancePerIteration(i))]);
              disp(' ');
         end
         % Check RecTempWaitbar
         RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         % Cancel button clicked
         if isempty(RecTempWaitbar)
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Process cancelled!');
                   disp(' ');
              end
              Rec3dReconstruction = [];
              NumOfIterationsProcessed = i;
              break;
         end
         % Aktualize RecTempWorkbar
         if ~isempty(RecTempWaitbar) 
              waitbar(i./Iterations, RecTempWaitbar); drawnow;
         end
    end
    % -----Feedback-----

    
% Check SIRT iterative loop
NumOfIterationsProcessed = i;

end

% DelMode
if strcmp(DelMode, 'on') == 1
    for k=1:NumOfProj
         workingdir = pwd;
         cd (ReconstructionPath);
         delete(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em']);
         delete(['SIRT_pProj_' ReconstructionName '_' num2str(k) '.em']);
         delete(['SIRT_dProj_' ReconstructionName '_' num2str(k) '.em']);
         cd (workingdir);
    end
end

% Check SIRT iterative loop
if isequal(NumOfIterationsProcessed, Iterations) ~= 1
    Rec3dReconstruction = [];
    return;
end


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
    if ~isempty(RecTempWaitbar)
         delete(RecTempWaitbar); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp('          SIRT iterations calculated!');
         disp(' ');
    end
end
% -----Feedback-----

% Apply 3d-Weighting
if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
              
    WeightingFunction = tom_calc_weight_function ...
                                        (double([SizeX SizeY SizeZ]), double([ProjDir; Tiltangles]'), double(ObjectThickness));
    
    Rec3dReconstruction = tom_apply_weight_function(Rec3dReconstruction, WeightingFunction);
         
end
         
% Save Reconstruction
for jjj=1:10
         try
              workingdir = pwd;
              cd (ReconstructionPath);
              tom_emwrite([ReconstructionName ReconstructionExt], Rec3dReconstruction);
              cd (workingdir);
              break;
         catch
              pause(0.5);
              disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
         end 
end

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(ProtocolMode, 'on') == 1
    % ProtocolMode
    disp(' ');
    disp('          Reconstruction saved!');
    disp(' ');
end
% -----Feedback-----


% -------------------------------------------------------------------------








% -------------------------------------------------------------------------
% ReconstructionEngine -> Rec3dProject
% -------------------------------------------------------------------------
% HEAD
Rec3dProject.ProjectStatus = 'reconstructed';
    
% SIRTRESIDUALS
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerProjection = EuclidianDistancePerProjection;
Rec3dProject.SIRTRESIDUALS.EuclidianDistancePerIteration = EuclidianDistancePerIteration;
Rec3dProject.SIRTRESIDUALS.dProjVal = dProjVal;
Rec3dProject.SIRTRESIDUALS.DifVolVal = DifVolVal;
Rec3dProject.SIRTRESIDUALS.RecVolVal = RecVolVal;
    
% MESSAGE
if strcmp(FeedbackStyle, 'guistyle') == 1
    uiwait(msgbox('Reconstruction done and saved!', ...
                  'Do Reconstruction', 'warn'));
end
% -------------------------------------------------------------------------










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
% SIRTmProjProcessor
% -------------------------------------------------------------------------
function [NumOfmProjProcessed] = SIRTmProjProcessor...
(FeedbackStyle, ProtocolMode, DemoMode, ...
NumOfProj, Imdim, ...
ReconstructionName, ReconstructionPath, ReconstructionExt, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsName, ProjectionsPath, ...
Tiltangles, Tiltaxis, ProjDir, tx, ty, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, AlignTiltaxis, ...
Normalization, SmoothBorders, Taper, FILTER)


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    RecTempWaitbar = waitbar(0, 'Calculate measured projections', ...
                                                 'CreateCancelBtn', 'closereq');
    set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
    set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
    % Open DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         open('tom_Rec3dDemoModeGUI.fig'); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          Calculate measured projections');
    end
end
% -----Feedback-----

% WeightingFunction for analytical weighting
if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, 'analytical') == 1
    WeightingFunction = tom_calc_weight_function([Imdim Imdim], 'analytical');
end

% Loop over projections
for k=1:NumOfProj
    
    % Load projection
    for jjj=1:10
         try
              workingdir = pwd; 
              cd (ProjectionsPath{k}); 
              Proj = tom_emread(ProjectionsName{k});  
              cd (workingdir);
              break;
         catch
              pause(0.5);
              disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
         end 
    end
    
    % Save original header
    ProjHeaderOriginal = Proj.Header;
    
    % Convert to matrix
    Proj = Proj.Value;
    
    % PreBinning
    if (PreBinning > 0)
         Proj = tom_bin(Proj, PreBinning);
    end
    
    % Convert to double
    Proj = double(Proj);
    
    % Normalization
    if strcmp(Normalization, 'phase') == 1
         ProjDev = mean(mean(Proj));
         Proj = (Proj-ProjDev)./ProjDev;
    end
    
    % Smooth borders
    Proj = tom_smooth(Proj, SmoothBorders);
    
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(DemoMode, 'on') == 1
         % DemoMode 'Projection'
         DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
         if ~isempty(DemoModeGUI)
              DemoModeFun('axesProjection', Proj); drawnow;
         end
    end
    % -----Feedback-----
    
    
    % Move projection
    Proj = tom_move(Proj, [-fix(tx(k)) -fix(ty(k))]);
    fProj = fftshift(fft2(Proj));
    fProj = shift_in_fs(fProj, [-(tx(k)-fix(tx(k))) -(ty(k)-fix(ty(k)))]);
    Proj = real(ifft2(ifftshift(fProj)));
    clear fProj;
    
    % Rotate projection / Handedness / Taper
    if strcmp(Taper, 'off') == 1
         Proj = tom_rotate(Proj, -Tiltaxis(k), 'linear', [Imdim./2 Imdim./2]);
    elseif strcmp(Taper, 'on') == 1
         Proj = tom_rotate(Proj, -Tiltaxis(k), 'linear', [Imdim./2 Imdim./2], 'taper');
    end
    
    % Weight projection
    if strcmp(ApplyWeighting, 'on') == 1
         
         if strcmp(WeightingMethod, 'exact') == 1
              WeightingFunction = tom_calc_weight_function...
              (double([Imdim Imdim]), double([ProjDir; Tiltangles]'), ...
              double(ObjectThickness), double([ProjDir(k) Tiltangles(k)]));
         end
         
         if strcmp(WeightingMethod, 'exact') == 1 || ...
            strcmp(WeightingMethod, 'analytical') == 1
              Proj = tom_apply_weight_function(Proj, WeightingFunction);
         end
         
         % ProtocolMode
         if strcmp(ProtocolMode, 'on') == 1
              if strcmp(WeightingMethod, 'exact') == 1 || ...
                 strcmp(WeightingMethod, 'analytical') == 1
                   disp(['                         Weighted']);
              end
         end
         
    end
    
    % Filter projection
    Proj = tom_apply_filter(Proj, FILTER);
    
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(DemoMode, 'on') == 1
         % DemoMode 'Processed'
         DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
         if ~isempty(DemoModeGUI)
              DemoModeFun('axesProcessed', Proj); drawnow;
         end
    end
    % -----Feedback-----
    
    
    % PostBinning
    if (PostBinning > 0)
         Proj = tom_bin(Proj, PostBinning);
    end
    
    % Convert to single
    Proj = single(Proj);
    
    % Add and aktualize header
    Proj = tom_emheader(Proj);
    Proj.Header = ProjHeaderOriginal;
    Proj.Header.Tiltangle = Tiltangles(k); 
    Proj.Header.Size = [size(Proj.Value,1),size(Proj.Value,2),1];
    Proj.Header.Magic(4) = 5;
    Proj.Header.Tiltaxis = Tiltaxis(k);
    
    % Save SIRT measured projections
    for jjj=1:10
         try
              workingdir = pwd;
              cd (ReconstructionPath);
              tom_emwrite(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em'], Proj);
              cd (workingdir);
              break;
         catch
              pause(0.5);
              disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
         end 
    end
    
    % Clear memory
    clear Proj;
    
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1
         % ProtocolMode
         if strcmp(ProtocolMode, 'on') == 1
              disp(['                    Processed projection: ' ProjectionsName{k}]);
         end
         % Check RecTempWaitbar
         RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         % Cancel button clicked
         if isempty(RecTempWaitbar)
              % DemoMode
              if strcmp(DemoMode, 'on') == 1
                   DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
                   if ~isempty(DemoModeGUI)
                        delete(DemoModeGUI); drawnow;
                   end
              end
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Process cancelled!');
                   disp(' ');
              end
              NumOfmProjProcessed = k;
              return;
         end
         % Aktualize RecTempWorkbar
         if ~isempty(RecTempWaitbar) 
              waitbar(k./NumOfProj, RecTempWaitbar(1)); drawnow;
         end
    end
    % -----Feedback-----


end


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
    if ~isempty(RecTempWaitbar)
         delete(RecTempWaitbar); drawnow;
    end
    % Delete DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
         if ~isempty(DemoModeGUI)
              delete(DemoModeGUI); drawnow;
         end
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp('          All measured projections calculated!');
         disp(' ');
    end
end
% -----Feedback-----


% Check SIRT Projections Preprocessor
NumOfmProjProcessed = k;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% SIRTCalculateStartVolume
% -------------------------------------------------------------------------
function [Rec3dReconstruction, NumOfProjBackprojected] = SIRTCalculateStartVolume...
(FeedbackStyle, ProtocolMode, DemoMode, ...
NumOfProj, Imdim, ...
ReconstructionName, ReconstructionPath, ReconstructionExt, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsName, ProjectionsPath, ...
Tiltangles, Tiltaxis, ProjDir, tx, ty, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, AlignTiltaxis, ...
Normalization, SmoothBorders, FILTER, SizeX, SizeY, SizeZ)


% Initialize Volume
Rec3dReconstruction = zeros(SizeX, SizeY, SizeZ, 'single');


% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    RecTempWaitbar = waitbar(0, 'Calculate start volume', ...
                                                 'CreateCancelBtn', 'closereq');
    set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
    set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
    % Open DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         open('tom_Rec3dDemoModeGUI.fig'); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          Calculate start volume');
    end
end
% -----Feedback-----

% Backprojection loop
for k=1:NumOfProj
    
    % Load processed projection
    for jjj=1:10
         try
              workingdir = pwd;
              cd (ReconstructionPath);
              Proj = tom_emread(['SIRT_mProj_' ReconstructionName '_' num2str(k) '.em']);
              cd (workingdir);
              break;
         catch
              pause(0.5);
              disp(['I/O Error !  ' datestr(now) '  try: ' num2str(jjj)]);
         end 
    end
    
    % Convert to matrix
    Proj = Proj.Value;
    
    % Convert to single
    Proj = single(Proj);
    
    % Backprojection / Align Tiltaxis
    tom_backproj3d(Rec3dReconstruction, Proj, ProjDir(k), Tiltangles(k), [0 0 0]);
    
    
    % -----Feedback-----
    % DemoMode 'Volume'
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(DemoMode, 'on') == 1
         ProjVolxy(:,:) = Rec3dReconstruction(1:SizeX, 1:SizeY, SizeZ/2);
         ProjVolyz(:,:) = Rec3dReconstruction(SizeX/2, 1:SizeY, 1:SizeZ);
         ProjVolxz(:,:) = Rec3dReconstruction(1:SizeX, SizeY/2, 1:SizeZ);
         DemoModeFun('axesVolxy', ProjVolxy); drawnow; clear ProjVolxy;
         DemoModeFun('axesVolyz', ProjVolyz); drawnow; clear ProjVolyz;
         DemoModeFun('axesVolxz', ProjVolxz); drawnow; clear ProjVolxz;
    end
    % -----Feedback-----
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1
         % ProtocolMode
         if strcmp(ProtocolMode, 'on') == 1
              disp(['                    Processed projection: ' ProjectionsName{k}]);
         end
         % Check RecTempWaitbar
         RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         % Cancel button clicked
         if isempty(RecTempWaitbar)
              % DemoMode
              if strcmp(DemoMode, 'on') == 1
                   DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
                   if ~isempty(DemoModeGUI)
                        delete(DemoModeGUI); drawnow;
                   end
              end
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Process cancelled!');
                   disp(' ');
              end
              Rec3dReconstruction = [];
              NumOfProjBackprojected = k;
              return;
         end
         % Aktualize RecTempWorkbar
         if ~isempty(RecTempWaitbar) 
              waitbar(k./NumOfProj, RecTempWaitbar); drawnow;
         end
    end
    % -----Feedback-----

end

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
    if ~isempty(RecTempWaitbar)
         delete(RecTempWaitbar); drawnow;
    end
    % Delete DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         DemoModeGUI = findall(0, 'Tag', 'Rec3dDemoModeGUI');
         if ~isempty(DemoModeGUI)
              delete(DemoModeGUI); drawnow;
         end
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp('          Start volume calculated!');
         disp(' ');
    end
end
% -----Feedback-----


% Check SIRT Calculate Start Volume
NumOfProjBackprojected = k;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Shift in Fourier Space
% -------------------------------------------------------------------------
function im_out=shift_in_fs(im,v)
[dimx,dimy]=size(im);
[x,y]=ndgrid( -floor(size(im,1)/2):-floor(size(im,1)/2)+(size(im,1)-1),...
-floor(size(im,2)/2):-floor(size(im,2)/2)+size(im,2)-1);
v = v./[dimx dimy];
x = v(1)*x + v(2)*y;clear y;
im_out=(im.*exp(-2*pi*i*x));
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% DemoModeFun
% -------------------------------------------------------------------------
function DemoModeFun(axesName, axesPicture)

axesHandle = findall(0, 'Tag', axesName);

if ~isempty(axesHandle)
    axes(axesHandle);
    tom_Rec3dImagesc(axesPicture, 'noinfo');
    colormap(gray);
end
% -------------------------------------------------------------------------

