function [Rec3dReconstruction, Rec3dProject] = tom_Rec3dWBPReconstructionEngine...
(Rec3dControl, Rec3dProject, RecMode, ProjNameConv, FeedbackStyle)
%TOM_REC3DWBPRECONSTRUCTIONENGINE is a module of TOM_REC3DGUI.
%
%   [Rec3dReconstruction,Rec3dProject] = tom_Rec3dWBPReconstructionEngine...
%   (Rec3dControl,Rec3dProject,RecMode,ProjNameConv,FeedbackStyle)
%
%   This programme operates a weighted-backprojection 3d-reconstruction
%   or creates temp-files alternatively.
%
%PARAMETERS
%
%  INPUT
%   Rec3dControl                   control-structure of Rec3dGUI
%   Rec3dProject                    main-structure of Rec3dGUI
%   RecMode                          'rec' or 'temp' 
%   ProjNameConv                 '001' or '1'
%   FeedbackStyle                  'guistyle' or 'funstyle'
%
%  OUTPUT
%   Rec3dReconstruction      3d-reconstruction
%   Rec3dProject                    aktualized main-structure of Rec3dGUI
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DALIGNMENTENGINE   TOM_BACKPROJ3D   TOM_CALC_WEIGHT_FUNCTION
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
WaitbarMode = Rec3dControl.WaitbarMode;
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
isoscale = Rec3dProject.PARAMETER.isoscale;
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
% RECONSTRUCTION / Filter
FILTER = Rec3dProject.FILTER;
% RECONSTRUCTION / Detail
DetailMode = Rec3dProject.DETAIL.DetailMode;
OverviewSizeZ = Rec3dProject.DETAIL.OverviewSizeZ;
OverviewPreBinning = Rec3dProject.DETAIL.OverviewPreBinning;
OverviewPostBinning = Rec3dProject.DETAIL.OverviewPostBinning;
NumberOfDetails = Rec3dProject.DETAIL.NumberOfDetails;
DetailCoordinates = Rec3dProject.DETAIL.DetailCoordinates;
Align = Rec3dProject.DETAIL.av3_alignstruct;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check ReconstructionMethod
if strcmp(ReconstructionMethod, 'WBP') ~= 1 && strcmp(RecMode, 'rec') == 1
    % ERROR
    msgbox('Reconstruction Method must be WBP!', 'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check ReconstructionName, ReconstructionPath, ReconstructionExt
if strcmp(RecMode, 'rec') == 1 && ...
   strcmp(ReconstructionName, '') == 1 && ...
   strcmp(ReconstructionPath, '') == 1 && ...
   strcmp(ReconstructionExt, '') == 1
    
    ReconstructionName = ProjectName;
    ReconstructionPath = pwd;
    ReconstructionExt = '.vol';
end

% Check ReconstructionPath
if strcmp(RecMode, 'rec') == 1 && strcmp(ReconstructionPath, '') == 1
    % ERROR
    msgbox('No Reconstruction Path!', 'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check TempFilesName, TempFilesPath, TempFilesExt
if strcmp(RecMode, 'temp') == 1 && ...
   strcmp(TempFilesName, '') == 1 && ...
   strcmp(TempFilesPath, '') == 1 && ...
   strcmp(TempFilesExt, '') == 1
    
    TempFilesName = [ProjectName '_'];
    TempFilesPath = pwd;
    TempFilesExt = '.em';
end

% Check TempFilesPath
if strcmp(RecMode, 'temp') == 1 && strcmp(TempFilesPath, '') == 1
    % ERROR
    msgbox('No TempFiles Path!', 'Create TempFiles', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check Imdim
if strcmp(RecMode, 'rec') == 1 && ...
   (SizeX > (Imdim/(2.^PreBinning)) || SizeY > (Imdim/(2.^PreBinning)))
    % ERROR
    msgbox('Reconstruction size is bigger than resulting image dimension.', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check Normalization
if strcmp(RecMode, 'rec') == 1 && strcmp(Normalization, 'off') == 1
    % ERROR
    msgbox('Normalization is set off. This is not possible for WBP reconstruction!', ...
           'Do Reconstruction', 'error');
    Rec3dReconstruction = [];
    return;
end

% Check SmoothBorders
if SmoothBorders > (Imdim/(2.^PreBinning))
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

% Check WeightingMethod
if strcmp(RecMode, 'temp') == 1
    if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
         % ERROR
         msgbox('3d-Weighting is not possible for TempFiles!', ...
                'Create TempFiles', 'error');
         Rec3dReconstruction = [];
         return;
    end
end

% Check DetailMode
if strcmp(RecMode, 'rec') == 1
    if strcmp(DetailMode, 'on') == 1 || strcmp(DetailMode, 'av3') == 1
         if isequal(NumberOfDetails, 0) == 1
              % ERROR
              msgbox('DetailMode is set on or av3 - but NumberOfDetails is zero!', ...
                     'Do Reconstruction', 'error');
              Rec3dReconstruction = [];
              return;
         end
    end
end
% -------------------------------------------------------------------------

% FeedbackStyle
if strcmp(FeedbackStyle, 'funstyle') == 1
    WaitbarMode = 'off';
    ProtocolMode = 'off';
    DemoMode = 'off';
end

% ProtocolMode
if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         if strcmp(RecMode, 'rec') == 1
              disp('Do Reconstruction');
         elseif strcmp(RecMode, 'temp') == 1
              disp('Create TempFiles');
         end
         disp(' ');
end


% -------------------------------------------------------------------------
% DetailOffset
% -------------------------------------------------------------------------
if strcmp(RecMode, 'rec') == 1
    
    % Normal Volume Reconstruction
    % DetailOffset = [0 0 0]
    if strcmp(DetailMode, 'off') == 1
         DetailOffset = [0 0 0];
    end
    
    % Multi Detail Reconstruction
    % DetailOffset = [dx1 dy1 dz1]
    %                         [dx2 dy2 dz2]
    %                         ...
    %                         [dxk dyk dzk]
    if (strcmp(DetailMode, 'on') == 1 || strcmp(DetailMode, 'av3') == 1) && ...
         NumberOfDetails >= 1
              
              for k=1:NumberOfDetails
                   
                   DetailOffsetPreXY = (-1)./(2.^(PreBinning+PostBinning)).* ...
                   ([Imdim./2 Imdim./2] - 2.^(OverviewPreBinning+OverviewPostBinning).* ...
                   [DetailCoordinates(k,1) DetailCoordinates(k,2)]);
                   
                   DetailOffsetPreZ = 1./(2.^(PreBinning+PostBinning)).* ...
                                                        2.^(OverviewPreBinning+OverviewPostBinning).* ...
                                                        (-OverviewSizeZ./2 + DetailCoordinates(k,3));
                   
                   DetailOffset(k,1) = DetailOffsetPreXY(1);
                   DetailOffset(k,2) = DetailOffsetPreXY(2);
                   DetailOffset(k,3) = DetailOffsetPreZ;
                   
              end
    
    end

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Prepare Parameter
% -------------------------------------------------------------------------
% Calculate shifts, image dimension, isoscale
txResult = tx./(2.^PreBinning);
tyResult = ty./(2.^PreBinning);
ImdimResult = Imdim./(2.^PreBinning); %Imdim is now binned projection dimension
isoscaleResult = 1./isoscale;

% AlignTiltaxis / Modify ProjDir
if strcmp(AlignTiltaxis, 'Y-axis') == 1
   ProjDirResult(1:NumOfProj) = 0;
elseif strcmp(AlignTiltaxis, 'ProjDir') == 1
    ProjDirResult = ProjDir;
end

% Handedness / Modify Tiltaxis
if isequal(Handedness, 180) == 1
    TiltaxisResult = Tiltaxis + 180;
elseif isequal(Handedness, 0) == 1
    TiltaxisResult = Tiltaxis;
end

% Prepare  ProjectionsName
if strcmp(ProjNameConv, '001') == 1
    ProjectionsNameResult = PrepareProjectionsName(NumOfProj, ProjectionsName);
elseif strcmp(ProjNameConv, '1') == 1
    ProjectionsNameResult = ProjectionsName;
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Prepare Projections
% -------------------------------------------------------------------------
% WBP Projections Preprocessor
[NumOfProjPreprocessed] = WBPProjectionsPreprocessor...
(RecMode, FeedbackStyle, WaitbarMode, ProtocolMode, DemoMode, ...
NumOfProj, ImdimResult, ...
ReconstructionName, ReconstructionPath, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsNameResult, ProjectionsPath, ...
Tiltangles, TiltaxisResult, ProjDirResult, txResult, tyResult, isoscaleResult, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, ...
Normalization, SmoothBorders, FILTER);

% Check WBP Projections Preprocessor
if isequal(NumOfProjPreprocessed, NumOfProj) ~= 1
    
    if strcmp(DelMode, 'on') == 1 && strcmp(RecMode, 'rec') == 1
         
         for k=1:NumOfProjPreprocessed
              workingdir = pwd;
              cd (ReconstructionPath);
              delete(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em']);
              cd (workingdir);
         end
              
    end
    
    if strcmp(DelMode, 'on') == 1 && strcmp(RecMode, 'temp') == 1
              
         for k=1:NumOfProjPreprocessed
              workingdir = pwd;
              cd (TempFilesPath);
              delete([TempFilesName num2str(k) TempFilesExt]);
              cd (workingdir);
         end
    
    end
    
    Rec3dReconstruction = [];
    return;

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Reconstruction
% -------------------------------------------------------------------------
if strcmp(RecMode, 'rec') == 1
    
    % Normal Volume Reconstruction
    % DetailOffset = [0 0 0]
    if strcmp(DetailMode, 'off') == 1
         
         % WBP Backprojection Operator
         [Rec3dReconstruction, NumOfProjBackprojected] = WBPBackprojectionOperator...
         (RecMode, FeedbackStyle, WaitbarMode, ProtocolMode, DemoMode, ...
         NumOfProj, ReconstructionName, ReconstructionPath, ...
         ProjectionsNameResult, Tiltangles, ProjDirResult, ...
         SizeX, SizeY, SizeZ, DetailOffset);
         
         % Check WBP Backprojection Operator
         if isequal(NumOfProjBackprojected, NumOfProj) ~= 1
         
              if strcmp(DelMode, 'on') == 1
              
                   for k=1:NumOfProjPreprocessed
                        workingdir = pwd;
                        cd (ReconstructionPath);
                        delete(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em']);
                        cd (workingdir);
                   end
              
              end
              
              Rec3dReconstruction = [];
              return;
         
         end
         
         % DelMode
         if strcmp(DelMode, 'on') == 1
              
                   for k=1:NumOfProjPreprocessed
                        workingdir = pwd;
                        cd (ReconstructionPath);
                        delete(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em']);
                        cd (workingdir);
                   end
              
         end
         
         % -----Feedback-----
         if strcmp(FeedbackStyle, 'guistyle') == 1 && ...
            strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Calculate and apply 3d-weighting');
                   disp(' ');
              end
         end
         % -----Feedback-----
         
         % Apply 3d-Weighting
         if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
              
              WeightingFunction = tom_calc_weight_function ...
              (double([SizeX SizeY SizeZ]), double([ProjDirResult; Tiltangles]'), double(ObjectThickness));
              
              Rec3dReconstruction = tom_apply_weight_function(Rec3dReconstruction, WeightingFunction);
         
         end
         
         % -----Feedback-----
         if strcmp(FeedbackStyle, 'guistyle') == 1
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Save reconstruction');
                   disp(' ');
              end
         end
         % -----Feedback-----
         
         % Save Reconstruction
         workingdir = pwd;
         cd (ReconstructionPath);
         tom_emwrite([ReconstructionName ReconstructionExt], Rec3dReconstruction);
         cd (workingdir);
    
    end
    
    % Multi Detail Reconstruction
    % DetailOffset = [dx1 dy1 dz1]
    %                         [dx2 dy2 dz2]
    %                         ...
    %                         [dxk dyk dzk]
    if (strcmp(DetailMode, 'on') == 1 || strcmp(DetailMode, 'av3') == 1) && ...
       NumberOfDetails >= 1
         
         % -----Feedback-----
         if strcmp(FeedbackStyle, 'guistyle') == 1 && ...
            strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
              % ProtocolMode
              if strcmp(ProtocolMode, 'on') == 1
                   disp('          Calculate 3d-weighting');
                   disp(' ');
              end
         end
         % -----Feedback-----
         
         % Calculate 3d-Weighting
         if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
              WeightingFunction = tom_calc_weight_function ...
              ([SizeX SizeY SizeZ], [ProjDirResult; Tiltangles]', ObjectThickness);
         end
              
         % Loop over Details
         for k=1:NumberOfDetails
              
              % WBP Backprojection Operator
              [Rec3dReconstruction, NumOfProjBackprojected] = WBPBackprojectionOperator...
              (RecMode, FeedbackStyle, WaitbarMode, ProtocolMode, DemoMode, ...
              NumOfProj, ReconstructionName, ReconstructionPath, ...
              ProjectionsNameResult, Tiltangles, ProjDirResult, ...
              SizeX, SizeY, SizeZ, [DetailOffset(k,1) DetailOffset(k,2) DetailOffset(k,3)]);
              
              % Check WBP Backprojection Operator
              if isequal(NumOfProjBackprojected, NumOfProj) ~= 1
                   
                   if strcmp(DelMode, 'on') == 1
                             
                        % Delete Reconstruction TempFiles
                        for i=1:NumOfProj
                             workingdir = pwd;
                             cd (ReconstructionPath);
                             delete(['WBP_TEMP_' ReconstructionName '_' num2str(i) '.em']);
                             cd (workingdir);
                        end
                             
                        % Delete DetailReconstructions
                        for j=1:(k-1)
                             workingdir = pwd;
                             cd (ReconstructionPath);
                             delete([ReconstructionName '_' num2str(j) ReconstructionExt]);
                             cd (workingdir);
                        end
                             
                   end
                        
                   Rec3dReconstruction = [];
                   return;
                        
              end

              % -----Feedback-----
              if strcmp(FeedbackStyle, 'guistyle') == 1 && ...
                 strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
                   % ProtocolMode
                   if strcmp(ProtocolMode, 'on') == 1
                        disp('          Apply 3d-weighting');
                        disp(' ');
                   end
              end
              % -----Feedback-----
              
              % Apply 3d-Weighting
              if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, '3d') == 1
                   Rec3dReconstruction = tom_apply_weight_function...
                   (Rec3dReconstruction, WeightingFunction);
              end
              
              % -----Feedback-----
              if strcmp(FeedbackStyle, 'guistyle') == 1
                   % ProtocolMode
                   if strcmp(ProtocolMode, 'on') == 1
                        disp(['          Save detail reconstruction ' num2str(k)]);
                        disp(' ');
                   end
              end
              % -----Feedback-----
              
              % Save DetailReconstruction
              workingdir = pwd;
              cd (ReconstructionPath);
              tom_emwrite([ReconstructionName '_' num2str(k) ReconstructionExt], Rec3dReconstruction);
              cd (workingdir);
              
              % Clear DetailReconstruction
              if isequal(NumberOfDetails, k) ~= 1
                   clear Rec3dReconstruction;
              end
              
              % Aktualize av3 structure
              if strcmp(DetailMode, 'av3') == 1
                   % Filename of DetailReconstruction             
                   Align(1,k).Filename = [ReconstructionName '_' num2str(k) ReconstructionExt];
                   % DetailCoordinates
                   Align(1,k).Tomogram.Position.X = DetailCoordinates(k,1);
                   Align(1,k).Tomogram.Position.Y = DetailCoordinates(k,2);
                   Align(1,k).Tomogram.Position.Z = DetailCoordinates(k,3);
                   % -------------------------------------------
                   % Additional information
                   % -------------------------------------------
                   Align(1,k).Tomogram.DetailCoordinates.X = DetailCoordinates(k,1);
                   Align(1,k).Tomogram.DetailCoordinates.Y = DetailCoordinates(k,2);
                   Align(1,k).Tomogram.DetailCoordinates.Z = DetailCoordinates(k,3);
                   Align(1,k).Tomogram.OverviewSizeZ = OverviewSizeZ;
                   Align(1,k).Tomogram.OverviewPreBinning = OverviewPreBinning;
                   Align(1,k).Tomogram.OverviewPostBinning = OverviewPostBinning;
                   Align(1,k).Tomogram.DetailOffset.X = DetailOffset(k,1);
                   Align(1,k).Tomogram.DetailOffset.Y = DetailOffset(k,2);
                   Align(1,k).Tomogram.DetailOffset.Z = DetailOffset(k,3);
                   Align(1,k).Tomogram.PreBinning = PreBinning;
                   Align(1,k).Tomogram.PostBinning = PostBinning;
                   Align(1,k).Tomogram.DetailReconstructionSize.X = SizeX;
                   Align(1,k).Tomogram.DetailReconstructionSize.Y = SizeY;
                   Align(1,k).Tomogram.DetailReconstructionSize.Z = SizeZ;
                   % -------------------------------------------
              end
                   
         end
         
         % DelMode
         if strcmp(DelMode, 'on') == 1
              for k=1:NumOfProj
                   workingdir = pwd;
                   cd (ReconstructionPath);
                   delete(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em']);
                   cd (workingdir);
              end
         end
              
         % Save av3 structure
         if strcmp(DetailMode, 'av3') == 1
              % Save Project
              workingdir = pwd;
              cd (ReconstructionPath);
              save(['av3_alignstruct_' ReconstructionName '.mat'], 'Align');
              cd (workingdir);
         end
         
    end
    
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ReconstructionEngine -> Rec3dProject
% -------------------------------------------------------------------------
% RecMode 'rec'
if strcmp(RecMode, 'rec') == 1
    
    % HEAD
    Rec3dProject.ProjectStatus = 'reconstructed';
    
    % MESSAGE
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(WaitbarMode, 'on') == 1
         uiwait(msgbox('Reconstruction done and saved!', ...
                                   'Do Reconstruction', 'warn'));
    end
    
end

% RecMode 'temp'
if strcmp(RecMode, 'temp') == 1
    
    % MESSAGE
    if strcmp(FeedbackStyle, 'guistyle') == 1 && strcmp(WaitbarMode, 'on') == 1 
         uiwait(msgbox('TempFiles created and saved!', ...
                                  'Create TempFiles', 'warn'));
    end
    
    % Rec3dReconstruction
    Rec3dReconstruction = [];
    
end
% -------------------------------------------------------------------------










% -------------------------------------------------------------------------
% PrepareProjectionsName
% -------------------------------------------------------------------------
function [ProjectionsNameResult] = PrepareProjectionsName(NumOfProj, ProjectionsNamePre)

for k=1:NumOfProj
    
    UnderlineIndex = findstr(ProjectionsNamePre{k}, '_');
    WhichUnderline = size(UnderlineIndex, 2);
    UnderlineIndex = UnderlineIndex(1,WhichUnderline);
    
    WordLength = size(ProjectionsNamePre{k}, 2);
    
    if strcmp(ProjectionsNamePre{k}(UnderlineIndex+1), '0') == 1 && ...
       strcmp(ProjectionsNamePre{k}(UnderlineIndex+2), '0') == 1
   
         for i=1:UnderlineIndex
              ProjectionsNameResult{k}(i) = ProjectionsNamePre{k}(i);
         end
         
         for i=(UnderlineIndex+3):WordLength
              ProjectionsNameResult{k}(i-2) = ProjectionsNamePre{k}(i);
         end
    
    elseif strcmp(ProjectionsNamePre{k}(UnderlineIndex+1), '0') == 1 && ...
       strcmp(ProjectionsNamePre{k}(UnderlineIndex+2), '0') ~= 1
         
         for i=1:UnderlineIndex
              ProjectionsNameResult{k}(i) = ProjectionsNamePre{k}(i);
         end
         
         for i=(UnderlineIndex+2):WordLength
              ProjectionsNameResult{k}(i-1) = ProjectionsNamePre{k}(i);
         end
    
    else
        
         ProjectionsNameResult{k} = ProjectionsNamePre{k};
    
    end

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% WBP Projections Preprocessor
% -------------------------------------------------------------------------
function [NumOfProjPreprocessed] = WBPProjectionsPreprocessor...
(RecMode, FeedbackStyle, WaitbarMode, ProtocolMode, DemoMode, ...
NumOfProj, ImdimResult, ...
ReconstructionName, ReconstructionPath, ...
TempFilesName, TempFilesPath, TempFilesExt, ...
ProjectionsNameResult, ProjectionsPath, ...
Tiltangles, TiltaxisResult, ProjDirResult, txResult, tyResult, isoscaleResult, PreBinning, PostBinning, ...
ApplyWeighting, WeightingMethod, ObjectThickness, ...
Normalization, SmoothBorders, FILTER)

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    if strcmp(WaitbarMode, 'on') == 1
         RecTempWaitbar = waitbar(0, 'Calculate projections', 'CreateCancelBtn', 'closereq');
         set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
         if strcmp(RecMode, 'rec') == 1
              set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
         elseif strcmp(RecMode, 'temp') == 1
              set(RecTempWaitbar, 'Name', 'Create TempFiles'); drawnow;
         end
    end
    % Open DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         open('tom_Rec3dDemoModeGUI.fig'); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          WBP Projections Preprocessor');
    end
end
% -----Feedback-----

% WeightingFunction for analytical weighting
if strcmp(ApplyWeighting, 'on') == 1 && strcmp(WeightingMethod, 'analytical') == 1
    WeightingFunction = tom_calc_weight_function([ImdimResult ImdimResult], 'analytical');
end

% Loop over projections
for k=1:NumOfProj
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(['                    Process projection: ' ProjectionsNameResult{k}]);
         disp('                    ---------------------------------------------------');
    end
    
    % Load projection
    workingdir = pwd;
    cd (ProjectionsPath{k});
    Proj = tom_emread(ProjectionsNameResult{k});     
    cd (workingdir);
    
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
    Proj = tom_move(Proj, [-fix(txResult(k)) -fix(tyResult(k))]);
    fProj = fftshift(fft2(Proj));
    fProj = shift_in_fs(fProj, [-(txResult(k)-fix(txResult(k))) -(tyResult(k)-fix(tyResult(k)))]);
    Proj = real(ifft2(ifftshift(fProj)));
    clear fProj;
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(['                         Moved']);
         disp(['                            in x-direction: ' num2str(txResult(k))]);
         disp(['                            in y-direction: ' num2str(tyResult(k))]);
    end
    
    % Rotate projection / Handedness
    Proj = tom_rotate(Proj, -TiltaxisResult(k), 'linear', [ImdimResult./2 ImdimResult./2]);
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(['                         Rotated with Tiltaxis: ' num2str(-TiltaxisResult(k))]);
    end
    
    % Weight projection
    if strcmp(ApplyWeighting, 'on') == 1
         
         if strcmp(WeightingMethod, 'exact') == 1
              WeightingFunction = tom_calc_weight_function...
              (double([ImdimResult ImdimResult]), double([ProjDirResult; Tiltangles]'), ...
              double(ObjectThickness), double([ProjDirResult(k) Tiltangles(k)]));
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
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         if isequal(FILTER.Apply, 1) == 1 || isequal(FILTER.Apply, 2) == 1
              disp(['                         Filtered']);
         end
    end
    
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
    %Proj.Header.Magic = [6 0 0 2]'; 
    Proj.Header.Magic(4) = 5;
    Proj.Header.Tiltaxis = TiltaxisResult(k);
    
    % Save ReconstructionTempFile
    if strcmp(RecMode, 'rec') == 1
         workingdir = pwd;
         cd (ReconstructionPath);
         tom_emwrite(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em'], Proj);
         cd (workingdir);
    end
    
    % Save TempFile
    if strcmp(RecMode, 'temp') == 1
         workingdir = pwd;
         cd (TempFilesPath);
         tom_emwrite([TempFilesName num2str(k) TempFilesExt], Proj);
         cd (workingdir);
    end
    
    % Clear memory
    clear Proj;
    
    % -----Feedback-----
    if strcmp(FeedbackStyle, 'guistyle') == 1
         % ProtocolMode
         if strcmp(ProtocolMode, 'on') == 1
              disp('                    ---------------------------------------------------');
              disp(['                    Processed projection: ' ProjectionsNameResult{k}]);
              disp(' ');
         end
         % Check RecTempWaitbar
         if strcmp(WaitbarMode, 'on') == 1
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
                   NumOfProjPreprocessed = k;
                   return;
              end
              % Aktualize RecTempWorkbar
              if ~isempty(RecTempWaitbar)
                   waitbar(k./NumOfProj, RecTempWaitbar); drawnow;
              end
         end
    end
    % -----Feedback-----

end

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    if strcmp(WaitbarMode, 'on') == 1
         RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         if ~isempty(RecTempWaitbar)
              delete(RecTempWaitbar); drawnow;
         end
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
         disp('          All projections calculated!');
         disp(' ');
    end
end
% -----Feedback-----

% Check WBP Projections Preprocessor
NumOfProjPreprocessed = k;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% WBP Backprojection Operator
% -------------------------------------------------------------------------
function [Rec3dReconstruction, NumOfProjBackprojected] = WBPBackprojectionOperator...
(RecMode, FeedbackStyle, WaitbarMode, ProtocolMode, DemoMode, ...
NumOfProj, ReconstructionName, ReconstructionPath, ...
ProjectionsNameResult, Tiltangles, ProjDirResult, ...
SizeX, SizeY, SizeZ, DetailOffset)

% Initialize Volume
Rec3dReconstruction = zeros(SizeX, SizeY, SizeZ, 'single');

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Initialize RecTempWaitbar
    if strcmp(WaitbarMode, 'on') == 1
         RecTempWaitbar = waitbar(0, 'Calculate backprojection', 'CreateCancelBtn', 'closereq');
         set(RecTempWaitbar, 'Tag', 'RecTempWaitbarTAG');
         if strcmp(RecMode, 'rec') == 1
              set(RecTempWaitbar, 'Name', 'Do Reconstruction'); drawnow;
         elseif strcmp(RecMode, 'temp') == 1
              set(RecTempWaitbar, 'Name', 'Create TempFiles'); drawnow;
         end
    end
    % Open DemoModeGUI
    if strcmp(DemoMode, 'on') == 1
         open('tom_Rec3dDemoModeGUI.fig'); drawnow;
    end
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(' ');
         disp('          WBP Backprojection Operator');
    end
end
% -----Feedback-----

% Backprojection loop
for k=1:NumOfProj
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(['                    Process projection: ' ProjectionsNameResult{k}]);
         disp('                    ---------------------------------------------------');
    end
    
    % Load processed projection
    workingdir = pwd;
    cd (ReconstructionPath);
    Proj = tom_emread(['WBP_TEMP_' ReconstructionName '_' num2str(k) '.em']);
    cd (workingdir);
    
    % Convert to matrix
    Proj = Proj.Value;
    
    % Convert to single
    Proj = single(Proj);
    
    % Backprojection / Align Tiltaxis
    tom_backproj3d(Rec3dReconstruction, Proj, ProjDirResult(k), Tiltangles(k), ...
                                [DetailOffset(1,1) DetailOffset(1,2) DetailOffset(1,3)]);
    
    % ProtocolMode
    if strcmp(ProtocolMode, 'on') == 1
         disp(['                         Backprojected with']);
         disp(['                            Tiltangle: ' num2str(Tiltangles(k))]);
         disp(['                            ProjDir: ' num2str(ProjDirResult(k))]);
    end
    
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
              disp('                    ---------------------------------------------------');
              disp(['                    Processed projection: ' ProjectionsNameResult{k}]);
              disp(' ');
         end
         % Check RecTempWaitbar
         if strcmp(WaitbarMode, 'on') == 1
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
    end
    % -----Feedback-----

end

% -----Feedback-----
if strcmp(FeedbackStyle, 'guistyle') == 1
    % Delete RecTempWaitbar
    if strcmp(WaitbarMode, 'on') == 1
         RecTempWaitbar = findall(0, 'Tag', 'RecTempWaitbarTAG');
         if ~isempty(RecTempWaitbar)
              delete(RecTempWaitbar); drawnow;
         end
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
         disp('          All projections backprojected!');
         disp(' ');
    end
end
% -----Feedback-----

% Check WBP Backprojection Operator
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

