function [Rec3dProject] = tom_Rec3dAlignmentEngine(Rec3dControl, Rec3dProject)
%TOM_REC3DALIGNMENTENGINE is a module of TOM_REC3DGUI.
%
%   Rec3dProject = tom_Rec3dAlignmentEngine(Rec3dControl,Rec3dProject)
%
%   It calculates reconstruction parameter for both singleaxis and dualaxis
%   reconstruction projects in the case of fiducial marker files.
%   In the case of feature tracking markerfiles it calculates
%   the reconstruction parameter for singleaxis tilting geometry.
%
%PARAMETERS
%
%  INPUT
%   Rec3dControl     control-structure of Rec3d
%   Rec3dProject     main-structure of Rec3d
%  
%  OUTPUT
%   Rec3dProject     aktualized main-structure of tom_Rec3dGUI
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DEULERANGLES   TOM_REC3DEULERANGLESRESIDUAL
%   TOM_REC3DNEWORIGIN2
%   TOM_REC3DRIGIDBODYALIGNMENT   TOM_REC3DRIGIDBODYRESIDUAL
%   TOM_REC3DALIGNMENTRESIDUAL
%
%   fiducial marker 01/01/07 ME / feature tracking 01/02/08 ME 
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
% Rec3dControl -> AlignmentEngine
% -------------------------------------------------------------------------
ProtocolMode     = Rec3dControl.ProtocolMode;
% -------------------------------------------------------------------------
% Rec3dProject -> AlignmentEngine
% -------------------------------------------------------------------------
% PROJECT
TiltingGeometry  = Rec3dProject.PROJECT.TiltingGeometry;
NumOfProj1       = Rec3dProject.PROJECT.NumOfProj1;
RefProj1         = Rec3dProject.PROJECT.RefProj1;
Tiltangles1      = Rec3dProject.PROJECT.Tiltangles1;
MarkersOnProj1   = Rec3dProject.PROJECT.MarkersOnProj1;
Markerfile1      = Rec3dProject.PROJECT.Markerfile1;
NumOfProj2       = Rec3dProject.PROJECT.NumOfProj2;
RefProj2         = Rec3dProject.PROJECT.RefProj2;
Tiltangles2      = Rec3dProject.PROJECT.Tiltangles2;
MarkersOnProj2   = Rec3dProject.PROJECT.MarkersOnProj2;
Markerfile2      = Rec3dProject.PROJECT.Markerfile2;
NumOfMarkers     = Rec3dProject.PROJECT.NumOfMarkers;
Imdim            = Rec3dProject.PROJECT.Imdim;
% ALIGNMENT
AlignmentMethod  = Rec3dProject.ALIGNMENT.AlignmentMethod;
ReferenceMarker  = Rec3dProject.ALIGNMENT.ReferenceMarker;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Which MarkerfileType
% -------------------------------------------------------------------------
MarkerfileType = questdlg('Please choose the type of your markerfile:', ...
                                             'Do Alignment', 'Fiducial marker', 'Feature tracking', ...
                                             'Cancel', 'Fiducial marker');

if strcmp(MarkerfileType, 'Cancel') == 1 || strcmp(MarkerfileType, '') == 1
    return;
end
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Fiducial marker
% -------------------------------------------------------------------------
if strcmp(MarkerfileType, 'Fiducial marker') == 1

% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check Origin
if strcmp(TiltingGeometry, 'singleaxis')
    if (Markerfile1.Value(2, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile1.Value(3, RefProj1, ReferenceMarker)) == (-1)
         % ERROR
         msgbox('Alignment Origin is not defined!', 'Do Alignment', 'error');
         return;
    end
end   
if strcmp(TiltingGeometry, 'dualaxis')
    if (Markerfile1.Value(2, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile1.Value(3, RefProj1, ReferenceMarker)) == (-1) || ...
       (Markerfile2.Value(2, RefProj2, ReferenceMarker)) == (-1) || ...
       (Markerfile2.Value(3, RefProj2, ReferenceMarker)) == (-1)
         % ERROR
         msgbox('Alignment Origin is not defined!', 'Do Alignment', 'error');
         return;
    end
end    
% -------------------------------------------------------------------------

% ProtocolMode
if strcmp(ProtocolMode, 'on')
    disp(' ');
    disp('Do Alignment');
    disp(' ');
end

% -------------------------------------------------------------------------
% ALIGNMENT
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'singleaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'singleaxis')

% -------------------------------------------------------------------------
% AlignmentMethod 'rigidbody'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'rigidbody')

% Define alignment origin of tiltseries1
Origin1 = [Markerfile1.Value(2, RefProj1, ReferenceMarker) ...
           Markerfile1.Value(3, RefProj1, ReferenceMarker) ...
           Imdim/2+1];

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Tiltaxis, translations, isoscale of tiltseries1
for k = 1:NumOfProj1
    Tiltaxis1(k) = (180/pi)*Rho1;
    tx1(k) = Markerfile1.Value(7, k, 1);
    ty1(k) = Markerfile1.Value(8, k, 1);
    isoscale1(k) = 1;
end

% Calculate ResidualMatrix and WarpAlignment of tiltseries1
[ResidualMatrix1, WarpAlignment1] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj1,NumOfProj1,NumOfMarkers,Origin1,Imdim,Tiltaxis1,Tiltangles1,m3d1,tx1,ty1);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% AlignmentMethod 'freetilt'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'freetilt')
    % ERROR
    msgbox('Alignment Method freetilt not implemented!', ...
                   'Do Alignment', 'error');
    return;
end
% -------------------------------------------------------------------------

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'dualaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'dualaxis')

% -------------------------------------------------------------------------
% AlignmentMethod 'rigidbody'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'rigidbody')

% Define alignment origin of tiltseries1
Origin1 = [Markerfile1.Value(2, RefProj1, ReferenceMarker) ...
           Markerfile1.Value(3, RefProj1, ReferenceMarker) ...
           Imdim/2+1];

% Define alignment origin of tiltseries2
Origin2 = [Markerfile2.Value(2, RefProj2, ReferenceMarker) ...
           Markerfile2.Value(3, RefProj2, ReferenceMarker) ...
           Imdim/2+1];

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries2
[Markerfile2.Value, Rho2, Sigma2, x2, y2, z2] = tom_Rec3dRigidBodyAlignment...
(Markerfile2.Value, ReferenceMarker, RefProj2, Origin2, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Combine 3d marker coordinates of tiltseries2
m3d2(1,1:NumOfMarkers) = x2(1:NumOfMarkers);
m3d2(2,1:NumOfMarkers) = y2(1:NumOfMarkers);
m3d2(3,1:NumOfMarkers) = z2(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Move tiltseries2 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d2(1,1:NumOfMarkers) = m3d2(1,1:NumOfMarkers) - Origin2(1) + (Imdim/2+1); 
m3d2(2,1:NumOfMarkers) = m3d2(2,1:NumOfMarkers) - Origin2(2) + (Imdim/2+1); 
m3d2(3,1:NumOfMarkers) = m3d2(3,1:NumOfMarkers) - Origin2(3) + (Imdim/2+1);

% Calculate 3d rotation between tiltseries1 and tiltseries2
[Psi, Theta, Phi, RotMatrix] = tom_Rec3dEulerAngles...
(m3d1, m3d2, NumOfMarkers, ProtocolMode);

% Calculate new alignment origin of tiltseries2
Origin2 = tom_Rec3dNewOrigin2...
(Origin1, Origin2, Imdim, RotMatrix, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries1
[Markerfile1.Value, Rho1, Sigma1, x1, y1, z1] = tom_Rec3dRigidBodyAlignment...
(Markerfile1.Value, ReferenceMarker, RefProj1, Origin1, Imdim, ProtocolMode);

% Calculate 'rigidbody' alignment of tiltseries2
[Markerfile2.Value, Rho2, Sigma2, x2, y2, z2] = tom_Rec3dRigidBodyAlignment...
(Markerfile2.Value, ReferenceMarker, RefProj2, Origin2, Imdim, ProtocolMode);

% Combine 3d marker coordinates of tiltseries1
m3d1(1,1:NumOfMarkers) = x1(1:NumOfMarkers);
m3d1(2,1:NumOfMarkers) = y1(1:NumOfMarkers);
m3d1(3,1:NumOfMarkers) = z1(1:NumOfMarkers);

% Combine 3d marker coordinates of tiltseries2
m3d2(1,1:NumOfMarkers) = x2(1:NumOfMarkers);
m3d2(2,1:NumOfMarkers) = y2(1:NumOfMarkers);
m3d2(3,1:NumOfMarkers) = z2(1:NumOfMarkers);

% Move tiltseries1 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d1(1,1:NumOfMarkers) = m3d1(1,1:NumOfMarkers) - Origin1(1) + (Imdim/2+1); 
m3d1(2,1:NumOfMarkers) = m3d1(2,1:NumOfMarkers) - Origin1(2) + (Imdim/2+1); 
m3d1(3,1:NumOfMarkers) = m3d1(3,1:NumOfMarkers) - Origin1(3) + (Imdim/2+1);

% Move tiltseries2 3d marker coordinates of ReferenceMarker to [0 0 0]
m3d2(1,1:NumOfMarkers) = m3d2(1,1:NumOfMarkers) - Origin2(1) + (Imdim/2+1); 
m3d2(2,1:NumOfMarkers) = m3d2(2,1:NumOfMarkers) - Origin2(2) + (Imdim/2+1); 
m3d2(3,1:NumOfMarkers) = m3d2(3,1:NumOfMarkers) - Origin2(3) + (Imdim/2+1);

% Tiltaxis, translations, isoscale of tiltseries1
for k = 1:NumOfProj1
    Tiltaxis1(k) = (180/pi)*Rho1;
    tx1(k) = Markerfile1.Value(7, k, 1);
    ty1(k) = Markerfile1.Value(8, k, 1);
    isoscale1(k) = 1;
end

% Tiltaxis, translations, isoscale of tiltseries2
for k = 1:NumOfProj2
    Tiltaxis2(k) = (180/pi)*Rho2;
    tx2(k) = Markerfile2.Value(7, k, 1);
    ty2(k) = Markerfile2.Value(8, k, 1);
    isoscale2(k) = 1;
end

% Calculate ResidualMatrix and WarpAlignment of tiltseries1
[ResidualMatrix1, WarpAlignment1] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj1,NumOfProj1,NumOfMarkers,Origin1,Imdim,Tiltaxis1,Tiltangles1,m3d1,tx1,ty1);

% Calculate ResidualMatrix and WarpAlignment of tiltseries2
[ResidualMatrix2, WarpAlignment2] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj2,NumOfProj2,NumOfMarkers,Origin2,Imdim,Tiltaxis2,Tiltangles2,m3d2,tx2,ty2);

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% AlignmentMethod 'freetilt'
% -------------------------------------------------------------------------
if strcmp(AlignmentMethod, 'freetilt')
    % ERROR
    msgbox('AlignmentMethod freetilt not implemented!', ...
                   'Do Alignment', 'error');
    return;
end
% -------------------------------------------------------------------------

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% RECONSTRUCTION / Parameter
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'singleaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'singleaxis')
    
    % Reconstruction parameter of Tiltseries1
    for k=1:NumOfProj1
    
    % Tiltangles
    Tiltangles(k) = Tiltangles1(k);
    % Tiltaxis
    Tiltaxis(k) = Tiltaxis1(k) + 90;
    % ProjDir
    ProjDir(k) =  Tiltaxis1(k) + 90;
    % tx
    tx(k) = tx1(k);
    % ty
    ty(k) = ty1(k);
    % isoscale
    isoscale(k) = isoscale1(k);
    
    end
     
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% TiltingGeometry 'dualaxis'
% -------------------------------------------------------------------------
if strcmp(TiltingGeometry, 'dualaxis')
    
    % Reconstruction parameter of Tiltseries1
    for k=1:NumOfProj1
    
    % Tiltangles
    Tiltangles(k) = Tiltangles1(k);
    % Tiltaxis
    Tiltaxis(k) = Tiltaxis1(k) + 90;
    % ProjDir
    ProjDir(k) =  Tiltaxis1(k) + 90;
    % tx
    tx(k) = tx1(k);
    % ty
    ty(k) = ty1(k);
    % isoscale
    isoscale(k) = isoscale1(k);

    end
    
    % Reconstruction parameter of Tiltseries2
    for k=1:NumOfProj2
    
    % Tiltangles
    a = cos((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k)) - ...
        sin((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sqrt(1-a*a);

    Tiltangles(k+NumOfProj1) = (180/pi)*atan2(b,a);
    
    % Tiltaxis
    a = cos((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k)) + ...
        sin((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sin((pi/180)*Theta)*sin((pi/180)*(-Tiltaxis2(k)+Psi));

    Tiltaxis(k+NumOfProj1) = (180/pi)*atan2(b,a) + Tiltaxis2(k) + 90;
    
    % ProjDir
    a = sin((pi/180)*Theta)*cos((pi/180)*Tiltangles2(k)) + ...
        cos((pi/180)*Theta)*sin((pi/180)*Tiltangles2(k))*cos((pi/180)*(-Tiltaxis2(k)+Psi));

    b = sin((pi/180)*Tiltangles2(k))*sin((pi/180)*(-Tiltaxis2(k)+Psi));

    ProjDir(k+NumOfProj1) = -(180/pi)*(atan2(b,a) + (pi/180)*Phi) + 90;
    
    % tx
    tx(k+NumOfProj1) = tx2(k);
    % ty
    ty(k+NumOfProj1) = ty2(k);
    % isoscale
    isoscale(k+NumOfProj1) = isoscale2(k);

    end

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ALIGNMENT / Residuals
% -------------------------------------------------------------------------
% Calculate AlignmentResidual of tiltseries1
[AveragePerProjection1, AveragePerMarker1] = ...
tom_Rec3dAlignmentResidual(ResidualMatrix1, NumOfMarkers, NumOfProj1);

% Calculate AlignmentResidual of tiltseries2
if strcmp(TiltingGeometry, 'dualaxis')
    [AveragePerProjection2, AveragePerMarker2] = ...
    tom_Rec3dAlignmentResidual(ResidualMatrix2, NumOfMarkers, NumOfProj2);
end

% Calculate EulerAnglesResidual
if strcmp(TiltingGeometry, 'dualaxis')
    [EulerAnglesResidual, ...
     MaximumResidual, ...
     ResidualSpheres, ...
     AverageResidualSphere] = tom_Rec3dEulerAnglesResidual...
    (m3d1, m3d2, NumOfMarkers, RotMatrix, ProtocolMode);
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Check
% -------------------------------------------------------------------------
% Check Translations
if strcmp(TiltingGeometry, 'singleaxis')
    for k=1:NumOfProj1
         if tx1(k) > Imdim || ty1(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end
    end
end
if strcmp(TiltingGeometry, 'dualaxis')
    for k=1:NumOfProj1
         if tx1(k) > Imdim || ty1(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end  
    end
    for k=1:NumOfProj2
         if tx2(k) > Imdim || ty2(k) > Imdim
         % ERROR
         msgbox('Alignment false! Click more marker or delete projections!', ...
                'Do Alignment', 'error');
         return;
         end
    end
end    
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% AlignmentEngine -> Rec3dProject
% -------------------------------------------------------------------------
% HEAD
Rec3dProject.ProjectStatus = 'aligned';
% ALIGNMENT
Rec3dProject.ALIGNMENT.Origin1 = Origin1;
Rec3dProject.ALIGNMENT.m3d1 = m3d1;
Rec3dProject.ALIGNMENT.Tiltaxis1 = Tiltaxis1;
Rec3dProject.ALIGNMENT.tx1 = tx1;
Rec3dProject.ALIGNMENT.ty1 = ty1;
Rec3dProject.ALIGNMENT.isoscale1 = isoscale1;
Rec3dProject.ALIGNMENT.WarpAlignment1 = WarpAlignment1;
if strcmp(TiltingGeometry, 'dualaxis')
    Rec3dProject.ALIGNMENT.Origin2 = Origin2;
    Rec3dProject.ALIGNMENT.m3d2 = m3d2;
    Rec3dProject.ALIGNMENT.Tiltaxis2 = Tiltaxis2;
    Rec3dProject.ALIGNMENT.tx2 = tx2;
    Rec3dProject.ALIGNMENT.ty2 = ty2;
    Rec3dProject.ALIGNMENT.isoscale2 = isoscale2;
    Rec3dProject.ALIGNMENT.WarpAlignment2 = WarpAlignment2;
    Rec3dProject.ALIGNMENT.RotMatrix = RotMatrix;
    Rec3dProject.ALIGNMENT.Psi = Psi;
    Rec3dProject.ALIGNMENT.Theta = Theta;
    Rec3dProject.ALIGNMENT.Phi = Phi;
end
% ALGRESIDUALS
Rec3dProject.ALGRESIDUALS.ResidualMatrix1 = ResidualMatrix1;
Rec3dProject.ALGRESIDUALS.Sigma1 = Sigma1;
Rec3dProject.ALGRESIDUALS.AveragePerProjection1 = AveragePerProjection1;
Rec3dProject.ALGRESIDUALS.AveragePerMarker1 = AveragePerMarker1;
if strcmp(TiltingGeometry, 'dualaxis')
    Rec3dProject.ALGRESIDUALS.ResidualMatrix2 = ResidualMatrix2;
    Rec3dProject.ALGRESIDUALS.Sigma2 = Sigma2;
    Rec3dProject.ALGRESIDUALS.AveragePerProjection2 = AveragePerProjection2;
    Rec3dProject.ALGRESIDUALS.AveragePerMarker2 = AveragePerMarker2;
    Rec3dProject.ALGRESIDUALS.EulerAnglesResidual = EulerAnglesResidual;
    Rec3dProject.ALGRESIDUALS.MaximumResidual = MaximumResidual;
    Rec3dProject.ALGRESIDUALS.ResidualSpheres = ResidualSpheres;
    Rec3dProject.ALGRESIDUALS.AverageResidualSphere = AverageResidualSphere;
end
% PARAMETER
Rec3dProject.PARAMETER.Tiltangles = Tiltangles;
Rec3dProject.PARAMETER.Tiltaxis = Tiltaxis;
Rec3dProject.PARAMETER.ProjDir = ProjDir;
Rec3dProject.PARAMETER.tx = tx;
Rec3dProject.PARAMETER.ty = ty;
Rec3dProject.PARAMETER.isoscale = isoscale;
% -------------------------------------------------------------------------

end
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Feature tracking
% -------------------------------------------------------------------------
if strcmp(MarkerfileType, 'Feature tracking') == 1

% Check TiltingGeometry
if strcmp(TiltingGeometry, 'dualaxis') == 1
    % ERROR
    msgbox('Feature tracking not implemented for dual axis tilting geometry yet!', ...
           'Do Alignment', 'error');
    return;
end

% Prepare markerfile
mfpre = Markerfile1.Value;
mf = mfpre(1:3,:,:);

% Add Tiltangles
[mf] = AddTiltangles2MarkerfileNEW(mf,Tiltangles1);

% Analyse FeatureChainLength
[FeatureChainLength] = AnalyseFeatureChainLengthNEW(mf);

% Filter Feature ChainLength
[mf] = FilterFeatureLengthNEW(mf,FeatureChainLength,round(NumOfProj1/2+3));%!!!!!!!!!!!!!
NumberOfFeaturesFullyTracked = size(mf,3); %!!!!!!!!!!!!!

% Which ProcessStep
ProcessStep = questdlg('Please choose process step:', ...
                       'Do Alignment', 'Fit Tiltaxis', 'Filter Features', 'Cancel', 'Fit Tiltaxis');

if strcmp(ProcessStep, 'Cancel') == 1 || strcmp(ProcessStep, '') == 1
    return;
end

% Fit Tiltaxis
if strcmp(ProcessStep, 'Fit Tiltaxis') == 1

    % Cluster Tiltaxis
    [TiltaxisPre] = ClusterTiltaxisNEW(mf);
    
    % Assignin variables to workspace
    assignin('base','TiltaxisPre',TiltaxisPre);
    
    % Calculate and display feature statistics
    disp(' ');
    disp('--------------------------------------');
    disp('Feature tracking alignment --- Fit Tiltaxis');
    disp('--------------------------------------');
    disp(['Mean tiltaxis orientation: ' num2str(mean(TiltaxisPre))]);
    disp(['Median tiltaxis orientation: ' num2str(median(TiltaxisPre))]);
    disp(['Variance tiltaxis orientation: ' num2str(var(TiltaxisPre))]);
    disp('--------------------------------------');
    
    % Display Histogram
    figure(1);hist(TiltaxisPre,NumberOfFeaturesFullyTracked);
    
    return;
    
end

% Filter Features
if strcmp(ProcessStep, 'Filter Features') == 1
    
    % Input Tiltaxis orientation
    prompt = {'Center','Width'};
    dlg_title = 'Enter Filter Features Values:';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines);
    if isempty(answer) == 1
        return;
    end
    TiltaxisCenter = str2num(answer{1});
    TiltaxisWidth = str2num(answer{2});
    
    % Cluster Tiltaxis
    [TiltaxisPre] = ClusterTiltaxisNEW(mf);
    
    % Filter Tiltaxis
    [mf] = FilterTiltaxisGaussNEW(mf,TiltaxisPre,TiltaxisCenter-TiltaxisWidth,TiltaxisCenter+TiltaxisWidth);
    NumberOfFeaturesFullyTrackedFiltered = size(mf,3);
    
    % Display
    disp(' ');
    disp('--------------------------------------');
    disp('Feature tracking alignment --- Filter Features');
    disp('--------------------------------------');
    disp(['Number of features fully tracked: ' num2str(NumberOfFeaturesFullyTracked)]);
    disp(['Number of features filtered: ' num2str(NumberOfFeaturesFullyTrackedFiltered)]);
    disp('--------------------------------------');
    disp(' ');
    
    % Calculate Alignment
    [Tiltaxis,TiltaxisMean,gmm,tx,txMean,txgmm,ty,tyMean,tygmm] = ClusterShiftsNEW(mf,Imdim);
    
    % -------------------------------------------------------------------------
    % AlignmentEngine -> Rec3dProject
    % -------------------------------------------------------------------------
    % HEAD
    Rec3dProject.ProjectStatus = 'aligned';
    % PARAMETER
    Rec3dProject.PARAMETER.Tiltangles = Tiltangles1;
    Rec3dProject.PARAMETER.Tiltaxis = TiltaxisMean;
    Rec3dProject.PARAMETER.ProjDir = TiltaxisMean;
    Rec3dProject.PARAMETER.tx = txgmm;%!!!!!!!!!!!!
    Rec3dProject.PARAMETER.ty = tygmm;%!!!!!!!!!!!!
    Rec3dProject.PARAMETER.isoscale = ones(1,NumOfProj1);
    % -------------------------------------------------------------------------
    
end

end
% -------------------------------------------------------------------------

% MESSAGE
uiwait(msgbox('Alignment done!', ...
              'Do Alignment', 'warn'));

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Feature tracking subfunctions
% -------------------------------------------------------------------------
function [mf_new] = AddTiltangles2MarkerfileNEW(mf,Tiltangles)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(mf,3);

% Add Tiltangles
mf_new = mf;
for k=1:NumOfFeatures
         mf_new(1,:,k) = Tiltangles(1,:);
end


function [FeatureChainLength] = AnalyseFeatureChainLengthNEW(mf)
%%%%%%%%%%
%%%%%%%%%%

% NumOfProj
NumOfProj = size(mf,2);

% NumOfFeatures
NumOfFeatures = size(mf,3);

% Build mask not measured data points
NotMeasuredMask = (mf == -1);

% Measure feature chain length
NumOfNotMeasuredDataPoints = sum(NotMeasuredMask,2);
NumOfNotMeasuredDataPoints = squeeze(NumOfNotMeasuredDataPoints);
NumOfNotMeasuredDataPointsRed = NumOfNotMeasuredDataPoints(2,1:NumOfFeatures);
FeatureChainLength = NumOfProj - NumOfNotMeasuredDataPointsRed;

% Display
Mask = (FeatureChainLength == NumOfProj);


function [mf_a_filt] = FilterTiltaxisGaussNEW(mf_a,Tiltaxis_a,TiltaxisLow_a,TiltaxisHigh_a)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(Tiltaxis_a,2);

% Mask_a
Mask_a = ((Tiltaxis_a >= TiltaxisLow_a) & (Tiltaxis_a <= TiltaxisHigh_a));

zaehler = 1;
for k=1:NumOfFeatures
        
        if (Mask_a(k) == 1)
            
            mf_a_filt(1:3,:,zaehler) = mf_a(1:3,:,k);
            
            zaehler = zaehler + 1;
        
        end

end


function [Origin,Tiltaxis,Sigma,tx,ty] = GetAlignmentParameterNEW(Alignment)
%%%%%%%%%%
%%%%%%%%%%

% NumberOfFeatures
NumberOfFeatures = size(Alignment,2);

% Get Alignment Loop
for k=1:NumberOfFeatures
    
    % Origin
    Origin(k,1:3) = Alignment{k}.Origin;
    
    % Tiltaxis
    Tiltaxis(k) = Alignment{k}.Tiltaxis;
    
    % Sigma
    Sigma(k) = Alignment{k}.Sigma;
    
    % tx
    tx(k,:) = Alignment{k}.tx;
    
    % ty
    ty(k,:) = Alignment{k}.ty;
    
end


function [mf_filtered] = FilterFeatureLengthNEW(mf,FeatureChainLength,MinChainLength)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(mf,3);

zaehler = 1;    
for j=1:NumOfFeatures
        
        if FeatureChainLength(j) >= MinChainLength
            
            mf_filtered(1:3,:,zaehler) = mf(1:3,:,j);
            
            zaehler = zaehler + 1;
        
        end

end


function [Tiltaxis] = ClusterTiltaxisNEW(mf)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(mf,3);

% NumOfProj
NumOfProj = size(mf,2);

% Allocate Tiltaxis
Tiltaxis = zeros(1,NumOfFeatures);

% Waitbar
h = waitbar(0,'Calculating tiltaxis orientation...');

% Calculation of Tiltaxis
for k=1:NumOfFeatures
    
    % Calculate means of difference vectors
    meanx = zeros(NumOfFeatures,1);
    meany = zeros(NumOfFeatures,1);
    numomdp = zeros(NumOfFeatures,1); %numomdp - number of measured data points
    
    for j=1:NumOfFeatures
         
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   meanx(j) = meanx(j) + mf(2,i,j) - mf(2,i,k);
                   meany(j) = meany(j) + mf(3,i,j) - mf(3,i,k);
                   numomdp(j) = numomdp(j) + 1;
              end
         
         end
         
    end
    
    meanx = meanx./numomdp;
    meany = meany./numomdp;
    
    % Calculate sums
    sumxx=0;
    sumyy=0;
    sumxy=0;
    
    for j=1:NumOfFeatures
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   sumxx = sumxx + (mf(2,i,j) - mf(2,i,k) - meanx(j)).^2;
                   sumyy = sumyy + (mf(3,i,j) - mf(3,i,k) - meany(j)).^2;
                   sumxy = sumxy + (mf(2,i,j) - mf(2,i,k) - meanx(j)).*...
                                   (mf(3,i,j) - mf(3,i,k) - meany(j));
              end
        
         end
    end

    % Calculate Tiltaxis with linear regression
    psi = (0.5).*atan((2.*sumxy)./(sumxx-sumyy));
    if (sumxx > sumyy)
         psi = psi - (0.5).*pi.*sign(psi);
    end
    Tiltaxis(k) = psi.*(180./pi);
    
    % Aktualize waitbar
    waitbar(k/NumOfFeatures);
    
end
close(h);


function [gmm] = ClusterMarkerModelsNEW(mf,Binning,Imdim)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(mf,3);

% NumOfProj
NumOfProj = size(mf,2);

% RefProj
[notneeded,RefProj] = min(abs(mf(1,:,1)));

% Tiltangles
theta = (pi./180).*mf(1,:,1);

% Initialize GlobalMarkerModel
globalMarkerModel_sum = zeros(3,NumOfFeatures);

% Loop over all reference features
for k=1:NumOfFeatures
    
    % ---------------------------------------------------------------------
    % Calculation of k'th Tiltaxis
    % ---------------------------------------------------------------------
    
    % Calculate means of difference vectors
    meanx = zeros(NumOfFeatures,1);
    meany = zeros(NumOfFeatures,1);
    numomdp = zeros(NumOfFeatures,1);
    
    for j=1:NumOfFeatures
         
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   meanx(j) = meanx(j) + mf(2,i,j) - mf(2,i,k);
                   meany(j) = meany(j) + mf(3,i,j) - mf(3,i,k);
                   numomdp(j) = numomdp(j) + 1;
              end
         
         end
         
    end
    
    meanx = meanx./numomdp;
    meany = meany./numomdp;
    
    % Calculate sums
    sumxx=0;
    sumyy=0;
    sumxy=0;
    
    for j=1:NumOfFeatures
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   sumxx = sumxx + (mf(2,i,j) - mf(2,i,k) - meanx(j)).^2;
                   sumyy = sumyy + (mf(3,i,j) - mf(3,i,k) - meany(j)).^2;
                   sumxy = sumxy + (mf(2,i,j) - mf(2,i,k) - meanx(j)).*...
                                   (mf(3,i,j) - mf(3,i,k) - meany(j));
              end
        
         end
    end

    % Calculate Tiltaxis with linear regression
    psi = (0.5).*atan((2.*sumxy)./(sumxx-sumyy));
    if (sumxx > sumyy)
         psi = psi - (0.5).*pi.*sign(psi);
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % Calculation of k'th marker model
    % ---------------------------------------------------------------------
    
    m3d = zeros(3,NumOfFeatures);
    numomdp = zeros(NumOfFeatures,1);
    
    for j=1:NumOfFeatures
    
         sumxx=0;
         sumyy=0;
         sumxy=0;
         sumyx=0;
         
         alpha_a=0;
         alpha_b=0;
         
         P = zeros(3,3);
         temp = zeros(3);
         
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   
                   sumxx = sumxx + (mf(2,i,j) - mf(2,i,k)).*cos(theta(i));
                   sumyy = sumyy + (mf(3,i,j) - mf(3,i,k)).*cos(theta(i));
                   sumxy = sumxy + (mf(2,i,j) - mf(2,i,k)).*sin(theta(i));
                   sumyx = sumyx + (mf(3,i,j) - mf(3,i,k)).*sin(theta(i));
                   
                   alpha_a = alpha_a + sin(theta(i)).^2;
                   alpha_b = alpha_b + sin(theta(i)).*cos(theta(i));
                   
                   numomdp(j) = numomdp(j) + 1;
                   
              end
              
         end
         
         P = [numomdp(j)-alpha_a.*sin(psi).^2   alpha_a.*cos(psi).*sin(psi)   alpha_b.*sin(psi);...
                alpha_a.*cos(psi).*sin(psi)   numomdp(j)-alpha_a.*cos(psi).^2   -alpha_b.*cos(psi);...
                    alpha_b.*sin(psi)               -alpha_b.*cos(psi)                alpha_a];
         
         temp(1) = (sumxx.*sin(psi)-sumyy.*cos(psi)).*sin(psi) + (cos(psi).*meanx(j)+sin(psi).*meany(j)).*cos(psi).*numomdp(j);
         temp(2) = -(sumxx.*sin(psi)-sumyy.*cos(psi)).*cos(psi) + (cos(psi).*meanx(j)+sin(psi).*meany(j)).*sin(psi)*numomdp(j);
         temp(3) = sumxy.*sin(psi) - sumyx.*cos(psi);
         
         dt = det(P);
         
         if (dt ~= 0)
              
              P_t = P;
              P_t(1,1) = temp(1);
              P_t(2,1) = temp(2);
              P_t(3,1) = temp(3);
              m3d(1,j) = det(P_t)./dt + mf(2,RefProj,k) - (Imdim./2+1);
              
              P_t = P;
              P_t(1,2) = temp(1);
              P_t(2,2) = temp(2);
              P_t(3,2) = temp(3);
              m3d(2,j) = det(P_t)./dt + mf(3,RefProj,k) - (Imdim./2+1);
              
              P_t = P;
              P_t(1,3) = temp(1);
              P_t(2,3) = temp(2);
              P_t(3,3) = temp(3);
              m3d(3,j) = det(P_t)./dt + (Imdim./2+1) - (Imdim./2+1);
                  
         else
              
              m3d(1,j) = 1000000;
              m3d(2,j) = 1000000;
              m3d(3,j) = 1000000;
              
              warning(['Feature ' num2str(j) ' : undefined! det = 0! Click more!']);
         end
    
    end
    % ---------------------------------------------------------------------
    
    % Save k'th marker model
    save(['m3d_' num2str(k) '.mat'],'m3d');
    
    % Build global marker model
    for j=1:NumOfFeatures
    
         % Move tiltseries a 3d marker coordinates of ReferenceFeature to [0 0 0]
         m3d(1,j) = m3d(1,j) - mf(2,RefProj,k) + mf(2,RefProj,k) + (Imdim./2+1);
         m3d(2,j) = m3d(2,j) - mf(3,RefProj,k) + mf(3,RefProj,k) + (Imdim./2+1);
         m3d(3,j) = m3d(3,j) - (Imdim./2+1) + (Imdim./2+1) + (Imdim./2+1);
         
    end
    
    % Sum GlobalMarkerModel
    globalMarkerModel_sum = globalMarkerModel_sum + m3d;
    
    disp(k);
    
end

% GlobalMarkerModel
gmm = globalMarkerModel_sum./NumOfFeatures;

% Binning
gmm(1,:) = gmm(1,:)./(2.^Binning);
gmm(2,:) = gmm(2,:)./(2.^Binning);
gmm(3,:) = gmm(3,:)./(2.^Binning);


function [Tiltaxis,TiltaxisMean,gmm,tx,txMean,txgmm,ty,tyMean,tygmm] = ClusterShiftsNEW(mf,Imdim)
%%%%%%%%%%
%%%%%%%%%%

% NumOfFeatures
NumOfFeatures = size(mf,3);

% NumOfProj
NumOfProj = size(mf,2);

% RefProj
[notneeded,RefProj] = min(abs(mf(1,:,1)));

% Tiltangles
theta = (pi./180).*mf(1,:,1);

% Initialize Tiltaxis
Tiltaxis = zeros(1,NumOfFeatures);

% Initialize GlobalMarkerModel
globalMarkerModel_sum = zeros(3,NumOfFeatures);

% Initialize Shifts
tx = zeros(NumOfFeatures,NumOfProj);
ty = zeros(NumOfFeatures,NumOfProj);

% Waitbar
h = waitbar(0,'Calculating statistical alignment...');

% Loop over all reference features
for k=1:NumOfFeatures
    
    % ---------------------------------------------------------------------
    % Calculation of k'th Tiltaxis
    % ---------------------------------------------------------------------
    % Calculate means of difference vectors
    meanx = zeros(NumOfFeatures,1);
    meany = zeros(NumOfFeatures,1);
    numomdp = zeros(NumOfFeatures,1);%numomdp - number of measured data points
    
    for j=1:NumOfFeatures
         
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   meanx(j) = meanx(j) + mf(2,i,j) - mf(2,i,k);
                   meany(j) = meany(j) + mf(3,i,j) - mf(3,i,k);
                   numomdp(j) = numomdp(j) + 1;
              end
         
         end
         
    end
    
    meanx = meanx./numomdp;
    meany = meany./numomdp;
    
    % Calculate sums
    sumxx=0;
    sumyy=0;
    sumxy=0;
    
    for j=1:NumOfFeatures
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   sumxx = sumxx + (mf(2,i,j) - mf(2,i,k) - meanx(j)).^2;
                   sumyy = sumyy + (mf(3,i,j) - mf(3,i,k) - meany(j)).^2;
                   sumxy = sumxy + (mf(2,i,j) - mf(2,i,k) - meanx(j)).*...
                                   (mf(3,i,j) - mf(3,i,k) - meany(j));
              end
        
         end
    end

    % Calculate Tiltaxis with linear regression
    psi = (0.5).*atan((2.*sumxy)./(sumxx-sumyy));
    if (sumxx > sumyy)
         psi = psi - (0.5).*pi.*sign(psi);
    end
    Tiltaxis(k) = (180./pi).*psi;
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % Calculation of k'th marker model
    % ---------------------------------------------------------------------
    m3d = zeros(3,NumOfFeatures);
    numomdp = zeros(NumOfFeatures,1);%numomdp - number of measured data points
    
    for j=1:NumOfFeatures
    
         sumxx=0;
         sumyy=0;
         sumxy=0;
         sumyx=0;
         
         alpha_a=0;
         alpha_b=0;
         
         P = zeros(3,3);
         temp = zeros(3);
         
         for i=1:NumOfProj
         
              if (mf(2,i,j) > -1) && (mf(2,i,k) > -1)
                   
                   sumxx = sumxx + (mf(2,i,j) - mf(2,i,k)).*cos(theta(i));
                   sumyy = sumyy + (mf(3,i,j) - mf(3,i,k)).*cos(theta(i));
                   sumxy = sumxy + (mf(2,i,j) - mf(2,i,k)).*sin(theta(i));
                   sumyx = sumyx + (mf(3,i,j) - mf(3,i,k)).*sin(theta(i));
                   
                   alpha_a = alpha_a + sin(theta(i)).^2;
                   alpha_b = alpha_b + sin(theta(i)).*cos(theta(i));
                   
                   numomdp(j) = numomdp(j) + 1;
                   
              end
              
         end
         
         P = [numomdp(j)-alpha_a.*sin(psi).^2   alpha_a.*cos(psi).*sin(psi)   alpha_b.*sin(psi);...
                alpha_a.*cos(psi).*sin(psi)   numomdp(j)-alpha_a.*cos(psi).^2   -alpha_b.*cos(psi);...
                    alpha_b.*sin(psi)               -alpha_b.*cos(psi)                alpha_a];
         
         temp(1) = (sumxx.*sin(psi)-sumyy.*cos(psi)).*sin(psi) + (cos(psi).*meanx(j)+sin(psi).*meany(j)).*cos(psi).*numomdp(j);
         temp(2) = -(sumxx.*sin(psi)-sumyy.*cos(psi)).*cos(psi) + (cos(psi).*meanx(j)+sin(psi).*meany(j)).*sin(psi)*numomdp(j);
         temp(3) = sumxy.*sin(psi) - sumyx.*cos(psi);
         
         dt = det(P);
         
         if (dt ~= 0)
              
              P_t = P;
              P_t(1,1) = temp(1);
              P_t(2,1) = temp(2);
              P_t(3,1) = temp(3);
              m3d(1,j) = det(P_t)./dt + mf(2,RefProj,k) - (Imdim./2+1);
              
              P_t = P;
              P_t(1,2) = temp(1);
              P_t(2,2) = temp(2);
              P_t(3,2) = temp(3);
              m3d(2,j) = det(P_t)./dt + mf(3,RefProj,k) - (Imdim./2+1);
              
              P_t = P;
              P_t(1,3) = temp(1);
              P_t(2,3) = temp(2);
              P_t(3,3) = temp(3);
              m3d(3,j) = det(P_t)./dt + (Imdim./2+1) - (Imdim./2+1);
                  
         else
              
              m3d(1,j) = 1000000;
              m3d(2,j) = 1000000;
              m3d(3,j) = 1000000;
              
              warning(['Feature ' num2str(j) ' : undefined! det = 0! Click more!']);
         end
    
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % Calculation of k'th shift
    % ---------------------------------------------------------------------
    shift_pre = zeros(2,NumOfProj);
    shift = zeros(2,NumOfProj,NumOfFeatures);
    
    for i = 1:NumOfProj
    
    sumxx = 0;
    sumyy = 0;
    ndif = 0;
    
    for j = 1:NumOfFeatures
        
        shift_pre(1,i) = m3d(1,j).*((sin(psi).^2).*cos(theta(i))+cos(psi).^2) + ...
                         m3d(2,j).*sin(psi).*cos(psi).*(1-cos(theta(i))) + ...
                         m3d(3,j).*sin(psi)*sin(theta(i)) + (Imdim./2+1);

        shift_pre(2,i) = m3d(1,j).*sin(psi).*cos(psi).*(1-cos(theta(i))) + ...
                         m3d(2,j).*((cos(psi).^2).*cos(theta(i))+sin(psi).^2) - ... 
                         m3d(3,j).*cos(psi).*sin(theta(i)) + (Imdim./2+1);
        
        if (mf(2,i,j) > -1) && (m3d(1,j) ~= 1000000)
            
            ndif = ndif + 1;
            shift(1,i,j) = mf(2,i,j) - shift_pre(1,i);
            shift(2,i,j) = mf(3,i,j) - shift_pre(2,i);
           
        else
            
            shift(1,i,j) = 0;
            shift(2,i,j) = 0;
            
        end
        
        sumxx = sumxx + shift(1,i,j);
        sumyy = sumyy + shift(2,i,j);
    
    end
    
    % Mean values of shifts
    if (ndif > 0)
        tx(k,i) = sumxx./ndif;
        ty(k,i) = sumyy./ndif;
    else
        tx(k,i) = 1000000;
        ty(k,i) = 1000000;
    end
    end
    % ---------------------------------------------------------------------
    
    
    % Build global marker model
    for j=1:NumOfFeatures
    
         % Move tiltseries a 3d marker coordinates of ReferenceFeature to [0 0 0]
         m3d(1,j) = m3d(1,j) - mf(2,RefProj,k) + mf(2,RefProj,k) + (Imdim./2+1);
         m3d(2,j) = m3d(2,j) - mf(3,RefProj,k) + mf(3,RefProj,k) + (Imdim./2+1);
         m3d(3,j) = m3d(3,j) - (Imdim./2+1) + (Imdim./2+1) + (Imdim./2+1);
         
    end
    
    % Sum GlobalMarkerModel
    globalMarkerModel_sum = globalMarkerModel_sum + m3d;
    
    % Aktualize waitbar
    waitbar(k/NumOfFeatures);
    
end
close(h);

% GlobalMarkerModel
gmm = globalMarkerModel_sum./NumOfFeatures;

% Mean tiltaxis
TiltaxisMean(1:NumOfProj) = mean(Tiltaxis) + 90;

% Mean translation
txMean = mean(tx);
tyMean = mean(ty);

% Global marker model translation
shift_pre = zeros(2,NumOfProj);
shift = zeros(2,NumOfProj,NumOfFeatures);

for i = 1:NumOfProj
    
    sumxx = 0;
    sumyy = 0;
    ndif = 0;
    
    for j = 1:NumOfFeatures
        
         shift_pre(1,i) = (gmm(1,j)-(Imdim./2+1)).*((sin(psi).^2).*cos(theta(i))+cos(psi).^2) + ...
                          (gmm(2,j)-(Imdim./2+1)).*sin(psi).*cos(psi).*(1-cos(theta(i))) + ...
                          (gmm(3,j)-(Imdim./2+1)).*sin(psi)*sin(theta(i)) + (Imdim./2+1);

         shift_pre(2,i) = (gmm(1,j)-(Imdim./2+1)).*sin(psi).*cos(psi).*(1-cos(theta(i))) + ...
                          (gmm(2,j)-(Imdim./2+1)).*((cos(psi).^2).*cos(theta(i))+sin(psi).^2) - ... 
                          (gmm(3,j)-(Imdim./2+1)).*cos(psi).*sin(theta(i)) + (Imdim./2+1);
        
         if (mf(2,i,j) > -1) && (gmm(1,j) ~= 1000000)
            
              ndif = ndif + 1;
              shift(1,i,j) = mf(2,i,j) - shift_pre(1,i);
              shift(2,i,j) = mf(3,i,j) - shift_pre(2,i);
           
         else
            
              shift(1,i,j) = 0;
              shift(2,i,j) = 0;
            
         end
        
         sumxx = sumxx + shift(1,i,j);
         sumyy = sumyy + shift(2,i,j);
    
    end
    
    % Mean values of shifts
    if (ndif > 0)
         txgmm(i) = sumxx./ndif;
         tygmm(i) = sumyy./ndif;
    else
         txgmm(i) = 1000000;
         tygmm(i) = 1000000;
    end

end

