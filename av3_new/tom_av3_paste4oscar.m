function Align = tom_av3_paste4oscar(filenames, volname, mask, Align, filter, normflag, waitbarflag)
%TOM_AV3_PASTE4OSCAR pastes big volumes with borders from subvolumes
%
%   Align = tom_av3_paste4oscar(filenames, volname, mask, Align, filter, normflag, waitbarflag)
%
%PARAMETERS
%
%  INPUT
%   filenames           2D char array with the full pathnames of the input volumes,
%                        can be generated with strvcat
%   volname             full pathname to output volume
%   mask                full pathname to mask
%   Align               particle alignment structure (optional)
%   filter              [low high] bandpass filter in pixel (optional)
%   normflag            norm particles to phase contrast (optional)
%   waitbarflag         1 = waitbar gui
%                       0 = text messages
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_paste4oscar(strvcat('test1.vol','test2.vol'),'outvol.vol',Align,0)
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI
%
%   created by AK 06/10/05
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


if nargin < 7
    waitbarflag = 0;
end
if nargin < 6
    normflag = 0;
end
if nargin < 5
    filter = [0 0];
end
if nargin < 4
    Align = [];
end

for j = 1:size(filenames,1)
    if ~exist(strtrim(filenames(j,:)),'file')
        error(['File ' strtrim(filenames(j,:)) ' does not exist!']);
    end;
end;

tempheader = tom_reademheader(strtrim(filenames(1,:)));
dimensions = [tempheader.Header.Size(1),tempheader.Header.Size(2),tempheader.Header.Size(3)];

mask = tom_emreadc(mask);

if ~(mask.Header.Size == tempheader.Header.Size') == [1 1 1]
    errordlg('Mask must have the same dimensions as subvolumes!');
    return;
end

xbuffer = floor(size(filenames,1)^(1./3));

sizex = xbuffer*dimensions(1)+dimensions(1);
sizey = xbuffer*dimensions(2)+dimensions(2);
zbuffer = ceil(size(filenames,1)./xbuffer^2);
sizez = zbuffer*dimensions(3)+dimensions(3);
if waitbarflag == 1
    h = waitbar(0,'Creating new volume');
else
    disp(['Creating new volume ', num2str(sizex), 'x',num2str(sizey), 'x',num2str(sizez)]);
end
tom_emwritec(volname, [sizex sizey sizez],'new');

borderx = floor(dimensions(1)./2);
bordery = floor(dimensions(2)./2);
borderz = floor(dimensions(3)./2);

header = tom_reademheader(volname);
header.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
%header.Header.Comment = strcat(num2str(borderx),':',num2str(bordery),':',num2str(borderz),':',num2str(xbuffer),':',num2str(zbuffer));
tom_writeemheader(volname, header.Header);

x = borderx; y = bordery; z = borderz;
j = size(filenames,1);

xcounter = 1;
ycounter = 1;

for i = 1:j
    
    f = tom_emreadc(strtrim(filenames(i,:)));

    %if alignment structure is provided, shift and rotate particles before
    %pasting
    if ~isempty(Align)

        outlist = tom_av3_align_sum(Align);
 %       outlist=Align;
        shiftx = outlist(i).Shift.X;
        shifty = outlist(i).Shift.Y;
        shiftz = outlist(i).Shift.Z;
        phi = outlist(i).Angle.Phi;
        psi = outlist(i).Angle.Psi;
        theta = outlist(i).Angle.Theta;
        rotmatrix = outlist(i).Angle.Rotmatrix;
        tshift = [shiftx shifty shiftz];
        %rotate
        if ~(phi == 0 && psi == 0 && theta == 0)
            f.Value = single(f.Value);
            %rotation matrix
%            f.Value = double(tom_shift(tom_rotate(f.Value,rotmatrix),tshift));
%            tom_dspcub(f.Value);drawnow;
            f.Value = tom_rotate(f.Value,rotmatrix,'linear'); % original version

%            f.Value = double(tom_rotate(tom_shift(f.Value,-[shiftx shifty shiftz]),[-psi -phi -theta]));
%            f.Value = double(tom_rotate(tom_shift(f.Value,[shiftx shifty shiftz]),[phi psi theta]));
            %f.Value = tom_rotate(f.Value,[phi psi theta],'linear');
        end

        %shift changed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if ~(shiftx == 0 & shifty == 0 & shiftz == 0)
             f.Value = tom_shift(double(f.Value),[shiftx shifty shiftz]);
         end
        
    end
    
    %filter image
    if filter(1) ~= 0 | filter(2) ~= 0
        f.Value = tom_bandpass(f.Value,filter(1),filter(2));
        f.Value = f.Value .* (dimensions(1) .* dimensions(2) .* dimensions(3));
        Align(size(Align,1),i).Filter = filter;
    end
    
    %norm image
    if normflag == 1
       f.Value = tom_norm(f.Value,'phase');
       Align(size(Align,1),i).NormFlag = 1;
    end
    
    %apply mask
    f.Value = f.Value.*mask.Value;
        
    tom_emwritec(volname,f.Value,'subregion',[x+1 y+1 z+1],[dimensions(1) dimensions(2) dimensions(3)]);
    
    if xcounter < xbuffer
        x = x + dimensions(1);
        xcounter = xcounter + 1;
    else
        x = borderx;
        y = y + dimensions(2);
        xcounter = 1;
        ycounter = ycounter + 1;
        if ycounter > xbuffer
            ycounter = 1;
            y = bordery;
            z = z + dimensions(3);
        end
    end
    
    if waitbarflag == 1
        waitbar(i./j,h,[num2str(i), ' of ', num2str(j), ' files done.']);
    elseif waitbarflag == 0
        disp([num2str(i), ' of ', num2str(j), ' files done.']);
    end
end;

if waitbarflag == 1
    close(h);
end