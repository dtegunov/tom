function average=tom_av3_average3(Align,method,threshold,iclass,alg_part_path,waitbarflag)
%  average=tom_av3_average3(Align,method,threshold,iclass,waitbarflag)
%  
%     tom_av3_average3(Align,method,threshold,iclass,waitbarflag)
%  
%  PARAMETERS
%  
%    INPUT
%     Align               align struct
%     method              ('inverse') or direct or direct_dd                            
%     threshold           threshold for Align.CCC
%     iclass              (use all) class number 
%     alg_part_path       (opt.) base path for aligned parts
%     waitbarflag         (0) flag for graphical waibar 
%    
%    OUTPUT
%     average           particle average
%        
%  
%  EXAMPLE
%    
%   %generate an avg
%   avg=tom_av3_average3(Align(end,:),'inverse',99,[0 1]);
%   
%   %generate an avg and write aligned particles
%   mkdir('my_alg_parts');
%   avg=tom_av3_average3(Align(end,:),'inverse',0,'','my_alg_parts/parts_');
%
%  NOTE:
% 
%   create folder 4 output !!
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by fb okt.2011
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom



if nargin < 6
    waitbarflag = 0;
end

if nargin < 5
    alg_part_path = '';
end

if nargin<4 || isempty(iclass)
    iclass = 0:2000;
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


try
    tsize = Align(1).Tomogram.Header.Size';
catch ME
    prompt = {'Enter templatesize x:','Enter templatesize y:','Enter templatesize z:'};
    dlg_title = 'Missing template size';
    num_lines = 1;
    def = {'','',''};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
    tsize = [str2num(answer{1}), str2num(answer{2}), str2num(answer{3})];
end



icount = 0;

for indpart = 1:size(Align,2)
    if (Align(indpart).CCC)>=threshold &&  (isempty(find(ismember(iclass,Align(indpart).Class )))==0)  
        icount = icount +1;
        itomo = Align(indpart).Filename;
        xshift = Align(indpart).Shift.X;
        yshift =Align(indpart).Shift.Y;
        zshift = Align(indpart).Shift.Z;
        tshift = [xshift yshift zshift];
        phi= Align(indpart).Angle.Phi;
        psi=Align(indpart).Angle.Psi;
        the=Align(indpart).Angle.Theta;
        ifile = indpart;
        
        name=Align(indpart).Filename;
        particle = tom_emreadc(name); particle = particle.Value;
        
        if icount == 1
            wei = zeros(size(particle,1),size(particle,2),size(particle,3));
            average = wei;
        end;
        
        maxangle = Align(indpart).Tomogram.AngleMin;
        minangle = Align(indpart).Tomogram.AngleMax;
        %FIXME
        if maxangle == 0
            maxangle = -65;
        end
        
        if minangle == 0
            minangle = 65;
        end
        
        wedge = tom_av3_wedge(particle,minangle,maxangle);
       
        if isequal(method,'inverse')
            particle = double(tom_rotate(tom_shift(particle,-tshift),[-psi -phi -the]));
        elseif isequal(method,'direct')
            particle = double(tom_rotate(tom_shift(particle,tshift),[phi psi the]));
        elseif isequal(method,'direct_dd')
            particle = double(tom_shift(tom_rotate(particle,[phi psi the]),tshift));
        else %sum
            particle = double(tom_shift(tom_rotate(particle,Align(indpart).Angle.Rotmatrix),tshift));
        end;
        
        
        particle = tom_norm(particle,1);
        average = average + particle;
        
         
        if isempty(alg_part_path)==0
            tom_emwritec([alg_part_path num2str(indpart) '.em'],particle);
        end;
        
        
        if isequal(method,'inverse')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi,-phi,-the])),0.5,1,'z'),0,0.5);
        elseif isequal(method,'direct')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[phi psi the])),0.5,1,'z'),0,0.5);
        elseif isequal(method,'direct_dd')
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[phi psi the])),0.5,1,'z'),0,0.5);
        else
            tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,Align(indpart).Angle.Rotmatrix)),0.5,1,'z'),0,0.5);
        end;
        wei = wei + tmpwei;
        if waitbarflag == 1
            waitbar(indpart./size(Align,2),h,[num2str(indpart), ' of ', num2str(size(Align,2)), ' files done.']);
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