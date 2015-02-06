function tom_av2_read_mrcstacks(flag,base_mrcfilename,base_emfilename,emfilenumberlow,def1,def2,obj_pixs,Cs,volt,mic)

%TOM_AV2_READ_MRCSTACKS reads a numbered list of mrc-stacks and extracts
%the .em-files (focal pairs or singles), modifying the em-header.
%
%   tom_av2_read_mrcstacks(flag,base_mrcfilename,base_emfilename,emfilenumberlow,def1,def2,obj_pixs,Cs,volt,mic)
%
%PARAMETERS
%
%  INPUT
%   flag                'pairs' or 'singles', depending on acquisition scheme
%   base_mrcfilename    wildcard for mrc files (e.g 100701_yys1278.*.mrc) ... like using ls
%   base_emfilename     output filename without number or extension
%   emfilenumberlow     starting number of output .em-file (singles or first image of pair)
%   def1                estimated defocus of first image of pair
%   def2                estimated defocus of second image of pair
%   obj_pixs            object pixelsize
%   Cs                  spherical abberation of microscope
%   volt                acceleration voltage
%   mic                 microscope type
%    
%
%EXAMPLE
%   tom_av2_read_mrcstacks('singles','100701_yys1278*.mrc','yys1278_',11433,2.5,2.5,2.21,2,200000,'F20');
%   tom_av2_read_mrcstacks('pairs','100701_yys1278*.mrc','yys1278_',11433,1.5,3.0,2.21,2,200000,'F20');
%
%REFERENCES
%
%SEE ALSO
%
%   created by SB 10/02/24
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

if nargin<10
    mic='F20';
end;
if nargin<9
    volt=200000;
end;
if nargin<8
    Cs=2;
end;
if nargin<7
    obj_pixs=2.21;
end;
if nargin<6
    def2=3.0;
end;

tilt_ang_flag=1;

dd=dir(base_mrcfilename);

base_path_mrc=fileparts(base_mrcfilename);

emfilenumberhigh=emfilenumberlow;

switch flag;

    case 'pairs'

        warning off;
        mkdir('low');
        mkdir('high');
        warning on;
        
    for k=1:length(dd)
        
        if (isempty(base_path_mrc)==0)
            tmp_filename_m=[base_path_mrc '/' dd(k).name];
        else
            tmp_filename_m=[dd(k).name];
        end;
        
        a=tom_mrcread(tmp_filename_m);
        szz=a.Header.Size;
        
        
        clear('a');
        s_tmp=dir(tmp_filename_m);
        s_tmp=s_tmp.bytes;
        sz_head=s_tmp-(szz(1).*szz(2).*szz(3)*2);
        disp(['Reading ' tmp_filename_m 'Sz: ' num2str(szz) '    Head: ' num2str(sz_head)]);
        a.Value=tom_rawread(tmp_filename_m,'int16','le',[szz(1) szz(2) szz(3)],0,sz_head);
        
        try
            setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH') ':/usr/local/apps/IMOD/lib']);
            [aa bb]=unix(['extracttilts ' tmp_filename_m ' tilts.txt']);
            tilts=importdata('tilts.txt');
            if (length(tilts) ~= size(a.Value,3))
                disp('warning: Tiltangle extraction not possible');
                tilt_ang_flag=0;
            else
                disp('Tiltangles extracted !');
            end;
        catch Me
            tilt_ang_flag=0;
            disp(' ');
            disp('warnin: Could not extract tilt angles');
            disp(Me.message);
        end;
        
        
        l=size(a.Value,3);
    
        for i=1:2:l; 
            tmp=tom_emheader(a.Value(:,:,i));
            tmp.Header.Defocus=-def1.*10000;
            tmp.Header.Objectpixelsize=obj_pixs;
            tmp.Header.Microscope=mic;
            tmp.Header.Voltage=volt;
            tmp.Header.Cs=Cs;
             if (tilt_ang_flag==1)
                tmp.Header.Tiltangle=tilts(i);
            end;
            tmp.Value=tom_xraycorrect2(tmp.Value);
            tom_emwritec(['low/' base_emfilename num2str(emfilenumberlow) '.em'],tmp); 
            disp(['low/' base_emfilename num2str(emfilenumberlow) '.em Tiltangle: ' num2str(tmp.Header.Tiltangle)]);  
            emfilenumberlow=emfilenumberlow+1; 
        end;          
       
        
        
        for i=2:2:l; 
            tmp=tom_emheader(a.Value(:,:,i));
            tmp.Header.Defocus=-def2.*10000;
            tmp.Header.Objectpixelsize=obj_pixs;
            tmp.Header.Microscope=mic;
            tmp.Header.Voltage=volt;
            tmp.Header.Cs=Cs;
            if (tilt_ang_flag==1)
                tmp.Header.Tiltangle=tilts(i);
            end;
            tmp.Value=tom_xraycorrect2(tmp.Value);
            tom_emwritec(['high/' base_emfilename num2str(emfilenumberhigh) '.em'],tmp); 
            disp(['high/' base_emfilename num2str(emfilenumberhigh) '.em  Tiltangle: ' num2str(tmp.Header.Tiltangle)]);  
            emfilenumberhigh=emfilenumberhigh+1;  
        end;  
        
    
    end;

    case 'singles'

        warning off;
        mkdir('low');
        warning on;
        
        
    for k=1:length(dd)
        
        if (isempty(base_path_mrc)==0)
            tmp_filename_m=[base_path_mrc '/' dd(k).name];
        else
            tmp_filename_m=[dd(k).name];
        end;
        
        a=tom_mrcread(tmp_filename_m);
        szz=a.Header.Size;
        clear('a');
        s_tmp=dir(tmp_filename_m);
        s_tmp=s_tmp.bytes;
        sz_head=s_tmp-(szz(1).*szz(2).*szz(3)*2);
        disp(['Reading ' tmp_filename_m '  Sz: ' num2str(szz) '     Head: ' num2str(sz_head)]);
        a.Value=tom_rawread(tmp_filename_m,'int16','le',[szz(1) szz(2) szz(3)],0,sz_head);
        l=size(a.Value,3);
    
        try
            setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH') ':/usr/local/apps/IMOD/lib']);
            [aa bb]=unix(['extracttilts ' tmp_filename_m ' tilts.txt']);
            tilts=importdata('tilts.txt');
            if (length(tilts) ~= size(a.Value,3))
                disp('warning: Tiltangle extraction not possible');
                tilt_ang_flag=0;
            else
                disp('Tiltangles extracted !');
            end;
        catch Me
            tilt_ang_flag=0;
            disp(' ');
            disp('warnin: Could not extract tilt angles');
            disp(Me.message);
        end;
        
        
        
        for i=1:l; 
            tmp=tom_emheader(a.Value(:,:,i));
            tmp.Header.Defocus=-def1.*10000;
            tmp.Header.Objectpixelsize=obj_pixs;
            tmp.Header.Microscope=mic;
            tmp.Header.Voltage=volt;
            tmp.Header.Cs=Cs;
            if (tilt_ang_flag==1)
                tmp.Header.Tiltangle=tilts(i);
            end;
            tmp.Value=tom_xraycorrect2(tmp.Value);
            tom_emwritec(['low/' base_emfilename num2str(emfilenumberlow) '.em'],tmp); 
            disp([base_emfilename num2str(emfilenumberlow) '.em Tiltangle ' num2str(tmp.Header.Tiltangle) ]);
            emfilenumberlow=emfilenumberlow+1; 
        end;          
      
    
    end;
    
    otherwise
            disp('Unknown method.')

end;