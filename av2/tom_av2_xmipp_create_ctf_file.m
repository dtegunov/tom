function tom_av2_xmipp_create_ctf_file(f_align2d,particle_sel,particle_binning,outpath,ctfDatFile,astig_flag)
%tom_av2_xmipp_create_ctf_file creates ctf-param and ctfdat files for xmipp
%
%   tom_av2_xmipp_create_ctf_file(f_align2d,particle_sel,particle_binning,outpath,ctfDatFile,astig_flag)
%
%  TOM_AV2_XMIPP_CREATE_CTF_FILE creates ctf-param and ctfdat files for xmipp
%  from the fitted micrographs                                
%  
%  
%
%PARAMETERS
%
%  INPUT
%   f_align2d           filename of the alignment file
%   particle_sel        xmipp sel file containing the particles
%   particle_binning    binning of the particles in the sel-file
%                       (refering to the org micrograph 
%                        ...pixelsize must be calculated for the ctf-param file)  
%   outpath             path and root-filename for the calculated ctf-param files 
%   ctfDatFile          path for the ctfDatFile      
%   astig_flag          (1) flag for taking mesured astig into account
%
%
%  OUTPUT
%
%EXAMPLE
%     
%  
% tom_av2_xmipp_create_ctf_file('14_corrf_low_128.mat','14_corrf_low_128.sel',1,'/fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_low_128_ctfDat/ctf_','14_corrf_low_128_ctfDat.txt');
% 
% creates a folder ctf_param_files where the xmipp ctf-param files are
% written and a text file outctf.txt where every particle is matched with a
% ctf-file
%
%NOTES
% 
% 1:
%   ctf-param file exist for every particle and carry the complete ctf
%   information 
%   
%
%   sampling_rate=        1
%   voltage=              200
%   defocusU=             -11600
%   defocusV=             -11600
%   azimuthal_angle=      0
%   spherical_aberration= 2
%   chromatic_aberration= 0
%   energy_loss=          0
%   lens_stability=       0
%   convergence_cone=     0
%   longitudinal_displace=0
%   transversal_displace= 0
%   Q0=                   -0.1
%   K=                    1
%
%
% 2:  ctfDatFiles are normal text files and describe the match between 
%     ctf-param files and particles 
%  
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_187794.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_187794.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_190033.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_190033.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_190215.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_190215.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_190449.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_190449.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_191343.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_191343.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_191414.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_191414.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_193621.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_193621.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_196693.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_196693.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_196816.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_196816.ctfparam 
%    /fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090805_p47f11/log/rec/parts_high_128/parts_196922.spi /fs/sun11/lv01/pool/pool-nickell2/scratch/test_fb_ctf_group/ctf_param_files/ctf_196922.ctfparam 
% 
%
%
% 3: 
%    Astigmatism from TOM 2 Xmipp (Gregorieff)
%    U(long Axis)  = Dz_det.*1e10 +  ((-Dz_delta_det.*1e10) ./2);
%    V(short Axis) = Dz_det.*1e10 -  ((-Dz_delta_det.*1e10) ./2); 
%    azimuthal_angle=Phi_0_det+135;
%
%
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by AK 09/09/09 by fb
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

if (isempty(particle_binning))
    particle_binning=0;
end;



if (nargin <6)
    astig_flag=1;
end;


try
    warning off;
    mkdir(fileparts(outpath)); 
    warning on;
catch
end;


particle_names=importdata(particle_sel);
particle_names=particle_names.textdata;


load(f_align2d);


fp=fopen(ctfDatFile,'wt');

[ctff_path ctff_base ctff_ext]=fileparts(ctfDatFile);

ctf_h_text_name=[ctff_path ctff_base '_allh' ctff_ext];

fp_ctfh=fopen(ctf_h_text_name,'wt');
fprintf(fp_ctfh,'%s %s %s %s %s %s %s %s %s %s\n','filename','file_idx','particlename','part_idx','pos_x','pos_y','pixelsize','defocus','delta defocus','astig angle');

disp(['Creating']);

file_name_old='';

for i=1:size(align2d,2)
    
    file_name=[strrep(align2d(1,i).filename,'_corr/','/') ,'.mat'];
    if (strcmp(file_name,file_name_old)==0)
        load(file_name);
        file_name_old=file_name;
    end;
    
    
    [a b c]=fileparts(particle_names{i});
    
    t_pos=strfind(b,'_');
    t_pos=min(t_pos);
    name=b(t_pos+1:end);
    %[name num]=strtok(b,'_');
    %num=strrep(strrep(num,'.spi',''),'_','');
    num=name;
    ctf_out_name=[outpath  num '.ctfparam'];
    
    %transfer parameters
    sampling_rate=st_out.Fit.EM.Objectpixelsize.*1e10;
    Voltage=st_out.Fit.EM.Voltage./1000;
    spherical_aberration=st_out.Fit.EM.Cs.*1000;
    chromatic_aberration=st_out.Fit.EM.Cc.*1000;
    
    if (astig_flag==1)
        dz_u=st_out.Fit.Dz_det.*1e10+((-st_out.Fit.Dz_delta_det.*1e10)./2);
        dz_v=st_out.Fit.Dz_det.*1e10-((-st_out.Fit.Dz_delta_det.*1e10)./2);
        ang=st_out.Fit.Phi_0_det+135;
    else
        dz_u=st_out.Fit.Dz_det.*1e10;
        dz_v=st_out.Fit.Dz_det.*1e10;
        ang=0;
    end;
    
    
    fp_ctf=fopen(ctf_out_name,'wt');
    fprintf(fp_ctf,['sampling_rate=        '  num2str(sampling_rate.*2^particle_binning) '\n']);
    fprintf(fp_ctf,['voltage=              '  num2str(Voltage) '\n']);
    fprintf(fp_ctf,['defocusU=             '  num2str(dz_u) '\n']);
    fprintf(fp_ctf,['defocusV=             '  num2str(dz_v) '\n']);
    fprintf(fp_ctf,['azimuthal_angle=      '  num2str(ang) '\n']);
    fprintf(fp_ctf,['spherical_aberration= '  num2str(spherical_aberration) '\n']);
    fprintf(fp_ctf,['chromatic_aberration= '  num2str(chromatic_aberration) '\n']);
    fprintf(fp_ctf,['energy_loss=          '  num2str(0) '\n']);
    fprintf(fp_ctf,['lens_stability=       '  num2str(0) '\n']);
    fprintf(fp_ctf,['convergence_cone=     '  num2str(0) '\n']);
    fprintf(fp_ctf,['longitudinal_displace='  num2str(0) '\n']);
    fprintf(fp_ctf,['transversal_displace= '  num2str(0) '\n']);
    fprintf(fp_ctf,['Q0=                   '  num2str(-0.1) '\n']);
    fprintf(fp_ctf,['K=                    '  num2str(1) '\n']);
    
    
    fclose(fp_ctf);
    
    
    fprintf(fp,[particle_names{i} ' ' ctf_out_name ' \n']);
 
    [dd tmp_part da]=fileparts(particle_names{i});
    [rest part_nr]=strtok(tmp_part,'_');
    part_nr=str2double(strrep(part_nr,'_',''));
   
    [dd tmp_file da]=fileparts(align2d(1,i).filename);
    [rest file_nr]=strtok(tmp_file,'_');
    file_nr=str2double(strrep(file_nr,'_',''));
    
    fprintf(fp_ctfh,'%s %d %s %d %d %d %f %f %f %f\n',align2d(1,i).filename,file_nr,particle_names{i},part_nr,round(align2d(1,i).position.x),round(align2d(1,i).position.y),sampling_rate.*2^particle_binning,st_out.Fit.Dz_det.*1e6,st_out.Fit.Dz_delta_det.*1e6,st_out.Fit.Phi_0_det);
        
    if (mod(i,1000)==0)
        disp([ num2str(i) ' particles done!']);
    end;

    
    
end;


fclose(fp);
fclose(fp_ctfh);

disp('end.');

