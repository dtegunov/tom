function num_of_ctf_gr=tom_av2_xmipp_ctf_group(outputfolder,particle_list,split_file,f_fit_st,pix_size,sz_im,wiener_const,first_max_corr,b_fun)
%TOM_AV2_XMIPP_CTF_GROUP creates a ctf-group folder for
% xmipp_protocol_projmatch
%
%   tom_av2_xmipp_ctf_group(outputfolder,particle_list,split_file,f_fit_st,pix_size,sz_im,wiener_const,first_max_corr,b_fun)
%
%  TOM_AV2_XMIPP_CTF_GROUP  creates a ctf-group folder for
%  xmipp_protocol_projmatch which has 2 reasons: memory leak of xmipp_ctf_group
%                                                ctf-models of xmipp and tom differ  
%          
%  PARAMETERS
%
%  INPUT
%   outputfolder      output foldername (use CtfGroups 4 xmipp_protocols)   
%   particle_list     2 column text-file containing particlename and defocus in Ang
%   split_file        xmipp split text-file containing the defocus (in Ang) interv for the groups
%   f_fit_st          filename of ctf-fit struct
%   pix_size          pixelsize in Ang of the padded filter (pixs is the same as the org images )
%   sz_im             size wiener filter should have (org img_size * padding of xmipp_protocols )
%   wiener_const      (0.1) 1./SNR            
%   first_max_corr    (0) use 1 for correcting before the first max
%   b_fun             (no b-factor applied) or b-factor
%  
%  OUTPUT
%   num_of_ctf_gr      number of resulting ctf groups
%
%
%EXAMPLE
%  %use tom functio 2 create particle vs defocus text-file
%  tom_av2_xmipp_ctf_gen_defocus_partlist('picklist.mat','all_part.sel','parts_defocus.txt');
%  
%  tom_av2_xmipp_ctf_group('CtfGroups','parts_defocus.txt','split.doc','low_1.em.mat',1,[256 256],0.1,1);
%  
%  %gerate particle-list (parts_defocus.txt) from tom_av2_xmipp_create_ctf_file text-file output
%  %cat 04_ctfDat_allh.txt | awk '{print $3 " " $8*10000}' > parts_defocus.txt 
%
%  % to use in with: xmipp_protocols_projmatch
%  % 1: create a directory: mkdir  ProjMatch/run1/
%  % 2: copy CtfGroups:  cp -R CtfGroups ProjMatch/run1/
%  % 3: create a empty ctfdat file: touch all_images.ctfdat
%  % 4: modify script (comment the call of xmipp_ctf_groups and set Reference correction 2 false)
%
%   # Make CTF groups
%
%   if (self._DoCtfCorrection):
%           self._NumberOfCtfGroups=execute_ctf_groups(self._mylog,
%                                                      self._SelFileName,
%                                                      self._CTFDatName,
%                                                      self._PaddingFactor,
%                                                      self._DataArePhaseFlipped,
%                                                      self._WienerConstant,
%                                                      self._DoAutoCtfGroup,
%                                                      self._CtfGroupMaxDiff,
%                                                      self._CtfGroupMaxResol,
%                                                      self._SplitDefocusDocFile)
%        else:
%           self._NumberOfCtfGroups=1
% 
%  to  
%
%  # Make CTF groups
%       
%   self._NumberOfCtfGroups=12
%
%
% and
%
%  # Initial reference is CTF-amplitude corrected?
%    if ( (_iteration_number == 1) and (_ReferenceIsCtfCorrected==False) ):
%       self._ReferenceIsCtfCorrected=False
%    else: 
%       self._ReferenceIsCtfCorrected=True
% 
% to
%
%   self._ReferenceIsCtfCorrected=False  
%
%
%REFERENCES
%
%  www.wadsworth.org/spider_doc/spider/docs/techs/ctf/ctf.html#ref
%
%NOTE:
% 
% 
% % head parts_defocus.txt 
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_1.spi -26000.00
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_2.spi -26000.00
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_3.spi -26000.00
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_4.spi -26000.00
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_5.spi -26000.00
% % /fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/autopick/result_slot3/all/org/part_fold1/part_6.spi -26000.00
%
% % head split.doc
% % ; Defocus values to split into 9 ctf groups
% %    1 1  -15000.00000
% %    2 1  -16000.00000
% %    3 1  -17000.00000
% %    4 1  -18000.00000
% %    5 1  -19000.00000
% %    6 1  -20000.00000
% %    7 1  -21000.00000
% %    8 1  -22000.00000
% %    9 1  -23000.00000
%
% %wiener defined as:
% wiener_f = ctf_theory ./ ((ctf_theory.*ctf_theory) + wiener_const);
%
%
%
%SEE ALSO
%   
% tom_av2_xmipp_ctf_gen_defocus_partlist, tom_calc_b_fact_weight_fun
%
% NOTE2:
% 

%
%
%   created by fb ...ole !!
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



disp('Reading input data!');

warning off; mkdir(outputfolder); warning on;

%parts_and_def=importdata(particle_list);

tmp_fid=fopen(particle_list,'r');
C = textscan(tmp_fid, '%s %f\n');
parts_and_def.textdata= C{1,1};
parts_and_def.data = C{1,2};
fclose(tmp_fid);

load(f_fit_st);
Fit=st_out.Fit;

pix_size=pix_size.*1e-10;
voltage=Fit.EM.Voltage;
img_size=sz_im;
sigma=Fit.decay_det;
q_0_r=Fit.decay_part_coh_ill_det;
Cs=Fit.EM.Cs;
Cc=Fit.EM.Cc;
deltaE=Fit.decay_energy_spread_det;
amplitude_contrast=Fit.amplitude_contrast_det;


split=importdata(split_file);


if (nargin < 9)
    b_fun=ones(sz_im);
else
    ttmp=max(size(b_fun));
    if(ttmp(1)==1)
        b_fun=tom_calc_b_fact_weight_fun(sz_im,b_fun,pix_size.*1e10,(pix_size.*1e10).*2+((pix_size.*1e10).*2).*0.1);
    end;
end;

tom_spiderwrite([outputfolder '/b_fun.spi'],b_fun);

u_bound=1; 

fid_log=fopen([outputfolder '/log_file.log' ],'wt');

fid_df=fopen([outputfolder '/ctf_groups.defocus' ],'wt');
fprintf(fid_df,'%s\n','# Defocus values for each group (avg, max, min & std) ');

fid_imgN=fopen([outputfolder '/ctf_groups.imgno' ],'wt');
fprintf(fid_imgN,'%s\n','# Number of images in each group');

zz_group=1;

disp('done!');

disp_param(fid_log,outputfolder,particle_list,split_file,f_fit_st,pix_size,sz_im,wiener_const,first_max_corr,b_fun,voltage,sigma,q_0_r,Cs,Cc,deltaE,amplitude_contrast);



for i=1:size(split.data,1)+1
    if (i > size(split.data,1))
        l_bound=-100000000000;
    else
        l_bound=split.data(i,3);
    end;
    
    gr_idx=find((parts_and_def.data >= l_bound) .* (parts_and_def.data < u_bound) );
    if (isempty(gr_idx)==0)
        me_df(zz_group)=mean(parts_and_def.data(gr_idx));
        std_df(zz_group)=std(parts_and_def.data(gr_idx));
        num(zz_group)=length(gr_idx);
        min_df(zz_group)=min(parts_and_def.data(gr_idx));
        max_df(zz_group)=max(parts_and_def.data(gr_idx));
        
        if (zz_group<10)
            base=[outputfolder '/ctf_group00000' num2str(zz_group)];
        end;
        if (zz_group>=10 && zz_group < 100)
            base=[outputfolder '/ctf_group0000' num2str(zz_group)];
        end;
        if (zz_group>=100)
            base=[outputfolder '/ctf_group000' num2str(zz_group)];
        end;
        
         fp_sel=fopen([base '.sel' ],'wt');
        
        for ii=1:length(gr_idx)
            fprintf(fp_sel,'%s 1\n',parts_and_def.textdata{gr_idx(ii)});
        end;
        fclose(fp_sel);
        fprintf(fid_df,'%d %f %f %f %f\n',zz_group,me_df(zz_group),max_df(zz_group),min_df(zz_group),std_df(zz_group));
        fprintf(fid_imgN,'%d %f\n',zz_group,num(zz_group));
        
        Dz=me_df(zz_group)*1e-10;
        [phase amplitude decay E_k E_i_k E_e_k]=tom_ctf2(Dz,0,0,pix_size,voltage,img_size,Cs,sigma,q_0_r,Cc,deltaE);
        ctf_theory=(sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude);
        %ctf_theory=E_e_k.*E_i_k.*decay.*(sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude);
        
        ctf_theory=abs(ctf_theory);
        
        
        %no correction up to 1. max
        if (first_max_corr==0)
            mid=floor(size(ctf_theory)./2)+1;
            inner_cut=min(tom_crossing(diff(ctf_theory(mid(1)+1:end,mid(1)) )))+1;
            inner_idx=find(tom_spheremask(ones(sz_im),inner_cut));
            ctf_theory(inner_idx)=1;
        end;
        
        wiener_f = ctf_theory ./ ((ctf_theory.*ctf_theory) + wiener_const);
        
        wiener_f = wiener_f .* b_fun;
        tom_spiderwrite([base '.wien'],fftshift(wiener_f));
        disp(['Group nr ' num2str(zz_group) ' from ' num2str(l_bound) ' to ' num2str(u_bound) ]);
        disp(['  nr parts: ' num2str(num(zz_group))]);
        disp(['  defocus mean: ' num2str(me_df(zz_group)) ' max ' num2str(max_df(zz_group))  ' min ' num2str(min_df(zz_group)) ' std ' num2str(std_df(zz_group)) ]);
        
        u_bound=l_bound;
        zz_group=zz_group+1;
        
    else
        disp(['warning: defoucus group: nr ' num2str(i) ' is empty']);
        disp(['No particles between ' num2str(l_bound) ' and ' num2str(u_bound) ]);
    end;
    
end;

disp(' ');
disp(' ');

disp(['length of input: ' num2str(length(parts_and_def.data ))])
disp(['Sum of all groups: ' num2str(sum(num))]);


fclose(fid_df);
fclose(fid_imgN);
[a b]=unix(['cp ' split_file ' ' outputfolder '/split_backup.doc' ]);

num_of_ctf_gr=zz_group-1;
disp('done !');

function disp_param(fid_log,outputfolder,particle_list,split_file,f_fit_st,pix_size,sz_im,wiener_const,first_max_corr,b_fun,voltage,sigma,q_0_r,Cs,Cc,deltaE,amplitude_contrast)

disp(' ');
disp('======================================>');
fprintf(fid_log,'%s\n','======================================>');
disp(['outputfolder: ' outputfolder]);
fprintf(fid_log,'%s\n',['outputfolder: ' outputfolder]);
disp(['particle_list: ' particle_list]);
fprintf(fid_log,'%s\n',['particle_list: ' particle_list]);
disp(['split_file: ' split_file]);
fprintf(fid_log,'%s\n',['split_file: ' split_file]);
disp(['f_fit_st: ' f_fit_st]);
fprintf(fid_log,'%s\n',['f_fit_st: ' f_fit_st]);
disp(['pix_size: ' num2str(pix_size)]);
fprintf(fid_log,'%s\n',['pix_size: ' num2str(pix_size)]);
disp(['sz_im: ' num2str(sz_im)]);
fprintf(fid_log,'%s\n',['sz_im: ' num2str(sz_im)]);
disp(['wiener_const: ' num2str(wiener_const)]);
fprintf(fid_log,'%s\n',['wiener_const: ' num2str(wiener_const)]);
disp(['first_max_corr: ' num2str(first_max_corr)]);
fprintf(fid_log,'%s\n',['first_max_corr: ' num2str(first_max_corr)]);
if (std(b_fun(:))==0 )
    disp(['b_fun: no b factor applied']);
    fprintf(fid_log,'%s\n',['b_fun: no b factor applied']);
else
    disp(['b_fun: b factor applied']);
    fprintf(fid_log,'%s\n',['b_fun: no b factor applied']);
end;
disp(['voltage: ' num2str(voltage)]);
fprintf(fid_log,'%s\n',['voltage: ' num2str(voltage)]);
disp(['sigma: ' num2str(sigma)]);
fprintf(fid_log,'%s\n',['sigma: ' num2str(sigma)]);
disp(['q_0_r: ' num2str(q_0_r)]);
fprintf(fid_log,'%s\n',['q_0_r: ' num2str(q_0_r)]);
disp(['Cs: ' num2str(Cs)]);
fprintf(fid_log,'%s\n',['Cs: ' num2str(Cs)]);
disp(['Cc: ' num2str(Cc)]);
fprintf(fid_log,'%s\n',['Cc: ' num2str(Cc)]);
disp(['deltaE: ' num2str(deltaE)]);
fprintf(fid_log,'%s\n',['deltaE: ' num2str(deltaE)]);
disp(['amplitude_contrast: ' num2str(amplitude_contrast)]);
fprintf(fid_log,'%s\n',['amplitude_contrast: ' num2str(amplitude_contrast)]);
disp('<======================================');
fprintf(fid_log,'%s\n','<======================================');
disp(' ');



