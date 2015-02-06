function tom_av2_xmipp_ctf_gen_defocus_partlist(f_sel_file,f_pickList,f_def_part_list,sumExt)
%tom_av2_xmipp_ctf_gen_defocus_partlist generates a particle vs defocus
%textfile 4 tom_av2_xmipp_ctf_group
%
%   tom_av2_xmipp_ctf_gen_defocus_partlist(ctf_fold,fit_st)
%
%  tom_av2_xmipp_ctf_gen_defocus_partlist 
%  
%
%PARAMETERS
%
%  INPUT
%   f_sel_file        name of sel file
%   f_pickList        name of picklist
%   f_def_part_list   name of the output defocus file
%   sumExt            extension of sums (only importend 4 frames)
%
%EXAMPLE
%  
%  tom_av2_xmipp_ctf_gen_defocus_partlist('all_part.sel','picklist.mat','parts_defocus.txt');
%  %example 4 frames
%  tom_av2_xmipp_ctf_gen_defocus_partlist('all_part_frames.sel','picklist_frmaes.mat','parts_defocus_frames.txt','.em');
%
%REFERENCES
%
%  
%
%SEE ALSO
%   
% tom_av2_xmipp_ctf_group
%
% NOTE:
% 
% sel-file and picklist have to have same order !! (...use tom_av2_xmipp_picklist2stack)
% for e.g. picklist(1,5) has to correspond 2 the paricle in row nr 5 in the sel-file
%
% convetion 4 frames is that sel and pickList correspond in Order !!!  
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


sel=importdata(f_sel_file);

load(f_pickList);

if (iscell(align2d(1,1).filename))
     if (size(align2d,2)*length(align2d(1,1).filename)~=length(sel.textdata))
        error(['Picklist and sel-file differ in size!']);
    end;
    frmaeList=1;
else
    if (size(align2d,2)~=length(sel.textdata))
        error(['Picklist and sel-file differ in size!']);
    end;
    frmaeList=0;
end;

if (frmaeList)
    filename_old=strrep(align2d(1,1).filename{1},'_corr/','/');
    filename_old=[fileparts(filename_old) sumExt];
    numFrames=length(align2d(1,1).filename);
else
    filename_old=strrep(align2d(1,1).filename,'_corr/','/');
    numFrames=1;
end;

load([filename_old '.mat']);
fid=fopen(f_def_part_list,'wt');
for i=1:length(sel.textdata)
    if (frmaeList)
        al_count=floor(i/(numFrames+0.0001))+1;
        fr_count=mod(i,numFrames);
        if (fr_count==0)
            fr_count=numFrames;
        end;
        filename=strrep(align2d(1,al_count).filename{fr_count},'_corr/','/');
        filename=[fileparts(filename) sumExt];
    else
        al_count=i;
        filename=strrep(align2d(1,al_count).filename,'_corr/','/');
    end;
    if (strcmp(filename,filename_old)==0)
        load([filename '.mat']);
        filename=filename_old;
    end;
    dz=st_out.Fit.Dz_det.*1e10;
    part_n=sel.textdata{i};
    fprintf(fid,'%s %f\n',part_n,dz);
    if (mod(i,5000)==0)
        disp([num2str(i) ' particles processed ']);
    end;
end;
fclose(fid);



