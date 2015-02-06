function tom_av2_align_hightolow_list_beate2(base_path1,ext1,base_path2,ext2,filename_list,flag)
%TOM_AV2_ALIGN_HIGHTOLOW_LIST_BEATE generates a text file for aligning 
% focal pairs
%
%   tom_av2_align_hightolow_list_beate(folder1,folder2,filename_list,flag)
%
%  tom_av2_align_hightolow_list_beate 
%  generates a text file for aligning focal pairs ..for tom_av2_align_hightolow_beate
%
%PARAMETERS
%
%  INPUT
%   base_path1             basepath of the images 
%   ext1                   extension of the images
%   base_path2             (opt.)basepath of the images
%                           only needed for flag   2folders  
%   ext2                   (opt.) extension of the images
%                            only needed for flag   2folders 
%   filename_list          filename of the created high 2 low list
%   flag                   flag describing the input data
%                          (high_odd,high_even, or 2folders )
%  OUTPUT
%
%EXAMPLE
%     
%   tom_av2_align_hightolow_list_beate('/fs/scratch/fbeck/test_htl_al/high_low/test_','.em','','','list.txt','high_odd');
%   %creates the the list.txt file for all micrographs in 1 folder (abs path output)
%   %the odd micrographs are the high-defocused images !!!!
%
%   tom_av2_align_hightolow_list_beate('high_low/test_','.em','','','list.txt','high_odd');
%   %creates the the list.txt file for all micrographs in 1 folder (rel path output)
%   %the odd micrographs are the high-defocused images !!!! 
%
%   tom_av2_align_hightolow_list_beate('high/test_','.em','low/test_','.em','list.txt','2folders');
%   %creates the the text file for micrographs in 2 folders 
%   %note the high-defocus images have to be the first input parameters!!!!!!!!!!!!! 
%
%   tom_av2_align_hightolow_list_beate('high_low/test_','.em','','','list.txt','high_odd');
%   %creates the the list.txt file for all micrographs in 1 folder
%   %the odd micrographs are the high-defocused images !!!!
%
% NOTE
% 
% format of the output text-file:
%     
% path of the high image  blank  path of the low image
%
%  !cat list.txt
% 
%  /fs/scratch/fbeck/test_htl_al/high_low/test_1.em /fs/scratch/fbeck/test_htl_al/high_low/test_2.em
%  /fs/scratch/fbeck/test_htl_al/high_low/test_3.em /fs/scratch/fbeck/test_htl_al/high_low/test_4.em
%  /fs/scratch/fbeck/test_htl_al/high_low/test_5.em /fs/scratch/fbeck/test_htl_al/high_low/test_6.em
%
%   
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_align_hightolow_beate
%
%   created  by fb
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


fp=fopen(filename_list,'wt');

if (isempty(base_path2) && isempty(ext2) && (strcmp(flag,'high_even') || strcmp(flag,'high_odd') ) )
    d1=dir([base_path1 '*' ext1]);
    
    vect=zeros(length(d1),1);
    for i=1:length(d1)
        [a b]=strtok(d1(i).name,'_');
        num=strrep(strtok(b,'.'),'_','');
        vect(i)=str2double(num);
    end;
    vect=sort(vect);
    
    if (strcmp(flag,'high_even'))
        zz=1;
        for i=1:length(d1)/2
            fprintf(fp,[base_path1 num2str(vect(zz)) ext1 ' ' base_path1 num2str(vect(zz+1)) ext1 '\n']);
            zz=zz+2;
        end;
    end;
    
    if (strcmp(flag,'high_odd'))
        zz=1;
        for i=1:length(d1)/2
            fprintf(fp,[base_path1 num2str(vect(zz+1)) ext1 ' ' base_path1 num2str(vect(zz)) ext1 '\n']);
            zz=zz+2;
        end;
    end;
    
else
    d1=dir([base_path1 '*' ext1]);
    d2=dir([base_path2 '*' ext2]);
    
    
    
    if (length(d1)==0)
        error(['no images with ' base_path1 '*' ext1 ' found']);
    end;
    
    if (length(d2)==0)
        error(['no images with ' base_path2 '*' ext2 ' found']);
    end;
    
    
    vect1=zeros(length(d1),1);
    for i=1:length(d1)
        [a b]=strtok(d1(i).name,'_');
        num=strrep(strtok(b,'.'),'_','');
        vect1(i)=str2double(num);
    end;
    vect1=sort(vect1);
    
    vect2=zeros(length(d2),1);
    for i=1:length(d2)
        [a b]=strtok(d2(i).name,'_');
        num=strrep(strtok(b,'.'),'_','');
        vect2(i)=str2double(num);
    end;
    vect2=sort(vect2);
    
    zz=1;
    for i=1:length(vect1)
        if (sum(ismember(vect2,vect1(i)))>0 )
            fprintf(fp,[base_path1 num2str(vect1(zz)) ext1 ' ' base_path2 num2str(vect1(zz)) ext2 '\n']);
            zz=zz+1;
        end;
    end;
    
    
end;


fclose(fp);




