function tom_av2_xmipp_ml3d2proj_match_scrt_adapt(proj_match_scrt,target_folder,use_classes)
%TOM_AV2_XMIPP_ML3D2PROJ_MATCH_SCRT_ADAPT adapts proj_match_scrt
%
%  tom_av2_xmipp_ml3d2proj_match_scrt_adapt(in_scrt_name,output_folder)
%
%  TOM_AV2_XMIPP_ML3D2PROJ_MATCH_SCRT_ADAPT  xmipp-projmatch
%  after ml3d classification
%  
%  
%
%PARAMETERS
%
%  INPUT
%   proj_match_scrt     filename of the xmipp proj_match script
%   target_folder       outputfolder of tom_av2_xmipp_ml3d2proj_match
%   use_classes         (opt) class numbers 2 adapt (default 'all')
%
%EXAMPLE
%    
% tom_av2_xmipp_ml3d2proj_match_scrt_adapt('xmipp_protocol_projmatch.py','/fs/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_low/');
%
%REFERENCES
%
%SEE ALSO
%  tom_xmippdocread,tom_av2_xmipp_ml3d2proj_match
%
%NOTE
% 
% !!! path for target_folder shoul be absolute !!!!
% !!!abs path for target_folder shoul be loner than 70 characters (>70 xmipp crashes with seg fault in Reconstruction Module) !!!! 
% % 
% % proj_match_scrt rules (start volume, ctf_dat and ctf-group files have to have the following names !)
% 
%       *************************************************  
%       *  start volume    ==>   start.spi              *
%       *  doc file name   ==>   my_doc.doc (expert o.) *      
%       *  ctf_dat         ==>   all_images.ctfdat      *
%       *  defoucus split  ==>   split_ctf.doc          *
%       *************************************************
%
%   created by fb 
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron ponography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

%get info from outputfolder

dd=dir([target_folder '/model_*']);
num_of_classes=length(dd);

disp([num2str(num_of_classes) ' classes in ' target_folder ' found!']);

if nargin<3
    use_classes=1:num_of_classes;
end;

for i=use_classes
    
    out_basepath=[target_folder '/model_' num2str(i) ];
    
    disp(' '); disp('Adapting xmipp script... ' ); 
    
    
    [ps_p ps_n ps_e]=fileparts(proj_match_scrt);
    ps_name=[out_basepath '/' ps_n ps_e];
    
    
    [a line]=unix(['grep "SelFileName=''" ' proj_match_scrt]);
    line=line(1:end-1);
    idx=strfind(line,'''');
    scrt_sel_name=line(idx(1)+1:idx(2)-1);
    new_sel=['model_' num2str(i) '.sel'];
    
    [a line]=unix(['grep "DocFileName=''" ' proj_match_scrt]);
    line=line(1:end-1);
    idx=strfind(line,'''');
    scrt_doc_name=line(idx(1)+1:idx(2)-1);
    new_doc=['model_' num2str(i) '.doc'];
    
    [a line]=unix(['grep "ProjectDir=''" ' proj_match_scrt]);
    if (isempty(line))
        [a line]=unix(['grep -m 1 "ProjectDir=" ' proj_match_scrt]);
    end;
    line=line(1:end-1);
    idx=strfind(line,'''');
    if (isempty(idx))
        idx=strfind(line,'"');
    end;
    scrt_proj_name=line(idx(1)+1:idx(2)-1);
    if (strcmp(out_basepath(1),'/'))
        new_proj_path=out_basepath;
    else
        new_proj_path=[pwd '/' out_basepath];
    end;
    
    if (length(new_proj_path) > 80)
        disp(['target_folder path: ' target_folder] );
       % error(['target_folder length: ' num2str(length(new_proj_path)) ' path too long !!']);
    end;
    
    call=['awk ' '''{' ...
        'gsub("' scrt_sel_name '","' new_sel '");' ...
        'gsub("' scrt_doc_name '","' new_doc '");' ...
        'gsub("' scrt_proj_name '","' new_proj_path '");' ...
        'print }''' ' > ' ps_name ' ' proj_match_scrt];
    [a b]=unix(call);
    unix(['chmod ugo+rwx ' ps_name]);
    disp([proj_match_scrt ' ==> ' ps_name]); 
    disp(' ');
    
    [a b]=unix(['head -30 ' ps_name ' | tail -20']);
    disp(b);
    [a b]=unix(['head -65 ' ps_name ' | tail -20']);
    disp(b);
    
end;
    
    