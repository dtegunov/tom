function tmpstack=tom_av2_xmipp_picklist2stack(picklist,particle_name,sel_file_name,radius,norm,invert,binning,start,mask,split_floder,align_flag,rad_fact,on_error,randperm,ctf_corr,field_corr,xray,taper_flag,log_file_name)
%TOM_AV2_XMIPP_PICKLIST2STACK transforms picklist 2 xmipp stack
%and creates a .sel file
%
%   tom_av2_xmipp_picklist2stack(picklist,particle_name,sel_file_name,norm,invert,bin,start,mask,split_folder,align_flag,rad_fact,on_error)
%
%PARAMETERS
%
%  INPUT
%   picklist           picklist filename 
%   particle_name      basename of the particles                   
%   sel_file_name      filename of sel file with extension!
%   radius             radius of the particle
%   norm               ('gradient&mean0+1std') gradient/statistics or use                                             
%                      'local_12&mean0+1std'    to subtract a background with a filter of 12
%                      'gradient&mean0+1std-bg' normalizes in respect 2 background spherical mask with radius of sz/2 is used    
%   invert             (1) 0 or 1
%   binning            (0) 0, 1, ... to do a rescal use a vect like [128 128]
%   start              (1) particle starting number
%   mask               ('no_mask') binary mask in em format (size has 2 be rad*2)           
%   split_folder       ('off') particles are written in multiple folders
%   align_flag         (0) write parts aligned
%   rad_fact           (1) cut out particles bigger 2 avoid crop edges
%   on_error           ('error') throws an error and stops
%                      'rand' go on and unse rand vals 4 the corresponding particle 
%   randperm           (0) permutes pixels which are outside of mask (works only if a mask is given)      
%   ctf_corr           ('none') flag for "on-the-fly" ctf correction (flipping has to be done!)
%                               use 'flip' to 4 phase flipping
%   field_corr         ('none') flag for  "on-the-fly" gain and dark correction  
%                               use field_corr{1}='gain.em' % is multiplied   
%                                   field_corr{2}='dark.em' % is subtrakted
%   xray               (0)  use 1 for "standard" xray (4.5 std) for K2 frames use value bigger than 1 as cut off    
%   taper_flag         (1) taper particles if 2 close 2 the micrograph edge 
%   log_file_name      (picklist '.log') 
%                       if picklist is given in memory and no log_file_name is specified no logfile is written   
%
%  OUTPUT
%    tmpstack          (opt.) stack in memory  
%
%EXAMPLE
%  
% tom_av2_xmipp_picklist2stack('picklist.mat','parts_','26S.sel',80,'gradient&mean0+1std',1,0,1)
%
% tom_av2_xmipp_picklist2stack('bad.mat','/fs/pool/pool-nickell3/fb/4Tueb/test_mult/pick_all/comb_results/stacks/bad/parts_','bad.sel',128,'gradient&mean0+1std',1,2,1,'no_mask',30000);
%
% %read picklist in memory ...use only for small pickLists 
%  
% stack=tom_av2_xmipp_picklist2stack('picklist.mat','','',128,'gradient&mean0+1std',1,2);
%
% %read MULTIPLE picklist in memory and rescale ...use only for small pickLists 
% 
% stack=tom_av2_xmipp_picklist2stack('pickLists/pick_high_*.mat','','',280,'gradient&mean0+1std',1,[64 64]);
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_bin_sel,tom_xmippsellread
%
%   created by SN 01/08/07
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


if nargin<5
    norm='gradient&mean0+1std';
end;

if nargin<6
    invert=1;
end;

if nargin<7
    binning=0;
end;

if nargin<8
    start=1;
end;

if nargin<9
   mask='no_mask'; 
end;

if nargin<10
   split_floder='off'; 
end;

if nargin<11
   align_flag=0; 
end;

if nargin<12
    rad_fact=1;
end;

if nargin<13
    on_error='error';
end;

if nargin<14
    randperm=0;
end;

if nargin<15
    ctf_corr='none';
end;

if nargin<16
    field_corr='none';
end;

if nargin<17
    xray=0;
end

if nargin < 18
   taper_flag=1;
end

if nargin < 19
    if (isstruct(picklist))
        log_file_name='';
    else
        log_file_name=[picklist '.log'];
    end;
end;
    
if (randperm==1 && strcmp(mask,'no_mask')==1)
    error('mask is needed for randperm falg !!!');
end;

if (split_floder==1)
    disp('sorry split flag not implemented so far try splitting pickList instead!');
    error('split flag not implemented so far!!');
end;


my_disp(' ');
my_disp(['Running on: ' num2str(matlabpool('size')) ' processors'],'info: ',1,log_file_name,1);

[binning,particle_name,sel_file_name]=transformInputs2Cell(binning,particle_name,sel_file_name);
[grad_filt,norm_str]=parseNormInput(norm);
mask_backGround=genBackgroundMask(norm,radius);
align2d=loadPickList(picklist,radius,log_file_name);
[align2dReOrg,lenAlign2dReorg]=reorgAlign2d(align2d,align_flag);
outfilestruct=buildOutFileStruct(particle_name,split_floder,radius,binning,nargout,lenAlign2dReorg,sel_file_name);
show_param(picklist,particle_name,sel_file_name,radius,norm,invert,binning,start,mask,split_floder,align_flag,rad_fact,outfilestruct,on_error,randperm,ctf_corr,field_corr,xray,taper_flag,log_file_name);
mynargout=nargout;

dispHeadTail('head',log_file_name,lenAlign2dReorg);
parfor i=1:length(align2dReOrg)
%for i=1:length(align2dReOrg)
   [tmpstack{i},tmppos{i}]=processPickList(align2dReOrg{i},outfilestruct,radius,norm_str,grad_filt,mask_backGround,invert,start,mask,split_floder,rad_fact,taper_flag,on_error,randperm,ctf_corr,field_corr,xray,mynargout,log_file_name);
end;

align2d2Sel(align2d,outfilestruct,'wt');
tmpstack=cell2PartStack(tmpstack,tmppos);

dispHeadTail('tail',log_file_name);



function stack=cell2PartStack(tmpstack,tmppos)

stack=[];
abspos=[];
for i=1:length(tmpstack)
    stack=cat(3,stack,tmpstack{i});
    abspos=cat(2,abspos,tmppos{i});
end;
tmpOrdered(:,:,abspos)=stack(:,:,1:end);
stack=tmpOrdered;
clear('tmpOrdered');




function [binning,particle_name,sel_file_name]=transformInputs2Cell(binning,particle_name,sel_file_name)

if (iscell(binning)==0)
    tmp=binning;
    clear('binning');
    binning{1}=tmp;
end;
if (iscell(particle_name)==0)
    tmp=particle_name;
    clear('particle_name');
    particle_name{1}=tmp;
end;

if (iscell(sel_file_name)==0)
    tmp=sel_file_name;
    clear('sel_file_name');
    sel_file_name{1}=tmp;
end;


function dispHeadTail(flag,log_file_name,numP)

if (strcmp(flag,'head'))
    my_disp('****************************************','info: ',0,log_file_name);
    my_disp('starting Extraction ','info: ',1,log_file_name,1);
    my_disp('****************************************','info: ',0,log_file_name);
    my_disp(['Total Num Of Particles: ' num2str(numP)],'',1);
    fprintf(1,'Number of Processed Particles:        ');
end;

if (strcmp(flag,'tail'))
    fprintf(1,'\n');
    my_disp('','',0,log_file_name);
    my_disp('****************************************','info: ',0,log_file_name);
    my_disp('finished Extraction ','info: ',1,log_file_name,1);
    my_disp('****************************************','info: ',0,log_file_name);
end;


function dispProgess(outfilestruct,log_file_name)

if (isempty(outfilestruct(1).fid)==0 )
     [a b]=unix(['wc -l ' outfilestruct(1).fid ' | awk ''{print $1} ''']);
     vv=sprintf('\b\b\b\b\b\b\b%6d',round(str2double(b)));
     my_disp([num2str(round(str2double(b))) ' particles processed '],'info: ',0,log_file_name,1);
     disp(vv);
     drawnow('update');
end;


function align2d2Sel(align2d,outfilestruct,openFlag)

if (isempty(outfilestruct(1).fid)==0)
    for is=1:length(outfilestruct)
        zz=0;
        for i=1:size(align2d,2)
            pName=[outfilestruct(is).filename num2str(i)];
            if (iscell(align2d(1,i).filename))
                for ii=1:length(align2d(1,i).filename)
                    zz=zz+1;
                    all_names{zz}=[outfilestruct(is).path '/' pName strrep(outfilestruct(is).ext,'.','-') '/' pName '_F_' num2str(ii) outfilestruct(is).ext];
                end;
            else
                zz=zz+1;
                all_names{zz}=[outfilestruct(is).path '/' pName outfilestruct(is).ext];
            end;
            
        end;
        writeSel(all_names,outfilestruct(is).fid,openFlag);
        clear('all_names');
    end;
end;

function [tmpstack,posVect]=processPickList(align2dReOrg,outfilestruct,radius,norm_str,grad_filt,mask_backGround,invert,start,mask,split_floder,rad_fact,taper_flag,on_error,randperm,ctf_corr,field_corr,xray,valback,log_file_name)

GainAndDark=loadGainAndDark(field_corr);

%Main loop over picklist
particleBuffer=-1;
img_buffer.name='';

for i=1:length(align2dReOrg)
    %transfer Variables
    imageName=align2dReOrg{i}{1};
    imageNameSum=align2dReOrg{i}{3};
    particlePosition=align2dReOrg{i}{2};
    particleTransform=align2dReOrg{i}{4};
    sumIdx=align2dReOrg{i}{5}+(start-1);
    frameIdx=align2dReOrg{i}{6};
    
    img_buffer=fillImageBuffer(imageName,img_buffer,GainAndDark,ctf_corr,[imageNameSum '.mat'],xray,log_file_name);
    particle=cutOutParticle(img_buffer,particlePosition,radius,rad_fact,taper_flag,on_error,log_file_name);
    particle=processParticle(particle,imageName,particlePosition,particleTransform,radius,rad_fact,grad_filt,xray,invert,on_error,log_file_name);
    particleBuffer=writeParticle(particle,particleBuffer,outfilestruct,sumIdx,frameIdx,length(align2dReOrg),mask,norm_str,randperm,mask_backGround,split_floder,valback,xray,log_file_name);
end;


drawnow('update');
dispProgess(outfilestruct,log_file_name);
drawnow('update');

if valback > 0 
    tmpstack=rescaleNormMask(particleBuffer.stack,outfilestruct(1).binning{1},mask,norm_str,randperm,mask_backGround,xray);
    posVect=particleBuffer.absPosInStack;
else
    tmpstack=[];
    posVect=[];
end;





function [list,lenTotal]=reorgAlign2d(align2d,alignFlag)
zz=1;
disp('ReOrganizing Picklist');
for i=1:size(align2d,2)
    if (iscell(align2d(1,i).filename))
        for is=1:length(align2d(1,i).filename)
            all_filenames{zz}=clean_path(align2d(1,i).filename{is});
            all_pos{zz}=[align2d(1,i).position.x align2d(1,i).position.y];
            all_Orgfilenames{zz}=clean_path(align2d(1,i).Orgfilename);
            if (alignFlag)
                all_AlignParam{zz}=[align2d(1,i).angle align2d(1,i).shift.x align2d(1,i).shift.y];
            else
                all_AlignParam{zz}=[0 0 0];
            end;
            
            all_PartIdx{zz}=i;
            all_FrameIdx{zz}=is;
            zz=zz+1;
        end;
    else
         all_filenames{zz}=clean_path(align2d(1,i).filename);
         all_pos{zz}=[align2d(1,i).position.x align2d(1,i).position.y];
         all_Orgfilenames{zz}=clean_path(align2d(1,i).filename);
         if (alignFlag)
            all_AlignParam{zz}=[align2d(1,i).angle align2d(1,i).shift.x align2d(1,i).shift.y];
         else
            all_AlignParam{zz}=[0 0 0];
         end;
         
         all_PartIdx{zz}=i;
         all_FrameIdx{zz}=-1;
         zz=zz+1;
    end;
end;


part_count=0;
u_filenames=unique(all_filenames);
for i=1:length(u_filenames)
    idx=find(ismember(all_filenames,u_filenames{i}));
    for ii=1:length(idx)
        tmpList{ii}={all_filenames{idx(ii)} ; all_pos{idx(ii)} ; all_Orgfilenames{idx(ii)}; all_AlignParam{idx(ii)} ; all_PartIdx{idx(ii)}; all_FrameIdx{idx(ii)}};
        part_count=part_count+1;
    end;
    list{i}=tmpList;
    clear('tmpList');
end;

lenTotal=part_count;
disp(['ReOrganized PickList==> total length of: ' num2str(lenTotal)]);
disp(['ReOrganized PickList==> number of micrographs: ' num2str(length(u_filenames))]);

function pathClean=clean_path(pathNotClean)

pathClean=regexprep(pathNotClean,'/fs[0-9][0-9]/lv[0-9][0-9]/','/');


function img_buffer=fillImageBuffer(img_name,img_buffer,GainAndDark,ctf_corr,cftFileName,xray,log_file_name)
  
if (strcmp(ctf_corr,'none')==0 || (isempty(GainAndDark{1})==0 && isempty(GainAndDark{2})) )
    if (strcmp(img_buffer.name,img_name)==0)
        try
            my_disp(['Reading: ' img_name],'info: ',0,log_file_name,1);
            img_buffer.img=tom_emreadc(img_name);
            [imgPath,imgFileName,imgExt]=fileparts(img_name);
        catch Me
            my_disp(['error reading file ' img_name ' ' Me.message],'error: ',1,log_file_name,1);
            error('Micrograph read error!!');
        end;
        img_buffer.img=img_buffer.img.Value;
        img_buffer.img=correctImage(img_buffer.img,img_name,GainAndDark,ctf_corr,cftFileName,xray,log_file_name);
        img_buffer.size=size(img_buffer.img);
        img_buffer.name=img_name;
    end;
else
    try
        header = tom_reademheader(img_name);
    catch Me
        
        disp(['error reading file ' img_name ' ' Me.message]);
        error('Micrograph read error!!');
    end;
    img_buffer.img=img_name;
    img_buffer.size=header.Header.Size;
    img_buffer.name=img_name;
end;



function img=correctImage(img,img_name,GainAndDark,ctf_corr,ctf_filename,xray,log_file_name)

img=single(img);
sz_img=size(img);
cutSz=(sz_img(1)./2)-(0.05.*sz_img(1));
mesg4log=['Correcting: ' img_name ];

if (length(GainAndDark)>1)
    if (isempty(GainAndDark{2})==0)
        mesg4log=[mesg4log ' -dark '];
        img=img-GainAndDark{2};
    end;
end;

if (isempty(GainAndDark{1})==0)
    mesg4log=[mesg4log ' .*gain '];
    img=img.*GainAndDark{1};
end;
if (xray==1)
    mesg4log=[mesg4log ' xray: 4.6std '];
    img=tom_xraycorrect2(img);
end;
if (xray > 1)
    mesg4log=[mesg4log ' xray: >' num2str(xray) ' '];
    img=tom_xraycorrect2(img,2,xray);
end;


if (strcmp(ctf_corr,'none')==0)
    try 
        load(ctf_filename);
        [cftPathtmp,ctfFilenametmp,ctfExttmp]=fileparts(ctf_filename);
    catch Me
        ctfErrorMsg=['error loading ' img_name '.mat ' Me.message];
        my_disp(ctfErrorMsg,'error: ',1,log_file_name);
        error(ctfErrorMsg);
    end;
    mesg4log=[mesg4log ' ' ctf_corr ' dz: ' num2str(st_out.Fit.Dz_det.*1e6) ' dzd: ' num2str(st_out.Fit.Dz_delta_det.*1e6) ' dzdAng: ' num2str(st_out.Fit.Phi_0_det) ' using: ' ctfFilenametmp '.mat'];
    img=tom_correct_for_ctf_and_mtf_new(img,st_out.Fit,ctf_corr,cutSz,'');    
end;

my_disp(mesg4log,'info: ',0,log_file_name,1);

function [particleBuffer]=writeParticle(particle,particleBuffer,outfilestruct,sumIdx,frameIdx,lenList,mask,norm_str,randperm,mask_backGround,split_flag,nargout,xray,log_file_name)

sz_part=size(particle);

if (isnumeric(particleBuffer))
    clear('particleBuffer');
    %initalize particle Buffer ...pseudo constructor :-)
    particleBuffer.stack=zeros(sz_part(1),sz_part(2),lenList);
    particleBuffer.count=0;
    for i=1:length(outfilestruct)
        particleBuffer.binning{i}=outfilestruct(i).binning;
    end;
end;

particleBuffer.count=particleBuffer.count+1;
particleBuffer.stack(:,:,particleBuffer.count)=particle;
particleBuffer.absPosInStack(particleBuffer.count)=sumIdx;
if (nargout < 1)
    
    for i=1:length(outfilestruct)
        if (frameIdx > 0)
            pName=[outfilestruct(i).filename num2str(sumIdx)];
            outputNames=[outfilestruct(i).path '/' pName strrep(outfilestruct(i).ext,'.','-') '/' pName '_F_' num2str(frameIdx) outfilestruct(i).ext];
        else
            outputNames=[outfilestruct(i).path '/' outfilestruct(i).filename num2str(sumIdx) outfilestruct(i).ext];
        end;
        particleBuffer.output{particleBuffer.count,i}=outputNames;
    end;
    
    if ((particleBuffer.count >= outfilestruct(1).fileBufferSize) || particleBuffer.count==lenList )
        my_disp(['flushing buffer sz: ' num2str(particleBuffer.count) ],'info: ',0,log_file_name,1);
        particleBuffer.stack=particleBuffer.stack(:,:,1:particleBuffer.count);
        for is=1:length(outfilestruct)
            tmpstackProc=rescaleNormMask(particleBuffer.stack,outfilestruct(is).binning,mask,norm_str,randperm,mask_backGround,xray);
            writeStack(tmpstackProc,particleBuffer.output(:,is));
            writeSel(particleBuffer.output(:,is),outfilestruct(is).fid);
        end;
        clear('particleBuffer');
        particleBuffer=-1;
    end;
   
end;



function writeSel(partNames,selName,openFlag)

if (nargin < 3)
    openFlag='a';
end;

for i=1:10
    fid = fopen(selName,openFlag);
    if (fid > 0)
        break;
    end;
    pause(0.3);
end;

for i=1:length(partNames)
    fprintf(fid,'%s 1\n',partNames{i});
end;

fclose(fid);


function show_param(picklist,particle_name,sel_file_name,radius,norm,invert,binning,start,mask,split_folder,align_flag,rad_fact,outfilestruct,on_error,randperm,ctf_corr,field_corr,xray,taper_flag,log_file_name)

%my_disp(message,status,out_screen,out_fileName

logStat='info: ';

my_disp(' ','',1,log_file_name);
my_disp('======================================>',logStat,1,log_file_name);
my_disp(['modus: ' outfilestruct(1).BufferMode],logStat,1,log_file_name);
if (isstruct(picklist))
    my_disp(['picklist: in memory'],logStat,1,log_file_name);
else
    my_disp(['picklist: ' picklist],logStat,1,log_file_name);
end;
if (iscell(particle_name))
    for is=1:length(particle_name)
         my_disp(['particle_name-' num2str(is)  ': '  particle_name{is}],logStat,1,log_file_name);
    end;
else
    my_disp(['particle_name: ' particle_name],logStat,1,log_file_name);
end;
if (iscell(particle_name))
    for is=1:length(particle_name)
        my_disp(['sel_file_name-' num2str(is) ': '   sel_file_name{is}],logStat,1,log_file_name);
    end;
else
    my_disp(['sel_file_name: ' sel_file_name],logStat,1,log_file_name);
end;

my_disp(['radius: ' num2str(radius)],logStat,1,log_file_name);
my_disp(['norm: ' norm],logStat,1,log_file_name);
my_disp(['invert: ' num2str(invert)],logStat,1,log_file_name);

if (iscell(binning))
    for is=1:length(binning)
        my_disp(['binning/rescale-' num2str(is) ': ' num2str(binning{is}) ],logStat,1,log_file_name);
    end;
else
    my_disp(['binning/rescale-' num2str(is) ': ' num2str(binning) ],logStat,1,log_file_name);
end;

my_disp(['start: ' num2str(start)],logStat,1,log_file_name);
if (isnumeric(mask)==0)
    my_disp(['mask: ' mask],logStat,1,log_file_name);
else
    my_disp('mask: user binary',logStat,1,log_file_name);
end;
if (isnumeric(split_folder)==0)
    my_disp(['split_folder: ' split_folder],logStat,1,log_file_name);
else
    my_disp(['split_folder: ' num2str(split_folder)],logStat,1,log_file_name);
end;
my_disp(['align_flag: ' num2str(align_flag)],logStat,1,log_file_name);
my_disp(['rad_fact: ' num2str(rad_fact)],logStat,1,log_file_name);
my_disp(['Buffer Size: ' num2str( outfilestruct(1).fileBufferSize)],logStat,1,log_file_name);
my_disp(['on error: ' on_error],logStat,1,log_file_name);
my_disp(['permute pixels outside mask: ' num2str(randperm)],logStat,1,log_file_name);
my_disp(['ctf_corr: ' ctf_corr],logStat,1,log_file_name);
if (iscell(field_corr))
    my_disp(['field_corr gain: ' num2str(field_corr{1})],logStat,1,log_file_name);
    if (length(field_corr)>1)
        my_disp(['field_corr dark: ' num2str(field_corr{2})],logStat,1,log_file_name);
    end;
else
    my_disp(['field_corr: ' num2str(field_corr)],logStat,1,log_file_name);
end;
my_disp(['xray: ' num2str(xray)],logStat,1,log_file_name);

my_disp(['taper_flag: ' num2str(taper_flag)],logStat,1,log_file_name);
my_disp(['log_file_name: ' num2str(log_file_name)],logStat,1,log_file_name);

my_disp('<======================================',logStat,1,log_file_name);
my_disp(' ','',1,log_file_name);


function particleStackProc=rescaleNormMask(particleStack,binning,mask,norm_str,randperm,mask_backGround,xray)

sz_org=size(particleStack);
new_sz=[sz_org(1) sz_org(2) size(particleStack,3)];

if (size(binning,2)>1)
    rescale=binning;
    binning=0;
else
    rescale=-1;
end;

if (rescale~=-1)
    new_sz=[rescale(1) rescale(2) size(particleStack,3)];
end;
    
if (binning > 0)
    new_sz=round([sz_org(1)./(2^binning) sz_org(2)./(2^binning) size(particleStack,3)]);
end;
    
particleStackProc=zeros(new_sz);


if (strcmp(mask,'no_mask')==0)
    if (rescale(1)~=-1)
        mask=imresize(mask,rescale)>0.08;
    else
        mask=tom_bin(mask,binning)>0.08;
    end;
end;
if (isempty(mask_backGround)==0)
    if (rescale(1)~=-1)
        mask_backGround=imresize(mask_backGround,rescale)>0.08;
    else
        mask_backGround=tom_bin(mask_backGround,binning)>0.08;
    end;
end;

for ip=1:size(particleStack,3)
    
    particle=particleStack(:,:,ip);
    %rescale particle and masks
    if (rescale(1)~=-1 || binning > 0)
        if (rescale(1)==-1)
            particle=tom_bin(particle,binning);
        else
            particle=tom_rescale(particle,rescale);
        end;
        if (xray==1)
            particle=tom_xraycorrect2(particle);
        end;
        if (xray > 1)
            particle=tom_xraycorrect2(particle,2,xray);
        end;
        
     end;
    
    
    %norm particle
    if strcmp(norm_str,'no_norm')==0
        if strcmp(mask,'no_mask')==0
            if (randperm==0)
                if (isempty(mask_backGround))
                    particle = tom_norm((particle+5).*2,norm_str,mask).*mask;
                else
                    particle = tom_norm(particle,norm_str,mask_backGround).*mask;
                end;
            else
                particle = tom_norm((particle+5).*2,norm_str);
            end;
        else
            if (isempty(mask_backGround))
                particle = tom_norm((particle+5).*2,norm_str);
            else
                particle = tom_norm((particle+5).*2,norm_str,mask_backGround);
            end;
        end;
    end;
    
    %mask is not set ones if no mask ...SPEED!
    if strcmp(mask,'no_mask')==0
        if (randperm==0)
            particle = particle.*mask;
        else
            particle = tom_permute_bg(particle,mask);
        end;
    end;
   
    particleStackProc(:,:,ip)=particle;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Process Particle                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function particle=processParticle(particle,img_name,pos,transForm,radius,rad_fact,grad_flag,xray,invert,on_error,log_file_name)

shift=transForm(1:2);
angle=transForm(3);
x=pos(1);
y=pos(2);
try
    
    %align particle
    if (sum(shift==[0 0])~=2) || (angle~=0)
        particle = tom_rotate(particle,angle);
        shift=round(shift./(2^binning));
        particle = tom_shift(particle,shift);
    end;
    
    %cut out particle to final size
    if (rad_fact > 1)
        sz_part=size(particle);
        try
            cut_start=round((sz_part-(2*radius))./2);
            cut_stop=cut_start+(radius)*2;
            particle = particle(cut_start(1):cut_stop(1)-1,cut_start(2):cut_stop(2)-1);
        catch Me
            particle = tom_taper(particle,[((radius)*rad_fact*2) ((radius)*rad_fact*2)]);
            particle = particle(radius:3*radius-1,radius:3*radius-1);
        end;
    end;
    
    
    if (size(particle,1) < 2*(radius) || size(particle,2) < 2*(radius) )
        particle=tom_taper(particle,[2*radius  2*radius]);
    end;
    
    %filter gradient
    if (grad_flag==1)
        particle = tom_xmipp_normalize(particle,'Ramp');
    end;
    if (grad_flag > 1)
        particle = particle-tom_filter(particle,grad_flag);
    end;
    
    if (invert==1)
        particle=-particle;
    end;
    
    if (xray==1)
        particle=tom_xraycorrect2(particle);
    end;
    if (xray > 1)
        particle=tom_xraycorrect2(particle,2,xray);
    end;
    
    particle=single(particle);
    
    
    
catch Me
    errorMsg=['error processing particle: ' img_name ' position: ' num2str(x) ' ' num2str(y) ' ' Me.message];
    if (strcmp(on_error,'error'))
        my_disp(errorMsg,'error: ',1,log_file_name,1);
        error(errorMsg);
    else
        my_disp(errorMsg,'error: ',1,log_file_name,1);
        my_disp(['using rand 4 particle'],'error: ',1,log_file_name,1);
        clear('particle');
        particle.Value = rand(2*radius,2*radius);
    end;
  
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  check cutout                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inside_flag,start_cut,size_cut,dir_flag]=check_cutout(x,y,radius,rad_fact,imsize)

inside_flag=1;



start_cut(1)=x-round(rad_fact*radius);
start_cut(2)=y-round(rad_fact*radius);
start_cut(3)=1;

size_cut(1)=2.*round(rad_fact*radius)-1;
size_cut(2)=2.*round(rad_fact*radius)-1;
size_cut(3)=0;

dir_flag=[0 0];


if (start_cut(1) < 1)
    inside_flag=0;
    size_cut(1)=size_cut(1)+start_cut(1);
    start_cut(1)=1;
    dir_flag(1)=-1;
end;

if (start_cut(2) < 1)
    inside_flag=0;
    size_cut(2)=size_cut(2)+start_cut(2);
    start_cut(2)=1;
    dir_flag(2)=-1;
end;

if (x+rad_fact*radius > imsize(1))
    inside_flag=0;
    size_cut(1)=round(rad_fact*radius)+(imsize(1)-x)-1;
    dir_flag(1)=1;
end;

if ((y+rad_fact*radius> imsize(2)) )
    inside_flag=0;
    size_cut(2)=round(rad_fact*radius)+(imsize(2)-y)-1;
    dir_flag(2)=1;
end;



function outfilestruct=buildOutFileStruct(particle_name,split_floder,radius,binning,nargout,lenAlign,sel_file_name)

[FileBufferSize,ReadMode]=calcFileBufferSize(radius,lenAlign,nargout);

if (nargout==1)
    outfilestruct(1).BufferMode='memory';
    outfilestruct(1).fileBufferSize=FileBufferSize;
    outfilestruct(1).binning=binning;
    outfilestruct(1).fid='';
    return;
end;

warning off;
if (iscell(particle_name))
    for is=1:length(particle_name)
        [a,b,c]=fileparts(particle_name{is});
        partDirName{is}=a;
        partFileName{is}=b;
        partExt{is}=c;
    end;
else
    [a,b,c]=fileparts(particle_name);
    partDirName{1}=a;
    partFileName{1}=b;
    partExt{1}=c;
end;
%mkdir(a);

for is=1:length(partDirName)
    superDir=fileparts(partDirName{is});
    if (exist(superDir,'dir')==0 )
        mkdir(superDir);
    end;
    if (exist(partDirName{is},'dir')==0 )
        mkdir(partDirName{is});
    end;
    outfilestruct(is).method='singleFiles';
    outfilestruct(is).path=partDirName{is};
    outfilestruct(is).filename=partFileName{is};
    outfilestruct(is).ext='.spi';
    outfilestruct(is).fileformat='spider';
    outfilestruct(is).size=[];
    if (strcmp(split_floder,'off'))
        outfilestruct(is).folder_num_offset=-1;
        outfilestruct(is).num_of_entries=5000000000000000;
    else
        outfilestruct(is).num_of_entries=split_floder;
        outfilestruct(is).folder_num_offset=-floor((start-1)./outfilestruct.num_of_entries)+1;
    end;
    outfilestruct(is).fileBufferSize=FileBufferSize;
    outfilestruct(is).targetSize=calcTargetSize(radius,binning{is});
    outfilestruct(is).binning=binning{is};
    outfilestruct(is).BufferMode=ReadMode;
    
    if (nargout<1)
        if (iscell(sel_file_name))
            fid=sel_file_name{is};
        else
            fid=sel_file_name;
        end;
        
        [a,b]=unix(['rm ' fid]);
        if (fid==-1)
             error(['error opening file for writing: ' sel_file_name{is}]);
        end; 
    end;
    outfilestruct(is).fid=fid;
end;
warning on;

if (lenAlign>50000 && outfilestruct.num_of_entries==5000000000000000)
    warning('picklist bigger > 50000');
    disp('consider using split_flag with 50000');
end;






function align2d=loadPickList(picklist,radius,log_file_name)

try 
    if (isstruct(picklist))
        align2d=picklist;
        picklist='in memory';
    else
        if (exist(picklist,'file'))
            load(picklist);
        else
            align2d=tom_av2_cat_pickLists(picklist);
        end;
    end;
catch Me
    errorM=['Cannot load picklist: ' picklist];
    my_disp(errorM,'error: ',1,log_file_name);
    error(errorM);
end;

if (isempty(radius))
    radius=align2d(1,1).radius;
end;

logStat='info';

my_disp(['Picklist contains: ' num2str(size(align2d,2)) ' particles'],logStat,1,log_file_name);
if (iscell(align2d(1,1).filename));
    my_disp(['Picklist contains: ' num2str(length(align2d(1,1).filename)) ' frames for each particle'],logStat,1,log_file_name);
end;

function [bufferSz,r_mod]=calcFileBufferSize(radius,lenAlign,nargout)

if (nargout<1)
    bufferSz=round((64./radius)^2.*20000);
    r_mod='standard';
else
    bufferSz=lenAlign;
    r_mod='memory';
end;

if (bufferSz > lenAlign)
    bufferSz=lenAlign;
end;


function targetSize=calcTargetSize(radius,binning)

if (length(binning) > 1)
    targetSize=binning;
else
    targetSize=[(radius.*2)./2^binning (radius.*2)./2^binning];
end;



function packages=update_packages(packages,idx_written,idx)

if isempty(idx_written)
    return;
end;

v=zeros(length(idx_written),1);

for i=1:length(idx_written)
    v(i)=find(idx==idx_written(i));
end;



for i=1:size(packages,1)
    new_start=max( (packages(i,1) < v(:)) .* (v(:) < packages(i,2)) .* (v(:)) );
    packages(i,1)=new_start;
    packages(i,3)=packages(i,2)-packages(i,1)+1;
end;



function [grad_filt,norm_str]=parseNormInput(norm)

if (findstr(norm,'gradient'))
    grad_filt=1;
end;

if (findstr(norm,'local_'))
    [a b]=strtok(norm,'&');
    [a b]=strtok(a,'_');
    grad_filt=str2num(strrep(b,'_',''));
end;

if (findstr(norm,'&'))
   [a b]=strtok(norm,'&');
   norm_str=strrep(b,'&','');
end;

if (findstr(norm_str,'-bg'))
   norm_str=strrep(norm_str,'-bg','');
end;

function mask_backGround=genBackgroundMask(norm,rad)
st_sz=[rad rad].*2;
if (findstr(norm,'-bg'))
    mask_backGround=tom_spheremask(ones(st_sz(1),st_sz(2)),round(st_sz(1)./2)-1)==0;
else
    mask_backGround='';
end;


function GainAndDark=loadGainAndDark(field_corr)

if (strcmp(field_corr,'none')==0)
    for iGAD=1:length(field_corr)
        try
            tmpImg=tom_emreadc(field_corr{iGAD});
        catch Me
            error(['Could not read: ' field_corr{iGAD} ' ' Me.message]);
        end;
        
        GainAndDark{iGAD}=tmpImg.Value;
    end;
    if (length(GainAndDark)==1)
        GainAndDark{2}='';
    end;
else
    GainAndDark{1}='';
    GainAndDark{2}='';
end;



function particle=cutOutParticle(img_buffer,pos,radius,rad_fact,taper_flag,on_error,log_file_name)

x=pos(1);
y=pos(2);
imagesize=img_buffer.size;
img_name=img_buffer.name;
part_sz=[(radius.*rad_fact).*2 (radius.*rad_fact).*2];

%check coordinates
[inside_flag,start_cut,size_cut,dir_flag]=check_cutout(x,y,radius,rad_fact,imagesize);
if ((inside_flag==0 && taper_flag==0) || x < 0 || y < 0 || x > imagesize(1) || y > imagesize(2))
    errorMsg=['error: particle coordinates are out of range: ' img_name ' position: ' num2str(x) ' ' num2str(y) ' ' ];
    if (strcmp(on_error,'error'))
        my_disp(errorMsg,'error: ',1,log_file_name,1);
        error(errorMsg);
    else
        my_disp([errorMsg ' unsing Rand Particle: '],'error: ',1,log_file_name,1);
    end;
    particle=rand(part_sz(1),part_sz(2),'single');
    return;
end;

%cut out particle
try
    clear('particle');
    if (isnumeric(img_buffer.img))
         particle.Value = tom_cut_out(img_buffer.img,start_cut,size_cut+1);
    else
         particle = tom_emreadc(img_name,'subregion',start_cut,size_cut);
    end;
catch Me
    msgTmp=['error reading subregion: ' img_name ' position: ' num2str(x) ' ' num2str(y) ' ' Me.message];
    if (strcmp(on_error,'error'))
        my_disp(msgTmp,'error: ',1,log_file_name);
        error(msgTmp);
    else
        my_disp(msgTmp,'error: ',1,log_file_name);
        my_disp(['using rand 4 particle'],'info: ',1,log_file_name);
        clear('particle');
        particle.Value = rand(part_sz(1),part_sz(2));
    end;
end;
particle = single(particle.Value);

if (sum(size(particle)==((round(radius*rad_fact)*2)))~=2)
    particle=tom_taper2(particle,[part_sz(1) part_sz(2)],dir_flag);
end;




function sucess_vector=writeStack(stack,names)

sucess_vector=ones(size(stack,3),1);
for i=1:size(stack,3)
    basename=fileparts(names{i});
    if (exist(basename,'dir')==0)
        mkdir(basename);
    end;
    sucess=save_write(names{i},stack(:,:,i),'spider');
    if (sucess==0)
       pause(1);
       sucess(i)=positions(i);
    end;
end;

   
function sucess=save_write(name_string,out,fileformat)

sucess=0;
for i=1:10
    try
        switch fileformat
            case 'em'
                tom_emwrite(name_string,out);
            case 'spider'
                tom_spiderwrite(name_string,out);
            otherwise
                error('Fileformat not implemented !!!');
        end;
        sucess=1;
        return;
    catch
        pause(3);
    end;
end;



function my_disp(message,status,out_screen,out_fileName,date)

if (nargin < 2)
    status='';
end;

if (nargin < 3)
    out_screen=1;
end;

if (nargin < 4)
    out_fileName='';
end;

if (nargin <5)
    date=0;
end;

if (date==0)
    mdate='';
else
    mdate=[datestr(now) ' '];
end;

Cmessage=[status mdate message];

if (out_screen==1)
    if (strfind(status,'error'))
        disp(' ');
    end;
    disp(Cmessage);
    if (strfind(status,'error'))
        disp('Number of Processed Particles:        ');
        pause(2);
    end;
 end;

if (isempty(out_fileName)==0)
    fid=fopen(out_fileName,'a');
    fprintf(fid,'%s\n',Cmessage);
    fclose(fid);
end;





