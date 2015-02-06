function  [tmpstack back_in_memory_idx]=tom_av2_createstack3(inalignname,outalign,outfilestruct,binning,radius,rad_fact,norm_str,grad_filt,align_flag,taper_flag,mask,start_nr,num_of_packages,file_buffer_size,gappy_write,restart_flag,verbose)
%TOM_AV2_CREATESTACK3 creates a particle stack from an alignment file
%
%    tom_av2_createstack3(inalignname,outalign,outfilestruct,binning,radius,rad_fact,norm_str,align_flag,taper_flag,mask,num_of_packages,file_buffer_size,gappy_write,restart_flag,verbose)
%
%PARAMETERS
%
%  INPUT
%   inalignname         full path to input alignment file (.mat)
%   outalign            output align struct
%   outfilestruct       structure for fileoutput
%                       outfilestruct.method             file structure '1stack' 'Nstacks' 'singleFiles'  
%                       outfilestruct.num_of_entries     number of entries
%                       outfilestruct.path               path for output             
%                       outfilestruct.filename           filename 
%                       outfilestruct.ext                extension     
%                       outfilestruct.fileformat         used fileformat (em,spider,imagic)
%                                               
%
%
%   binning             binning
%   radius              radius the particles should be cut out
%   rad_fact            factor for cutting the particles ...for writing
%                       aligned particles without cropping
%   norm_str            string for norming the particle ('mean0+1std','phase'...) 
%   grad_filt           filters gradient
%   align_flag          1: save particles aligned, 0: save unaligned (optional)
%   taper_flag          flag for tapering particles
%   mask                mask in memory use 'no_mask' not ones() for no mask ...speed!
%   start_nr            start_nr for particle index
%   num_of_packages     nuber of packages for paraell computing 
%   file_buffer_size    size of the filebuffer hackster use back_in_memory to get stack in memory instead of writiting 2 HD      
%   gappy_write         0: no jumps in the numbering sceme 
%                       1: position in the inalignment file is equal position                             
%                       2: position in the inalignment file is equal
%                          position gap are filled with zeros 
%   restart_flag        restarts writing procedure 
% 
%  OUTPUT
%  tmpstack           stack output in memory ...set  file_buffer_size 2
%                     'back_in_memory' hackster
%
%EXAMPLE
%  
%
% outfilestruct.method='singleFiles'; 
% outfilestruct.num_of_entries=10000; 
% outfilestruct.path='part_fold';
% outfilestruct.folder_num_offset=1;
% outfilestruct.filename='part_';
% outfilestruct.ext='.em';
% outfilestruct.fileformat='em';
% 
% [st outfilestruct]=tom_os3_calc_meanCC(align2d,outfilestruct,0.46);
%
% tom_av2_createstack3('pick_small.mat','',outfilestruct,2,160,1,'mean0+1std',1,0,1,'no_mask',4);
%
%
%
%
%
%
%
%
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER TOM_AV2_CHECK_WRITTEN
%
%   cre0ted by FB 06/16/08
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





if nargin < 15 || isempty(verbose)
    verbose=0;
end;

 

h=waitbar_wrapper(1,1,'reading alignment file ...',verbose,'');

if ~isstruct(inalignname)
    try
        s = load(inalignname);
    catch
        error('Could not load alignment file.');
    end
    alignstruct = s.align2d;
else
    alignstruct = inalignname;
end;
clear('s');




if nargin < 3 || isempty(outfilestruct)
    outfilestruct.method='singleFiles'; 
    outfilestruct.num_of_entries=10000;
    outfilestruct.path='part_fold';
    outfilestruct.filename='part_';           
    outfilestruct.ext='.em';
    outfilestruct.fileformat='em';
end;

if nargin < 2 || isempty(outalign)
    outalign_filename=[outfilestruct.path '_align'];
end;



if nargin < 4 || isempty(binning)
    binning=0;
end;

if nargin < 5 || isempty(radius)
    radius = alignstruct(1,1).radius;
end;

if nargin < 6 || isempty(rad_fact)
    rad_fact=1;
end;

if nargin < 7 || isempty(norm_str)
    norm_str='mean0+1std';
end;

if nargin < 8 || isempty(grad_filt)
    grad_filt=0;
end;

if nargin < 9 || isempty(align_flag)
    align_flag=0;
end;

if nargin < 10 || isempty(taper_flag)
    taper_flag=1;
end;

if nargin < 11 || isempty(mask)
    mask='no_mask';
end;

if nargin < 12 || isempty(start_nr)
    start_nr=1;
end;

if nargin < 13 || isempty(num_of_packages)
    num_of_packages=1;
end;

if nargin < 14 || isempty(file_buffer_size)
    file_buffer_size=750;
end;

if nargin < 15 || isempty(gappy_write)
    gappy_write=1;
end;

if nargin < 16  || isempty(restart_flag)
    restart_flag=0;
end;



%transfer often used variable
if (isfield(outfilestruct,'idx')==0) 
    %warning('Index missing ...syc index to align st !');
    outfilestruct.idx.written_index=1:size(alignstruct,2);
end;

idx=outfilestruct.idx.written_index;

if (num_of_packages > 1 && gappy_write==0)
    disp('set num_of_packages to 1');
    error('Paralell and no gappy write is not implemented ! ');
end;



radius_bin=round(radius./2^binning);
packages=tom_calc_packages(num_of_packages,length(idx));

if (gappy_write==0)
    new_align=alignstruct;
end;


if (restart_flag==1)
    st_idx=tom_av2_check_written(outfilestruct);
    packages=update_packages(packages,st_idx.written_index,idx);
end;



for ii=1:num_of_packages
%parfor ii=1:num_of_packages

    partcounter = 1;
    stackcounter = 1;
    al_shift=[0 0];
    al_rot=0;
    used_part=packages(ii,1);
    
    if strcmp(file_buffer_size,'back_in_memory')
        tmpstack = zeros(round((radius./2^binning)*2),round((radius./2^binning)*2),length(alignstruct),'single');
        back_in_memory_idx=zeros(length(alignstruct),1);
    else
        tmpstack = zeros(radius_bin*2,radius_bin*2,file_buffer_size,'single');
    end

    for i=packages(ii,1):packages(ii,2)
        
       

        %h=waitbar_wrapper(packages(ii,2),i,'cutting particles ',verbose,h);

        %check if image is readable
         try
            header = tom_reademheader(alignstruct(1,idx(i)).filename);
        catch
            disp(['error reading file ' alignstruct(1,idx(i)).filename]);
            continue;
        end
        imagesize = header.Header.Size;

        x = alignstruct(1,idx(i)).position.x;
        y = alignstruct(1,idx(i)).position.y;

        
        use_all_fl=1;
        if (use_all_fl==1)
            if (x < 0)
                x=1;
            end;
            if (y < 0)
                y=1;
            end;
            if x > imagesize(1)
                x=imagesize(1) -1;
            end;
            if y > imagesize(2)
                y=imagesize(2)-1;
            end;
         end;
        
        
        %check coordinates
        [inside_flag,start_cut,size_cut,dir_flag]=check_cutout(x,y,radius,rad_fact,imagesize);
        if ((inside_flag==0 && taper_flag==0) || x < 0 || y < 0 || x > imagesize(1) || y > imagesize(2))
            disp([alignstruct(1,idx(i)).filename ' skipped. Coordinates out of range']);
            disp(['x: ' num2str(x) '  y: ' num2str(y) ]);
            continue;
        end;

        %cut out particle
        try
            particle = tom_emreadc(alignstruct(1,idx(i)).filename,'subregion',start_cut,size_cut,'binning',binning);
        catch
            disp(['error cutting particle: ' alignstruct(1,idx(i)).filename ' position: ' num2str(x) ' ' num2str(y)]);
            continue;
        end;
        particle = -single(particle.Value);


        %process particle
        if (align_flag==1)
            al_shift=[alignstruct(1,idx(i)).shift.x alignstruct(1,idx(i)).shift.y];
            al_rot=alignstruct(1,idx(i)).angle;
        end;

        try
            particle=process_particle(particle,al_shift,al_rot,radius,rad_fact,binning,mask,norm_str,grad_filt,dir_flag);
        catch
            disp(['error processing particle: ' alignstruct(1,idx(i)).filename ' position: ' num2str(x) ' ' num2str(y)]);
            continue;
        end;

        tmpstack(:,:,partcounter)=particle;

        if (gappy_write==1)
            used_part_idx(partcounter)=idx(i);
        else
            used_part_idx(partcounter)=packages(ii,1)+used_part-1;
        end;


        partcounter = partcounter + 1;
        used_part=used_part+1;

        if strcmp(file_buffer_size,'back_in_memory')
            back_in_memory_idx(i)=1;
        end;
        
        
        %partcounter
        if (strcmp(file_buffer_size,'back_in_memory')==0)
            if (partcounter > file_buffer_size || (i==packages(ii,2) && partcounter > 2) )
                if (i==packages(ii,2) && partcounter > 2)
                    tmpstack = tmpstack(:,:,1:partcounter-1);
                end;
                %h=waitbar_wrapper(packages(ii,2),i,'writing buffered files...',verbose,h);
                error_pos=tom_av2_stack_write(tmpstack,outfilestruct,used_part_idx+(start_nr-1));
                partcounter = 1;
                used_part_idx=[];
                stackcounter = stackcounter + file_buffer_size;
            end
        end;
 
        if(mod(i,2000)==0)
               disp(num2str(i));
         end;
    
    end; % end of for i

    
    
end; % end of parfor


if (strcmp(file_buffer_size,'back_in_memory')==0)
    
    %save align
    if (gappy_write==1)
        align2d=alignstruct; clear('alignstruct');
        [idx align2d]=tom_av2_check_written(outfilestruct,align2d);
        save(outalign_filename,'align2d');
    else
        idx=tom_av2_check_written(outfilestruct);
        new_align=new_align(1,1:used_part);
        save(outalign_filename,'new_align');
    end;
    
    
    outfilestruct.idx=idx;
    outfilestruct.size=[radius_bin*2 radius_bin*2 length(outfilestruct.idx.written_index)];
    save(outfilestruct.path,'outfilestruct');
    
    waitbar_wrapper(packages(1,2),packages(1,2)+2,'writing buffered files...',verbose,h);
else
    tmpstack=tmpstack(:,:,1:(partcounter-1));
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Process Particle                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function particle=process_particle(particle,shift,angle,radius,rad_fact,binning,mask,norm_str,grad_flag,dir_flag)


if (sum(size(particle)==((round(radius*rad_fact)*2)./2^binning))~=2)
    rad_t=round(radius./2^binning);
    particle=tom_taper2(particle,[round(rad_t*rad_fact)*2 round(rad_t*rad_fact)*2],dir_flag);
end;

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
        cut_start=round((sz_part-(2*radius/2^binning))./2);
        cut_stop=cut_start+(radius/2^binning)*2;
        particle = particle(cut_start(1):cut_stop(1)-1,cut_start(2):cut_stop(2)-1);
    catch
        particle = tom_taper(particle,[((radius/2^binning)*rad_fact*2) ((radius/2^binning)*rad_fact*2)]);
        particle = particle(radius/2^binning:3*radius/2^binning-1,radius/2^binning:3*radius/2^binning-1);
    end;
end;

%mask is not set ones if no mask ...SPEED!
if strcmp(mask,'no_mask')==0
    particle = particle.*mask;
end

if (size(particle,1) < 2*(radius./2^binning) || size(particle,2) < 2*(radius./2^binning) )
    particle=tom_taper(particle,[2*radius  2*radius]);
end;


%filter gradient
if (grad_flag==1)
     particle = tom_xmipp_normalize(particle,'Ramp');
end;

%norm particle
if strcmp(norm_str,'no_norm')==0
    if strcmp(mask,'no_mask')==0
        particle = tom_norm((particle+5).*2,norm_str,mask);
    else
        particle = tom_norm((particle+5).*2,norm_str);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Waitbar wrapper%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=waitbar_wrapper(max,progress,message,verbose,h)

if (verbose==0)
    return;
end;

if (progress > max)
    close(h);
    return;
end;

if (isempty(h))
    h = waitbar(progress./max,message);
end;


h = waitbar(progress./max,h,message);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update_h_array%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_array=update_h_array(h_array,count,iswritten,isaligned,normstr)
    

h_array(1,count)=iswritten;


if nargin < 4
    return; %hack for speed !
end;

h_array(2,count)=isaligned;

switch normstr
        case 'phase'
            normed_num=1;
        case '2std'
            normed_num=2;
        case '3std'
            normed_num=3;
        case 'mean0+1std'
            normed_num=4;
        otherwise
            error('norm method not implemeted');
end;


h_array(3,count)=normed_num;



function [align2d idx]=update_align2d(align2d,restart_flag,outfilestruct,align_flag,norm_str,gappy_write)

for i=1:size(align2d,2)
    
    if (restart_flag==0)
        align2d(1,i).iswritten=h_array(1,i);
        align2d(1,i).isaligned=h_array(2,i);
        switch h_array(3,i)
            case 1
                hstr='phase';
            case 2
                hstr='2std';
            case 3
                hstr='3std';
            case 4
                hstr='mean0+1std';
            otherwise
                error('norm method not implemeted');
        end;
        align2d(1,i).normed=hstr;
    else
        align2d(1,i).isaligned=align_flag;
        align2d(1,i).normed=norm_str;
    end;
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


