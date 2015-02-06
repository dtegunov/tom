function tom_av2_match_lib(f_stack,lib_path,output_doc,psi_incre,max_sh,filt_ker,num_of_nodes,chunk_size)
%TOM_AV2_MATCH_LIB matches a particle stack against a proj-lib 
% matching is done in cart coordinates ...useful 4 pre centering
% inplane rotation has 2 be sampled
%
%   
%
%PARAMETERS
%
%  INPUT
%   f_stack               *.sel (xmipp stack)
%   lib_path              input proj lib
%   output_doc            output_doc name of the output .doc
%   psi_incre             psi increment
%   max_sh                max shift 4 == radius 4 ccf mask 
%   filt_ker              value for tom_filter    
%   num_of_nodes          number of processor  
%   chunk_size            number of parts that can be n times loaded + lib 
%  
%  
%  OUTPUT
%
%
%EXAMPLE
%
% tom_av2_match_lib('parts.sel','/fs/pool/pool-titan1/CWT/com_data/pre_center/lib_low_10deg_2/ref','doc_cent.doc',10,100,2,8,100);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

%load sel

show_param(f_stack,lib_path,output_doc,psi_incre,max_sh,filt_ker,num_of_nodes,chunk_size);

sel=importdata(f_stack);

packages=tom_calc_packages(num_of_nodes,length(sel.textdata));

parfor i=1:size(packages,1)
    tmp_sel=sel.textdata(packages(i,1):packages(i,2));
    result{i}=proc_pack(tmp_sel,lib_path,chunk_size,max_sh,filt_ker,i);
end;

all_results=[];
for i=1:length(result)
    all_results=cat(1,all_results,result{i});
end;

%writ into doc-file
doc_st=tom_av2_xmipp_empty_doc(f_stack);
load([lib_path '_ang.mat']);


for i=1:length(doc_st)
    org_ref=floor(all_results(i,2)./psi_incre);
    doc_st(i).ref=org_ref;
    doc_st(i).xoff=all_results(i,3);
    doc_st(i).yoff=all_results(i,4);
    doc_st(i).rot=ref_ang(all_results(i,2),3);
    doc_st(i).tilt=ref_ang(all_results(i,2),4);
    doc_st(i).psi=ref_ang(all_results(i,2),5);
    doc_st(i).maxCC=all_results(i,1);
end;

disp(' ');
disp(['writing doc: ' output_doc]);
tom_xmippdocwrite(output_doc,doc_st);
disp('finished');

function result=proc_pack(partlist,lib_path,chunk_size,max_sh,ker,run_num)


pack_str=['pack Nr :' num2str(run_num) ' '];

%load lib
disp([pack_str 'loading ' lib_path]);
proj_lib=load_lib(lib_path,ker);

packages=tom_calc_packages(ceil(length(partlist)./chunk_size),length(partlist));

mid=floor([size(proj_lib,1) size(proj_lib,2)]./2)+1;

%build mask
num_of_pix=(size(proj_lib,1).*size(proj_lib,1));
mask_cc=tom_spheremask(ones(size(proj_lib,1),size(proj_lib,2)),max_sh);
idx_mask_cc=find(mask_cc==0);

result=zeros(length(partlist),4);
zz=1;
flag_outp=0;
for pack_nr=1:size(packages,1)
    
    %load parts in buffer
    idx=packages(pack_nr,1):packages(pack_nr,2);
    %part_buffer=zeros(size(proj_lib,1),size(proj_lib,1),length(idx));
    part_buffer=zeros(size(proj_lib,1),size(proj_lib,1),length(idx))+ j*ones(size(proj_lib,1),size(proj_lib,1),length(idx));
    disp([pack_str 'loading particle buffer ' num2str(packages(pack_nr,1)) ' to ' num2str(packages(pack_nr,2)) ]);
    for i=1:length(idx)
        tmp_part=tom_spiderread(partlist{idx(i)});
        %part_buffer(:,:,i)=tom_norm(tom_filter(tmp_part.Value,ker),'mean0+1std');
        part_buffer(:,:,i)=tom_fourier(tom_norm(tom_filter(tmp_part.Value,ker),'mean0+1std'));
    end;
   
   
    for i=1:size(part_buffer,3)
         
         all_sh=zeros(2,size(proj_lib,3));
         all_val=zeros(1,size(proj_lib,3));
        for ii=1:size(proj_lib,3)
            %cc=tom_corr(proj_lib(:,:,ii),part_buffer(:,:,i));
            cc=real(ifftshift(tom_ifourier(part_buffer(:,:,i).*proj_lib(:,:,ii) ))) / num_of_pix;
            cc(idx_mask_cc)=-1;
            [pos val]=tom_peakc(cc);
            all_sh(:,ii)=mid-pos(1,1:2);
            all_val(ii)=val(1);
        end;
        [mv mp]=max(all_val);
        result(zz,1)=mv;
        result(zz,2)=mp;
        result(zz,3:4)=all_sh(1:2,mp);
        tmp_co = double(int16(zz*10/length(partlist)));
        if tmp_co > flag_outp
            flag_outp = tmp_co;
            disp([' pack nr: ' num2str(run_num) ' ' num2str(flag_outp*10) ' % done']);
        end;
         zz=zz+1;
    end;
    
end;


function proj_lib=load_lib(lib_path,ker)

load([lib_path '_ang.mat']);

im_tmp=tom_spiderread([lib_path lead_zer(1) '.xmp']);
sz=size(im_tmp.Value);
%proj_lib_real=zeros(sz(1),sz(2),size(ref_ang,1));
proj_lib=zeros(sz(1),sz(2),size(ref_ang,1))+j*zeros(sz(1),sz(2),size(ref_ang,1));
for i=1:size(ref_ang,1)
    im_tmp=tom_spiderread([lib_path lead_zer(i) '.xmp']);
    proj_lib(:,:,i)=conj(tom_fourier( tom_norm(tom_filter(im_tmp.Value,ker),'mean0+1std')));
    proj_lib_real(:,:,i)=im_tmp.Value;
end;




function zer=lead_zer(num)

if (ischar(num))
    num=str2double(num);
end;

if (num<10)
    zer='00000';
end;

if (num>=10 && num<100)
    zer='0000';
end;

if (num>=100 && num<1000)
    zer='000';
end;

if (num>=1000 && num<10000)
    zer='00';
end;

if (num>=10000 && num<100000)
    zer='0';
end;

zer=[zer num2str(num)];


function show_param(f_stack,lib_path,output_doc,psi_incre,max_sh,filt_ker,num_of_nodes,chunk_size)


disp(' ');
disp('======================================>');
disp(['f_stack: ' f_stack]);
disp(['lib_path: ' lib_path]);
disp(['output_doc: ' output_doc]);
disp(['psi_incre: ' num2str(psi_incre)]);
disp(['max_sh: ' num2str(max_sh)]);
disp(['filt_ker: ' num2str(filt_ker)]);
disp(['num_of_nodes: ' num2str(num_of_nodes)]);
disp(['chunk_size: (in parts) ' num2str(chunk_size)]);
disp('<======================================');
disp(' ');


