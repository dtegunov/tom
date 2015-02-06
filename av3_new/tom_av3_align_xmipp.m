function aligned_vol=tom_av3_align_xmipp(ref,vol,rot_v,tilt_v,psi_v,trans_v,mask,opt,tmp_dat)
%tom_av3_align_xmipp is a interface 2 xmipp_align_volumes
% 
%   aligned_vol=tom_av3_align_xmipp(ref,vol,rot_v,tilt_v,psi_v,trans_v,opt,tmp_dat)
%
%PARAMETERS
%
%  INPUT
%   ref                reference
%   vol                vol 2 be aligned
%   rot_v              vector of Euler-angles to be scanned (start stop step)
%   tilt_v             vector of Euler-angles to be scanned (start stop step)
%   psi_v              vector of Euler-angles to be scanned (start stop step)
%   trans_v            (-2 2 2) mask volume inside the ccf is calculated (start stop step)
%   mask               mask in memory 
%   opt                ('local') use local optimizer istead of exhaustive search  
%   tmp_dat            mask volume applied to the ccf
%
%  
%  OUTPUT
%   aligned_vol       aligned volume 
% 
%
%EXAMPLE
% 
%   test_vol=(tom_spheremask(ones(32,32,32),4,0,[22 22 16])+tom_cylindermask(ones(32,32,32),6))>0;
%   test_vol(:,:,1:6)=0; test_vol(:,:,end-6:end)=0;
%   test_vol_trans=tom_move(tom_rotate(test_vol,[2 3 6]),[-2 3 1]);
%   mask=tom_spheremask(ones(32,32,32),15);
%   
%   cc=tom_corr(test_vol,test_vol_trans,'norm');
%   [a b]=tom_peak(cc)  
%
%   aligned_vol=tom_av3_align_xmipp(test_vol,test_vol_trans,[-8 8 2],[-8 8 2],[-8 8 2],[-6 6 2],mask);
%   
%   cc=tom_corr(test_vol,aligned_vol,'norm');
%   [a b]=tom_peak(cc)
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av3_align
%
%   created by FB 11/21/12
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


if (nargin < 6)
    trans_v=[-2 2 2];
end;

if (nargin < 7)
    mask='';
end;

if (nargin < 8)
    opt='local';
end;

if (nargin < 9)
    tmp_dat='xxx_4alg_xxx';
end;


refHD=[tmp_dat '_ref.vol'];
volHD=[tmp_dat '_vol.vol'];
maskHD=[tmp_dat '_mask.vol'];

try
tom_spiderwrite(refHD,ref);
tom_spiderwrite(volHD,vol);
if (isempty(mask)==0)
    tom_spiderwrite(maskHD,mask);
end
catch Me
    disp(Me.message);
    error('cannot write tmp files ... check permissions!');
end;


call_base=['xmipp_align_volumes -i1 ' refHD ' -i2 ' volHD ' ' ];
call_param='';
call_param=[call_param '-rot ' my_num2str(rot_v)  ' -tilt ' my_num2str(tilt_v) ' -psi ' my_num2str(psi_v) ' '];
call_param=[call_param '-x ' my_num2str(trans_v) ' -x ' my_num2str(trans_v) ' -z ' my_num2str(trans_v) ' '];
call_param=[call_param ' -apply '];

if (isempty(opt)==0)
    call_param=[call_param '-' opt ' '] ;
end;

if (isempty(mask)==0)
    call_param=[call_param '-mask ' maskHD];
end;

call=[call_base call_param];

disp(call);
unix(call);

aligned_vol=tom_spiderread(volHD);
aligned_vol=aligned_vol.Value;

call='rm xxx_4alg_xxx*.vol';
disp(call);
unix(call);


function str_out=my_num2str(v)

str_out='';
for i=1:length(v)
    str_out=[str_out ' ' num2str(v(i))];
end;

str_out=str_out(2:end);
