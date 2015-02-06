function [res,mat_out]=tom_compare_deluxe(vol1,vol2,pixs,crit,shells,mask,comp_map,use_vox,samp,disp_flag)
% TOM_COMPARE_DELUXE performs FSC 
%  
%     [res mat_out]=tom_compare_deluxe(vol1,vol2,pixs,crit,shells,mask,comp_map,use_vox,samp,disp_flag)
%  
%  TOM_COMPARE_DELUXE performs FSC overall or with a cernel mask as 3D
%  resolution map yes!!
%
%
%  PARAMETERS
%  
%    INPUT
%     vol1                volume1 or filename
%     vol1                volume2 or filename
%     pixs                pixelsize in Ang
%     critria             (0.5) threshold for resolution usually 0.5 or 0.3
%     shells              (size of vol ./2 ) number of shells
%     mask                mask to mesure resolution in a certain area ...use
%                         default to get spherical mask with radis= 0.125 * size of vol 
%     comp_map            (0) compute 3d resolution map flag 0/1
%     use_vox             tigh hard mask to reduce calculations  
%     samp                sampling 2 use every second voxel 
%     disp_flag           (1) flag for display
%                          
%    
%    OUTPUT
%     res                  mesured resolution in Ang or 3d res map
%     mat_out              matrix with mesured values x-mipp style !
%                           
%         mat_out(:,1) Resolution in 1./Ang 
%         mat_out(:,2) FRC 
%         mat_out(:,3) FRC-noise ...not implemented (-1) !!!!
%         mat_out(:,4) Resolution in Ang
%
%
%  EXAMPLE
%
%   [res data]=tom_compare_deluxe(im_w.Value,im_wo.Value,3.6,0.5);
%   %calc normal fsc ole
%
%   
%   matlabpool open local 8;
%   use_vox=tom_cylindermask(ones(128,128,128),40);
%   mask=tom_spheremask(ones(128,128,128),10,3); 
%  
%   [res data]=tom_compare_deluxe(vol1,vol2,4.42,0.5,64,mask,1,use_vox,2);
%   %rock 3d resolution map !! 
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_analyse_3d_res
%  
%     created by SN/FB 01/24/06
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom



if (ischar(vol1) || ischar(vol2))
    try
        vol1=tom_emreadc(vol1);
    catch
        vol1=tom_spiderread(vol1);
    end;
    vol1=vol1.Value;
    try
        vol2=tom_emreadc(vol2);
    catch
        vol2=tom_spiderread(vol2);
    end;
    vol2=vol2.Value;
end;

sz=size(vol1);


if (nargin<3)
    pixs=1;
end;

if (nargin<4)
    crit=0.5;
end;

if (nargin<5)
    shells=round(sz(1)./2);
end;

if (nargin<6)
    mask='';
end;

if (nargin<7)
    comp_map=0;
end;


if (nargin<8)
    disp_flag=1;
end;

if (strcmp(mask,'default'))
    mask=tom_spheremask(ones(size(vol1)),round(size(vol1,1)./8),3);
end;

if (isempty(mask))
    mask=ones(size(vol1));
end;


ny=pixs.*2;

if (comp_map==1)
    zz=1;
    list=zeros(size(vol1,1)*size(vol1,2)*size(vol1,3),3);
    test=zeros(size(vol1));
    disp(['Building list of Coordinates... ']);
    for ix=1:samp:size(vol1,1)
        for iy=1:samp:size(vol1,2)
            for iz=1:samp:size(vol1,3)
                if (use_vox(ix,iy,iz) > 0 )
                    list(zz,:)=[ix iy iz];
                    test(list(zz,1),list(zz,2),list(zz,3))=1;
                    zz=zz+1;
                end;
            end;
        end;
    end;
    l_tmp=list(1:zz-1,:); clear('list'); list=l_tmp;  clear('l_tmp');
    mid=floor(size(vol1)./2)+1;
    out_tmp=zeros(size(list,1),1);
    disp(['number of calculations: ' num2str(size(list,1)) ]);
    tic;
    parfor ii=1:size(list,1)
    % for ii=1:size(list,1)
    
        pos=list(ii,:);
        sh=pos-mid;
        
        mask_tmp=tom_move(mask,sh);
        mat=tom_compare(vol1.*mask_tmp,vol2.*mask_tmp,shells);
        
        num_of_sh=size(mat,1);
        
        out_tmp(ii)=calc_res(mat,ny,num_of_sh,crit);
        
%         if (out_tmp(ii)>50)
%             disp(num2str(ii));
%         end;

        
        if (mod(ii,1000)==0)
            disp([num2str(ii) ' done!']);
        end;
        
    end;
    vol_out=ones(size(vol1));
    mat=tom_compare(vol1.*mask,vol2.*mask,shells);
    num_of_sh=size(mat,1);
    [r mat_out]=calc_res(mat,ny,num_of_sh,crit);
    vol_out=vol_out.*mat_out(1,4); %initialize volume by min res
    for i=1:size(list,1)
        vol_out(list(i,1),list(i,2),list(i,3))=out_tmp(i);
    end;
    toc;
    res=vol_out;
    disp_flag=0;
    
else
    mat=tom_compare(vol1.*mask,vol2.*mask,shells);
    [res data]=calc_res(mat,ny,shells,crit);
end;




if (disp_flag)
    
    h=figure;
    
    p1=plot(data(:,1),data(:,2),'Linewidth',2); hold on;
    p2=plot(data(:,1),data(:,2),'ro','MarkerSize',3,'Linewidth',2);
    p=get(p1,'Parent');
    %set(fig,'Position',[80 12 60 20]);
    set(p,'Ylim',[0 1.1]);
    %set(t,'Color','white');
    xlabel('Resolution [1/Ang]');
    set(p,'Ytick',[0:0.1:1.1]);
    set(p,'Xtick',[0:data(end,1)./5:data(end,1)]);
    ylabel('FSC');
    grid on;
    p=get(p1,'Parent');
    xtick=get(p,'Xtick');
    xtick_nm=zeros(size(xtick));
    xtick_nm(2:end)=(1./xtick(2:end));
    
    for i=1:size(data(:,2),1)
        v_05=NaN;
        if data(i,2)<0.5
            v1=data(i-1,1);
            v2=data(i,1);
            v=(v1+v2)./2;
            %        text(v,0.5,sprintf('%0.3g A',1./v));
            text(0,0.5,sprintf('  %0.3g A >>>',1./v));
            v_05=v;
            break;
        end;
    end;
    
    for i=1:size(data(:,2),1)
        v_03=NaN;
        if data(i,2)<0.3
            v1=data(i-1,1);
            v2=data(i,1);
            v=(v1+v2)./2;
            %        text(v,0.3,sprintf('%0.3g A',1./v));
            text(0,0.3,sprintf('  %0.3g A >>>',1./v));
            v_03=v;
            break;
        end;
    end;
    
    
    for i=2:size(xtick,2)
        t=text(xtick(i),0.05,sprintf('%0.3g A',xtick_nm(i)));
    end;
    
    title(['Fourier Shell Correlation, Nyquist @ ' num2str(data(end,4)) ' Ang. FSC: 0.5 @ ' sprintf('%0.3g',1./v_05) ' A, 0.3 @ ' sprintf('%0.3g',1./v_03) ' A.']);
   
    
end;

if (comp_map==0)
    
    
    for i=1:size(data(:,2),1)
        if data(i,2)<crit
            v1=data(i-1,1);
            v2=data(i,1);
            v=(v1+v2)./2;
            % text(v,0.5,sprintf('%0.3g A',1./v));
            res=1./v;
            break;
        end;
    end;
    
    
    mat_out=data;
    
    
end;



function [v data]=calc_res(mat,ny,num_of_sh,crit)


for i=1:num_of_sh
    res=(num_of_sh./i).*ny;
    data(i,1)=1./res;
    if (isnan(mat(i,9)) )
        data(i,2)=0.999;
    else
        data(i,2)=mat(i,9);
    end
    data(i,3)=-1;
    data(i,4)=res;
end;

v=data(1,4); %ini with lowest res
for i=1:size(data(:,2),1)
    if data(i,2)<crit
        if (i>1)
            v1=data(i-1,4);
        else
            v1=data(i,4);
        end;
        v2=data(i,4);
        v=(v1+v2)./2;
        break;
    end;
end;







