function tom_av2_trace_particle(align2d,iter_num,interv,samp_phi,samp_the,filt_num)
%TOM_AV2_TRACE_PARTICLE creates ...
%
%   tom_av2_trace_particle(align2d,iter_num,interv,samp_phi,samp_the,filt_num)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   iter_num            ...
%   interv              ...
%   samp_phi            ...
%   samp_the            ...
%   filt_num            ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_trace_particle(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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


h=tom_reademheader(align2d(2,1).rec.file.Stack_Path);

size_stack=h.Header.Size;

stack=tom_emread(align2d(iter_num,1).rec.file.Stack_Path);

angular_scan=align2d(iter_num,1).rec.project.angular_scan;

avg_st_im=zeros(size_stack(1),size_stack(2),size(angular_scan,2));
avg_st_num=zeros(iter_num,size(angular_scan,2));

num2=zeros(size(samp_phi(1):samp_phi(2):samp_phi(3),2),size(samp_the(1):samp_the(2):samp_the(3),2));

avg_deluxe=zeros(size_stack(1),size_stack(1),size(samp_phi(1):samp_phi(2):samp_phi(3),2).*size(samp_the(1):samp_the(2):samp_the(3),2));
avg_deluxe_num=zeros(1,size(samp_phi(1):samp_phi(2):samp_phi(3),2).*size(samp_the(1):samp_the(2):samp_the(3),2) );


for i=interv(1):interv(2)
    %num(align2d(1,i).angleclass.angle(2)+181)=num(align2d(1,i).angleclass.angle(2)+181)+1;
    proj_nr=align2d(iter_num,i).rec.classify.angleclass.proj_nr;
    phi(i)=align2d(iter_num,i).rec.classify.angleclass.angle(1);
    theta(i)=align2d(iter_num,i).rec.classify.angleclass.angle(2);
    %num(align2d(1,i).angleclass.angle(2)+181)=num(align2d(1,i).angleclass.angle(2)+181)+1;
    ind1=( (phi(i)-samp_phi(1))./samp_phi(2) )+1;
    
    shift=align2d(iter_num,i).rec.classify.angleclass.shiftxy;
    rot=align2d(iter_num,i).rec.classify.angleclass.rot;
    
    ind2=( (theta(i)-samp_the(1))./samp_the(2) )+1;
    
    if (ind1==7 & ind2==1)
       
    end;
    
    ind_stack=(ind2-1).*size(samp_phi(1):samp_phi(2):samp_phi(3),2)+ind1;
    
    avg_deluxe(:,:,ind_stack)=avg_deluxe(:,:,ind_stack) + tom_shift(tom_rotate(stack.Value(:,:,i),rot),[shift]);
    avg_deluxe_num(ind_stack)=avg_deluxe_num(ind_stack)+1;
    
    if (ind2==0)
        disp('k');
    end;
    num2(ind1,ind2)=num2(ind1,ind2)+1;
    
    avg_st_im(:,:,proj_nr)=avg_st_im(:,:,proj_nr)+stack.Value(:,:,i);
    avg_st_num(proj_nr)=avg_st_num(proj_nr)+1;
    
    disp(i);
end;
avg=sum(stack.Value(:,:,interv(1):interv(2)),3);

%norm it norman

ttt=1;
for i=1:size(avg_st_im,3)
    if (avg_st_num(i)~=0)
        avg_st_im(:,:,i)=avg_st_im(:,:,i)./avg_st_num(i);
    end;
   
    
end;


h1=figure; tom_imagesc(num2); %set(h1,'title',[num2str(interv)]); 
%h2=figure; surfc(num2); shading interp; colormap hot; axis image;  %set(h2,'title',[num2str(interv)]); 
for i=1:size(avg_deluxe,3)
    if (avg_deluxe_num(i)> filt_num)
        avg_deluxe(:,:,i)=avg_deluxe(:,:,i)./avg_deluxe_num(i);
        avg_deluxe2(:,:,i)=tom_bin(avg_deluxe(:,:,i),0);
    end;
    

    
         

end;

immm=tom_dspcub2(avg_deluxe2,[size(samp_the(1):samp_the(2):samp_the(3),2) size(samp_phi(1):samp_phi(2):samp_phi(3),2)],[2 2],[2 2],0);
disp('tt'); tom_imagesc(immm);
%h3=figure; plot(sum(num2,1),'-');
%h4=figure; plot(sum(num2,2),'-');
%h5=figure; tom_imagesc(avg);
%h6=figure; tom_dspcub(avg_st_im);

