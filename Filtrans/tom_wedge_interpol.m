function vol_out=tom_wedge_interpol(vol,wedge_angle,particle_angle)


%try it in real space
vol=tom_rotate(vol,[270 90 90]);
sampling=10;

num=360./sampling;
ang=0;

for i=1:num
    ang=ang+sampling;
    
    vol_rot=tom_rotate(vol,[ang 0 0]);
    proj(:,:,i)=sum(vol_rot,1);
    figure; tom_dspcub(vol_rot);
    
end;

disp('test');


return;


% builde wedge
wedge=tom_wedge(ones(size(vol)),wedge_angle);
vol_out=zeros(size(vol));

%rotate to a more confortable position
wedge=tom_rotate(wedge,[0 0 90]);
vol=tom_rotate(vol,[0 0 90]);

% hike through all z-slices

for i=1:size(vol,3)
    %go to polar coordinate system
    im_tmp=fftshift(fftn(double(vol(:,:,i))));
    im_tmp=tom_cart2polar(im_tmp);
    im_tmp_wedge=tom_cart2polar(wedge(:,:,i));
    im_tmp_wedge_bin=(im_tmp_wedge>0);
   
    for ii=1:size(im_tmp,1) %loop over all slices
         %find start stop
         line=strrep(num2str(im_tmp_wedge_bin(ii,:)),' ','');
         line_org=im_tmp(ii,:);
         if (isempty(findstr(line,'0'))==0)  %loop over all lines
            tmp=findstr(line,'01');
            stop(1)=tmp(1);
            stop(2)=tmp(2);
            stop_interp(1)=stop(1)+30;
            stop_interp(2)=stop(2)+30;
            if (stop_interp(1) > size(im_tmp,2))
                stop_interp(1)=size(im_tmp,2);
            end;
             if (stop_interp(2) > size(im_tmp,2))
                stop_interp(2)=size(im_tmp,2);
            end;
            
            tmp=findstr(line,'10');  
            start(1)=tmp(1)-0;
            start(2)=tmp(2)-0;
            
            start_interp(1)=start(1)-30;
            start_interp(2)=start(2)-30;
            if (start_interp(1)  < 1)
                start_interp(1)=1;
            end;
             if (start_interp(2) < 1)
                start_interp(2)=1;
            end;
            
            %interpolate
            new_line=line_org;
            
            start_mean(1)=mean(im_tmp(ii,start_interp(1):start(1)));
            start_mean(2)=mean(im_tmp(ii,start_interp(2):start(2)));
            
            stop_mean(1)=mean(im_tmp(ii,stop(1):stop_interp(1)));
            stop_mean(2)=mean(im_tmp(ii,stop(2):stop_interp(2)));
            
            vals=interpol(start_mean(1),stop_mean(1),stop(1)-start(1));
            new_line(start(1):stop(1))=vals(1,:);
            vals=interpol(start_mean(2),stop_mean(2),stop(2)-start(2));
            new_line(start(2):stop(2))=vals(1,:);
            im_tmp_out(ii,:)=new_line; 
         else
             im_tmp_out(ii,:)=line_org; 
         end;
        
         disp('murat');
    end;
    im_tmp_out=tom_polar2cart(im_tmp_out);
    vol_out(:,:,i)=real(ifftn(ifftshift(im_tmp_out)));
    clear('im_tmp_out');
end;


function out=interpol(val1,val2,length)



point1=[0 val1];
point2=[length val2];

m=(point2(2)-point1(2))./(point2(1)-point1(2));
t=val1;
x=0:length;

out=m.*x+t;










