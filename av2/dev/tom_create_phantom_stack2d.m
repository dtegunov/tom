function [stack align2d log]=tom_create_phantom_stack2d(filename,template,angles,sh,num,noise,mix,disp)

%parse inputs
narg=nargin;
switch narg
    case 3
        sh=[0 0];
        num=1;
        noise=[0 0];
        mix=0;
        disp='off';
    case 4
        num=1;
        noise=[0 0];
        mix=0;
        disp='off';
    case 5
        noise=[0 0];
        mix=0;
        disp='off';
    case 6
        mix=0;
        disp='off';
   case 7
        disp='off';
   case 8
        disp=disp;
    otherwise
        error('Wrong number of input arguments');
end;


if (strcmp(disp,'on'))
    if (isempty(findobj('tag','proj')==1))
      figure; set(gcf,'tag','proj');
    end;
end;

if (isnumeric(num)==0)
    sz=5;
     for ii=1:size(angles,2)
         number(ii)=round(rand(1).*sz);
    end;
else
    if (size(num,2)==1)
        for ii=1:size(angles,2)
            number(ii)=num;
        end;
    else
        number=num;
    
    end;
end;

if (isnumeric(sh)==0)
    sz=size(template,1).*0.05;
    for ii=1:size(angles,2)
         shift(ii,:)=[round(sz.*rand(1)) round(sz .*rand(1))];
    end;
else
    if (size(sh,1)==1)
        for ii=1:size(angles,2)
            shift(ii,:)=sh;
        end;
    end
end;

z_stack=1;
mask=tom_spheremask(ones(size(template,1)),(size(template,1)./2)-15,14);
mask=ones(size(template,1));
for i=1:size(angles,2)
    
    %first approach
    rot_tmp=tom_rotate(template,[0 0 angles(2,i)]);
    rot_tmp=tom_rotate(rot_tmp,[270 90 angles(1,i)]);
    
    %second approach
    eu=tom_sum_rotation([0 0 angles(2,i); 270 90 angles(1,i)],[0 0 0 ; 0 0 0]);
    rot_tmp2=tom_rotate(template,[eu]);
    rot_tmp2=single(rot_tmp2);
    proj2=double(sum(rot_tmp2,3));
    % end second approach
    
    
    rot_tmp=single(rot_tmp);
    proj=double(sum(rot_tmp,3));
    
    
    if (strcmp(disp,'on'))
        figure(findobj('tag','proj'));
        subplot(3,1,1);
        tom_imagesc(double(proj));
        subplot(3,1,2);
        tom_imagesc(double(proj2));
        drawnow;
    end;
    
    
    
    m=(mask==0).*mean2(proj);
    
    proj_org=proj;
    for j=1:number(i);
        proj=proj_org;
        if (noise(1)~=0) %| noise(2)~=0)
            %proj=imnoise(proj,'gaussian',noise(1),noise(2));
            proj=(proj+0.0005*(rand(size(proj))) );
        end;
        if (isempty(find(shift==0))==0) 
            angg=rand(1).*10;
            %proj=tom_rotate(proj,angg);
            proj=tom_shift(proj,shift(i,:));
        end;

        proj=(proj.*mask)+m;
        stack(:,:,z_stack)=(proj.*mask)+m;
        log(z_stack,1)= angles(1,i);
        log(z_stack,2)=angles(2,i);
        log(z_stack,3)=eu(1);
        log(z_stack,4)=eu(2);
        log(z_stack,5)=eu(3);
        log(z_stack,6)=shift(i,1);
        log(z_stack,7)=shift(i,2);;
        align2d(1,z_stack).filename=filename;
        align2d(1,z_stack).dataset=['phantom_' strrep(strrep(strrep(datestr(now),' ','_'),':','_'),'-','_')];
        
        z_stack=z_stack+1;    
       
        
    end;
    
   if (strcmp(disp,'on'))
    figure(findobj('tag','proj')); 
    subplot(3,1,3);
    tom_imagesc(double(proj)); 
    set(gcf,'Name',['axis  '  num2str(angles(1,i)) '  angle ' num2str(angles(2,i)) ]);
    drawnow;
   end;
   
end;


for i=1:size(stack,3)
    align2d(1,i).stack_size=size(stack);
end;

tom_emwrite(filename,stack);


