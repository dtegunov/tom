function tom_build_phantom_templates2(num_of_particles,filenames,wedge_size,binning)


if nargin==3
    binning=0;
end;



h=tom_reademheader(filenames{1});

wedge=tom_wedge(ones(h.Header.Size'),wedge_size);
zz=1;

filen(1,:)='part_1000.em';


for ii=1:size(num_of_particles,2)
    im=tom_emread(filenames{ii});
    im=im.Value;
    im=tom_filter(im,2);
    im_org=im./1;
    for i=1:num_of_particles(ii)
        im=im_org;
        %ang=rand(1,3).*360;
                  %ang=[270 90 (rand.*360)];
        ang=round((rand(1).*360)./30).*30;
        
        ang=tom_sum_rotation([ang 0 0; 270 90 90],[0 0 0; 0 0 0]);
        ang_out(:,zz)=round(ang);
        im=im+(rand(size(im)).*4);
        
        %wedge_rot=tom_rotate(wedge,round(ang));
        
        
         %add noise
        im=im.*(rand(size(im)).*0.001);
        
        
        im=tom_apply_weight_function(tom_rotate(im,ang),wedge);
        mask=tom_spheremask(ones(h.Header.Size'),h.Header.Size(1)./2);
        %im=tom_apply_weight_function(tom_rotate(im,ang),wedge);
        
        %im=tom_rotate(im,round([-ang(2) -ang(1) -ang(3)]));
       
        %norm it norman 
        im=tom_norm(im+1,'phase').*mask;
        
        if zz < 10
            tom_emwrite(['part_' num2str(zz) '.em'],im);
            filen(zz,:)=['part_' num2str(zz) '.em   '];
        end;
        
        if zz >= 10 && zz < 100
            tom_emwrite(['part_' num2str(zz) '.em'],im);
            filen(zz,:)=['part_' num2str(zz) '.em  '];
        end;
            
        
        if zz > 99
            try
                tom_emwrite(['part_' num2str(zz) '.em'],im);
                filen(zz,:)=['part_' num2str(zz) '.em '];
            catch
                disp('tt');
            end;
            
            
        end;
        
        disp(['murat  part_' num2str(zz) '.em']);
        zz=zz+1;
    end;
end;

disp('');
Align = tom_av3_create_alignlist(filen,[],1);

for i=1:size(Align,2)
   Align(1,i).Angle.Psi=-ang_out(1,i);
   Align(1,i).Angle.Phi=-ang_out(2,i);
   Align(1,i).Angle.Theta=-ang_out(3,i);
   Align(1,i).Tomogram.AngleMax=-(90-wedge_size);
   Align(1,i).Tomogram.AngleMin=(90-wedge_size); 
end;

Align=tom_av3_align_sum(Align);
save('align.mat','Align');

% save('angina.mat','ang_out');
% 
% 
% for i=1:zz-1 
%     if (i < 10)
%         ccc_cell{i}=['part_00' num2str(i) '.em']; 
%         ccc(i,:)=['part_00' num2str(i) '.em'];
%    if (i>10 && i < 100)
%         ccc_cell{i}=['part_0' num2str(i) '.em']; 
%         ccc(i,:)=['part_0' num2str(i) '.em'];
%     
%    end;
%    
%    if (i > 99)
%         ccc_cell{i}=['part_' num2str(i) '.em']; 
%         ccc(i,:)=['part_' num2str(i) '.em'];
%    end;
%    
%    end;
% end;
% 
% st.rotmask.Apply=0;
% % align2d=tom_av3_create_alignlist(ccc,[],1);
% % st_out=tom_reshape_stack(ccc_cell,align2d,1,'','',st);
% % 
% 
% 
% disp('end');
% 
% 
% save('files.mat','ang_out');
%tom_emwrite('stack_reshape.em',st_out);
%stack_reshape=