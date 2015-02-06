function tom_av2_em_apply_refinement(class_doc,ref_doc_sym,output_doc)
% TOM_AV2_EM_APPLY_REFINEMENT 
%     tom_av2_em_apply_refinement(class_doc,ref_doc_sym,output_doc)
%  
%  PARAMETERS
%  
%    INPUT
%     ref                 reference image
%
%    
%    OUTPUT
%     
%  
%  EXAMPLE
%
%   tom_av2_em_apply_refinement('testcl5.doc','sym_Iter73.doc','out.doc')
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
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


 
 
%read corresponding docs
disp(['Reading ' class_doc]);
doc_cl_org=tom_xmippdocread(class_doc);

%sort doc_ref according 2 particle index
[a new_i]=sort([doc_cl_org(:).part_idx]);
doc_cl=doc_cl_org(new_i);
[a i_f c]=unique([doc_cl(:).part_idx],'first');
[a i_l c]=unique([doc_cl(:).part_idx],'last');
diff=i_l-i_f;

disp(['Reading ' ref_doc_sym]);
doc_ref_org=tom_xmippdocread(ref_doc_sym);
 
%sort doc_ref according 2 particle index
[a new_i]=sort([doc_ref_org(:).part_idx]);
doc_ref=doc_ref_org(new_i);

disp(['Sync index ']);
%get corresponding particles
idx_sm=find(ismember([doc_ref(:).part_idx],[doc_cl(:).part_idx]));

idx=zeros(length([doc_cl(:).part_idx]),1);

zz=0;
for i=1:length(idx_sm)
    for ii=1:diff(i)+1
        zz=zz+1;
        idx(zz)=idx_sm(i);
    end;
end;


disp(['Getting Alg Data ']);
doc_cl_ref=doc_ref(idx);
zz=0;

use_euler=0;

if (use_euler==1)
    
    for i=1:length(doc_cl_ref)
        flip_cl=doc_cl(i).flip;
        ang_cl=[doc_cl(i).rot doc_cl(i).tilt doc_cl(i).psi];
        flip_ref=doc_ref(idx(i)).flip;
        ang_ref=[doc_ref(idx(i)).rot doc_ref(idx(i)).tilt doc_ref(idx(i)).psi];
        ang_ref_trans=tom_sum_rotation_zyz([180 90 270;ang_ref]);
        
        
        if (mod(i,10000)==0)
            disp(num2str(i));
        end;
        
        if (flip_cl==flip_ref) && max(abs(ang_ref_trans-ang_cl))<45
            zz=zz+1;
            doc_cl_ref(zz)=doc_ref(idx(i));
            doc_cl_ref(zz).rot=ang_ref_trans(1);
            doc_cl_ref(zz).tilt=ang_ref_trans(2);
            doc_cl_ref(zz).psi=ang_ref_trans(3);
        end;
        
        if (flip_cl==flip_ref) && max(abs((ang_ref_trans+[180 0 0])-ang_cl))<45
            zz=zz+1;
            doc_cl_ref(zz)=doc_ref(idx(i));
            doc_cl_ref(zz).rot=ang_ref_trans(1)+180;
            doc_cl_ref(zz).tilt=ang_ref_trans(2);
            doc_cl_ref(zz).psi=ang_ref_trans(3);
        end;
        
        if (flip_cl~=flip_ref) && max(abs(ang_ref_trans-ang_cl))<45
            zz=zz+1;
            doc_cl_ref(zz)=doc_ref(idx(i));
            doc_cl_ref(zz).rot=ang_ref_trans(1);
            doc_cl_ref(zz).tilt=ang_ref_trans(2);
            doc_cl_ref(zz).psi=ang_ref_trans(3);
            doc_cl_ref(zz).flip=flip_cl;
        end;
        
        %     if (flip_cl==flip_ref) && max(abs(ang_ref_trans-ang_cl))>45
        %         zz=zz+1;
        %         doc_cl_ref(zz)=doc_ref(idx(i));
        %         doc_cl_ref(zz).rot=ang_ref_trans(1);
        %         doc_cl_ref(zz).tilt=ang_ref_trans(2);
        %         doc_cl_ref(zz).psi=ang_ref_trans(3);
        %         doc_cl_ref(zz).flip=flip_cl;
        %     end;
        
        
    end;
    
    
else
    
    for i=1:length(doc_cl_ref)
        flip_cl=doc_cl(i).flip;
        ang_cl=[doc_cl(i).rot doc_cl(i).tilt doc_cl(i).psi];
        flip_ref=doc_ref(idx(i)).flip;
        ang_ref=[doc_ref(idx(i)).rot doc_ref(idx(i)).tilt doc_ref(idx(i)).psi];
        ang_ref_trans=tom_sum_rotation_zyz([180 90 270;ang_ref]);
        dist1= tom_calc_euler_distance(ang_ref_trans,ang_cl);
        dist2= tom_calc_euler_distance(ang_ref_trans+[180 0 0],ang_cl);
%         dist3= tom_calc_euler_distance(tom_sum_rotation(ang_ref_trans,[0 180 0]),ang_cl);
%         dist4= tom_calc_euler_distance(tom_sum_rotation(ang_ref_trans,[180 180 0]),ang_cl);
        % disp([ num2str(dist1) ' ' num2str(dist2)  ' ' num2str(flip_cl) ' '  num2str(flip_ref)]);
        
        if (dist1<dist2 && dist1 < 20 && flip_cl==flip_ref)
            zz=zz+1;
            doc_cl_ref(zz)=doc_ref(idx(i));
            doc_cl_ref(zz).rot=ang_ref_trans(1);
            doc_cl_ref(zz).tilt=ang_ref_trans(2);
            doc_cl_ref(zz).psi=ang_ref_trans(3);
        end;
        
        if (dist2<dist1 && dist2 < 20 && flip_cl==flip_ref)
            zz=zz+1;
            doc_cl_ref(zz)=doc_ref(idx(i)); 
            doc_cl_ref(zz).rot=ang_ref_trans(1)+180;
            doc_cl_ref(zz).tilt=ang_ref_trans(2);
            doc_cl_ref(zz).psi=ang_ref_trans(3);
        end;
        
%         if (dist3<dist4 && dist3 < 20 && flip_cl~=flip_ref)
%             zz=zz+1;
%             doc_cl_ref(zz)=doc_ref(idx(i));
%             doc_cl_ref(zz).rot=ang_ref_trans(1)+180;
%             doc_cl_ref(zz).tilt=ang_ref_trans(2);
%             doc_cl_ref(zz).psi=ang_ref_trans(3);
%         end;
%     
%         if (dist4<dist3 && dist4 < 20 && flip_cl~=flip_ref)
%             zz=zz+1;
%             doc_cl_ref(zz)=doc_ref(idx(i));
%             doc_cl_ref(zz).rot=ang_ref_trans(1)+180;
%             doc_cl_ref(zz).tilt=ang_ref_trans(2);
%             doc_cl_ref(zz).psi=ang_ref_trans(3);
%         end;
        
        
        
    end;
end;


cut=doc_cl_ref(1:zz);
cut(1).header=doc_cl_org(1).header;
cut(1).part_idx_unique=0;

disp(['Writing doc file']);
tom_xmippdocwrite(output_doc,cut);
disp(['Used: ' num2str(zz) ' of ' num2str(length(doc_cl_ref)) ]);
 
 

%do high 2 low

% idx=zeros(length([doc_cl(:).part_idx]),1);
% tic;
% %for i=1:length([doc_cl(:).part_idx])
% tmpp=ismember([doc_ref(:).part_idx],[doc_cl(:).part_idx]);
% for i=1:100
%     idx(i)= find(([doc_ref(:).part_idx]==doc_cl(i).part_idx));
%     
%     if (mod(i,50)==0)
%         toc;
%         disp(num2str(i));
%         tic;
%     end;
% end;


%idx= find((ismember([doc_ref(:).part_idx],[doc_cl(:).part_idx]))==1);








 