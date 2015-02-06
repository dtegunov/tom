function model=tom_av2_reconstruction(startmodel,size_particle,align2d,angles,iterations,sym,corrflag,filter,num_of_tasks,demo,thresh)

moviemode = 1;

%**************************************************************************
% Welcome to main program
%**************************************************************************


%load ws_test_bapro;
model=startmodel;
error_m=0;
error_count=0;
iter_num=1;

for i=1:iterations

    if demo >= 1
        disp_reconstruction2d('init',size(angles,2),size(angles,2),size_particle,5);
        if nargin == 11
            disp_reconstruction2d('init_surface',startmodel, thresh,angles);
        end
        if moviemode == 1
            disp_reconstruction2d('movieopen','/fs/sally/pool-baumeister/demo/demo.avi',15,100);
            disp_reconstruction2d('moviepicture');
        end
    end

    fprintf('\n%s \n \n',['2d reconstruction start: ' num2str(iter_num) ', start time: ' datestr(clock)]);
    
    align2d=feed_align2d_structure(align2d,iter_num,size_particle,size(model),sym,filter,corrflag,num_of_tasks);
    
    if (error_m==0)
        [align2d error_m]=tom_av2_create_projections(model,angles,align2d,iter_num,demo);
    end;
    
    if (error_m==0)
        [align2d error_m]=tom_av2_angular_classification(align2d,iter_num,demo); % correlate projections and aligned particles
    end;
    
    if (error_m==0)
        [model errror_m]=tom_av2_backproj(align2d,iter_num,0);
        %        [model errror_m]=tom_av2_backproj(align2d,iter_num,demo);
    end;
    
    if (error_m==0)
        [model error_m]=tom_av2_post_processing(align2d,iter_num,model);
    end;
    
    if (error_m==0)
        iter_num=iter_num+1;
        error_count=0;
    else
        switch error_count
            case 2
                tom_purgequeue(align2d(iter_num,1).paraell.jobmanager);
            case 5
               % unix('sudo /etc/init.d/matlab');
            case 100
        end;
        error_m=0;
        error_count=error_count+1;
    end;
    mod = tom_emreadc(['model/model_' num2str(i)]);
    disp_reconstruction2d('disp_reconst',mod.Value(:,:,40:120));
    
    tom_emwrite(['model/model_bin_' num2str(i)],model);
    save(['align/align2d_save'],'align2d');
    disp_reconstruction2d('destroy');
end;



function align2d=feed_align2d_structure(align2d,iternum,size_particle,size_model,sym,filter,corrflag,num_of_tasks)

%h=tom_reademheader(align2d(1,1).filename); 
%sz=h.Header.Size(3);
sz=align2d(1,1).stack_size(3);


% for i=1:sz
    if (iternum==1)
        
        
        %cyro settings heiliger stand !
        align2d(iternum,1).model.size=size_model;
        align2d(iternum,1).model.size_part=size_particle;
        align2d(iternum,1).model.sym=sym;
        align2d(iternum,1).model.sym_angle=[180 0 161];
        align2d(iternum,1).model.mass=2700;
        align2d(iternum,1).model.pixel_size=3.6;
        align2d(iternum,1).reconstruction.weight_flag='volume';    
        align2d(iternum,1).filter.mask.model_mask_sp.Apply=1;
        align2d(iternum,1).filter.mask.model_mask_sp.Method='sphere3d';
        align2d(iternum,1).filter.mask.model_mask_sp.Value=[160 160 160 70 5];
        align2d(iternum,1).filter.mask.model_mask_cy.Apply=1;
        align2d(iternum,1).filter.mask.model_mask_cy.Method='cylinder3d';
        align2d(iternum,1).filter.mask.model_mask_cy.Value=[160 160 160 30 5];
        align2d(iternum,1).filter.mask.classify1.Apply=1;
        align2d(iternum,1).filter.mask.classify1.Value=[160 160 70 5 ];
        align2d(iternum,1).filter.mask.classify1.Method='sphere';
        align2d(iternum,1).filter.mask.classify2.Apply=1;
        align2d(iternum,1).filter.mask.classify2.Value=[160 160 20 5];
        align2d(iternum,1).filter.mask.classify2.Inverse=1;
        align2d(iternum,1).filter.mask.classify2.Method='sphere';
        align2d(iternum,1).filter.mask.align.Apply=1;
        align2d(iternum,1).filter.mask.align.Value=[160 160 35 5];
        align2d(iternum,1).filter.mask.align.Method='sphere';
        align2d(iternum,1).filter.mask.ccf_rot.Apply=1;
        align2d(iternum,1).filter.mask.ccf_rot.Value=[80 320 5 3];
        align2d(iternum,1).filter.mask.ccf_rot.Method='sphere';
        align2d(iternum,1).filter.mask.ccf_trans.Apply=1;
        align2d(iternum,1).filter.mask.ccf_trans.Value=[160 160 77 3];
        align2d(iternum,1).filter.mask.ccf_trans.Method='sphere';
        align2d(iternum,1).filter.filter.classify.Apply=1;
        align2d(iternum,1).filter.filter.classify.Value=[1 15 10];
        align2d(iternum,1).filter.filter.classify.Method='bandpass';
        align2d(iternum,1).filter.filter.align.Apply=1;
        align2d(iternum,1).filter.filter.align.Value=[1 32 0];
        align2d(iternum,1).filter.filter.align.Method='bandpass';
        align2d(iternum,1).corr_flag=corrflag;
        align2d(iternum,1).paraell.number_of_tasks=num_of_tasks;
        align2d(iternum,1).paraell.time_out=15000;
        align2d(iternum,1).paraell.packageloss=0.05;
        align2d(iternum,1).paraell.jobmanager='cluster02';
        align2d(iternum,1).paraell.workers.min=1;
        align2d(iternum,1).paraell.workers.max=8;
% ende cyro settings heiliger stand
        


%neagitve stain settings hl
%         align2d(iternum,1).model.size=size_model;
%         align2d(iternum,1).model.size_part=size_particle;
%         align2d(iternum,1).model.sym=sym;
%         align2d(iternum,1).model.sym_angle=[180 0 161];
%         align2d(iternum,1).model.mass=2700;
%         align2d(iternum,1).model.pixel_size=5.0;
%         align2d(iternum,1).model.size_part=size_particle;
%         align2d(iternum,1).model.sym=sym;
%         align2d(iternum,1).model.sym_angle=[180 0 161];
%         align2d(iternum,1).model.mass=2700;
%         align2d(iternum,1).model.pixel_size=3.6;
%         align2d(iternum,1).reconstruction.weight_flag='volume';    
%         align2d(iternum,1).filter.mask.model_mask_sp.Apply=1;
%         align2d(iternum,1).filter.mask.model_mask_sp.Method='sphere3d';
%         align2d(iternum,1).filter.mask.model_mask_sp.Value=[128 128 128 70 5];
%         align2d(iternum,1).filter.mask.model_mask_cy.Apply=1;
%         align2d(iternum,1).filter.mask.model_mask_cy.Method='cylinder3d';
%         align2d(iternum,1).filter.mask.model_mask_cy.Value=[128 128 128 32 5];
%         align2d(iternum,1).filter.mask.classify1.Apply=1;
%         align2d(iternum,1).filter.mask.classify1.Value=[128 128 70 5 ];
%         align2d(iternum,1).filter.mask.classify1.Method='sphere';
%         align2d(iternum,1).filter.mask.classify2.Apply=1;
%         align2d(iternum,1).filter.mask.classify2.Value=[128 128 16 5];
%         align2d(iternum,1).filter.mask.classify2.Inverse=1;
%         align2d(iternum,1).filter.mask.classify2.Method='sphere';
%         align2d(iternum,1).filter.mask.align.Apply=1;
%         align2d(iternum,1).filter.mask.align.Value=[128 128 20 5];
%         align2d(iternum,1).filter.mask.align.Method='sphere';
%         align2d(iternum,1).filter.mask.ccf_rot.Apply=1;
%         align2d(iternum,1).filter.mask.ccf_rot.Value=[64 256 5 3];
%         align2d(iternum,1).filter.mask.ccf_rot.Method='sphere';
%         align2d(iternum,1).filter.mask.ccf_trans.Apply=1;
%         align2d(iternum,1).filter.mask.ccf_trans.Value=[128 128 5 3];
%         align2d(iternum,1).filter.mask.ccf_trans.Method='sphere';
%         align2d(iternum,1).filter.filter.classify.Apply=1;
%         align2d(iternum,1).filter.filter.classify.Value=[1 32 0];
%         align2d(iternum,1).filter.filter.classify.Method='bandpass';
%         align2d(iternum,1).filter.filter.align.Apply=1;
%         align2d(iternum,1).filter.filter.align.Value=[1 32 0];
%         align2d(iternum,1).filter.filter.align.Method='bandpass';
%         align2d(iternum,1).corr_flag=corrflag;
%         align2d(iternum,1).paraell.number_of_tasks=num_of_tasks;
%         align2d(iternum,1).paraell.time_out=15000;
%         align2d(iternum,1).paraell.packageloss=0.05;
%         align2d(iternum,1).paraell.jobmanager='cluster02';
%         align2d(iternum,1).paraell.workers.min=1;
%         align2d(iternum,1).paraell.workers.max=8;
% 
%         align2d(iternum,1).filter.mask.ccf_trans.Apply=1;
%         align2d(iternum,1).filter.mask.ccf_trans.Value=[128 128 5 3];
%         align2d(iternum,1).filter.mask.ccf_trans.Method='sphere';
%         align2d(iternum,1).filter.filter.classify.Apply=1;
%         align2d(iternum,1).filter.filter.classify.Value=[1 32 0];
%         align2d(iternum,1).filter.filter.classify.Method='bandpass';
%         align2d(iternum,1).filter.filter.align.Apply=1;
%         align2d(iternum,1).filter.filter.align.Value=[1 32 0];
%         align2d(iternum,1).filter.filter.align.Method='bandpass';
%         align2d(iternum,1).corr_flag=corrflag;
%         align2d(iternum,1).paraell.number_of_tasks=num_of_tasks;
%         align2d(iternum,1).paraell.time_out=15000;
%         align2d(iternum,1).paraell.packageloss=0.05;
%         align2d(iternum,1).paraell.jobmanager='cluster02';
%         align2d(iternum,1).paraell.workers.min=1;
%         align2d(iternum,1).paraell.workers.max=8;

%end negative stain hl














%  align2d(iternum,1).filter.filter.classify.Apply=1;
%         align2d(iternum,1).filter.filter.classify.Value=[1 32 0];
%         align2d(iternum,1).filter.filter.classify.Method='bandpass';
% 
% 

        
%         align2d(iternum,1).model.mass=5000;
%         align2d(iternum,1).model.pixel_size=3.4;
%         align2d(iternum,1).filter.mask=[51 2 0]; %used in classification
%         align2d(iternum,1).filter.mask_ccf=[64 5 0]; %used in classification
%         align2d(iternum,1).filter.mask_ccf_rot=[15 5 0];
%         align2d(iternum,1).filter.mask_align=[64 64 0];
%         align2d(iternum,1).filter.spmask_model=[64 5 0]; %proj and backproj 
%         align2d(iternum,1).filter.cymask_model=[64 5 0]; %proj and backproj     
%         align2d(iternum,1).filter.bandpass=[2 50 2];
%         
%         align2d(iternum,1).filter.mask=[70 5 0]; %used in classification
%         align2d(iternum,1).filter.mask_ccf=[46 5 0]; %used in classification
%         align2d(iternum,1).filter.mask_ccf_rot=[20 5 0];
%         align2d(iternum,1).filter.mask_align=[50 30 4];        
%        align2d(iternum,1).filter.bandpass=[1 32 0]; 

%        align2d(iternum,1).filter.spmask_model=[70 5 0]; %proj and backproj 
%        align2d(iternum,1).filter.cymask_model=[24 5 0]; %proj and backproj     


       
        
       
        if (isfield(align2d(iternum,1),'file_path')==0)
            align2d(iternum,1).file_path=pwd;
            if (isempty(align2d(iternum,1).file_path)==1)
                align2d(iternum,1).file_path=pwd;
            end;
        end;
    else
        %updat new iteration
        align2d(iternum,1).model=align2d(iternum-1,1).model;
        align2d(iternum,1).reconstruction=align2d(iternum-1,1).reconstruction;
        align2d(iternum,1).filter=align2d(iternum-1,1).filter;
        align2d(iternum,1).filename=align2d(iternum-1,1).filename;
        align2d(iternum,1).dataset=align2d(iternum-1,1).dataset;
        align2d(iternum,1).stack_size=align2d(iternum-1,1).stack_size;
        align2d(iternum,1).corr_flag=align2d(iternum-1,1).corr_flag;
        align2d(iternum,1).paraell=align2d(iternum-1,1).paraell;
        align2d(iternum,1).file_path=align2d(iternum-1,1).file_path;
        align2d(iternum,1).model.sym=align2d(iternum-1,1).model.sym;
    end;
% end;









