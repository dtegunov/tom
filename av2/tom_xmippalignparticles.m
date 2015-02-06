function tom_xmippalignparticles(doc_filename)


fprintf('\n%s ', ['Reading  ' doc_filename ':']);  
st=tom_xmippdocread(doc_filename);
fprintf('%s \n', ['...done!']); 
disp('');


fprintf('%s ', ['Converting Angles: ' ]); 
for i=1:size(st,2)
    tmp(i)=st(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
    proj_angles(i,:) = angles';
    tmp_angles(i,:)=[st(i).rot st(i).tilt st(i).psi];
    if (mod(i,2000)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n','done!');

fprintf('%s ', ['Allocating Memory: ' ]); 
num_of_ref=max(tmp);
im_tmp=tom_spiderread(['./' st(1).name]);
count=zeros(num_of_ref);
count2=zeros(num_of_ref);
class=zeros(size(st,2),1);
fprintf('%s \n', ['...done! ' ]); 

stack_all = zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),length(st),'single');
jj =1;
fprintf('%s ', ['Building Classes: ' ]); 
for i=1:size(st,2)
    im=tom_spiderread(['./' st(i).name]);
    if (st(i).flip==0)
        im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(i).xoff st(i).yoff]),tmp_angles(i,3));
        stack_all(:,:,jj) = im_tmp_alg;
        count(st(i).ref)=count(st(i).ref)+1;
        [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
        proj_ang(:,jj)=tttemp;
        class(jj)=st(jj).ref;
    else
        im.Value=tom_mirror(im.Value,'x');
        im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
        stack_all(:,:,jj) = im_tmp_alg;
        count2(st(i).ref)=count2(st(i).ref)+1;
        [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
        proj_ang(:,jj)=tttemp;
        class(jj)=st(jj).ref;
    end;
    jj = jj +1;
    if (mod(i,2000)==0)    
        fprintf('%s','.');
    end;
end;

save('alignmentdata.mat','class','proj_ang');
tom_emwrite('stack_aligned.em',stack_all);
tom_emstack2spiderseries('stack_aligned.em','parts_aligned/26S_', 'parts_aligned/26S.sel');


fprintf('%s \n', ['...done! ' ]);



