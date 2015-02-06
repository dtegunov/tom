function tom_av2_xmipp_lst2stacks(filename,path_out)
%
%
%
%
%


[path_in filen ext]=fileparts(filename);

lst_st=lst_file2struct(filename);

%build up Filestruct
%im=tom_imagicread([path_in '/start.img']);

%im=tom_emheader(zeros(80,80,10000));

zzz=1;
for i=1:size(lst_st,2)

   
   %parse eman file
%    filename=[path_in 'cls000' num2str(i-1) '.lst'];
%    if i-1 >=10
%         filename=[path_in 'cls00' num2str(i-1) '.lst'];
%    end; 
%    if i-1 >= 100
%         filename=[path_in 'cls0' num2str(i-1) '.lst'];
%    end; 
%    if i-1 >=1000
%          filename=[path_in 'cls' num2str(i-1) '.lst'];
%    end;
%    lst_st=lst_file2struct(filename);
%    
%    if (isempty(lst_st))
%         continue;
%    end;
   
   %create folder
   folder_name=[path_out  '/class_' num2str(zzz)]; mkdir(folder_name);
   disp(folder_name);
   zzz=zzz+1;
   
   % write temporary stack for alignment
   write_tmp_stack(lst_st,folder_name,i,path_in);
   
   if (isempty(lst_st{i}))
       continue;
   end;
   
   %align_stack
 %  align_stack(lst_st,folder_name,i);
   
   %clear('lst_st');
    

end;

function align_stack(alignstruct,folder_name,counter)

stack_path=[folder_name '/tmp_st.em'];
ref_path=[folder_name ,'/proj.em'];
ref_out=[folder_name ,'/avg.em'];
stack_out=[folder_name ,'/class_' num2str(counter) '.em'];
align2d=tom_av2_create_alignfromstack(stack_path);

[align2d,b]=tom_av2_align_stack(stack_path,ref_path,align2d,stack_out,ref_out,alignstruct.filter_param,alignstruct.parallel_param,alignstruct.iterations,alignstruct.demo,'no_sub');

copyfile([folder_name ,'/avg_' num2str(alignstruct.iterations(2)) '.em'],[folder_name ,'/avg.em']);
save([folder_name ,'/class_' num2str(counter) '.mat'],'align2d');

%calculate Variance
h=tom_reademheader(stack_out);
dm=tom_reshape_stack(stack_out);

var_img=tom_calc_variance(dm);
var_img=reshape(var_img,[h.Header.Size(1) h.Header.Size(2)]);
tom_emwrite([folder_name '/variance.em'],var_img);
disp('emd');

function write_tmp_stack(lst_st,folder_name,count,path_in)


%sz_im=size(im.Value);
%tmp_st=zeros(sz_im(1),sz_im(2),1);


if (isempty(lst_st{count})==0)
    zz=1;
    tmp=lst_st{count};
    for i=1:size(tmp,2)
        try
            imm=tom_spiderread(tmp(i).filename);
            tmp_st(:,:,zz)=imm.Value;
%             tmp_st_align(:,:,zz)=tom_shift(tom_rotate(imm.Value,tmp(i).theta),[tmp(i).shiftx tmp(i).shifty]);
%             tmp_st_align2(:,:,zz)=tom_shift(tom_rotate(imm.Value,tmp(i).theta),[-tmp(i).shiftx -tmp(i).shifty]);
            tmp_st_align(:,:,zz)=tom_rotate(tom_shift(imm.Value,[tmp(i).shiftx tmp(i).shifty]),tmp(i).theta);
%             tmp_st_align4(:,:,zz)=tom_rotate(tom_shift(imm.Value,[-tmp(i).shiftx -tmp(i).shifty]),tmp(i).theta);
                
            zz=zz+1;
        catch
            disp('murat');
        end;

    end;
    tom_emwritec([folder_name '/tmp_st.em'],tmp_st);
    
   
    
    tom_emwritec([folder_name '/class_' num2str(count) '.em'],tmp_st_align);
    align2d=tom_av2_create_alignfromstack([folder_name '/class_' num2str(count) '.em']);
    save([folder_name '/class_' num2str(count) '.mat'],'align2d');

    h=tom_reademheader([folder_name '/class_' num2str(count) '.em']);
    dm=tom_reshape_stack([folder_name '/class_' num2str(count) '.em']);

    var_img=tom_calc_variance(dm);
    var_img=reshape(var_img,[h.Header.Size(1) h.Header.Size(2)]);
    tom_emwrite([folder_name '/variance.em'],var_img);
    tom_emwrite([folder_name '/avg.em'],sum(tmp_st_align,3));
end;

% path_in='';

    filename=[path_in '/' path_in '_lib0000' num2str(count) '.proj'];
   if count >=10
        filename=[path_in '/' path_in '_lib000' num2str(count) '.proj'];
   end; 
   if count >= 100
        filename=[path_in '/' path_in '_lib00' num2str(count) '.proj'];
   end; 
   if count >=1000
         filename=[path_in '/' path_in '_lib0' num2str(count) '.proj'];
   end;
  % lst_st=lst_file2struct(filename);
   
   if (isempty(lst_st))
       %continue;
   end;

projeska=tom_spiderread(filename);
tom_emwritec([folder_name '/proj.em'],projeska.Value);


function lst_st=lst_file2struct(filename)

disp(filename);
fid=fopen(filename);

if (fid==-1)
    lst_st='';
    return;
end;
zz=1;

%skip first line
line=fgetl(fid);


while 1
    line=fgetl(fid);
    if (line==-1)
        break;
    end;
    tmp_st(zz).filename=strtrim(strrep(line,';',' '));
    line=fgetl(fid);  
    line=sscanf(line,'%f');
    tmp_st(zz).num=line(1);
    tmp_st(zz).val=line(2);
    tmp_st(zz).phi=line(3);
    tmp_st(zz).psi=line(4);
    tmp_st(zz).theta=line(5);
    tmp_st(zz).shiftx=line(6);
    tmp_st(zz).shifty=line(7);
    tmp_st(zz).ref=line(8);
    tmp_st(zz).ccc=line(9);
    
    zz=zz+1;
end;

%build up final structure

for i=1:size(tmp_st,2)
    try
        kk=lst_st{tmp_st(i).ref};
        kk(size(kk,2)+1)=tmp_st(i); 
    catch
        kk=tmp_st(i);
    end;
    lst_st{tmp_st(i).ref}=kk;
end;

disp('end');

% 
% [val rest]=strtok(line,',');
%     lst_st(zz).part_nr=str2num(val)+1;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).stack=val;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).angle1=val;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).angle2=val;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).shift1=val;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).shift2=val;
%     [val rest]=strtok(rest,',');
%     lst_st(zz).val=val;
%     lst_st(zz).ref_nr=str2num(ref_nr);
%     lst_st(zz).ref_file=ref_file;



fclose(fid);