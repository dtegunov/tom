function tom_av2_eman_lst2stacks(path_in,path_out,alignstruct)
%
%
%
%
%



list_of_files=dir('*.lst');


%build up Filestruct
im=tom_imagicread([path_in '/start.img']);

%im=tom_emheader(zeros(80,80,10000));

zzz=1;
for i=1:size(list_of_files,1)

   
   %parse eman file
   filename=[path_in 'cls000' num2str(i-1) '.lst'];
   if i-1 >=10
        filename=[path_in 'cls00' num2str(i-1) '.lst'];
   end; 
   if i-1 >= 100
        filename=[path_in 'cls0' num2str(i-1) '.lst'];
   end; 
   if i-1 >=1000
         filename=[path_in 'cls' num2str(i-1) '.lst'];
   end;
   lst_st=lst_file2struct(filename);
   
   if (isempty(lst_st))
        continue;
   end;
   
   %create folder
   folder_name=[path_out  '/class_' num2str(zzz)]; mkdir(folder_name);
   zzz=zzz+1;
   
   % write temporary stack for alignment
   write_tmp_stack(im,lst_st,folder_name);
   
   if (isfield(lst_st,'part_nr')==0)
       continue;
   end;
   
   %align_stack
   align_stack(alignstruct,folder_name,i);
   
   clear('lst_st');
    

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

function write_tmp_stack(im,lst_st,folder_name)


sz_im=size(im.Value);
tmp_st=zeros(sz_im(1),sz_im(2),1);
if (isfield(lst_st,'part_nr')==1)

    zz=1;
    for i=1:size(lst_st,2)
        try
            if (str2num(lst_st(i).val)==0)
                tmp_st(:,:,zz)=im.Value(:,:,lst_st(i).part_nr);
                zz=zz+1;
            end;
        catch
            disp('murat');
        end;

    end;
    tom_emwritec([folder_name '/tmp_st.em'],tmp_st);

end;

all_proj=tom_imagicread('proj.img');
tom_emwritec([folder_name '/proj.em'],all_proj.Value(:,:,lst_st(1).ref_nr+1));


function lst_st=lst_file2struct(filename)

disp(filename);
fid=fopen(filename);

if (fid==-1)
    lst_st='';
    return;
end;
zz=1;

%get reference
line=fgetl(fid);
line=fgetl(fid);
line=strrep(char(line),char(9),',');
[val rest]=strtok(line,',');
ref_nr=val;
[val rest]=strtok(rest,',');
ref_file=val;
lst_st(1).ref_nr=str2num(ref_nr);
lst_st(1).ref_file=ref_file;


while 1
    line=fgetl(fid);
    if (line==-1)
        break;
    end;
    
    line=strrep(char(line),char(9),',');
    line=strrep(char(line),'  ','');
    [val rest]=strtok(line,',');
    lst_st(zz).part_nr=str2num(val)+1;
    [val rest]=strtok(rest,',');
    lst_st(zz).stack=val;
    [val rest]=strtok(rest,',');
    lst_st(zz).angle1=val;
    [val rest]=strtok(rest,',');
    lst_st(zz).angle2=val;
    [val rest]=strtok(rest,',');
    lst_st(zz).shift1=val;
    [val rest]=strtok(rest,',');
    lst_st(zz).shift2=val;
    [val rest]=strtok(rest,',');
    lst_st(zz).val=val;
    lst_st(zz).ref_nr=str2num(ref_nr);
    lst_st(zz).ref_file=ref_file;
    
    zz=zz+1;
end;


fclose(fid);