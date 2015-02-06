function tom_robust_start_engine

command=['xmipp_mpi_MLrefine3D -i ./data/26S.sel -vol ./three_seeds/three_seeds.sel -iter 3 -o ./MLsplit/MLsplit'];
logfile='./MLsplit/MLsplit';
command_restart=['xmipp_mpi_MLrefine3D -i ./data/26S.sel -vol ./three_seeds/three_seeds.sel -iter 3 -o ./MLsplit/MLsplit'];

%cluster={'cluster01'; 'cluster02'; 'cluster03'; 'cluster04'; 'cluster05'; 'cluster06'; 'cluster08'};
cluster={'calypso'; 'atlas'; 'prometheus'};
nr_of_cpus=2;
repeat_lamboot_max=30;

%%
cluster_present=check_for_present_cluster(cluster);

%%
for i=1:size(cluster_present,2)
    [status,result] = system(['ssh ' cluster_present{i} ' lamnodes']);
    if status==0
        [status,laminfo] = system(['ssh ' cluster_present{i} ' lamnodes']);break;
    end;    
end;

lam='';

fp=fopen(['~/lamboot.robust'], 'wt');


for i=1:size(cluster_present,2)
    lam=[lam cluster_present{i} ' cpu=' num2str(nr_of_cpus) '\n'];
    fprintf(fp,'%s\n',[cluster_present{i} ' cpu=' num2str(nr_of_cpus)]);
end;
fclose(fp);
cluster_present=check_for_present_cluster(cluster);
st=try_lamboot(cluster_present, repeat_lamboot_max);

error_count=0;
for i=1:100000

    if error_count==0
        [status,result] = system(['mpirun -np ' num2str(st(1).sum_num_procs) commmand]);
    end;

    if status==0 break; end;

    error_count=error_count+1;

    if error_count==1
        [status,result] = system(['mpirun -np ' num2str(st(1).sum_num_procs) commmand_restart]);
        disp(['tom: try first restart, i=' num2str(i)]);
    end;
    
    if error_count > 1
        cluster_present=check_for_present_cluster(cluster);
        st=try_lamboot(cluster_present, repeat_lamboot_max);
        [status,result] = system(['mpirun -np ' num2str(st(1).sum_num_procs) commmand_restart]);
        disp(['tom: try next restart, i=' num2str(i)]);
        error_count=1;
    end;
    
    if status==0 break; end;



end;

%%
function cluster_present=check_for_present_cluster(cluster)

ii=1;
disp(['tom: check for nodes:']);
for i=1:size(cluster,1)
    [status,message]=unix(['ping -w 1 -c 1 ' cluster{i}]);
    if status==0
        cluster_present{ii}=cluster{i};
        ii=ii+1;
    end;    
end;
disp(['tom: present nodes:']);
for i=1:size(cluster_present,2)
    disp([cluster_present{i} '']);
end;

function st=parse_lamnodes(str_lamnodes)

num_of_hosts=size(findstr(str_lamnodes,char(10)),2);

idx=findstr(str_lamnodes,':');

 if ( (size(idx,2)./2) ~= num_of_hosts)
    st=-1;
    return;
 end;
 zz=1;
 st(zz).num_of_procs=0;
 
 for i=1:2:(size(idx,2))
    st(zz).num_of_procs=str2num(str_lamnodes(idx(i)+1:idx(i+1)-1));
    
    zz=zz+1;
 end;
 
 idx2=findstr(str_lamnodes,char(9));
 
 zz=1;
 k=1;
 for i=1:size(idx2,2)
   tmp=str_lamnodes(idx2(i)+1:idx(k)-1);
   st(zz).host=strtok(tmp,'.');
   zz=zz+1;
   k=k+2;
 end;
 
 sum_num_procs=0;
 for i=1:size(st,2)
     sum_num_procs=sum_num_procs+st(i).num_of_procs;
 end;
 for i=1:size(st,2)
     st(i).sum_num_procs=sum_num_procs;
 end; 
 

function st=try_lamboot(cluster_present, repeat_lamboot_max)

disp(['tom: try lamboot ...']);
for repeat_lamboot=1:repeat_lamboot_max
    [status,laminfo] = system(['ssh ' cluster_present{1} ' lamhalt']);
    [status,laminfo] = system(['ssh ' cluster_present{1} ' lamboot ~/lamboot.robust &']);
    pause(4);
    [status,laminfo] = system(['ssh ' cluster_present{1} ' lamnodes']);
    st=parse_lamnodes(laminfo);
    if (isstruct(st))
        disp(['tom: successful']);
        break;
    end;    
end;





