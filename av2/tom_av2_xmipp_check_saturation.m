function tom_av2_xmipp_check_saturation(script_file,part_nr_vect)



%make backup copy of scrt
copyfile(script_file,[script_file '_tom_backup']);


%parse xmipp_scrt
file = textread('xmipp_protocol_projmatch.py','%s','delimiter','\n');

[line idx_sel]=find_in_file(file,'SelFileName=');
idx=findstr(line,'''');

selfilename=line(idx(1)+1:idx(2)-1);

org_cell=importdata(selfilename);

[a b c]=fileparts(selfilename);
%main loop

for i=1:length(part_nr_vect)
    
    %creat new cell
    cell_filenmae=[b '_' num2str(i) c];
    write_cell(org_cell.textdata,part_nr_vect(i),cell_filenmae);
    
    %modify script
    modify_in_file(script_file,'SelFileName=',['SelFileName=' cell_filenmae]);
    
    modify_in_file(script_file,['WorkingDir=ProjMatch/' num2str(i) org_cell],['SelFileName=' cell_filenmae]);
    
    %call new script
    %assume lamboot is done!
    unix('./script_file');
    

end;






function line=modify_in_file(file,str,new_str)

fid=fopen(file,'r+');

for i=1:1000000
    line=fgetl(fid);
    pos=strfind(line,str);
    if (isempty(pos)==0)
          fprintf(fid,[new_str ' \n']);
    end;
    
end;





function line=find_in_file(file,str)

for i=1:length(file)
    if (isempty(strfind(file{i},str) )==0)
        line=file{i};
        break;
    end;
end;




function write_cell(part_str,num_of_parts,cell_filenmae)

fid=fopen(cell_filenmae,'w');

for i=1:num_of_parts
    fprintf(fid,[part_str{i} '.spi 1 \n']);
end;

fclose(fid);



