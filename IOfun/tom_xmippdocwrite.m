function tom_xmippdocwrite(doc_name,doc_str)
% TOM_XMIPPDOCWRITE writes xmipp doc file
%  
%     tom_xmippdocwrite(doc_name,doc_str)
%  
%  PARAMETERS
%  
%    INPUT
%     doc_name      filename of xmipp doc file
%     doc_str       matlab doc-str
%    
%  
%  EXAMPLE
%      st=tom_xmippdocread('mydoc.doc',st);
%
%  REFERENCES
%  
%  NOTE:
% 
%  
%  
%
%  SEE ALSO
%     tom_xmippdocread,tom_spiderread
%  
%     created by FB 15/02/10
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


head_st=get_header_info(doc_str(1).header);

if (doc_str(1).cols > 9)
    base_out_st='%5d %d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f';
else
    base_out_st='%5d %d  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f'; 
end;


for i=1:length(head_st)
    base_out_st=[base_out_st ' %10.5f'];
end;
base_out_st=[base_out_st '\n'];

fid=fopen(doc_name,'wt');
fprintf(fid,' %s\n',doc_str(1).header);
for i=1:length(doc_str)
    fprintf(fid,' ; %s\n',doc_str(i).name);
    %fprintf(fid,'%5d %d  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
     %                                                                                            ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,doc_str(i).maxCC);
    
     if (length(head_st)==0)
         fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
             ,doc_str(i).yoff, doc_str(i).weight,doc_str(i).flip);
     end;
     
     if (length(head_st)==1)
         fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
             ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_')));
     end;
    
    if (length(head_st)==2)
        fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
                               ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_')),getfield(doc_str(i),strrep(head_st{2},'/','_')));
    end;
    if (length(head_st)==3)
        fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
                               ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_')) ...
                               ,getfield(doc_str(i),strrep(head_st{2},'/','_')),getfield(doc_str(i),strrep(head_st{3},'/','_')));
    end;
    if (length(head_st)==4)
        fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
                               ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_')) ...
                               ,getfield(doc_str(i),strrep(head_st{2},'/','_')),getfield(doc_str(i),strrep(head_st{3},'/','_')),getfield(doc_str(i),strrep(head_st{4},'/','_')));
    end;
    if (length(head_st)==5)
        fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
                               ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_')) ...
                               ,getfield(doc_str(i),strrep(head_st{2},'/','_')),getfield(doc_str(i),strrep(head_st{3},'/','_')),getfield(doc_str(i),strrep(head_st{4},'/','_'))...
                               ,getfield(doc_str(i),strrep(head_st{5},'/','_')));
    end;
    if (length(head_st)==6)
        fprintf(fid,base_out_st,doc_str(i).run_num,doc_str(i).cols,doc_str(i).rot,doc_str(i).tilt,doc_str(i).psi,doc_str(i).xoff ...
                               ,doc_str(i).yoff, doc_str(i).ref,doc_str(i).flip,getfield(doc_str(i),strrep(head_st{1},'/','_'))...
                               ,getfield(doc_str(i),strrep(head_st{2},'/','_')),getfield(doc_str(i),strrep(head_st{3},'/','_')),getfield(doc_str(i),strrep(head_st{4},'/','_'))...
                               ,getfield(doc_str(i),strrep(head_st{5},'/','_')),getfield(doc_str(i),strrep(head_st{5},'/','_')) );
    end;
     
end;

fclose(fid);


function head_st=get_header_info(line_in)


idx=strfind(line_in,',');

if (length(idx) > 6 )
    if (strfind(line_in(idx(6):idx(7)),'Flip')==0)
        disp('strange Header in doc !!');
        return;
    end;
else
    if (strfind(line_in(idx(5):idx(6)),'Flip')==0)
        disp('strange Header in doc !!');
        return;
    end;
    head_st='';
end;

if length(idx) > 6
    zz=0;
    for i=7:length(idx)
        if (i==length(idx))
            tmp=line_in(idx(i):end);
        else
            tmp=line_in(idx(i):idx(i+1));
        end;
        idx_tmp=strfind(tmp,' ');
        zz=zz+1;
        head_st{zz}=tmp(idx_tmp(1)+1:idx_tmp(2)-1);
    end;
end;





% function head_st=get_header_info(line_in)
% 
% 
% idx=strfind(line_in,',');
% 
% if (strfind(line_in(idx(6):idx(7)),'Flip')==0)
%     disp('strange Header in doc !!');
%     return;
% end;
% zz=0;
% for i=7:length(idx)
%     if (i==length(idx))
%         tmp=line_in(idx(i):end);
%     else
%         tmp=line_in(idx(i):idx(i+1));
%     end;
%     idx_tmp=strfind(tmp,' ');
%     zz=zz+1;
%     head_st{zz}=tmp(idx_tmp(1)+1:idx_tmp(2)-1);
% end;