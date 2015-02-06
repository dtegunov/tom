function align_out=tom_av2_change_align2d_basedir(align2d,new_basedir)


align_out=align2d;
for i=1:size(align2d,2)
     [old_path filen ext]=fileparts(align2d(1,i).filename);
     
     align_out(1,i).filename=[new_basedir '/' filen ext];
     align_out(1,i).isaligned=0; 
end;