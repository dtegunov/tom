function tom_av2_create_angleclasslist(align2d,iter_num,datname)

 fid = fopen(datname,'W');
 fprintf(fid,'Stackname    %s Number of Particles  %d  \n',align2d(iter_num,1).filename,size(align2d(iter_num,:),2));
 fclose(fid);


for i=1:size(align2d(iter_num,:),2)
    proj_nr=align2d(iter_num,i).angleclass.proj_nr;
    angle_rot=align2d(iter_num,1).angular_scan(2,proj_nr);
    angle_nut=align2d(iter_num,1).angular_scan(1,proj_nr);
    fid = fopen(datname,'a');
    fprintf(fid,'part nr: %d   angle1: %4.2f   angle2 %4.2f \n',i,angle_rot,angle_nut);
    fclose(fid);
end;