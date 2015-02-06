function tom_av2_xmipp_backup_rec(ori,dest,start_iter,final_iter)

% ori='ProjMatch/r3_m1m2chilo256_testmask_pad2/';
% dest='/fs/pool/pool-nickell3/26S/em/data/bohn/2d/rec/11__25_chilo256_5A_c2/ProjMatch/r3_m1m2chilo256_testmask_pad2/';
% 
% start_iter=12;
% final_iter=49;

for i=start_iter:final_iter;
    call=['cp ' ori 'Iter_' num2str(i) '/Iter_' num2str(i) '_resolution.fsc ' dest];
    disp(call);
    unix(call);
end;

call=['cp ' ori 'Iter_' num2str(final_iter) '/Iter_' num2str(final_iter) '_current_angles.doc ' dest];
disp(call);
unix(call);
call=['cp ' ori 'Iter_' num2str(final_iter) '/Iter_' num2str(final_iter) '_reconstruction.vol ' dest];
disp(call);
unix(call);
call=['cp ' ori '*.py ' dest 'xmipp_protocol_backup_Iter_' num2str(final_iter) '.py'];
disp(call);
unix(call);
