function [all_peaks] = tom_peak_paralell(A,R,num_of_peaks,num_of_nodes,tmp_dir)

%to do change from file interface to d-Array

%split image
packages_im=tom_calc_packages(num_of_nodes,size(A,1));
for i=1:num_of_nodes
    im_tmp=A(packages_im(i,1):packages_im(i,2),:);
    tom_emwritec([tmp_dir '/tmppeakim_' num2str(i) '.em'],im_tmp);
end;

packages_peaks=tom_calc_packages(num_of_nodes,num_of_peaks);

parfor i=1:num_of_nodes
   
   M=tom_emreadc([tmp_dir '/tmppeakim_' num2str(i) '.em']);
   M=M.Value;
   for ii=1:packages_peaks(i,3)
        [c(ii,:) val(ii,:) M] = tom_peakc(M,R); 
   end; 
   tmp_peaks{i}.c=c;
   tmp_peaks{i}.val=val;
end;

zz=1;
for i=1:num_of_nodes
    for ii=1:length(tmp_peaks{i}.val)
        all_peaks(zz,1)=tmp_peaks{i}.c(ii,1)+packages_im(i,1)-1;
        all_peaks(zz,2)=tmp_peaks{i}.c(ii,2);
        all_peaks(zz,3)=tmp_peaks{i}.val(ii);
        zz=zz+1;
    end;
end;

