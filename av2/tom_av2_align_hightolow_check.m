function tom_av2_align_hightolow_check(text_file)


st=importdata(text_file);

for i=1:size(st.data,1)
    shift=st.data(i,1:2);
    filename_high= st.textdata{i+1};
    filename_low=strrep(st.textdata{i+1},'/high/','/low/');
    %im1=tom_emread(filename_high);
    [a b c]=fileparts(filename_low);
    im2=tom_emreadc(filename_low);
    tom_emwrite(['shifted_low/' b c],tom_shift(single(im2.Value),shift));
    %diff1=(tom_norm(im1.Value,'mean0+1std')-tom_norm(tom_shift(im2.Value,shift),'mean0+1std') );
    %diff2=(tom_norm(im1.Value,'mean0+1std')-tom_norm(tom_shift(im2.Value,[0 0]),'mean0+1std') );
%     figure; tom_imagesc(diff1);
%     figure; tom_imagesc(diff2);
    disp(['shifted_low/' b c]);
end;