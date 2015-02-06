function tom_av2_xmipp_process_stack(input_sel,output_parts_dir,sel_file_name)


%
%example
%

fp=fopen(sel_file_name,'wt');

st=importdata(input_sel);


for i=1:length(st.data)
    im=tom_spiderread(st.textdata{i});
    im=tom_bin(im.Value,1);
    [a b c]=fileparts(st.textdata{i});
    tom_spiderwrite([output_parts_dir '' b c],im);
    if (mod(i,1000)==0)
        disp([b c]);
    end;

    fprintf(fp,[output_parts_dir '/' b c ' 1 \n']);
    
end;

fclose(fp);