function tom_convert_TietzTiff2EM(basedir,ext,offset,form,xray)

if nargin < 2
    ext='em';
end;

if nargin < 3
    %offset=129115;
    offset=0;
end;

if nargin < 4
    form='int16';
end;

if nargin < 5
    xray=1;
end;



files=dir([basedir '*.mat']);
num_of_files=length(files);

error_conf_c=0;
error_c_header=0;
root=fileparts(basedir);

disp(' ');
disp(['Starting Conversion of ' num2str(num_of_files) ' files' ]);
disp(' ');
disp(' ');



for i=1:num_of_files

    for ii=1:5
        try
            fname=[root '/' files(i).name];
            
            disp(['Processing: ' fname]);
            [b1 b2 dext]=fileparts(fname);
            fname=[b1 '/' b2 '.em'];
            %inf_st=imfinfo(fname);
            im_sz=[8192 8192];
            %if strcmp(inf_st.Format,'tif')
                if (xray==1)
                    im=tom_xraycorrect(tom_rawread(fname,form,'le',[im_sz(1) im_sz(2) 1],0,offset));
                else
                    im=tom_rawread(fname,form,'le',[im_sz(1) im_sz(2) 1],0,offset);
                end;
            %end;
            im=tom_emheader(im);
            try
                load([b1 '/' b2 '.mat']);
                im.Header=Header_EM;
            catch
                disp(['ERROR: Header Structure missing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!']);
                error_c_header=error_c_header+1;
            end; % try
            tom_emwrite([b1 '/' b2 '.em'],im);

            %delete(fname);
            delete([b1 '/' b2 '.mat']);
            break;   
        catch
            disp(['ERROR: Converting File: ' fname]);
            error_conf_c=error_conf_c+1;
            error_files{error_conf_c}=fname;
        end; %try
    end; % for ii=1:5
    disp(' ');
end; % i=1:num_of_files

disp(' ');
disp(' ');
disp('--------------ERRORS-------------------------------------------------------------------------------');

if (error_conf_c>0)
    for i=1:error_conf_c
        disp(error_files{i});
    end;
else
    disp('NO ERRORS =========> YES!')
end;