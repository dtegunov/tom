function tom_av2_check_xray(input)

if (tom_isemfile(input))
    header=tom_reademheader(input);
    for i=1:size(header.Header,3)
        im=tom_emreadc(input,'subregion',[1 1 i],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
        im.Value=tom_xraycorrect(im.Value);
        tom_emwritec(input,im.Value,'subregion',[1 1 i],[header.Header.Size(1)-1 header.Header.Size(1)-1 0]);
    end;
    
    
else
    
    st=importdata(input);
    
    for i=1:length(st.textdata)
        try
            im=tom_spiderread(st.textdata{i});
            d=std2(im.Value);
            meanimg=mean(mean(im.Value));
            J=im.Value;
            
            for ii=5:0.5:10
                [x y]=find(J>(meanimg+ii*d) | J<(meanimg-ii*d)); % factor 10 added by SN
                if isempty(x)==0
                    disp(st.textdata{i});
                    disp(['values outside range found: ' num2str(size(x,1)) '  '  num2str(ii) '*std = ']);
                    %pause(1);
                end;
            end;
           
            disp('*****************************************');
            
            %pause(2);
            
            
        catch ME
            disp(st.textdata{i});
            disp( ME.stack);
            disp( ME.message);
        end;
    end;
    
end;