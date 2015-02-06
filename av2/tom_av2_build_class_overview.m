function st=tom_av2_build_class_overview(folder,sz)


dd=dir(folder);

zz=1;

for i=1:size(dd,1)
    b_path=[folder '/class_' num2str(i) ];
    if (isdir(b_path)==1)
        [a num]=strtok(dd(3).name,'_');
        num=strrep(num,'_','');
        
        
        try
            disp([b_path '/avg.em']);
            avg=tom_emread([b_path '/avg.em']);
            avg=avg.Value;
            avg=tom_norm(avg,1);
        catch
            avg=ones(sz(1),sz(2)).*0.5;
        end;
        st(:,:,zz)=avg;
        zz=zz+1;
        
       
        
        
        try
             disp([b_path '/proj.em']);
             proj=tom_emread([b_path '/proj.em']);
             proj=proj.Value;
             proj=tom_norm(proj,1);
             %proj=zeros(sz(1),sz(2));
        catch
            proj=ones(sz(1),sz(2)).*0.5;
        end;
        
        st(:,:,zz)=proj;
        zz=zz+1;
       
        
         try
            disp([b_path '/variance.em']);
             var=tom_emread([b_path '/variance.em']);
            var=var.Value;
            var=tom_norm(var,1);
            if (mean2(var)==0)
                var=var+0.5;
            end;
            % var=zeros(sz(1),sz(2));
        catch
            var=ones(sz(1),sz(2)).*0.5;
        end;
        st(:,:,zz)=var;
        zz=zz+1;
        
        disp('');
        
        
        
    end;
        
end;