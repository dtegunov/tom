function tom_av2_index_analyse_tree(align2d,flag,lookup,im_st,param)


if strcmp(flag,'ref_and_min_ratio') ||   strcmp(flag,'all')
        
    zz=1;
    for i=1:size(align2d,2)
        try
            rat(i)=align2d(1,i).ratio_min;
        catch
            rat(i)=10;
        end;
        
        cl(i)=align2d(1,i).ref_class;
        if (zz~=cl(i))
            errorp(i)=rat(i);
        end;
        if (mod(i,10)==0)
            zz=zz+1;
        end;
    end;
    mean_rat(1:i)=mean(rat); 
    std_rat_up(1:i)=mean(rat)+0.7*std(rat);
    std_rat_low(1:i)=mean(rat)-0.7*std(rat);
    max_rat=max(rat);
    figure; set(gcf,'Name','min');
    plot(rat);
   % hold on; plot(tom_norm(cl,max_rat),'r-'); hold off;
    hold on; plot(errorp,'ro'); hold off;
    hold on; plot(mean_rat,'b-'); hold off;
    hold on; plot(std_rat_low,'g-'); hold off;
    hold on; plot(std_rat_up,'g-'); hold off;
end;



if strcmp(flag,'all_ratios') ||   strcmp(flag,'all')
    figure;
    ff=zeros(size(align2d,2),6); 
    for i=1:size(align2d,2) 
        ff(i,1:length(align2d(1,i).ratio))=align2d(1,i).ratio; 
    end; 
    surfc(ff); shading interp; colormap hot;
end;


if strcmp(flag,'corr_all_ratios') || strcmp(flag,'all')
    ff=zeros(size(align2d,2),6); 
    for i=1:size(align2d,2) 
        ff(i,1:length(align2d(1,i).ratio))=align2d(1,i).ratio; 
    end; 
    avg=sum(ff,1);
    for i=1:size(align2d,2) 
         cc(i)=tom_ccc(avg,ff(i,:),'norm'); 
    end; 
    for i=1:size(align2d,2)
        rat(i)=align2d(1,i).ratio_min;
        cl(i)=align2d(1,i).ref_class;
    end;
    figure; plot(cc); 
     hold on; plot(tom_norm(cl,1),'g-'); hold off;
    %hold on; plot(tom_norm(cl.*rat,1),'bl-'); hold off;
    legend('correlation','class');
    figure; plot(tom_norm(cl,1));
    hold on; plot(tom_norm(rat,1),'r-'); hold off;
    legend('min ratio','class');
    figure; plot(tom_norm(cl,1));
    hold on; plot(tom_norm(rat.*cc,1),'k-'); hold off;
    legend('min ratio .* cl','class');
    figure; plot(avg);
end;


if strcmp(flag,'calc_erros') || strcmp(flag,'all')

    zz=1;    
    error_c=0;
    for i=1:size(align2d,2)
        if (align2d(1,i).ref_class~=zz)
            error_c=error_c+1;
        end;
        if (mod(i,1)==0)
            zz=zz+1;
        end;
    end;
    error_c(2)=size(align2d,2);
    figure; bar(error_c);
    disp(['Wrong classified particles: ' num2str((error_c(1)./error_c(2)).*100 )  ' %  (' num2str(error_c(1)) ' of '  num2str(error_c(2)) ' )']);
    title(['Wrong classified particles: ' num2str((error_c(1)./error_c(2)).*100 ) ' %  (' num2str(error_c(1)) ' of ' num2str(error_c(2))  ' )' ]);
end;

if strcmp(flag,'track_ref') || strcmp(flag,'all')
    
    num=0;
    base=max([lookup{:,3}]);
    for i=1:length(lookup)
        if (find(lookup{i,2}==param))
            disp(num2str(lookup{i,1}));
            num=num+1;
            idx=tom_av2_index_bintree_not_2_index(num2str(lookup{i,1}),base );
            refs(:,:,num)=im_st(:,:,idx);
            last_i=i;
        end;
   end;
    figure; tom_imagesc(tom_gallery(refs,[num 1]),'noinfo'); title(num2str(lookup{last_i,1}));
end;









