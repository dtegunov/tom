function [nr shift]=tom_av3_search_index(index_path,part)



level=[];



for i=1:10000
    
    %load left and right
    left_name=strrep([index_path 'idx_ ' num2str(level) num2str(0) '.em'],' ','');
    right_name=strrep([index_path 'idx_' num2str(level) num2str(1) '.em'],' ','');
    
    if (exist(left_name,'file')==0 | exist(right_name,'file')==0 )
        break;
    end;
        
    
    idx_left=tom_emread(left_name); 
    idx_right=tom_emread(right_name); 
    
    disp(left_name);
    disp(right_name);
    
    
    
    cc=tom_corr(part,idx_left.Value,'norm');
    [a b]=tom_peak(cc);
    vals(1)=b;
    cc=tom_corr(part,idx_right.Value,'norm');
    [a b]=tom_peak(cc);
    vals(2)=b;
    [a b]=max(vals);
    
%     subplot(2,3,1); tom_imagesc(idx_left); title('left');
%     subplot(2,3,2); tom_imagesc(idx_right); title('right');
%     subplot(2,3,3); plot(vals); title('ccs');
%     subplot(2,3,5); tom_imagesc(part); title('particle');
    
    
    
    level=[level num2str(b-1)];
    
end;

if (b==1)
    nr=str2num(char(idx_left.Header.Comment)');
else
    nr=str2num(char(idx_right.Header.Comment)');
end;


figure;

 subplot(2,3,1); tom_imagesc(idx_left); title('left');
    subplot(2,3,2); tom_imagesc(idx_right); title('right');
    subplot(2,3,3); plot(vals); title('ccs');
    subplot(2,3,5); tom_imagesc(part); title('particle');
    





