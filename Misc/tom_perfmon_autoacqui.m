function tom_perfmon_autoacqui(filenmane,filename2,refresh,data_window,plot_type)

if nargin < 2
    refresh=0;
end;


if nargin < 3
    data=tom_perfmon_importfile(filenmane);
    data_window=size(data,1);
end;



figure;

if (refresh==0)
    data=tom_perfmon_importfile(filenmane);
    subplot(3,1,1); plot(data(:,1)); title('Available Harddisk'); 
    subplot(3,1,2); plot(data(:,2),'r-'); title('Available Memoery');
    subplot(3,1,3); plot(data(:,3),'g-'); title('CPU in %');

else
    
    for i=1:20000
        disp('reading data');
        data=tom_perfmon_importfile(filenmane);
        neg_tilt=tom_emread([filename2 '/log_data/neg_tilt.em']);
        pos_tilt=tom_emread([filename2 '/log_data/pos_tilt.em']);
        ccf=tom_emread([filename2 '/log_data/ccf.em']);
        [a b]=unix('ls -1 /fs/titan/test_neg_stain2/high/ | wc -l');
        num=str2num(strrep(b,' ',''));
        try
            im_high=tom_emreadc([filename2 'high/26S_' num2str(num)  '.em'],'resample',[8 8 1]);

            im_low=tom_emreadc([filename2 'low/26S_' num2str(num)  '.em'],'resample',[8 8 1]);
        catch
            disp(['Error cannot read:' filename2  'high/26S_  ' num2str(num)  '.em']);
        end;
        
        sz=size(data,1);
        
        if (sz-data_window+1) > 0
            start=sz-data_window+1;
        else
           start=1;     
        end;
        
        
        stop=sz;
        subplot(4,3,1); %hold on; 
        plot(data(start:stop,1),'g-'); %hold off;  
        title('Available Harddisk');
        subplot(4,3,2); %hold on; 
        plot(data(start:stop,2),'r-'); %hold off;  
        title('Available Memoery');
        subplot(4,3,3); %hold on; 
        plot(data(start:stop,3),'b-'); %hold off; 
        title('Processor');
        pause(0.1);
        drawnow;
        pause(refresh);
        
        subplot(4,3,4); %hold on; 
        imagesc(neg_tilt.Value'); colormap gray;  axis image;%hold off;  
        title('neg tilt');
        
        subplot(4,3,5); %hold on; 
        imagesc(pos_tilt.Value'); colormap gray;  axis image;%hold off;  
        title('pos tilt');
        
        subplot(4,3,6); %hold on; 
        imagesc(ccf.Value'); colormap gray;  axis image;%hold off;  
        title('ccf');
        
        
        subplot(4,3,7); %hold on; 
        imagesc(im_high.Value'); colormap gray; axis image; %hold off;  
        title('image low');
        
        subplot(4,3,8); %hold on; 
        imagesc(im_low.Value'); colormap gray;  axis image;%hold off;  
        title('image high');
        
%         subplot(4,3,9); %hold on; 
%         tom_imagesc(ccf,'noinfo'); %hold off;  
%         title('ccf');
        
        %superuser#1
        
        
        
        
        
        
        
        
    end;

end;