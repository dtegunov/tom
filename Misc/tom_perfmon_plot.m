function tom_perfmon_plot(filenmane,refresh,data_window,plot_type)

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
    
    for i=1:200
        disp('reading data');
        data=tom_perfmon_importfile(filenmane);
        sz=size(data,1);
        start=sz-data_window+1;
        stop=sz;
        subplot(3,1,1); %hold on; 
        plot(data(start:stop,1),'g-'); %hold off;  
        title('Available Harddisk');
        subplot(3,1,2); %hold on; 
        plot(data(start:stop,2),'r-'); %hold off;  
        title('Available Memoery');
        subplot(3,1,3); %hold on; 
        plot(data(start:stop,3),'b-'); %hold off; 
        title('Processor');
        pause(0.1);
        drawnow;
        pause(refresh);
        
    end;

end;