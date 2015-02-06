function [res mat_out]=tom_compare_deluxe2(fsc)

crit = 0.5;

sz = 96;

ny = (2.1).*2;

mat(:,9) = fsc';

num_of_sh = 48;

for i=1:num_of_sh
    
    res=(num_of_sh./i).*ny;
    
    data(i,1)=1./res;
    
    if (isnan(mat(i,9)) )
        data(i,2)=0.999;
    else
        data(i,2)=mat(i,9);
    end
    
    data(i,3)=-1;
    
    data(i,4)=res;
    
end;



    
    h=figure;
    
    p1=plot(data(:,1),data(:,2),'Linewidth',2); hold on;
    p2=plot(data(:,1),data(:,2),'ro','MarkerSize',3,'Linewidth',2);
    p=get(p1,'Parent');
    %set(fig,'Position',[80 12 60 20]);
    set(p,'Ylim',[0 1.1]);
    %set(t,'Color','white');
    xlabel('Resolution [1/Ang]');
    set(p,'Ytick',[0:0.1:1.1]);
    set(p,'Xtick',[0:data(end,1)./5:data(end,1)]);
    ylabel('FSC');
    grid on;
    p=get(p1,'Parent');
    xtick=get(p,'Xtick');
    xtick_nm=zeros(size(xtick));
    xtick_nm(2:end)=(1./xtick(2:end));
    
    for i=1:size(data(:,2),1)
        if data(i,2)<0.5
            v1=data(i-1,1);
            v2=data(i,1);
            v=(v1+v2)./2;
            %        text(v,0.5,sprintf('%0.3g A',1./v));
            text(0,0.5,sprintf('  %0.3g A >>>',1./v));
            v_05=v;
            break;
        end;
    end;
    
    
    
    
    for i=1:size(data(:,2),1)
        if data(i,2)<0.3
            v1=data(i-1,1);
            v2=data(i,1);
            v=(v1+v2)./2;
            %        text(v,0.3,sprintf('%0.3g A',1./v));
            text(0,0.3,sprintf('  %0.3g A >>>',1./v));
            v_03=v;
            break;
        end;
    end;
    
    
    for i=2:size(xtick,2)
        t=text(xtick(i),0.05,sprintf('%0.3g A',xtick_nm(i)));
    end;
    
    title(['Fourier Shell Correlation, Nyquist @ ' num2str(data(end,4)) ' Ang. FSC: 0.5 @ ' sprintf('%0.3g',1./v_05) ' A, 0.3 @ ' sprintf('%0.3g',1./v_03) ' A.']);

    


for i=1:size(data(:,2),1)
    if data(i,2)<crit
        v1=data(i-1,1);
        v2=data(i,1);
        v=(v1+v2)./2;
        %        text(v,0.5,sprintf('%0.3g A',1./v));
        
        res=1./v;
        break;
    end;
end;


mat_out=data;










% p1=plot(data(:,1),data(:,2),'Linewidth',2); hold on;
% ax=get(h,'children');
%
% p2=plot(data(:,1),data(:,2),'ro','MarkerSize',2,'Linewidth',2); hold off;
% %set(handles.display,'Position',[80 12 60 20]);
% set(ax,'Ylim',[0 1.1]);
% t=title(ax,['Fourier Shell Correlation, Nyquist @ ' num2str(data(end,4)) ' Ang']);
% set(t,'Color','white');
% xlabel(ax,'Resolution [1/Ang]');
% set(ax,'Ytick',[0:0.1:1.1]);
% set(ax,'Xtick',[0:data(end,1)./5:data(end,1)]);
% ylabel(ax,'FSC');








