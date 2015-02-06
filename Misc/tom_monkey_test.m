function tom_monkey_test(nr,tim)

%TOM_MONKEY_TEST is the monkey test for particle picking.
%
%   tom_monkey_test(nr,tim)
%
%PARAMETERS
%
%  INPUT
%   nr                  nr of numbers
%   tim                 time before the numbers are hidden
%                       in seconds
%
%EXAMPLE
%   tom_monkey_test(10,1)
%   % 10 numbers, 1 second
%
%REFERENCES
%
%SEE ALSO
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


f=figure;
%set(gcf,'Renderer','Opengl');
b=zeros(2000,1500);
bl=zeros(2000,1500);
bo=ones(120,120);
bo_b=zeros(120,120);
image(b');colormap gray;
set(f,'Position',[350           5        1350        2005]);
set(f,'Units','points');
axis off;
set(f,'DoubleBuffer','on');
for i=1:nr
    m=2;
    while m>1
    x=randi(size(b,1)-120)+60;
    y=randi(size(b,2)-120)+60;
    xc(i)=x;
    yc(i)=y;
    bc=tom_paste(bl,bo,[xc(i)-20 yc(i)-60]);
    bc=bc+b;
    m=max(max(bc));
    end;

    c(i,1)=x;
    c(i,2)=y;
    t(i)=text(x,y,num2str(i-1));
    set(t(i),'FontSize',64);set(t(i),'Color','White');
    b=tom_paste(b,bo,[xc(i)-20 yc(i)-60]);
end;
%g=getimage;

pause(tim);
imagesc(b');colormap gray;axis off;
tic;
for i=1:nr
    bl=0;
    while bl==0
%         keydown = waitforbuttonpress;
%         if (keydown == 0)
%             disp('Mouse button was pressed');
%         else
%             disp('Key was pressed');
%         end
%         p=get(gcf,'CurrentPoint')
                p=ginput(1);
        
        [idx, coords] = tom_nearestpoint(p, c);
        bl=b(round(p(1)),round(p(2)));
    end;
    if i~=idx
        disp('wrong!');
        break;
    end;
    b=tom_paste(b,bo_b,[xc(i)-20 yc(i)-60]);
    imagesc(b');colormap gray;axis off;
end;
toc

