function tom_ctfmap(filename,start,stop,clusterno,matfile,quality)

s = load(matfile);
align2d = s.align2d;
clear('s');

mapstruct = struct();

min_x = 999999;
max_x = -999999;
min_y = 999999;
max_y = -999999;
k = 1;
for i=start:stop
   disp(num2str(i));
    %read em header and extract position of image on grid
    h = tom_reademheader([filename num2str(i) '.em']);
    [xx,pos] = strtok(char(h.Header.Comment'), ';');
    [pos,remain] = strtok(pos,';');
    [xx,nominal] = strtok(remain,';');
    [nominal,cluster] = strtok(nominal,'!');
    nominal = nominal(2:end);
    cluster = cluster(2:end);
    [pos_x,remain] = strtok(pos,' ');
    [pos_y,pos_z] = strtok(remain,' ');

    %find filename in structure
    for kk=1:length(align2d)
        [pathstr, name, ext, versn] = fileparts([filename num2str(i) '.em']);
        if strcmp(align2d(kk).filename, [name ext]) ==1
            matidx = kk;
            break;
        end
    end   
    
    if clusterno == str2num(cluster') && align2d(matidx).quality == quality
        mapstruct(k).filename = [filename num2str(i) '.em'];
        mapstruct(k).pos.x = str2num(pos_x');
        mapstruct(k).pos.y = str2num(pos_y');
        mapstruct(k).pos.z = str2num(pos_z');
        mapstruct(k).nominal = str2num(nominal');
        mapstruct(k).fitted = h.Header.FocusIncrement./10000;
        mapstruct(k).cluster = str2num(cluster');
        %fit DZ
    %    [Dz,success] = tom_ctffitter2(mapstruct(k).filename,1,0,7,150);
    %    if success == 1
    %        mapstruct(k).fitted = Dz;
    %    else
    %        mapstruct(k).fitted = NaN;
    %    end

        %find position extrema
        if mapstruct(k).pos.x < min_x
            min_x = mapstruct(k).pos.x;
        elseif mapstruct(k).pos.x > max_x
            max_x = mapstruct(k).pos.x;
        end

        if mapstruct(k).pos.y < min_y
            min_y = mapstruct(k).pos.y;
        elseif mapstruct(k).pos.y > max_y
            max_y = mapstruct(k).pos.y;
        end
        map_x(k) = mapstruct(k).pos.x;
        map_y(k) = mapstruct(k).pos.y;
        try
            map_z(k) = mapstruct(k).nominal+mapstruct(k).pos.z;
            tmp_z(k)=mapstruct(k).pos.z;
            map_def_only(k)=mapstruct(k).nominal;
            
        catch
            disp(['skipping image ' filename num2str(i) '.em']);
        end
        map_z2(k) = mapstruct(k).fitted+mapstruct(k).pos.z;
        map_fit_only(k)=mapstruct(k).fitted;
        error_def(k)=mapstruct(k).fitted-(-4);
        cor_def(k)=(mapstruct(k).fitted-(-4))+mapstruct(k).nominal+mapstruct(k).pos.z;
        %map_z2(k) = mapstruct(k).pos.z;
        
        
        k = k + 1;
    end
end

l = 1;
for i=1:length(map_z2)
    if ~isnan(map_z2(i))
        map_z3(l) = map_z2(i);
        map_x3(l) = map_x(i);
        map_y3(l) = map_y(i);
        tmp_z3(l)= tmp_z(i);
        l = l + 1;
    end
end

[XI,YI] = meshgrid([min_x:2.24:max_x],[min_y:2.24:max_y]);
ZI = griddata(map_x,map_y,map_z,XI,YI);
ZI3 = griddata(map_x3,map_y3,map_z3,XI,YI);
ZZZ= griddata(map_x3,map_y3,tmp_z3,XI,YI);
def=griddata(map_x,map_y,map_def_only,XI,YI);
def_fitted=griddata(map_x,map_y,map_fit_only,XI,YI);
def_corrected=griddata(map_x,map_y,cor_def,XI,YI);
%plot nominal Dz
figure; set(gcf,'Name','nominal+Z');
mesh(XI,YI,ZI);
%surfl(XI,YI,ZI,'light');shading interp;

%hold on;plot3(map_x,map_y,map_z,'ro');grid on; colormap cool;
%contour3(XI,YI,ZI,40,'-k');

%plot fitted Dz
figure; set(gcf,'Name','fitted + Z');
% surfl(XI,YI,ZI3,'light');shading interp;
mesh(XI,YI,ZI3);
%hold on;plot3(map_x3,map_y3,map_z3,'ro');grid on; colormap cool;
%contour3(XI,YI,ZI3,40,'-k');

%plot z

figure; set(gcf,'Name','Z');
%surfl(XI,YI,ZZZ,'light');shading interp;
meshc(XI,YI,ZZZ);
%hold on;plot3(map_x3,map_y3,map_z3,'ro');grid on; colormap cool;


%plot defocus
figure; set(gcf,'Name','defocus');
%surfl(XI,YI,ZZZ,'light');shading interp;
surfc(XI,YI,def);
%hold on;plot3(map_x3,map_y3,map_z3,'ro');grid on; colormap cool;


%plot defocus fitted

figure; set(gcf,'Name','defocus corrected+Z');
%surfl(XI,YI,ZZZ,'light');shading interp;
meshc(XI,YI,def_corrected);
%hold on;plot3(map_x3,map_y3,map_z3,'ro');grid on; colormap cool;





















