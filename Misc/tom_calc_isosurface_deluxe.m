function [thresh_all mass vol vol_rest vol_and_rest mass2 thresh]=tom_calc_isosurface_deluxe(volume,molmass,pixelsize,accuracy,start_coordinate,radius,add_molmass,symmetry)




mask=tom_cylindermask(ones(size(volume)),radius);

vol_org=volume;
volume=volume.*mask;

[mass thresh]=iterate_mass(vol_org,molmass,pixelsize,'threshold',accuracy);

[mass stop]=iterate_mass(volume,add_molmass,pixelsize,'z++',start_coordinate,thresh);

thresh_all=thresh;



start=start_coordinate;


for i=1:1
     disp(['start: ' num2str(start) '  stop: ' num2str(stop) ]);  
    
    [mass thresh]=iterate_mass(tom_symref(volume(:,:,start:stop),symmetry),add_molmass,pixelsize,'threshold',accuracy);   
    %[mass start stop]=iterate_mass(tom_symref(volume(:,:,:),6),add_molmass,pixelsize,'+m+',round(start+stop)./2,thresh);

    
end;

vol_rest=zeros(size(volume));


vol_rest(:,:,start:stop)=(tom_symref(volume(:,:,start:stop),symmetry)>thresh).*vol_org(:,:,start:stop);
vol=(vol_org>thresh_all).*vol_org;
vol_and_rest=vol_rest + vol(:,:,1:start_coordinate);
mass2=mass;

disp('end');




function [mass out1 out2]=iterate_mass(volume,target_mass,pixelsize,flag,param1,param2)


z=1;

mass=max(max(max(volume)));
lauf=max(max(max(volume)));

if (strcmp(flag,'+m+'))
    start=param1;
    stop=param1;
    zsta=1;
    zsto=1;
end;



while mass<target_mass
    
    if (strcmp(flag,'threshold'))
        lauf=lauf-param1;
        start=1;
        stop=size(volume,3);
    end;
    
     if (strcmp(flag,'z++'))
        lauf=param2;
        start=param1;
        stop=start+z;
        z=z+1;
    end;
    
     if (strcmp(flag,'+m+'))
        lauf=param2;
        if mod(z,2)==1
            start=param1-zsta;
            zsta=zsta+1;
        else
            stop=param1+zsto;
            zsto=zsto+1;
        end;
        z=z+1;
    end;
    
    
    nrmass=sum(sum(sum(volume(:,:,start:stop)>=lauf)));
    rho=1.3; % g/cm^3 for protein
    pixelsize_cm3=(pixelsize.*1.0e-10).^3./(0.01.^3);
    nr=nrmass;
    mass_g=nr.*rho.*pixelsize_cm3;
    mass_kg=mass_g./(10.^3);
    mass=mass_kg./(1.66e-27)./1000;
end;




if (strcmp(flag,'threshold') )
    out2=0;
    out1=lauf;
end;


if (strcmp(flag,'z++') )
    out2=0;
    out1=start+z;
end;

if (strcmp(flag,'+m+') )
    out2=stop;
    out1=start;
end;





