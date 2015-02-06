function tom_av2_gen_amira_script(filename,vols,pos_flag,analysis_st)
%TOM_AV2_GEN_AMIRA_SCRIPT is a hardcoded config ... to be loaded as mat file yes !
%
%   tom_av2_gen_amira_script(filename,vols,pos_flag,analysis_st)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%   vols                ...
%   pos_flag            ...
%   analysis_st         ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_gen_amira_script(filename,vols,pos_flag,analysis_st)
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

 
const_st='';
analysis_st='';
loader_st='';
destructor_st.filename='test.hx';
analizer_st(1).type='isosurface';
analizer_st(1).thresh=-0.0328601;
analizer_st(2).type='orthoslice';

pos_st.flag='line';
%analizer_st(2).type='isosurface';









% loop over all volumes

out=constructor(const_st);

for i=1:size(vols,2)
    
    % initialize structures hurz !
    h=tom_reademheader(vols{i});
    loader_st.filename=vols{i};
    loader_st.size=[h.Header.Size'];
    loader_st.counter=i;
    loader_st.num_of_analysis=size(analizer_st,2);
    [pathinger label ext]=fileparts(vols{i});
    pos_st.label=label;
    pos_st.size=[h.Header.Size'];
    pos_st.num_of_vols=size(vols,2);
    pos_st.counter=i;
    
    for ii=1:size(analizer_st,2)
        analizer_st(ii).label=label;
        analizer_st(ii).counter=i;
    end;
    
    out=loader(loader_st,out);
    out=analizer(analizer_st,out);
    out=positioner(pos_st,out);



end;

destructor(destructor_st,out);






function out=constructor(const_st)

out{1}='# Amira Script';
out{2}='remove -all';
out{3}='remove model_1 Isosurface';
out{4}='';
out{5}='# Create viewers';
out{6}='viewer setVertical 0';
out{7}='';
out{8}='viewer 0 setBackgroundMode 1';
out{9}='viewer 0 setBackgroundColor 0.06 0.13 0.24';
out{10}='viewer 0 setBackgroundColor2 0.72 0.72 0.78';
out{11}='viewer 0 setTransparencyType 5';
out{12}='viewer 0 setAutoRedraw 0';
out{13}='viewer 0 show';
out{14}='mainWindow show';
out{15}='';



function out=loader(loader_st,out)

ind=size(out,2);

label_posy=(loader_st.counter.*(20.*loader_st.num_of_analysis));

[pathinger label ext]=fileparts(loader_st.filename);

out{ind+1}='set hideNewModules 0';
out{ind+3}=['[load -raw ' loader_st.filename ' little xfastest float 1 ' num2str(loader_st.size) ' ' num2str(0) ' ' num2str(loader_st.size(1)-1) ' ' num2str(0) ' ' num2str(loader_st.size(2)-1) ' ' num2str(0) ' ' num2str(loader_st.size(3)-1) ' -header 512 ] setLabel ' label];
out{ind+4}=[label ' fire'];
out{ind+5}=[label ' setIconPosition ' num2str(20) ' ' num2str(label_posy)];








function out=positioner(pos_st,out)

ind=size(out,2);


switch (pos_st.flag)

    case 'line'
        posx=(-((pos_st.size(1).*pos_st.num_of_vols)./2)) + (pos_st.counter.*pos_st.size(1));
        out{ind+1}=[pos_st.label ' setTransform 1 0 0 0 0 1 0 0 0 0 1 0 ' num2str(posx) ' 0 0 1'];
        
    case 'square'
        
    case 'cubic'

end;



function out=analizer(analizer_st,out)

ind=size(out,2);

%label_posy=analizer_st(1).counter.*20;
label_posy=((analizer_st(1).counter-1).*(20.*size(analizer_st,2) ));

for i=1:size(analizer_st,2)

   switch (analizer_st(i).type)
        case 'isosurface'
            label_posy=label_posy+20;
            out{ind+1}=['set hideNewModules 0'];
            out{ind+2}=['create HxIsosurface {Isosurface' num2str(analizer_st(i).counter) '}'];
            out{ind+3}=['Isosurface' num2str(analizer_st(i).counter) ' setIconPosition 359 ' num2str(label_posy)];
            out{ind+4}=['Isosurface' num2str(analizer_st(i).counter) ' data connect ' analizer_st(i).label];
            out{ind+5}=['Isosurface' num2str(analizer_st(i).counter) ' colormap setDefaultColor 1 0.8 0.4'];
            out{ind+6}=['Isosurface' num2str(analizer_st(i).counter) ' colormap setDefaultAlpha 0.500000'];
            out{ind+7}=['Isosurface' num2str(analizer_st(i).counter) ' fire'];
            out{ind+8}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setValue 1'];
            out{ind+9}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setSpecularLighting 1'];
            out{ind+10}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setTexture 0'];
            out{ind+11}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setAlphaMode 1'];
            out{ind+12}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setNormalBinding 1'];
            out{ind+13}=['Isosurface' num2str(analizer_st(i).counter) ' drawStyle setCullingMode 0'];
            out{ind+14}=['Isosurface' num2str(analizer_st(i).counter) ' threshold setMinMax -107.460914611816 47.4665031433105'];
            out{ind+15}=['Isosurface' num2str(analizer_st(i).counter) ' threshold setButtons 0'];
            out{ind+16}=['Isosurface' num2str(analizer_st(i).counter) ' threshold setIncrement 10.3285'];
            out{ind+17}=['Isosurface' num2str(analizer_st(i).counter) ' threshold setValue ' num2str(analizer_st(i).thresh)];
            out{ind+18}=['Isosurface' num2str(analizer_st(i).counter) ' threshold setSubMinMax -107.460914611816 47.4665031433105'];
            out{ind+19}=['Isosurface' num2str(analizer_st(i).counter) ' options setValue 0 1'];
            out{ind+20}=['Isosurface' num2str(analizer_st(i).counter) ' options setValue 1 0'];
            out{ind+21}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setMinMax 0 -1.00000001384843e+24 1.00000001384843e+24'];
            out{ind+22}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setValue 0 2'];
            out{ind+23}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setMinMax 1 -1.00000001384843e+24 1.00000001384843e+24'];
            out{ind+24}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setValue 1 2'];
            out{ind+25}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setMinMax 2 -1.00000001384843e+24 1.00000001384843e+24'];
            out{ind+26}=['Isosurface' num2str(analizer_st(i).counter) ' resolution setValue 2 2'];
            out{ind+27}=['Isosurface' num2str(analizer_st(i).counter) ' fire'];
            out{ind+28}=['Isosurface' num2str(analizer_st(i).counter) ' setViewerMask 65535'];
            out{ind+29}=['Isosurface' num2str(analizer_st(i).counter) ' select'];
            out{ind+30}=['Isosurface' num2str(analizer_st(i).counter) ' doIt snap 0 1'];
            out{ind+31}=['{Isosurface' num2str(analizer_st(i).counter) '} doIt hit'];
            out{ind+32}=['Isosurface'  num2str(analizer_st(i).counter) ' fire'];
        
       case 'orthoslice'
            label_posy=label_posy+20; 
           out{ind+1}=['set hideNewModules 0'];
           out{ind+2}=['create HxOrthoSlice {OrthoSlice' num2str(analizer_st(i).counter) '}'];
           out{ind+3}=['OrthoSlice' num2str(analizer_st(i).counter) ' setIconPosition 359 ' num2str(label_posy)];
           out{ind+4}=['OrthoSlice' num2str(analizer_st(i).counter) ' data connect ' analizer_st(i).label];
           out{ind+5}=['{OrthoSlice' num2str(analizer_st(i).counter) '} fire'];
           out{ind+6}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceOrientation setValue 0'];
           out{ind+7}=['{OrthoSlice' num2str(analizer_st(i).counter) '} fire'];
           out{ind+8}=['OrthoSlice' num2str(analizer_st(i).counter)  ' options setValue 0 1'];
           out{ind+9}=['OrthoSlice' num2str(analizer_st(i).counter) ' options setValue 1 0'];
           out{ind+10}=['OrthoSlice' num2str(analizer_st(i).counter) ' options setValue 2 0'];
           out{ind+11}=['OrthoSlice' num2str(analizer_st(i).counter) ' mappingType setValue 0 0'];
           out{ind+12}=['OrthoSlice' num2str(analizer_st(i).counter) ' linearRange setMinMax 0 -1.00000001384843e+24 1.00000001384843e+24'];
           out{ind+13}=['OrthoSlice' num2str(analizer_st(i).counter) ' linearRange setValue 0 -0.141395404934883'];
           out{ind+14}=['OrthoSlice' num2str(analizer_st(i).counter) ' linearRange setMinMax 1 -1.00000001384843e+24 1.00000001384843e+24'];
           out{ind+15}=['OrthoSlice' num2str(analizer_st(i).counter) ' linearRange setValue 1 0.0578559637069702'];
           out{ind+16}=['OrthoSlice' num2str(analizer_st(i).counter) ' contrastLimit setMinMax 0 -1.00000001384843e+24 1.00000001384843e+24'];
           out{ind+17}=['OrthoSlice' num2str(analizer_st(i).counter) ' contrastLimit setValue 0 7'];
           out{ind+18}=['OrthoSlice' num2str(analizer_st(i).counter) ' colormap setDefaultColor 1 0.8 0.5'];
           out{ind+19}=['OrthoSlice' num2str(analizer_st(i).counter) ' colormap setDefaultAlpha 1.000000'];
           out{ind+20}=['OrthoSlice' num2str(analizer_st(i).counter) ' colormap connect glow.col'];
           out{ind+21}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceNumber setMinMax 0 159'];
           out{ind+22}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceNumber setButtons 1'];
           out{ind+23}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceNumber setIncrement 1'];
           out{ind+24}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceNumber setValue 80'];
           out{ind+25}=['OrthoSlice' num2str(analizer_st(i).counter) ' sliceNumber setSubMinMax 0 159'];
           out{ind+26}=['OrthoSlice' num2str(analizer_st(i).counter) ' transparency setValue 0'];
           out{ind+27}=['OrthoSlice' num2str(analizer_st(i).counter) ' setFrameWidth 0'];
           out{ind+28}=['OrthoSlice' num2str(analizer_st(i).counter) ' setFrameColor 1 0.5 0'];
           out{ind+29}=['OrthoSlice' num2str(analizer_st(i).counter) ' frame 1'];
           out{ind+30}=['OrthoSlice' num2str(analizer_st(i).counter) ' fire'];

           out{ind+31}=['OrthoSlice' num2str(analizer_st(i).counter) ' fire'];
           out{ind+32}=['OrthoSlice' num2str(analizer_st(i).counter) ' setViewerMask 65535'];
           out{ind+33}=['OrthoSlice' num2str(analizer_st(i).counter) ' select'];

           out{ind+33}=['set hideNewModules 0'];
           
           
        case 'voltex'
    
    end;
    ind=size(out,2);
    

end;


out{ind+1}=['viewer 0 setCameraPosition 94.993 -247.417 16.7983'];
out{ind+2}=['viewer 0 setCameraOrientation 0.999463 0.0291039 0.0150786 1.76036'];
out{ind+3}=['viewer 0 setCameraFocalDistance 333.236'];
out{ind+4}=['viewer 0 setAutoRedraw 1'];
out{ind+5}=['viewer 0 redraw'];
        




function destructor(destructor_st,out)

fid = fopen(destructor_st.filename,'w');
for i=1:size(out,2)
    fprintf(fid, '%s \n',out{i});
end;
fclose(fid);









