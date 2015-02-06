function syn_tomogram = tom_subtomo_template_position_to_tomogram(template,alignFile, pickList )
    templ = tom_emread(template); templ = templ.Value;
    pickList = tom_emread(pickList); pickList = pickList.Value;
    load(alignFile);
    syn_tomogram = single(zeros(pickList(2,1),pickList(3,1),pickList(4,1)));

    for i = 1:20%length(Align) 
        shift     = Align(i).Shift;
        rot       = Align(i).Angle;

        rotTempl = single(tom_rotate(templ,[ -rot.Phi, -rot.Psi, -rot.Theta]));
        tom_shift(rotTempl,[-shift.X -shift.Y -shift.Z]);
        rotTemplSmall = tom_bin(rotTempl,pickList(5,1));
        splitF = regexp(Align(i).Filename,'/','split');
        pname=splitF(length(splitF));
        partNum = str2double(regexprep(pname{1}, '[a-zA-Z_\.]+', ''));
        syn_tomogram = tom_paste(syn_tomogram, rotTemplSmall, [pickList(6,partNum),pickList(7,partNum),pickList(8,partNum)]);
        
        disp(['Particle: ' pname{1} ' positioned in tomogram ' ...
              num2str(pickList(2,1)) 'x' num2str(pickList(3,1)) 'x' num2str(pickList(4,1)) ...
              ' at X:' num2str(pickList(6,partNum)) ' Y:' num2str(pickList(7,partNum)) ' Z:' num2str(pickList(8,partNum)) ]);
    end;
    tom_emwrite(['syn_tomogram.em'], syn_tomogram  );
end