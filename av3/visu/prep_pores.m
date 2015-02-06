for ii=1:size(indx)
    ifile = indx(ii);
    name = ['_' num2str(ifile) '.em'];
    part = tom_emread(name);part = tom_filter(part.Value,3,'circ');
    name = ['_' num2str(ifile) '_filt.em'];
    tom_emwrite(name, part);
    phi_opt = motl(17,indx(ii));
    psi_opt = motl(18,indx(ii));
    the_opt = motl(19,indx(ii));
    tshift = motl(11:13,indx(ii))';
    im = double(tom_rotate3d(tom_shift(part,-tshift),-psi_opt,-phi_opt,-the_opt));
    name = ['_' num2str(ifile) '_alig.em'];
    tom_emwrite(name, im);
end;