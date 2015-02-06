function tom_tilt_lines(alig,ref,pkt);
%TOM_TILT_LINES plots marker points
%
%   tom_tilt_lines(alig,ref,pkt)
%
%   TOM_TILT_LINES plots the difference vectors (2D) of a marker point PKT
%       in respect to a reference marker point REF. The positions of the
%       marker point are in ALIG (for format see tom_alignment3d). The
%       position vectors should be located on a straight line. Incorrectly
%       clicked markers can be identified by locating outliers.
%
%PARAMETERS
%
%  INPUT
%   alig                Alignment array (from marker file)
%   ref                 Reference marker point
%   pkt                 Which Point to be plotted against ref
%                        if == ref: all marker points are plotted
%  
%  OUTPUT
%
%EXAMPLE
%   tom_tilt_lines(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_ALIGNMENT3D, TOM_SETMARK, TOM_REC3D
%
%   created by ... (author date)
%   updated by ...
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


txt=1;
all=0;
if ref==pkt
    all=1;
    txt=0;
end
[s1,s2,s3]=size(alig);
if all==0
  plot(0,0,'r+');hold on;zoom on;
  inds=find((alig(2,:,ref) > -1) & (alig(3,:,ref) > -1) & (alig(2,:,pkt) > -1) & (alig(3,:,pkt) > -1));
  alig_x=alig(2,inds,ref)-alig(2,inds,pkt);
  alig_y=alig(3,inds,ref)-alig(3,inds,pkt);
  for kk=1:size(inds,2)
    plot(alig_x(kk),alig_y(kk),'r+');
    if txt==1
      text('Position',[alig_x(kk)+1 alig_y(kk)+1],'String',alig(11,inds(kk)),1); 
    end
  end
grid off;
hold on;

end

if all==1
  clf; plot(0,0,'r+');hold on;
  for ref=1:s3
    for pkt=ref+1:s3
        inds=find((alig(2,:,ref) > -1) & (alig(3,:,ref) > -1) & (alig(2,:,pkt) > -1) & (alig(3,:,pkt) > -1));
        alig_x=alig(2,inds,ref)-alig(2,inds,pkt);
        alig_y=alig(3,inds,ref)-alig(3,inds,pkt);
        plot(alig_x,alig_y,'+');
        text('Position',[alig_x(1) alig_y(1)],'String',ref,'Fontsize',20);   
        text('Position',[alig_x(1)+40 alig_y(1)-40],'String',pkt,'Fontsize',20);   
    end
  end
end
xlabel(' X-difference');
ylabel(' Y-difference');

