function out=tom_bfactor_filter3d(in,objectpixelsize,bfactor,apply_range,step)

out=in;

for ix=1:step:size(in,1)-step
    for iy=1:step:size(in,2)-step
        for iz=1:step:size(in,3)-step
        inb=in(ix:ix+step-1,iy:iy+step-1,iz:iz+step-1);
        [corrected decay_restore decay_restore_3d]=tom_apply_bfactor(inb,objectpixelsize,bfactor,1,apply_range);
        out(ix:ix+step-1,iy:iy+step-1,iz:iz+step-1)=corrected;
        end;
    end;
end;
