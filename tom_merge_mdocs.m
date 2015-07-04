function [ ] = tom_merge_mdocs( filename1, filename2, mergedname )

mdoc1 = tom_read_mdoc(filename1);
mdoc2 = tom_read_mdoc(filename2);

wrong = false;
wrong = wrong || mdoc1.pixelsize ~= mdoc2.pixelsize;
wrong = wrong || mdoc1.size(1) ~= mdoc2.size(1);
wrong = wrong || mdoc1.size(2) ~= mdoc2.size(2);
if wrong
    error('Image dimensions do not match');
end;

merged = mdoc1;
mdoc2.zvalue = mdoc2.zvalue + max(mdoc1.zvalue(:)) + 1;
merged.zvalue = [merged.zvalue mdoc2.zvalue];
merged.numframes = [merged.numframes mdoc2.numframes];
merged.rot = [merged.rot mdoc2.rot];
merged.theta = [merged.theta mdoc2.theta];
merged.defocus = [merged.defocus mdoc2.defocus];
merged.exptime = [merged.exptime mdoc2.exptime];
merged.expdose = [merged.expdose mdoc2.expdose];
merged.shift = [merged.shift mdoc2.shift];
merged.timestamp = [merged.timestamp mdoc2.timestamp];

mergedstackname = mergedname(1:end-5);
merged.imagefile = mergedstackname;

tom_write_mdoc(merged, mergedname);

end