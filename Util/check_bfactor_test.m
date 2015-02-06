function b=check_bfactor_test(p,range)

% calculates the bfactor for all EM files in folder 'p'.
% parallel support.
% example:
% fit range is 3.5 to 7 Angstrom:
%
% b=check_bfactor_test('/fs/sun16/lv01/pool/pool-nickell3/26S/em/data/Titan/CWT/14022011/15022011/low',[3.5 7])

dt=dir(p);
d=dt(3:end);

for i=1:size(d,1)-500
%parfor i=1:size(d,1)-500
%    tic
    b(i).name=[p filesep d(i).name];
%    try
        in=tom_emreadc3([p filesep d(i).name]);
        b(i).valid=1;
        [bfactor]=tom_fit_bfactor2D(in.Value, in.Header.Objectpixelsize, range,256,0);
        b(i).bfactor=bfactor;
        [result quality]=tom_check_image_quality(in);
        b(i).result=result;
        b(i).quality=quality;
%    catch
%        b(i).bfactor=0;
%        b(i).valid=0;
%    end;
%    toc
end;
