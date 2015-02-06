function b=check_quality(p)

% calculates the bfactor for all EM files in folder 'p'.
% parallel support.
% example:
% fit range is 3.5 to 7 Angstrom:
%
% b=check_quality('/fs/sun16/lv01/pool/pool-nickell3/26S/em/data/Titan/CWT/14022011/15022011/low')

d=dir([p '/*.em'])


%parfor i=1:size(d,1)
parfor i=1:30
%    tic
    disp(['processing: ' [p filesep d(i).name] ', #' num2str(i) ' out of ' num2str(30) '.']); 
    try
        in=tom_emreadc3([p filesep d(i).name]);
        [me,ma,mi,st,va]=tom_dev(in.Value,'noinfo');
        b(i).values.mean=me;
        b(i).values.max=ma;
        b(i).values.min=mi;
        b(i).values.std=st;
        b(i).values.var=va;
        [result quality]=tom_check_image_quality(in);
        b(i).result=result;
        b(i).quality=quality;
        b(i).valid=1;
    catch
        b(i).bfactor=0;
        b(i).valid=0;
    end;
%    toc
end;
