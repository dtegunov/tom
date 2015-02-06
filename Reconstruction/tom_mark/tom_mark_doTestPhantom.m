




markerset = tom_emread('/fs/sandy01/lv03/pool/bmsan/apps/tom_dev/data/Tiltseries/mark_pyrodictium.em');
markerset = markerset.Value;



tiltangles = markerset(1,:,1);
markerset = markerset(2:3,:,:);
markerset(1:2, any(markerset == -1, 1)) = nan;

irefmark = 1;
imdim = 1024;
irefmarkX = [];

volname = '/fs/sandy01/lv03/pool/bmsan/apps/tom_dev/data/Tiltseries/test.vol';
volname = '/fs/sandy01/lv03/pool/bmsan/apps/tom_dev/data/Tiltseries/q.vol';
vol = tom_emreadc(volname, 'binning', 0);

X = [[200, 111, 63];
     [158, 192, 75];
     [100, 234, 81];
     [222, 135, 68];
     [167,  62, 58];
     [ 84,  65, 59];
     [126, 207, 74];
     ]' * 1024/vol.Header.Size(1);


[psi, tx, ty, sigma, align_X]  = tom_mark_cvaf_alignment3d(markerset, tiltangles, irefmark, X(:,irefmark), imdim);


[P, x_proj, X_trian, x_reproj] = tom_mark_cvaf_alignment3d_reproj(tiltangles, psi, tx, ty, imdim, X, markerset);


