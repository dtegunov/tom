function track(filename,numproj)

thresh = 5;

matrix = zeros(12,numproj,1);

markercounter = 1;

for i=1:numproj-1
    
    baseImage = tom_emreadc([filename num2str(i) '.em']);
    registeredImage = tom_emreadc([filename num2str(i+1) '.em']);
    
    disp(sprintf('Tracking between %g deg and %g deg',baseImage.Header.Tiltangle,registeredImage.Header.Tiltangle));
    
    baseImage.Value = single(baseImage.Value);
    registeredImage.Value = single(registeredImage.Value);
    
    baseImage.Value = tom_norm(baseImage.Value,255);
    registeredImage.Value = tom_norm(registeredImage.Value,255);
    
    [fa,da] = vl_sift(baseImage.Value);
    [fb,db] = vl_sift(registeredImage.Value);
    
    [matches, scores] = vl_ubcmatch(da,db,thresh);
    
    [unused, perm] = sort(scores, 'descend');
    matches = matches(:, perm);
    
    controlPointsBase = [fa(1,matches(1,:))', fa(2,matches(1,:))'];
    controlPointsRegistered = [fb(1,matches(2,:))', fb(2,matches(2,:))'];
    
    disp(sprintf('%g control point pairs found.',size(controlPointsBase,1)));
    
    for j=1:size(controlPointsBase,1)
        matrix(1,i,markercounter) = baseImage.Header.Tiltangle;
        matrix(2,i,markercounter) = controlPointsBase(j,1);
        matrix(3,i,markercounter) = controlPointsBase(j,2);
        
        matrix(1,i+1,markercounter) = registeredImage.Header.Tiltangle;
        matrix(2,i+1,markercounter) = controlPointsRegistered(j,1);
        matrix(3,i+1,markercounter) = controlPointsRegistered(j,2);
        
        markercounter = markercounter + 1;
    end
    
    
    
end

tom_emwrite('markerfile_track.em',matrix);

