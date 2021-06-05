function [ R_Error, R2_Error, R4_Error ] = foveaLocalizationRefuge1(retDir, dataFile)

% retDir: Path to the folder with Refuge images (available for participants at Refuge Challenge, https://refuge.grand-challenge.org/)
% groundTruthData: Path to the file with the groundtruth x, y coords of the center of the fovea (available for participants at Refuge Challenge)

warning('OFF', 'all')

imgFiles = dir(fullfile(retDir, '*.jpg'));
groundTruthData = importdata(dataFile);
nImages = 400;
a = 2;
b = 1;
threshold = 2;
euclideanError = zeros(1, nImages);
R_Error = zeros(1, nImages);
R2_Error = zeros(1, nImages);
R4_Error = zeros(1, nImages);
    
for i = 1:nImages
    i
    retinalImage = imread(fullfile(retDir, imgFiles(i).name));
    nrows = size(retinalImage, 1);
    ncols = size(retinalImage, 2);
    retinalImageG = retinalImage(:, :, 2);

    % This piece of code corresponds to the method for optic disc localization and segmentation published in:
    % "Contrast based circular approximation for accurate and robust optic disc segmentation in retinal images"
    scaleFactor = 0.3455;
    retinalImage1 = imresize(retinalImage, scaleFactor);
    [ xCenter, yCenter ] = OpticDiscLocUsingFilters(retinalImage1); % % Get estimated optic disc center
    nrows1 = size(retinalImage1, 1);
    ncols1 = size(retinalImage1, 2);
    xCenter = round(xCenter/(ncols1/ncols));
    yCenter = round(yCenter/(nrows1/nrows));
    FOD = FOV_Diameter( retinalImage );
    odRadius = 0.15*((FOD)/2);
    
    % The previous anatomical information allows to estimate the center of the fovea localization
    % The original image is cropped around the estimated center
    minYFoveaBox = yCenter - (a/4) * odRadius;
    maxYFoveaBox = yCenter + a * odRadius;
    if xCenter <= round(ncols/2)
        minXFoveaBox = xCenter + (2.5 * 2 * odRadius) - (0.75* b * odRadius);
        maxXFoveaBox = xCenter + (2.5 * 2 * odRadius) + (0.75 * b * odRadius);
    else
        minXFoveaBox = xCenter - (2.5 * 2 * odRadius) - (0.75 * b * odRadius);
        maxXFoveaBox = xCenter - (2.5 * 2 * odRadius) + (0.75 * b * odRadius);
    end
    foveaImage = retinalImageG(minYFoveaBox:maxYFoveaBox, minXFoveaBox:maxXFoveaBox);
    foveaImage = imgaussfilt(foveaImage, 1);
    
    % The spatial-color histogram for x, G is computed
    schxG = histc(foveaImage, 0:255);
    schxG1 = schxG > threshold; % The histogram is thresholded
    CC = bwconncomp(schxG1);
    L = labelmatrix(CC);
    areaProps = regionprops(L, 'Area');
    areas = [ areaProps.Area ];
    [ maxAreas, indMaxAreas ] = max(areas);
    schxG2 = ismember(L, indMaxAreas);  % The largest connected component is selected
    [ yCoord, xCoord ] = find(schxG2);
    [ minyCoord, indminyCoord ] = min(yCoord);  % The x coordinate of the the pixel (x, G) belonging to the connected component with G = minimum(G) 
    xFovea = xCoord(indminyCoord) + (minXFoveaBox - 1); % The obtained value for x is rescaled to the dimensions of the original image
    
    % The spatial-color histogram for y, G is computed
    schyG = histc(foveaImage', 0:255);
    schyG1 = schyG > threshold; % The histogram is thresholded
    CC = bwconncomp(schyG1);
    L = labelmatrix(CC);
    areaProps = regionprops(L, 'Area');
    areas = [ areaProps.Area ];
    [ maxAreas, indMaxAreas ] = max(areas);
    schyG2 = ismember(L, indMaxAreas);  % The largest connected component is selected
    [ yCoord, xCoord ] = find(schyG2);
    [ minyCoord, indminyCoord ] = min(yCoord);  % The y coordinate of the the pixel (y, G) belonging to the connected component with G = minimum(G)
    yFovea = xCoord(indminyCoord) + (minYFoveaBox - 1); % The obtained value for y is rescaled to the dimensions of the original image    
    
    euclideanError(i) = sqrt((xFovea - groundTruthData.data(i, 4))^2 + (yFovea - groundTruthData.data(i, 5))^2);
    R = 0.15*((1563.5)/2);
    if euclideanError(i) <= R
        R_Error(i) = 1;
    end
    if euclideanError(i) <= R/2
        R2_Error(i) = 1;
    end
    if euclideanError(i) <= R/4
        R4_Error(i) = 1;
    end
end
   
end