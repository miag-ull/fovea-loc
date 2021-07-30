function [ R_Error, R2_Error, R4_Error, R8_Error, Dnorm ] = foveaLocalizationMessidor(retDir, groundTruthData)

% retDir: Path to the folder with Messidor images (available at https://www.adcis.net/es/descargas-de-software-de-terceros/messidor-es/)
% groundTruthData: Path to the file with the groundtruth x, y coords of the center of the fovea (available at http://www.uhu.es/retinopathy/eng/bd.php)
% For a detailed explanation of the output variables and the program, itself, see the associated paper

warning('OFF', 'all')

groundTruthData = importdata(groundTruthData);
groundTruthData1 = groundTruthData.data.Sheet1;
nImages = 1136;
a = 2;
b = 1;
threshold = 2;
euclideanError = zeros(1, nImages);
R_Error = zeros(1, nImages);
R2_Error = zeros(1, nImages);
R4_Error = zeros(1, nImages);
R8_Error = zeros(1, nImages);
Dnorm = zeros(1, nImages);

for i = 1:nImages
    i
    retinalImage = imread(fullfile(retDir, groundTruthData.textdata.Sheet1{i+1}));
    nrows = size(retinalImage, 1);
    ncols = size(retinalImage, 2);
    retinalImageG = retinalImage(:, :, 2);  % the G color channel is selected
        
    % This piece of code corresponds to the method for optic disc localization and segmentation published in:
    % "Contrast based circular approximation for accurate and robust optic disc segmentation in retinal images"
    if nrows == 960
        scaleFactor = 0.5960;              
    elseif nrows == 1488
        scaleFactor = 0.3930;            
    elseif nrows == 1536
        scaleFactor = 0.3730;            
    else
        error('Unexpected image dimensions')
    end
    retinalImage1 = imresize(retinalImage, scaleFactor);
    [ xCenter, yCenter ] = OpticDiscLocUsingFilters(retinalImage1); % Get estimated optic disc center
    nrows1 = size(retinalImage1, 1);
    ncols1 = size(retinalImage1, 2);
    xCenter = round(xCenter/(ncols1/ncols));
    yCenter = round(yCenter/(nrows1/nrows));
    FOD = FOV_Diameter( retinalImage ); % Diameter of the retina
    odRadius = 0.15*((FOD)/2);  % Optic disc radius
    
    % The previous anatomical information is used to obtain an initial estimate of the center of the fovea
    % The original image is cropped around this approximate location
    minYFoveaBox = yCenter - (a/4) * odRadius;
    maxYFoveaBox = yCenter + a * odRadius;
    if xCenter <= round(ncols/2)
        minXFoveaBox = max(1, xCenter + (2.5 * 2 * odRadius) - (0.75 * b * odRadius));
        maxXFoveaBox = min(ncols, xCenter + (2.5 * 2 * odRadius) + (0.75 * b * odRadius));
    else
        minXFoveaBox = max(1, xCenter - (2.5 * 2 * odRadius) - (0.75 * b * odRadius));
        maxXFoveaBox = min(ncols, xCenter - (2.5 * 2 * odRadius) + (0.75 * b * odRadius));
    end
    foveaImage = retinalImageG(minYFoveaBox:maxYFoveaBox, minXFoveaBox:maxXFoveaBox);
    foveaImage = imgaussfilt(foveaImage, 1);    % The cropped image is smoothed by applying a gaussian filter
        
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
    
    euclideanError(i) = sqrt((xFovea - groundTruthData1(i, 2))^2 + (yFovea - groundTruthData1(i, 3))^2);
    Dnorm(i) = (euclideanError(i)*100)/FOD;
    R = groundTruthData1(i, 1)/2;
    if euclideanError(i) <= R
        R_Error(i) = 1;
    end
    if euclideanError(i) <= R/2
        R2_Error(i) = 1;
    end
    if euclideanError(i) <= R/4
        R4_Error(i) = 1;
    end
    if euclideanError(i) <= R/8
        R8_Error(i) = 1;
    end
end
   
end