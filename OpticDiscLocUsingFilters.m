function [ xCenter, yCenter ] = OpticDiscLocUsingFilters(retinalImage)
% Finds the estimated optic disc center of a retinal image.

% xCenter: x coordinate of the optic disc center
% yCenter: y coordinate of the optic disc center

% retinalImage: input retinal image. It must be resized before calling this
% method so that its FOV diameter is equal to 540 pixels

warning('OFF', 'all')

nrows = size(retinalImage, 1);
ncols = size(retinalImage, 2);
retinalImageR = retinalImage(:, :, 1);
retinalImageG = retinalImage(:, :, 2);
retinalImageB = retinalImage(:, :, 3);

% FOV mask extraction (this piece of code may only be applicable to certain images)
level = graythresh(retinalImageR);
blackMask = retinalImageR < level*255;
ccBlackMask = bwconncomp(blackMask);
blackMaskLabels = labelmatrix(ccBlackMask);
blackMaskProps = regionprops(blackMaskLabels, 'Area');
blackMaskAreas = [ blackMaskProps.Area ];
[ ~, indmaxblackMaskAreas ] = max(blackMaskAreas);
FOVmask = ~(blackMaskLabels == indmaxblackMaskAreas);
se = strel('disk', 1);
FOVmask = imerode(FOVmask, se);
retinalImageR(~FOVmask) = mean(retinalImageR(FOVmask));
retinalImageG(~FOVmask) = mean(retinalImageG(FOVmask));
retinalImageB(~FOVmask) = mean(retinalImageB(FOVmask));

% A vessel-enhanced image is obtained by applying morphological operations to the green channel
retinalImageGc = 255 - retinalImageG;
se = strel('disk', 8);
vesselEnhancedImage = imtophat(retinalImageGc, se);
vesselEnhancedImage(~FOVmask) = mean(vesselEnhancedImage(FOVmask));
% imtool(vesselEnhancedImage, [])
vesselMask = vesselEnhancedImage > prctile(vesselEnhancedImage(:), 98.5);
%imtool(vesselMask)

% The density of vessels mask is obtained by thresholding the subtraction of the output of two average filters computed on the vessel-enhanced image
opticDiscRadius = 40;
h1 = fspecial('average', [ opticDiscRadius*2 opticDiscRadius ]);
h2 = fspecial('average', [ opticDiscRadius*2 3*opticDiscRadius ]);
vesselEnhancedImage = double(vesselEnhancedImage);
vesselEnhancedImageh1 = imfilter(vesselEnhancedImage, h1, 'replicate');
vesselEnhancedImageh2 = imfilter(vesselEnhancedImage, h2, 'replicate');
vesselDensity = double(vesselEnhancedImageh1) - double(vesselEnhancedImageh2);
vesselDensity(~FOVmask) = 0;
% imtool(vesselDensity, [])
vesselDensityMask = vesselDensity > 0.3*max(vesselDensity(:));
% imtool(vesselDensityMask, [])

% The convergence of vessels mask is computed by thresholding the number of intersections of lines detected using the Hough transform
vesselEdges = edge(vesselEnhancedImage, 'canny', 0.25);
[H, theta, rho] = hough(vesselEdges);
P = houghpeaks(H, 500, 'Threshold', 0);
lines = houghlines(vesselEdges, theta, rho, P, 'FillGap', 3, 'MinLength', 15);
% figure, imshow(vesselEdges)
% hold on
% for i = 1:length(lines)
%    xy = [lines(i).point1; lines(i).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end
% The detected lines are sorted according to their y coordinates
vec1 = zeros(length(lines), 2);
vec2 = zeros(length(lines), 2);
for i = 1:length(lines)
    xy1 = lines(i).point1;
    xy2 = lines(i).point2;
    vec1(i, 1) = xy1(1);
    vec1(i, 2) = xy2(1);
    vec2(i, 1) = xy1(2);
    vec2(i, 2) = xy2(2);
end
[ ~, indvec2sort ] = sort(vec2, 1, 'ascend');
indvec2sort = indvec2sort(:, 1);
vec21 = vec2(:, 1);
vec22 = vec2(:, 2);
vec2sort = [ vec21(indvec2sort) vec22(indvec2sort) ];
vec11 = vec1(:, 1);
vec12 = vec1(:, 2);
vec1sort = [ vec11(indvec2sort) vec12(indvec2sort) ];
vesselEdgeIntersections = zeros(nrows, ncols);
xi = 1:nrows;
for i = 1:size(vec2sort, 1)
    mi = (vec1sort(i, 2) - vec1sort(i, 1))/(vec2sort(i, 2) - vec2sort(i, 1));
    bi = vec1sort(i, 1) - mi*vec2sort(i, 1);
    yi = mi*xi + bi;
    for j = 1:size(vec2sort, 1)        
        mj = (vec1sort(j, 2) - vec1sort(j, 1))/(vec2sort(j, 2) - vec2sort(j, 1));
        bj = vec1sort(j, 1) - mj*vec2sort(j, 1);
        yj = mj*xi + bj;
        [ miny, indminy ] = min(abs(yi - yj));
        yindminy = mi*indminy + bi;
        if miny < 0.5 && sum(abs(yi - yj) == miny) == 1 && yindminy >=1 && yindminy <=ncols
            vesselEdgeIntersections(indminy, round(yindminy)) = vesselEdgeIntersections(indminy, round(yindminy)) + 1;
        end
    end
end
h3 = fspecial('disk', opticDiscRadius);
vesselEdgeIntersectionDensity = imfilter(double(vesselEdgeIntersections), h3, 'replicate');
vesselEdgeIntersectionDensity(~FOVmask) = 0;
% imtool(vesselEdgeIntersectionDensity, [])
vesselEdgeIntersectionDensityMask = vesselEdgeIntersectionDensity > 0.3 * max(vesselEdgeIntersectionDensity(:));
% imtool(vesselEdgeIntersectionDensityMask)

% A final constraint mask based on the retinal vessel information is obtained as a logical And between the two previous masks
vesselConstraintMask = vesselDensityMask & vesselEdgeIntersectionDensityMask;
%imtool(vesselConstraintMask)

% An intensity image is obtained by adding the three color channels
retinalIntensityImage = double(retinalImageR) + double(retinalImageG) + double(retinalImageB);
retinalIntensityImage = medfilt2(retinalIntensityImage);
%imtool(retinalIntensityImage, [])

% The intensity and vessel-enhanced images are combined and filtered to build the optic disc center detector
intensityVesselEnhancedImage = double(retinalIntensityImage);
intensityVesselEnhancedImage(vesselMask) = max(retinalIntensityImage(:));
intensityVesselEnhancedImage(~FOVmask) = mean(intensityVesselEnhancedImage(FOVmask));
%imtool(intensityVesselEnhancedImage, [])
h4 = fspecial('average', [ opticDiscRadius*2 opticDiscRadius*4 ]);
opticDiscDetector = imfilter(double(intensityVesselEnhancedImage), h3, 'replicate') - imfilter(double(intensityVesselEnhancedImage), h4, 'replicate');
%imtool(opticDiscDetector, [])
opticDiscDetector(~vesselConstraintMask) = 0;
%imtool(opticDiscDetector, [])

% The x and y coordinates where the detector output reaches the maximum value are taken as the best approximation to the optic disc center  
[ ~, indmaxopticDiscDetector ] = max(opticDiscDetector(:));
[ yCenter, xCenter ] = ind2sub(size(opticDiscDetector), indmaxopticDiscDetector);

end