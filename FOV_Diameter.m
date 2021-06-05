function D = FOV_Diameter( retinalImage )

retinalImageR = retinalImage(:, :, 1);

% FOV mask extraction (this piece of code may only be applicable to certain images)
level = graythresh(retinalImageR);
blackMask = retinalImageR < level*255;
ccBlackMask = bwconncomp(blackMask);
blackMaskLabels = labelmatrix(ccBlackMask);
blackMaskProps = regionprops(blackMaskLabels, 'Area');
blackMaskAreas = [ blackMaskProps.Area ];
[ ~, indmaxblackMaskAreas ] = max(blackMaskAreas);
FOVmask = ~(blackMaskLabels == indmaxblackMaskAreas);
% imtool(FOVmask)
[ y, x ] = find(FOVmask);
miny = min(y(:));
maxy = max(y(:));
D = maxy - miny;

end
