function [imageRct, tform, scaling, xMesh, yMesh] = rectifyImagePinhole(img, calibStruct, dz, theta)
% rectifyImagePinhole - Applies a perspective transformation that rectifies
%       an image to a fronto-parallel view based on pinhole camera 
%       calibration parameters.
%
% Usage: outputImage = rectifyImagePinhole(img, calibStruct);
%         - Transforms the image "img" using the pinhole calibration model
%           in calibStruct to appear as if taken by a camera placed above
%           the world coordinate system origin, with its optical axis
%           perpendicular to the XY-plane and its image axes aligned with 
%           the world X and Y axes. The object plane is Z = 0 by default.
%        outputImage = rectifyImagePinhole(img, calibStruct, 1.5, 10);
%         - Same as above, but the object plane is at Z = 1.5 mm and the
%           camera is rotated by 10 degrees about the Z-axis, with positive
%           rotation defined by the right-hand rule.
%
% Inputs:
%   img         - Numeric array representing the image.
%   calibStruct - Calibration struct containing intrinsics and extrinsics.
%                 This struct must contain the fields
%                    - intrinsicsMatlab:    Matlab cameraIntrinsics object.
%                    - extrinsics:          Matlab rigidtform3d object.
%   dz          - Optional scalar specifying the offset of the object plane
%                 along the world Z-axis (in mm). When dz = 0 (default), 
%                 the object plane coincides with the world XY-plane 
%                 (Z = 0). Positive values move the plane towards +Z.
%   theta       - Optional rotation angle (in degrees) about the Z-axis of 
%                 the fronto-parallel camera.
%
% Outputs:
%   imageRct - Fronto-parallel view of the input image.
%   tform    - Object of type projtform2d containing a homography matrix
%              that transforms coordinates from the camera view to the 
%              fronto-parallel view. For efficiency, use this object with
%              the imwarp function when rectifying multiple images. 
%   scaling  - Structure array containing the X and Y coordinates 
%              corresponding to pixel (1,1) in the dewarped image under the
%              "Offset" field and the scaling (in mm/pixel) for each axis. 
%   xMesh    - Array of X coordinates.
%   yMesh    - Array of Y coordinates. 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%
% Author: Ali Semati
% May 2025; Last revision: 26-June-2025

%------------- BEGIN CODE --------------

arguments
    img
    calibStruct
    dz (1,1) double = 0
    theta (1,1) double = 0
end

n = [0;0;1];    % plane normal vector in world coordinates

imageSize = calibStruct.intrinsicsMatlab.ImageSize;     % (Ny, Nx)
K = calibStruct.intrinsicsMatlab.K;   % intrinsic matrix
R1 = calibStruct.extrinsics.R;  % rotation matrix
% vector from camera origin to world origin in camera coordinates
t1 = calibStruct.extrinsics.Translation';
% vector from world origin to camera origin in world coordinates
t1_w = -R1'*t1;

R2 = [1 0 0; ...
      0 -1 0; ...
      0 0 -1];      % fronto-parallel rotation matrix
if theta~=0
    R2 = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1] * R2;
end
t2 = -R2*[0; 0; abs(t1_w(3))];   % fronto-parallel translation vector

n_c = R1*n;         % plane normal vector in camera coordinates
d = abs(n_c'*t1 + dz);   % distance from camera center to homography plane

% compute homography matrix
% adapted from Hartley & Zisserman (2003), Eq. 13.2
H = K*(R2*R1' - (R1*R2'*t2 - t1)*n_c'/d)*inv(K);
H = H./H(3,3);  % normalize

% determine optimum image resolution
% dewarped image should be large enough so that downscaling does not occur,
% but no larger

% corner points, starting from top left and moving counter-clockwise
% first row is u, second is v
cornerPoints = [1 1 1; 
                1, imageSize(1), 1;
                imageSize(2), imageSize(1), 1;
                imageSize(2), 1, 1
               ]';
tform_corner= H*cornerPoints;
% coordinates of corners in dewarped image
tform_corner = tform_corner./tform_corner(3,:);


edgeLengths = abs(tform_corner(1:2,:) - tform_corner(1:2,[2 3 4 1]));
edgeLengths = sort(edgeLengths, 2);

% find smallest edges
xEdgeMin = edgeLengths(1,3);
yEdgeMin = edgeLengths(2,3);

% scale so that the pixel density of the smallest edge is equal to the
% original image
cDiff = diff(tform_corner, 1, 2);
if abs(cDiff(2,1)) > abs(cDiff(1,1))    % world x axis aligned with u
    scale = max(imageSize(2)/xEdgeMin, imageSize(1)/yEdgeMin);    
else                                    % world x axis aligned with v
    scale = max(imageSize(1)/xEdgeMin, imageSize(2)/yEdgeMin);
end

H(1:2,:) = H(1:2,:)*scale;
tform = projtform2d(H);

% homography cannot correct for non-linear distortions
% undistort first
imgUndist = undistortImage(img, calibStruct.intrinsicsMatlab);

% apply perspective transformation
[imageRct, RB]= imwarp(imgUndist, tform);

% calculate image scaling factor
% find world vectors (1, 0, 0) mm and (0, 1, 0) mm in the original image (cam1)
sPointsCam1 = K*[R1 t1]*[0 0 dz 1; 1 0 dz 1; 0 1 dz 1]';
sPointsCam1 = sPointsCam1 ./ sPointsCam1(3,:);
% apply transformation to obtain vectors in the fronto-parallel image (cam2)
sPointsCam2 = transformPointsForward(tform, sPointsCam1(1:2,:)');

% find coordinates of pixel (1,1) in the coordinate system of the original image
% RB.worldLimits are pixel edge to pixel edge so +0.5 to get center coordinates
firstPixelX = RB.XWorldLimits(1) + 0.5;
firstPixelY = RB.YWorldLimits(1) + 0.5;

sPointsCam2(:,1) = sPointsCam2(:,1) - firstPixelX + 1;
sPointsCam2(:,2) = sPointsCam2(:,2) - firstPixelY + 1;

M_image = [sPointsCam2, ones(3,1)];
P_world = [ 0, 0;    
            1, 0;   
            0, 1 ];

T = (M_image \ [P_world(:,1), P_world(:,2)])';  

Width = size(imageRct, 2);
Height = size(imageRct, 1);
[U, V] = meshgrid(1:Width, 1:Height); 
xMesh = T(1,1)*U + T(1,2)*V + T(1,3); 
yMesh = T(2,1)*U + T(2,2)*V + T(2,3);

scaling = struct();
scaling.X.Unit = 'mm';
scaling.Y.Unit = 'mm';

scaling.imageRB = RB;
scaling.thetaDeg = theta;

vecX = sPointsCam2(2,:) - sPointsCam2(1,:);   % pixel displacement for +1 mm world X
scaling.mmPerPixel = 1/(norm(vecX) + eps);

if theta == 0
    slopeX = 1/( sPointsCam2(2,1) - sPointsCam2(1,1) + eps);    % mm per pixel
    slopeY = 1/( sPointsCam2(3,2) - sPointsCam2(1,2) + eps);    % mm per pixel
    scaling.X.Slope = slopeX;   
    scaling.Y.Slope = slopeY;

    % physical x coordinate at center of pixel (1,1)
    scaling.X.Offset = xMesh(1,1); 
    % y coordinate at center of pixel (1,1)
    scaling.Y.Offset = yMesh(1,1); 
end

end