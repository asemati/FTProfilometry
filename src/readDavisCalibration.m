function calibCell = readDavisCalibration(xmlAddr)
% READDAVISCALIBRATION - Produces a cell containing calibration data read
%                        from a Davis .xml file.
%
% Usage: calib = readDavisCalibration('D:/path/to/Calibration.xml');
%
% Inputs:
%   xmlAddr - address of calibration .xml file.
%
% Outputs:
%   calibCell - (N,2) cell structure, where N is the number of camera
%   calibrations contained in the file. The first column contains the data
%   and the second contains a string denoting the type of the calibration:
%   "Pinhole" or "Polynomial". For polynomial calibrations, the cell value
%   is itself a cell (M, 1), where M is the number of planes used in the
%   polynomial calibration. For pinhole calibrations, it is a struct which 
%   contains the intrinsic and extrinsic parameters. 
% 
%
% Other m-files required: xml2struct
% Subfunctions: readDavisPinholeCalibration, readDavisPolyCalibration
% MAT-files required: none
%
%
% Author: Ali Semati
% May 2025; Last revision: 21-August-2026
%
%------------- BEGIN CODE --------------

calibStruct = xml2struct(xmlAddr);

if iscell(calibStruct.Calibration.CoordinateSystemsForEachView.CoordinateSystem.CoordinateMapper)
    cMap = calibStruct.Calibration.CoordinateSystemsForEachView.CoordinateSystem.CoordinateMapper;
else
    cMap{1} = calibStruct.Calibration.CoordinateSystemsForEachView.CoordinateSystem.CoordinateMapper;
end

calibCell = cell(length(cMap), 2);

for i = 1:length(cMap)
    if isfield(cMap{i}, 'PinholeParameters')
        calibCell{i,1} = readDavisPinholeCalibration(cMap{i});
        calibCell{i,2} = "Pinhole"; 
    elseif isfield(cMap{i}, 'PolynomialParameters')
        calibCell{i,1} = readDavisPolyCalibration(cMap{i});
        calibCell{i,2} = "Polynomial"; 
    end
end

end

function calibStruct = readDavisPinholeCalibration(cMap)
%   Extracts pinhole parameters from Davis pinhole calibration struct. 
%   The field "intrinsics" in calibStruct contains intrinsics in OpenCV 
%   format (used by Davis). This format has zero based indexing, making the
%   first pixel (0,0). The field "intrinsicsMatlab" contains the same
%   information, but with the principal point shifted by one to account for
%   the fact that in Matlab, the first pixel is (1,1).

    if isfield(cMap.Attributes, 'CameraIdentifier')
        camNum = str2double(cMap.Attributes.CameraIdentifier);
    else
        camNum = [];
    end

    % radial distortion coefficients
    radialDist = zeros(1,3);
    % tangential distortion coefficients
    tangentialDist = zeros(1,3);
    % size (in pixels) of the camera frame (rows,columns)
    imageSize = zeros(1,2);
    
    % cannot handle multiview calibrations
    pinholeParams = cMap.PinholeParameters;
    internalParams = pinholeParams.InternalCameraParameters;
    commonParams = pinholeParams.CommonParameters;
    
    % Davis stores data as strings, convert to double
    radialDist(1) = str2double(internalParams.RadialDistortion.Attributes.radialDistortionCoefficient1);
    radialDist(2) = str2double(internalParams.RadialDistortion.Attributes.radialDistortionCoefficient2);
    radialDist(3) = 0;
    if isfield(internalParams, 'TangentialDistortion')
        tangentialDist(1) = str2double(internalParams.TangentialDistortion.Attributes.tangentialDistortionCoefficient1);
        tangentialDist(2) = str2double(internalParams.TangentialDistortion.Attributes.tangentialDistortionCoefficient2);
    else
        
        tangentialDist(1:2) = 0;
    end
    
    % pixel at optical center
    principlePointX = str2double(internalParams.PrincipalPoint.Attributes.x);
    principlePointY = str2double(internalParams.PrincipalPoint.Attributes.y);
    
    % length of pixel in camera sensor
    pixelSize = str2double(internalParams.SensorPixelSizeMm.Attributes.Value);
    
    focalLengthX = str2double(internalParams.FocalLengthPixel.Attributes.x)/pixelSize;
    focalLengthY = str2double(internalParams.FocalLengthPixel.Attributes.y)/pixelSize;
    
    imageSize(1) = str2double(commonParams.OriginalImageSize.Attributes.Height);
    imageSize(2) = str2double(commonParams.OriginalImageSize.Attributes.Width);
    
    distortionCoefficients = [radialDist(1) radialDist(2) tangentialDist(1) tangentialDist(2) radialDist(3)];
    
    intrinsics = cameraIntrinsics([focalLengthX, focalLengthY], ...
                            [principlePointX principlePointY], ...
                            imageSize, ...
                            "RadialDistortion",distortionCoefficients(1:2), ...
                            "TangentialDistortion", distortionCoefficients(3:4));
    
    intrinsicsMatlab = cameraIntrinsics([focalLengthX, focalLengthY], ...
                            [principlePointX + 1, principlePointY + 1], ...
                            imageSize, ...
                            "RadialDistortion",distortionCoefficients(1:2), ...
                            "TangentialDistortion", distortionCoefficients(3:4));
    
    translation = pinholeParams.ExternalCameraParameters.TranslationMm.Attributes;
    translation = [str2double(translation.Tx); str2double(translation.Ty); str2double(translation.Tz)];
    
    % Euler angles (x,y,z)
    rotation = pinholeParams.ExternalCameraParameters.RotationAngles.Attributes;
    rotation = [str2double(rotation.Rx); str2double(rotation.Ry); str2double(rotation.Rz)];
    
    scaleStructDavis = struct();
    Scales = pinholeParams.Scales;
    
    scaleStructDavis.X.Slope = ...
        str2double(pinholeParams.Scales.LinearScaleX.Attributes.FactorMmPerPixel);
    scaleStructDavis.X.SlopeOrg = scaleStructDavis.X.Slope;
    scaleStructDavis.X.Offset = ...
        str2double(pinholeParams.Scales.LinearScaleX.Attributes.OffsetMm);
    scaleStructDavis.X.Unit = Scales.LinearScaleX.Attributes.Unit;
    
    scaleStructDavis.Y.Slope = ...
        str2double(Scales.LinearScaleY.Attributes.FactorMmPerPixel);
    scaleStructDavis.Y.SlopeOrg = scaleStructDavis.Y.Slope;
    scaleStructDavis.Y.Offset = ...
        str2double(Scales.LinearScaleY.Attributes.OffsetMm);
    scaleStructDavis.Y.Unit = Scales.LinearScaleY.Attributes.Unit;
    
    
    
    extrinsics = rigidtform3d(rad2deg(rotation'), translation);
    
    calibStruct = struct();
    calibStruct.intrinsics = intrinsics;
    calibStruct.intrinsicsMatlab = intrinsicsMatlab;
    calibStruct.extrinsics = extrinsics;
    calibStruct.EulerAngles = rotation';    % X,Y,Z in radians
    calibStruct.scaling = scaleStructDavis;
    calibStruct.camNum = camNum;
end

function calibCell = readDavisPolyCalibration(cMap)
%   Extracts polynomial coefficients from Davis polynomial calibration struct. 
    if isfield(cMap.Attributes, 'CameraIdentifier')
        camNum = str2double(cMap.Attributes.CameraIdentifier);
    else
        camNum = [];
    end
    calibType = cMap.Attributes.Type;
    
    params = cMap.PolynomialParameters;
    if iscell(params.PolynomialMapping)
        mappingCell = params.PolynomialMapping;
    else
        mappingCell{1} = params.PolynomialMapping;
    end
    
    calibCell = cell(1, length(mappingCell));

    for planeNum = 1:length(mappingCell)
        mapping = mappingCell{planeNum};
        A_txt = mapping.Polynomial3rdOrder.CoefficientsA.Attributes;
        B_txt = mapping.Polynomial3rdOrder.CoefficientsB.Attributes;
        S0 = str2double(mapping.Origin.Attributes.s_o);
        T0 = str2double(mapping.Origin.Attributes.t_o);
        dwSize = params.CommonParameters.CorrectedImageSize.Attributes;
        dwSize = [str2double(dwSize.Width), str2double(dwSize.Height)];
        orgSize = params.CommonParameters.OriginalImageSize.Attributes;
        orgSize = [str2double(orgSize.Width), str2double(orgSize.Height)];
        pixelsPerMM = params.CommonParameters.PixelPerMmFactor.Attributes.Value;
        pixelsPerMM = str2double(pixelsPerMM);
        mmPerPixel_x = params.Scales.LinearScaleX.Attributes.FactorMmPerPixel;
        mmPerPixel_x = str2double(mmPerPixel_x);
        mmPerPixel_y = params.Scales.LinearScaleY.Attributes.FactorMmPerPixel;
        mmPerPixel_y = str2double(mmPerPixel_y);
        
        offsetX = params.Scales.LinearScaleX.Attributes.OffsetMm;
        offsetX = str2double(offsetX);
        offsetY = params.Scales.LinearScaleY.Attributes.OffsetMm;
        offsetY = str2double(offsetY);
        
        fieldsA = fieldnames(A_txt);
        fieldsB = fieldnames(B_txt);
        
        correctFieldOrd = [1 2 3 5 8 9 10 6 4 7];
        fieldsA = fieldsA(correctFieldOrd);
        fieldsB = fieldsB(correctFieldOrd);
        
        A = zeros(10,1);
        B = zeros(10,1);
        
        for i = 1:10
            A(i) = str2double(A_txt.(fieldsA{i}));
            B(i) = str2double(B_txt.(fieldsB{i}));
        end
                
        calibStruct = struct();
        calibStruct.A = A;
        calibStruct.B = B;
        calibStruct.origin = [S0 T0];
        calibStruct.orgSize = orgSize;
        calibStruct.dwSize = dwSize;
        calibStruct.offset = [offsetX offsetY];
        calibStruct.pixelsPerMM = pixelsPerMM;
        
        [x, y] = meshgrid(0:dwSize(1)-1, 0:dwSize(2)-1);
        
        GSx = x - (A(1)+A(2)*(2*(x - S0)/orgSize(1))+...
                   A(3)*(2*(x - S0)/orgSize(1)).^2+...
                   A(4)*(2*(x - S0)/orgSize(1)).^3+...
                   A(5)*(2*(y - T0)/orgSize(2))+...
                   A(6)*(2*(y - T0)/orgSize(2)).^2+...
                   A(7)*(2*(y - T0)/orgSize(2)).^3+...
                   A(8)*(2*(x - S0)/orgSize(1)).*(2*(y - T0)/orgSize(2))+...
                   A(9)*(2*(x - S0)/orgSize(1)).^2.*(2*(y - T0)/orgSize(2))+...
                   A(10)*(2*(y - T0)/orgSize(2)).^2.*(2*(x - S0)/orgSize(1)));
        
        GSy = y - (B(1)+B(2)*(2*(x - S0)/orgSize(1))+...
                   B(3)*(2*(x - S0)/orgSize(1)).^2+...
                   B(4)*(2*(x - S0)/orgSize(1)).^3+...
                   B(5)*(2*(y - T0)/orgSize(2))+...
                   B(6)*(2*(y - T0)/orgSize(2)).^2+...
                   B(7)*(2*(y - T0)/orgSize(2)).^3+...
                   B(8)*(2*(x - S0)/orgSize(1)).*(2*(y - T0)/orgSize(2))+...
                   B(9)*(2*(x - S0)/orgSize(1)).^2.*(2*(y - T0)/orgSize(2))+...
                   B(10)*(2*(y - T0)/orgSize(2)).^2.*(2*(x - S0)/orgSize(1)));
        calibStruct.GSx = GSx;
        calibStruct.GSy = GSy;
        
        scaling = struct();
        scaling.X.Slope = mmPerPixel_x;
        scaling.X.SlopeOrg = mmPerPixel_x;
        scaling.X.Unit = 'mm';
        scaling.Y.Slope = mmPerPixel_y;
        scaling.Y.SlopeOrg = mmPerPixel_y;
        scaling.Y.Unit = 'mm';
        scaling.X.Offset = offsetX;
        scaling.Y.Offset = offsetY;
        
        calibStruct.scaling = scaling;
        zLevelPixel = str2double(mapping.ZPosition.Attributes.Value);
        calibStruct.zLevel = zLevelPixel*abs(mmPerPixel_x);
        calibStruct.camNum = camNum;
        calibStruct.calibType = calibType;
    
        calibCell{planeNum} = calibStruct;
    end

end

