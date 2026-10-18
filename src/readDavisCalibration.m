function calibCell = readDavisCalibration(xmlAddr)
% Produces a cell containing calibration data read from a Davis .xml file.
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
%   "pinhole" or "polynomial". For polynomial calibrations, the cell value
%   is itself a cell (M, 1), where M is the number of planes used in the
%   polynomial calibration. For pinhole calibrations, it is a struct which 
%   contains the intrinsic and extrinsic parameters. 
% 
%
% Other m-files required: none
% Subfunctions: readDavisPinholeCalibration, readDavisPolyCalibration,
%               hasNode
% MAT-files required: none
%
%
% Author: Ali Semati
% May 2025
%
%------------- BEGIN CODE --------------

calibStruct = readstruct(xmlAddr);

cMap = calibStruct.CoordinateSystemsForEachView.CoordinateSystem.CoordinateMapper;

calibCell = cell(numel(cMap), 2);

for i = 1:numel(cMap)
    if hasNode(cMap(i), 'PinholeParameters')
        calibCell{i,1} = readDavisPinholeCalibration(cMap(i));
        calibCell{i,2} = "pinhole"; 
    elseif hasNode(cMap(i), 'PolynomialParameters')
        calibCell{i,1} = readDavisPolyCalibration(cMap(i));
        calibCell{i,2} = "polynomial"; 
    end
end

end

function calibStruct = readDavisPinholeCalibration(cMap)
% Extracts pinhole parameters from Davis pinhole calibration struct. 
%   Davis uses the OpenCV format for its pinhole model. This format uses
%   zero based indexing, making the first pixel (0,0). The intrinsics
%   object produced here contains the same information as the Davis
%   calibration, but with the principal point shifted by one to account for
%   the fact that in Matlab, the first pixel is (1,1).

    if hasNode(cMap, 'CameraIdentifierAttribute')
        camNum = double(cMap.CameraIdentifierAttribute);
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
    
    radialDist(1) = double(internalParams.RadialDistortion.radialDistortionCoefficient1Attribute);
    radialDist(2) = double(internalParams.RadialDistortion.radialDistortionCoefficient2Attribute);
    radialDist(3) = 0;
    if hasNode(internalParams, 'TangentialDistortion')
        tangentialDist(1) = double(internalParams.TangentialDistortion.tangentialDistortionCoefficient1Attribute);
        tangentialDist(2) = double(internalParams.TangentialDistortion.tangentialDistortionCoefficient2Attribute);
    else
        
        tangentialDist(1:2) = 0;
    end
    
    % pixel at optical center
    principlePointX = double(internalParams.PrincipalPoint.xAttribute);
    principlePointY = double(internalParams.PrincipalPoint.yAttribute);
    
    % length of pixel in camera sensor
    pixelSize = double(internalParams.SensorPixelSizeMm.ValueAttribute);
    
    focalLengthX = double(internalParams.FocalLengthPixel.xAttribute)/pixelSize;
    focalLengthY = double(internalParams.FocalLengthPixel.yAttribute)/pixelSize;
    
    imageSize(1) = double(commonParams.OriginalImageSize.HeightAttribute);
    imageSize(2) = double(commonParams.OriginalImageSize.WidthAttribute);
    
    distortionCoefficients = [radialDist(1) radialDist(2) tangentialDist(1) tangentialDist(2) radialDist(3)];
    
    intrinsics = cameraIntrinsics([focalLengthX, focalLengthY], ...
                            [principlePointX + 1, principlePointY + 1], ...
                            imageSize, ...
                            "RadialDistortion",distortionCoefficients(1:2), ...
                            "TangentialDistortion", distortionCoefficients(3:4));
    
    translation = pinholeParams.ExternalCameraParameters.TranslationMm;
    translation = [double(translation.TxAttribute); double(translation.TyAttribute); double(translation.TzAttribute)];
    
    % Euler angles (x,y,z)
    rotation = pinholeParams.ExternalCameraParameters.RotationAngles;
    rotation = [double(rotation.RxAttribute); double(rotation.RyAttribute); double(rotation.RzAttribute)];
    
    scaleStructDavis = struct();
    Scales = pinholeParams.Scales;
    
    scaleStructDavis.X.Slope = ...
        double(Scales.LinearScaleX.FactorMmPerPixelAttribute);
    scaleStructDavis.X.SlopeOrg = scaleStructDavis.X.Slope;
    scaleStructDavis.X.Offset = ...
        double(Scales.LinearScaleX.OffsetMmAttribute);
    scaleStructDavis.X.Unit = char(Scales.LinearScaleX.UnitAttribute);
    
    scaleStructDavis.Y.Slope = ...
        double(Scales.LinearScaleY.FactorMmPerPixelAttribute);
    scaleStructDavis.Y.SlopeOrg = scaleStructDavis.Y.Slope;
    scaleStructDavis.Y.Offset = ...
        double(Scales.LinearScaleY.OffsetMmAttribute);
    scaleStructDavis.Y.Unit = char(Scales.LinearScaleY.UnitAttribute);
    
    
    
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
% Extracts polynomial coefficients from Davis polynomial calibration struct. 
    if hasNode(cMap, 'CameraIdentifierAttribute')
        camNum = double(cMap.CameraIdentifierAttribute);
    else
        camNum = [];
    end
    calibType = char(cMap.TypeAttribute);
    
    params = cMap.PolynomialParameters;
    mappings = params.PolynomialMapping;    % struct array, one per plane
    
    calibCell = cell(1, numel(mappings));

    % Coefficient name suffixes in the order the polynomial below uses
    % them: 1, s, s^2, s^3, t, t^2, t^3, s*t, s^2*t, s*t^2. Davis names
    % the attributes a_<suffix> and b_<suffix>.
    coeffTerms = ["o" "s" "s2" "s3" "t" "t2" "t3" "st" "s2t" "st2"];

    for planeNum = 1:numel(mappings)
        mapping = mappings(planeNum);
        coeffNodeA = mapping.Polynomial3rdOrder.CoefficientsA;
        coeffNodeB = mapping.Polynomial3rdOrder.CoefficientsB;
        S0 = double(mapping.Origin.s_oAttribute);
        T0 = double(mapping.Origin.t_oAttribute);
        dwSize = params.CommonParameters.CorrectedImageSize;
        dwWidth = double(dwSize.WidthAttribute);
        dwHeight = double(dwSize.HeightAttribute);
        dwSize = [dwWidth, dwHeight];
        orgSize = params.CommonParameters.OriginalImageSize;
        orgWidth = double(orgSize.WidthAttribute);
        orgHeight = double(orgSize.HeightAttribute);
        orgSize = [orgWidth, orgHeight];
        pixelsPerMM = double(params.CommonParameters.PixelPerMmFactor.ValueAttribute);
        mmPerPixel_x = double(params.Scales.LinearScaleX.FactorMmPerPixelAttribute);
        mmPerPixel_y = double(params.Scales.LinearScaleY.FactorMmPerPixelAttribute);
        
        offsetX = double(params.Scales.LinearScaleX.OffsetMmAttribute);
        offsetY = double(params.Scales.LinearScaleY.OffsetMmAttribute);
        
        A = zeros(10,1);
        B = zeros(10,1);

        for i = 1:10
            A(i) = double(coeffNodeA.("a_" + coeffTerms(i) + "Attribute"));
            B(i) = double(coeffNodeB.("b_" + coeffTerms(i) + "Attribute"));
        end
                
        calibStruct = struct();
        calibStruct.A = A;
        calibStruct.B = B;
        calibStruct.origin = [S0 T0];
        calibStruct.orgSize = orgSize;
        calibStruct.dwSize = dwSize;
        calibStruct.offset = [offsetX offsetY];
        calibStruct.pixelsPerMM = pixelsPerMM;
        
        [x, y] = meshgrid(0:dwWidth-1, 0:dwHeight-1);
        
        GSx = x - (A(1)+A(2)*(2*(x - S0)/orgWidth)+...
                   A(3)*(2*(x - S0)/orgWidth).^2+...
                   A(4)*(2*(x - S0)/orgWidth).^3+...
                   A(5)*(2*(y - T0)/orgHeight)+...
                   A(6)*(2*(y - T0)/orgHeight).^2+...
                   A(7)*(2*(y - T0)/orgHeight).^3+...
                   A(8)*(2*(x - S0)/orgWidth).*(2*(y - T0)/orgHeight)+...
                   A(9)*(2*(x - S0)/orgWidth).^2.*(2*(y - T0)/orgHeight)+...
                   A(10)*(2*(y - T0)/orgHeight).^2.*(2*(x - S0)/orgWidth));
        
        GSy = y - (B(1)+B(2)*(2*(x - S0)/orgWidth)+...
                   B(3)*(2*(x - S0)/orgWidth).^2+...
                   B(4)*(2*(x - S0)/orgWidth).^3+...
                   B(5)*(2*(y - T0)/orgHeight)+...
                   B(6)*(2*(y - T0)/orgHeight).^2+...
                   B(7)*(2*(y - T0)/orgHeight).^3+...
                   B(8)*(2*(x - S0)/orgWidth).*(2*(y - T0)/orgHeight)+...
                   B(9)*(2*(x - S0)/orgWidth).^2.*(2*(y - T0)/orgHeight)+...
                   B(10)*(2*(y - T0)/orgHeight).^2.*(2*(x - S0)/orgWidth));
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
        zLevelPixel = double(mapping.ZPosition.ValueAttribute);
        calibStruct.zLevel = zLevelPixel*abs(mmPerPixel_x);
        calibStruct.camNum = camNum;
        calibStruct.calibType = calibType;
    
        calibCell{planeNum} = calibStruct;
    end

end

function tf = hasNode(s, name)
% True if scalar struct s has field name holding real content. readstruct
%   gives every element of a struct array the union of its siblings'
%   fields and fills the ones an element lacks with missing.
%   Inputs:     s       - scalar struct
%               name    - field name
%   Output:     tf      - logical
    tf = isfield(s, name) && ~isa(s.(name), 'missing');
end