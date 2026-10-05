% data used here is available on Dataverse.no at https://doi.org/10.18710/MWEHEM
% download ProfilometryCalibration.zip and unzip it in the working directory
workDir = "path/to/workingDirectory";
imgSaveFolder = "calibImageSet";
%% read NetCDF files and convert
filePath = fullfile(workDir, "imageSet_ProfCalibration.nc");
imageSet = ncread(filePath, 'imageSet');
heightVec = ncread(filePath, 'height'); 
nPlanes = numel(heightVec);

if ~isfolder(fullfile(workDir, imgSaveFolder))
    mkdir(fullfile(workDir, imgSaveFolder))
end

for i = 1:size(imageSet, 3)
    imName = sprintf("image_%04d.tiff", i);
    imwrite(imageSet(:,:,i), fullfile(workDir, imgSaveFolder, imName));
end

%% set up FtpSolver object and solve
% The reference (zero-height plane) is also frame 6 of the stack, so
% FtpSolver warns that it is in the data sequence. The calibration module 
% discards the zero-height plane itself.

calibCase = FtpSolver("calibCase", ...
    refPath       = fullfile(workDir, imgSaveFolder, "image_0006.tiff"), ...    % flat reference fringe image
    dataPath      = fullfile(workDir, imgSaveFolder, "image_0001.tiff"), ...    % fringe image sequence
    camCalibPath  = fullfile(workDir, "CameraPinholeCalibration.xml"), ...      % camera calibration path
    cropRect      = [584 685 1420 1150] ...       % computational domain (X, Y, Width, Height)
    );    

calibCase.calibModeOn();
calibCase.setOutputRectFromMargin(0);
calibCase.showROI();

% use the 15th fringe peak for phase correction
% set tracked peak after setOutputRectFromMargin(), which resets it
calibCase.setPhaseCorrection(method="spatial", trackedPeak=15)
calibCase.setDemodulation(method="wavelet", waveletWavelengths=[28 36], ...
                                    waveletDivisions=10, waveletAngles=-90);

calibCase.solve(); 

calibCase.animate(1, nPlanes, 2, quantity='phase'); % (first frame, last frame, animation frames per second)

%% calibrate and save
% heights +5 and +10 mm appear to be outliers, exclude them from the calibration
outlierHeights = ismember(heightVec, [5 10]);       
calibCase.calibratePoly(2, heightVec, excludeHeights=outlierHeights) % fit a quadratic polynomial to each pixel

calibCase.writePolyCalibration(fullfile(workDir, "polynomialCalibration.mat"));

%% rerun to calculate surface elevation field
calibCase.calibModeOff();
calibCase.solve();
calibCase.animate(1, nPlanes, 2)    
