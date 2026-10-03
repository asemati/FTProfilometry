% data used here is available on Dataverse.no at https://doi.org/10.18710/MWEHEM
workDir = "path/to/workingDirectory";
%% read NetCDF files and convert
ncFile = fullfile(workDir, "imageSet_ProfCalibration.nc");
imageSet = ncread(ncFile, 'imageSet');
heightVec = ncread(ncFile, 'height'); 
nPlanes = numel(heightVec);

if ~isfolder(fullfile(workDir, "imageSet"))
    mkdir(fullfile(workDir, "imageSet"))
end

for i = 1:size(imageSet, 3)
    imName = sprintf("image_%04d.tiff", i);
    imwrite(imageSet(:,:,i), fullfile(workDir, "imageSet", imName));
end

%% set up FtpSolver object and solve
% The reference (zero-height plane) is also frame 6 of the stack, so
% FtpSolver warns that it is in the data sequence. calibratePoly 
% drops the zero-height plane itself.

calibCase = FtpSolver("calibCase", ...
    refPath       = fullfile(workDir, "imageSet", "image_0006.tiff"), ...       % flat reference fringe image
    dataPath      = fullfile(workDir, "imageSet", "image_0001.tiff"), ...       % fringe image sequence
    camCalibPath  = fullfile(workDir, "CameraPinholeCalibration.xml"), ...  % camera calibration path
    cropRect      = [584 685 1420 1150] ...       % computational domain (X, Y, Width, Height)
    );    

calibCase.calibModeOn();
calibCase.setOutputRectFromMargin(0);
calibCase.showROI();

% use the 15th fringe peak for phase correction
% set tracked peak after setOutputRectFromMargin, which resets it
calibCase.setPhaseCorrection(method="spatial", trackedPeak=15)
calibCase.setDemodulation(method="wavelet", waveletWavelengths=[28 36], ...
                                    waveletDivisions=10, waveletAngles=-90);

calibCase.solve(); 

calibCase.animate(1, nPlanes, 1, quantity='phase'); % animate the calculated phase

%% calibrate and save
% heights +5 and +10 mm appear to be outliers, exclude them from the calibration
outlierHeights = false(size(heightVec));
outlierHeights(8:9) = true;                 
calibCase.calibratePoly(2, heightVec, excludeHeights=outlierHeights)

calibCase.writePolyCalibration(fullfile(workDir, "polyCalibrationParams.mat"));

%% rerun to calculate surface elevation field
calibCase.calibModeOff();
calibCase.solve();
calibCase.animate(1, nPlanes, 1)
