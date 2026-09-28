%% read NetCDF files and convert
filename = "path\to\imageSet_ProfCalibration.nc";
imageSet = ncread(filename, 'imageSet', [1 1 1], [inf, inf, inf]);
for i = 1:size(imageSet, 3)
    imName = sprintf("image_%04d.png", i);
    imwrite(imageSet(:,:,i), fullfile("path\to\workingDirectory\imageSet", imName));
end
heightVec = ncread(filename, 'height', 1, inf); 

%% set up FtpSolver object and solve
calibCase = FtpSolver("calibCase", ...              % case ID must be a valid MATLAB variable name
    refAddr       = "path\to\workingDirectory\imageSet\image_0006.png", ...       % flat reference fringe image
    dataAddr      = "path\to\workingDirectory\imageSet\image_0001.png", ...       % fringe image sequence
    camCalibAddr  = "path\to\workingDirectory\CameraPinholeCalibration.xml", ...  % camera calibration path
    cropRect      = [584 685 1411 1149] ...       % computational domain (X, Y, Width, Height)
    );    

calibCase.clbModeOn();
calibCase.setDisplayToMargin();
calibCase.pCorrOpts.peakInd = 15;   % use the 15th fringe peak for phase correction
calibCase.solve(); 
calibCase.animate(1, 12, 1, target='phase');

%% calibrate and save
% heights +5 and +10 mm appear to be outliers, exclude them from the calibration
outlierHeights = false(size(heightVec));
outlierHeights(8:9) = true;                 
calibCase.calibratePoly(2, heightVec, excludeHeights=outlierHeights)

calibCase.writePolyCalibration("path\to\workingDirectory\polyCalibrationParams.mat");

%% rerun to calculate surface elevation field
demoCase.clbModeOff();
demoCase.solve();
demoCase.animate(1, 12, 1)