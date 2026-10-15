%% FtpSolver tutorial: calibrating the phase-to-height conversion
%
% FtpSolver measures how the surface shifts a projected fringe pattern.
% This shift is expressed as a phase difference from a reference image of
% the same fringe pattern on a flat surface. Turning that phase into an
% elevation needs a model of the optical geometry. This script shows how to
% build and check one from a calibration experiment: a flat plane imaged at
% twelve known heights, from -25 mm to +25 mm, one of which (0 mm) is the
% reference.
%
% Two models are compared here. Takeda's model describes an ideal
% camera-projector geometry with two parameters. The polynomial model fits,
% at every pixel separately, a polynomial from phase to height, which
% absorbs lens distortion and a non-uniform fringe period and grants full
% freedom in the placement of the camera and projector.
%
% Data: Dataverse.no, https://doi.org/10.18710/MWEHEM
%   Download ProfilometryCalibration.zip and unzip it into workDir below.


%% 1. Paths and data
% The zip file holds the calibration images and the heights of the planes
% as NetCDF, together with the camera calibration .xml file.

workDir       = "path/to/workingDirectory"; % folder containing the .nc and .xml files
ncFile        = fullfile(workDir, "imageSet_ProfCalibration.nc");
camCalibFile  = fullfile(workDir, "CameraPinholeCalibration.xml");
imageFolder   = fullfile(workDir, "calibImageSet");
polyCalibFile = fullfile(workDir, "polynomialCalibration.mat");     % written in Sec. 10


%% 2. Converting the calibration images
% FtpSolver reads an image sequence from a folder, ordering the frames by
% file name, so the images are written out as zero-padded tiff files.
% heightVec holds the height of the plane in each image, in mm.
%
% The zero-height image doubles as the reference. It stays in the image
% folder, so it is also solved as one of the planes. FtpSolver warns about
% this when the object is created and its phase should come out as zero
% everywhere. The calibration functions leave the zero-height plane out of
% the fit.

heightVec = ncread(ncFile, "height");
nPlanes = numel(heightVec);
refIndex = find(heightVec == 0);
fprintf("%d planes at heights (mm): %s\n", nPlanes, mat2str(heightVec'))

if ~isfolder(imageFolder)
    mkdir(imageFolder)
end
if numel(dir(fullfile(imageFolder, "*.tiff"))) ~= nPlanes
    imageSet = ncread(ncFile, "imageSet");      % 12-bit images stored as uint16
    for k = 1:nPlanes
        imwrite(imageSet(:,:,k), fullfile(imageFolder, sprintf("image_%04d.tiff", k)));
    end
end


%% 3. Creating the solver object for calibration
% A pixelwise calibration has to be made on the full-resolution grid, so
% that every pixel of a later measurement has its own coefficients,
% whatever resize factor that measurement uses. prepForCalibration() sets
% both resize factors to 1, and switches storage to phase only, since there
% is no elevation model yet.
%
% compRect restricts the work to the part of the image covered by the
% fringe. Only pixels inside it, minus the unwrap margin described next,
% get coefficients. Elsewhere, the calibration coefficients are set to NaN.

calibCase = FtpSolver("calibCase", ...
    refPath      = fullfile(imageFolder, sprintf("image_%04d.tiff", refIndex)), ...
    dataPath     = fullfile(imageFolder, "image_0001.tiff"), ...
    camCalibPath = camCalibFile, ...
    compRect     = [584 685 1420 1140]);

calibCase.prepForCalibration();


%% 4. The calibration region
% Near the edges of the computational domain the demodulated phase is
% unreliable, so a band along the edges, the unwrap margin, is set to NaN.
% The calibration is fitted over everything inside that margin.
% setOutputRectFromMargin(0) places the output window exactly there, so
% the phase and elevation read back in sections 7 and 12 cover the same
% pixels as the fit. In showDomains(), the red rectangle indicates the
% computational domain, the blue the unwrap margin and the green the output
% window, which here coincides with the margin.
% The computational domain (compRect) was chosen so that it does not exceed
% the fringe boundaries. A sudden truncation of the fringe produces edge
% artifacts in the phase and should be avoided.

calibCase.setOutputRectFromMargin(0, 0);
% overlay the domains on the reference image
calibCase.showDomains();
% overlay the domains on an animation of the image sequence
calibCase.animateRectifiedImages(1,nPlanes,2) 


%% 5. Phase correction for a calibration
% Unwrapping fixes the phase of each pixel relative to its neighbors, but
% leaves the phase map as a whole uncertain by 2*pi*n. Temporal phase
% correction, the default for measurements, resolves this by assuming the
% phase at each pixel changes by less than pi radians from one frame to the
% next. Calibration planes jump by several millimeters between images,
% which leads to a phase shift greater than pi between frames, so that
% assumption fails here. Spatial correction is used instead. It tracks the
% same fringe peak in every frame and measures its shift with respect to the
% reference, obtaining the true phase shift at that point. This fixes n in
% each image independently.
%
% Peaks are counted from one edge of the uncropped image. On construction,
% the solver chooses the tracked peak automatically. This peak is the first
% one in the reference image that is (peakClearance*fringe period) pixels
% away from the edge of the output window. peakClearance is 2.5 by default
% and can be modified. Below, the automatic choice is shown first. However,
% this peak falls outside of the unwrap margin in images 10 to 12. To
% prevent this, peakClearance is set to 6, after which the solver
% recalculates the tracked peak index. The tracked peak index can also be
% set manually. Note that changing the domain geometry (by calling
% setOutputRectFromMargin for example) resets the tracked peak based on the
% peakClearance value and the size of the output domain.

calibCase.setPhaseCorrection(method="spatial");
fprintf("Automatically chosen peak with default edge safety factor: %d\n", calibCase.phaseCorrOpts.trackedPeakInd)

calibCase.setPhaseCorrection(peakClearance=6);
fprintf("Automatically chosen peak: %d\n", calibCase.phaseCorrOpts.trackedPeakInd)

% showTrackedPeaks shows the fringe peaks identified on a line that goes
% across the fringe, with the tracked one highlighted in green. This is
% done for the reference image on the left and a selection of data images
% on the right. Check that the highlighted peak is the same physical
% fringe in every panel. If it is not, the planes where it differs will be
% off by whole multiples of 2*pi. By default the middle correction line is
% shown, but a fourth argument to the function picks another.

calibCase.showTrackedPeaks(1, nPlanes, 2);


%% 6. Demodulation
% The wavelet method analyzes the images with Morlet wavelets over a range
% of fringe wavelengths and picks, at every pixel, the one that fits best.
% The range is given in pixels of the rectified full-resolution image and
% should bracket the fringe period of the reference, printed here. The
% VortexSheddingAndWaves example compares this method with Fourier
% demodulation.

fprintf("Reference fringe period: %.1f pixels\n", calibCase.fringe.periodOrg)
calibCase.setDemodulation(method="wavelet", waveletWavelengths=[18 25], ...
                          waveletDivisions=10, waveletAngles=-90);

%% 7. Solving and inspecting the phase
% solve() processes all twelve planes and stores their phase. Before
% calibrating, it is often helpful to check the phase for errors. It should
% vary smoothly in space and monotonically with height.

calibCase.solve();
calibCase.animate(1, nPlanes, 2, quantity="phase", DataAspectRatio=[1 1 2]);    % 2 planes per second

phase = calibCase.probe([0 0], quantity="phase");
p = polyfit(phase, heightVec, 1);   % fit a line to the phase

% The phase-height relation at the center pixel. A straight line would mean
% a linear model suffices. The curvature, and the way it differs from pixel
% to pixel, is what the models below must capture.
figure
hold on
plot(phase, heightVec, "o-", LineWidth=1.5)
plot(phase, phase*p(1) + p(2), "k--", LineWidth=1.5)
hold off

xlabel("\Delta\phi (rad)"); ylabel("Plane height (mm)")
legend("h(\Delta\phi) at center pixel", "Linear fit", Location="northwest")
grid on
set(gca, 'fontsize', 15)


%% 8. A first polynomial fit, and finding outliers
% calibratePoly(order, heightVec) fits, at every pixel, a polynomial
%
%   h = a1*dPhi + a2*dPhi^2 + ...
%
% by least squares over the planes. Be default, there is no constant a0
% term since the reference plane has zero phase by definition. Use
% constOffset=true when calling calibratePoly() to add a constant term.
%
% The fit opens two figures showing the error of the model. One for the
% pixelwise coefficients and one for coefficients obtained by fitting a
% second-order surface to the pixelwise coefficients (see Sec. 10). Each
% shows the mean absolute error for every plane, a height vs phase plot for
% the pixel located at the center of the domain, and spatial maps of the
% error at four of the heights.
%
% A first fit with all planes shows that the +5 mm and +10 mm planes sit
% off the trend of the others in the error-versus-height panel. Whatever
% the cause, a plane that disagrees with all its neighbors distorts the fit
% for every pixel, so they are excluded from the calibrations in Sec. 10.
calibCase.calibratePoly(2, heightVec);

outlierHeights = ismember(heightVec, [5 10]);


%% 9. Takeda's model
% For a camera and a projector at height L above the reference plane and
% a distance d apart, the phase difference caused by a surface at height h
% gives
%
%   h = L*dPhi / (dPhi - 2*pi*f0*d),
%
% where f0 is the fringe frequency on the reference plane. calibrateTakeda
% finds the L and d that minimize the elevation error in a least squares
% sense. Bounds can be given when the geometry is roughly known, which
% helps the fit converge:
%   calibCase.calibrateTakeda(heightVec, boundsL=[Lmin Lmax], boundsD=[dmin dmax]);

calibCase.calibrateTakeda(heightVec, excludeHeights=outlierHeights);
L = calibCase.elevModel.L;
d = calibCase.elevModel.d;

% The fitted d is negative here, which is expected. Defining the image axes
% in the conventional sense, that is, the horizontal axis points to the
% right and the vertical down, then d < 0 means the projector sits on the
% negative x or y side of the camera. In the setup here, the projector sits
% on the negative-y side, i.e. beyond the top edge of the image. For the
% polynomial calibration, this sign convention is absorbed by the
% coefficients. 
%
% The sign of d follows from the phase convention used throughout the
% solver: phaseDiff is positive when the fringe moves (relative to the
% reference) in the negative axis direction. A surface raised toward the
% camera intercepts the projected light early, so its fringes shift towards
% the projector, here placed in the negative y direction. Therefore the
% measured phase shift is positive for positive elevation. The height model
%           h = L*dPhi./(dPhi - 2*pi*f0*d) 
% turns that positive phase shift into a positive height only when d
% is negative.
%
% Had the projector been placed on the other side of the camera, the same
% calibration would return the same L and a positive d with the same
% magnitude. |d| should approximate the physical camera-projector
% separation, which is a useful check on the calibration.


%% 10. The polynomial calibration
% The calibratePoly() function takes as its first argument the order of the
% polynomial. If the order is too low, it can fail to model the optical
% geometry. If it is too high, it can fit to noise instead. In our
% experiments, going beyond a second degree polynomial did not confer
% substantially increased accuracy for the intended measurement range of
% +-15 mm, though a fourth-degree polynomial does considerably reduce error
% in the (15, 25) mm range, as seen in the comparison below. The last fit
% made is the one kept by the solver, so the second-order fit comes last.
%
% Planes can also be weighted rather than dropped, for example to favor
% the heights a measurement will mostly see:
%   calibCase.calibratePoly(2, heightVec, weights=w);
%
% Every call to calibratePoly produces two things: pixelwise coefficients
% and the surfaces that were fitted to these pixelwise coefficients. The
% fitted surfaces can smooth over noise and erroneous phase results that
% affect a small number of pixels, say from a spec of dirt on the
% calibration surface. These surfaces are defined beyond the calibrated
% region too, although there they are an extrapolation. The default setting
% for a polynomial model is pixelwise, but this can be toggled on and off
% by calling setPolyPixelwise().

outlierHeightsAndZero = outlierHeights;
outlierHeightsAndZero(refIndex) = true;
h = heightVec(~outlierHeightsAndZero);

linearErr = calibCase.calibratePoly(1, heightVec, ...
            excludeHeights=outlierHeights, plot=false);
fourthOrderErr = calibCase.calibratePoly(4, heightVec, ...
            excludeHeights=outlierHeights, plot=false);
secondOrderErr = calibCase.calibratePoly(2, heightVec, ...
            excludeHeights=outlierHeights);

figure; 
plot(h, linearErr, '-o', h, secondOrderErr, '-^', ...
                         h, fourthOrderErr, '-sq', 'linewidth', 1.5)

legend("Linear fit", "Second-order fit", "Fourth-order fit", 'Location', 'north')
xlabel("Height (mm)")
ylabel("Mean absolute error (mm)")
set(gca, 'fontsize', 15)

% writePolyCalibration saves both pixelwise coefficients and the six
% coefficients of the second-order surface fitted to each coefficient plane
% to a .mat file. The file is loaded into an FtpSolver object with
% setElevModelPoly(filePath). It stays valid only for the same camera
% calibration, rectified image size, and projected fringe period. The
% Takeda model, however, is only affected by the optical setup and is valid
% as long as this setup is unchanged.
calibCase.writePolyCalibration(polyCalibFile);

%% 11. Inspecting the polynomial calibration
% The three coefficient arrays produced by the second-order polynomial fit
% are plotted here. The constant term is zero everywhere by default unless
% requested in the call to calibratePoly(). The first and second order
% coefficients show only a minor variation in space for our setup, despite
% the fact that the camera was tilted by 12.5 degrees. This is because the
% camera and projector were placed roughly at the same height (about 1 m
% above the surface), and the images were rectified to a fronto-parallel
% view before calibration, approximating the geometry of Takeda. As a
% result, extrapolating the coefficients beyond the calibrated region can
% be done with minimal error. In the next section, the accuracy of the
% polynomial and Takeda models are compared.

polyCalibStruct = load(polyCalibFile);
figTitles = {"Constant term (zero)", "First-order coefficients", "Second-order coefficients"};
[Ny, Nx, ~] = size(polyCalibStruct.coeffMatOrg);
figure
tiledlayout(1,3, 'Padding', 'compact')
for i = 1:3
    nexttile
    imagesc(polyCalibStruct.coeffMatOrg(:,:,i))
    axis equal
    xlabel('x (pixel)')
    ylabel('y (pixel)')
    colorbar
    title(figTitles{i})
    xlim([1 Nx])
    ylim([1 Ny])
    set(gca, 'fontsize', 15)
end

%% 12. Checking the models on the calibration planes
% The calibration images are now solved again, this time converting phase
% to elevation with each model in turn. The resize factor is set to 0.5 for
% speed. For each model and plane, the RMS of the elevation error over the
% region is plotted. The excluded planes are shown with open markers.

checkCase = FtpSolver("checkCase", ...
    refPath      = fullfile(imageFolder, sprintf("image_%04d.tiff", refIndex)), ...
    dataPath     = fullfile(imageFolder, "image_0001.tiff"), ...
    camCalibPath = camCalibFile, ...
    compRect     = [584 685 1420 1140], ...
    resizeFactor = 0.5);

checkCase.setOutputRectFromMargin(0,0);
checkCase.setPhaseCorrection(method="spatial", peakClearance=6);
checkCase.setDemodulation(method="wavelet", waveletWavelengths=[18 25], ...
                          waveletDivisions=10, waveletAngles=-90);

modelNames = ["Takeda", "Polynomial, pixelwise", "Polynomial, smoothed"];
meanError = zeros(nPlanes, numel(modelNames));
rmsError  = zeros(nPlanes, numel(modelNames));
symbols   = strings(1, numel(modelNames));
for m = 1:numel(modelNames)
    switch m
        case 1
            checkCase.setElevModelTakeda(L, d);
            symbols(m) = "o";
        case 2
            checkCase.setElevModelPoly(polyCalibFile);
            checkCase.setPolyPixelwise(true);
            symbols(m) = "^";
        case 3
            checkCase.setPolyPixelwise(false);
            symbols(m) = "sq";
    end
    checkCase.solve();
    elev = checkCase.outputSubarray("elevation", 1, nPlanes);
    err = elev - reshape(heightVec, 1, 1, []);
    rmsError(:, m)  = squeeze(sqrt(mean(err.^2, [1 2], "omitmissing")));
end

outlierHeightsAndZero = outlierHeights;
outlierHeightsAndZero(refIndex) = true;
figure
hold on
for m = 1:numel(modelNames)
    h = plot(heightVec(~outlierHeightsAndZero), rmsError(~outlierHeightsAndZero, m), symbols(m) + "-", LineWidth=1.5);
    plot(heightVec(outlierHeights), rmsError(outlierHeights, m), symbols(m), ...
         Color=h.Color, MarkerFaceColor="none", HandleVisibility="off", LineWidth=2)
end

hold off
set(gca, 'fontsize', 15)
xlabel("Plane height (mm)");
ylabel("RMS error (mm)");
grid on
legend(modelNames, Location="northwest")

% The calibration file written in Sec. 10 can now be used in a measurement:
%   measurement.setElevModelPoly(polyCalibFile);
%   measurement.setPolyPixelwise(false);    % or true, inside the calibrated region

%% 13. Full calibration workflow
% The script below shows what a typical calibration workflow looks like.
% calibCase = FtpSolver("calibCase", ...
%     refPath      = fullfile(imageFolder, sprintf("image_%04d.tiff", refIndex)), ...
%     dataPath     = fullfile(imageFolder, "image_0001.tiff"), ...
%     camCalibPath = camCalibFile, ...
%     compRect     = [584 685 1420 1140]);
% 
% calibCase.prepForCalibration();
% calibCase.setPhaseCorrection(method="spatial", peakClearance=6);
% calibCase.setDemodulation(method="wavelet", waveletWavelengths=[18 25], ...
%     waveletDivisions=10, waveletAngles=-90);
% 
% calibCase.solve();
% 
% calibCase.calibratePoly(2, heightVec);
% 
% calibCase.writePolyCalibration(polyCalibFile);
