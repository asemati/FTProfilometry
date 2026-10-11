%% FtpSolver tutorial: running a full case
%
% This script walks through the capabilities of FtpSolver on a recording of
% a free surface carrying both surface waves and the dimples left by
% vortices shed in the wake of a cylinder. The recording has 600 frames
% taken at 15 frames per second, together with a still-water reference
% image.
%
% The script is meant to be run one section at a time. Sections 1 to 3 must
% be run first. After that, most sections build on the one before, but the
% full run in section 14 restates every setting it uses, so it can be run
% directly after section 3 if you only want the result.
%
% The early sections run the solver on a short stretch of frames to
% compare options quickly. Only section 14 processes the whole
% recording.
%
% Data: Dataverse.no, https://doi.org/10.18710/MWEHEM
%   Download imageSet_C10_vortexWaves.nc (3.7 GB) and place it in workDir
%   below. The camera calibration and the phase-to-height calibration used
%   here are in the sampleData folder of the project.


%% 1. Paths and dependencies
% If FtpSolver is not found, add the repository to the path first with
% addpath(genpath("path/to/FTProfilometry")).

workDir = "path/to/workingDirectory";     % folder containing the .nc file

% Locate the repository from the solver itself, so that the files in
% sampleData are found wherever the repository was cloned.
repoRoot        = fileparts(fileparts(which("FtpSolver")));
sampleData      = fullfile(repoRoot, "sampleData");
camCalibFile    = fullfile(sampleData, "CameraPinholeCalibration.xml");
polyCalibFile   = fullfile(sampleData, "profCalibrationVtxWaves.mat");

ncFile          = fullfile(workDir, "imageSet_C10_vortexWaves.nc");
imageFolder     = fullfile(workDir, "imageSet");
refFile         = fullfile(workDir, "vtxWaves_refImage.tiff");
resultsFolder   = fullfile(workDir, "Results");


%% 2. Converting the dataset to an image sequence
% FtpSolver reads a recording as a folder of image files (or a LaVision
% DaVis .set/.im7 file). The dataset is distributed as NetCDF, so the
% frames are written out once as tiff files. Two details matter here.
%
% First, the solver orders frames by file name, so the frame numbers are
% zero-padded. Without padding, text order would put frame 10 before
% frame 2. Second, the reference image is written outside the image folder.
% Every file with the same data extension in that folder is treated as a
% frame, so a reference stored there would be solved as if it were part of
% the recording.
%
% The image set takes 6.8 GB of space when decompressed.

info = ncinfo(ncFile, "imageSet");
nFrames = info.Size(3);
frameTimes = ncread(ncFile, "time");
frameRate = 1/diff(frameTimes(1:2));
fprintf("%d frames at %.2g frames per second\n", nFrames, frameRate)

if ~isfolder(imageFolder)
    mkdir(imageFolder)
end

% skip the conversion if it was done on an earlier run
if numel(dir(fullfile(imageFolder, "*.tiff"))) ~= nFrames
    blockSize = 50;
    for k0 = 1:blockSize:nFrames
        nRead    = min(blockSize, nFrames - k0 + 1);
        imgBlock = ncread(ncFile, "imageSet", [1 1 k0], [Inf Inf nRead]);
        for j = 1:nRead
            imwrite(imgBlock(:, :, j), ...
                fullfile(imageFolder, sprintf("image_%04d.tiff", k0 + j - 1)));
        end
    end
end

if ~isfile(refFile)
    imwrite(ncread(ncFile, "refImage"), refFile)
end


%% 3. Creating the solver object
% A solver object holds one measurement case: its images, its settings
% and its results. The constructor takes a case ID, which names the output
% files and must be a valid MATLAB variable name, and a set of name-value
% options.
%
%   refPath       The reference image: the fringe pattern on the flat,
%                 undisturbed surface. Every frame is compared with it,
%                 and elevation is measured relative to this surface.
%   dataPath      Any one file of the image sequence. All files with the
%                 same extension in its folder are taken as frames.
%   camCalibPath  A pinhole or polynomial type camera calibration. Each
%                 image is undistorted and rectified to a view straight
%                 down onto the reference plane, and pixel positions are
%                 given world coordinates based on the units used in the
%                 calibration.
%   compRect      The computational domain, [x y width height] in pixels
%                 of the rectified image. Everything outside it is ignored,
%                 with the exception of fringe peaks used for spatial phase
%                 correction.
%   resizeFactor  Scales the images before they are processed. Halving the
%                 image dimensions quarters the number of pixels and so
%                 reduces the computation time considerably. The fringes
%                 must still be well resolved after resizing. A good rule
%                 of thumb is to have a minimum of 10 pixels per fringe
%                 period.
%
% Other options include mmPerPixel (a plain scale factor when there is no
% camera calibration), period and normalAxis (to give the fringe period
% and direction instead of detecting them), imgRotAngle (to rotate the
% rectified view), and burntPixelMask (see section 12).
%
% During construction the reference image is rectified and the fringe
% period and orientation are detected from its Fourier spectrum. The
% result is stored in vtxWaves.fringe and printed below. periodOrg is the
% period in pixels of the rectified image, period is the same quantity
% after resizing, and normalAxis is the image axis closest to the fringe
% normal, i.e. the direction along which image intensity oscillates.

vtxWaves = FtpSolver("vtxWaves", ...
    refPath      = refFile, ...
    dataPath     = fullfile(imageFolder, "image_0001.tiff"), ...
    camCalibPath = camCalibFile, ...
    resizeFactor = 0.5, ...
    compRect     = [120 320 2500 2080]);

disp(vtxWaves.fringe)

% All settings are readable as properties, grouped by topic. For example,
% procOpts holds the processing settings and demodOpts the demodulation
% settings. They are changed through the set...() methods, never by
% assigning to them directly.


%% 4. Coordinate frames and the computational domain
% A run passes through three coordinate frames.
%
% The RAW frame is the rectified camera image at full resolution. All
% user inputs given in pixels (compRect, unwrap margins, wavelet
% wavelengths) are in this frame, so they do not change meaning when the
% resize factors change.
%
% The computational frame (COMP) is the RAW image after resizing by
% resizeFactor (default = 1) and cropping to compRect. Demodulation and
% unwrapping happen here.
%
% The OUTPUT frame is what you get back. The unwrapped phase is resized
% once more by resizeFactorOutput and cropped to the output window.
% Demodulation always leads to loss of information at length-scales smaller
% than the fringe period and therefore, nothing is gained by storing the
% results at full resolution. We found the default output downsampling
% factor of 0.5 to strike a good balance between resolution and storage
% size.
%
% showDomains() draws rectangles that define the COMP (red) and OUTPUT
% (green) frames on the reference image. The unwrap margin is shown in
% blue. Near the edges of the computational domain, the phase field is
% unreliable and so a band along the border is left out of unwrapping and
% set to NaN. By default the margin is one fringe period wide across the
% fringes and a fraction of a period (0.2) along them. The output window
% sits just inside the unwrap margin by default.

vtxWaves.showDomains();

% The computational domain can also be drawn by hand, or set numerically:
%   vtxWaves.drawCompDomain();              % draw on the image
%   vtxWaves.setGeometry(compRect=[120 320 2500 2080]);
%
% The margin and the output window can be adjusted independently. After
% changing the margin, call setOutputRectFromMargin to set the output
% window relative to the new margin. Its arguments are insets in RAW
% pixels. Zero places the window exactly at the margin. The command below
% sets the unwrap margin to 10 pixels in the X (horizontal) direction and
% 40 pixels in the Y (vertical) direction, then sets the output window
% relative to the margin.
%   vtxWaves.setUnwrapMargin(10, 40);
%   vtxWaves.setOutputRectFromMargin(5,5);  % (pixels X, pixels Y)

% set the output window to lie just inside the unwrap margin (the default)
vtxWaves.setOutputRectFromMargin(); 
% this is equivalent to 
vtxWaves.setOutputRectFromMargin(0,0);

% get coordinate meshes in the OUTPUT frame
[xOut, yOut] = vtxWaves.outputCoordinates();
fprintf("Output arrays: %d x %d pixels, covering %.0f x %.0f mm\n", ...
    size(xOut, 1), size(xOut, 2), ...
    max(xOut, [], "all") - min(xOut, [], "all"), ...
    max(yOut, [], "all") - min(yOut, [], "all"))


%% 5. Looking at the images
% Before solving anything, it can be useful to observe the images. Problems
% such as reflections, bubbles, and debris are much easier to spot here
% than in the results.
%
% animateRectifiedImages shows the rectified images. animateImages shows
% them after preprocessing, exactly as the demodulation sees them:
% resized, cropped to the computational domain, and with the slowly
% varying background intensity removed (see section 6).

vtxWaves.animateRectifiedImages(400, 500, 20);    % frames 400 to 500 at 20 fps
% set the showDomains option to false to remove the domain overlays
% vtxWaves.animateRectifiedImages(400, 500, 20, showDomains=false);
vtxWaves.animateImages(400, 500, 20);


%% 6. How demodulation works: the fringe spectrum
% A fringe image can be written as
%
%   I(x, y) = a(x, y) + b(x, y) cos(2*pi*f0*x + phi(x, y))
%
% where a is the background intensity, b the fringe contrast, f0 the
% carrier frequency of the projected pattern and phi the phase that the
% surface adds by shifting the fringes. Its Fourier transform has three
% main lobes: one at zero frequency from a, and two at +f0 and -f0 carrying
% (b/2)exp(+i*phi) and (b/2)exp(-i*phi). As long as a, b and phi vary
% slowly compared with the fringes, the lobes are separate.
%
% Fourier demodulation keeps only the lobe at +f0 with a band-pass filter
% and transforms back, which leaves (b/2)exp(i*(2*pi*f0*x + phi)). Doing
% the same with the reference image and multiplying one by the complex
% conjugate of the other cancels the carrier, and the angle of the
% product is the phase difference phi - phi_ref. 
%
% The plot below shows this spectrum for one preprocessed frame. The
% markers are at +f0 and -f0, and the circle is the band-pass filter: the
% filter has fallen to 1/e at its edge, with a radius set by
% setDemodulation(filterWidthFrac=...) as a fraction of f0. Background
% removal (dcSubtraction, on by default) suppresses the central lobe and
% reduces its contribution to the filtered region. Note the presence of
% harmonics in the spectrum. Any overlap between the filter and the
% harmonics (or the DC lobe) will introduce errors into the demodulation
% process. On the other hand, a filter set too tight around the carrier
% frequency will not encompass the full lobe which will, again, introduce
% errors.

testFrame = 300;
img = vtxWaves.getImages(testFrame, testFrame);
img = img{1};
refImg = vtxWaves.getRefImage();

f0 = 1/vtxWaves.fringe.period;        % carrier, cycles per computational pixel
nVec = vtxWaves.fringe.normVec;
filterRadius = f0*vtxWaves.demodOpts.filterWidthFrac;
[ny, nx] = size(img);
fx = ((0:nx-1) - floor(nx/2))/nx;     % cycles per pixel
fy = ((0:ny-1) - floor(ny/2))/ny;

refSpectrum = abs(fftshift(fft2(refImg .* (hannWindow(ny) * hannWindow(nx)'))));   % Hann window limits spectral leakage
frameSpectrum = abs(fftshift(fft2(img .* (hannWindow(ny) * hannWindow(nx)'))));

figure
tiledlayout(1,2)
nexttile

imagesc(fx, fy, log10(refSpectrum))
axis image
axis xy
hold on
plot(f0*nVec(1)*[1 -1], f0*nVec(2)*[1 -1], "r+", MarkerSize=8, LineWidth=2)
theta = linspace(0, 2*pi, 200);
plot(f0*nVec(1) + filterRadius*cos(theta), f0*nVec(2) + filterRadius*sin(theta), "r--")
hold off
xlim(4*f0*[-1 1]); ylim(4*f0*[-1 1])
xlabel("f_x (cycles per pixel)"); ylabel("f_y (cycles per pixel)")
title("log_{10} spectrum of reference image")
colorbar

nexttile
imagesc(fx, fy, log10(frameSpectrum))
axis image
axis xy
hold on
plot(f0*nVec(1)*[1 -1], f0*nVec(2)*[1 -1], "r+", MarkerSize=8, LineWidth=2)
plot(f0*nVec(1) + filterRadius*cos(theta), f0*nVec(2) + filterRadius*sin(theta), "r--")
hold off
xlim(4*f0*[-1 1]); ylim(4*f0*[-1 1])
xlabel("f_x (cycles per pixel)"); ylabel("f_y (cycles per pixel)")
title(sprintf("log_{10} spectrum of frame %d", testFrame))
colorbar



%% 7. Fourier versus wavelet demodulation
% The Fourier method uses a single filter for the whole image. When the
% surface is steep, the local fringe frequency departs from f0 and part
% of the signal can fall outside the filter. If the filter is widened to
% catch it, noise and the other lobes start to leak in.
%
% The wavelet method instead analyses the image with Morlet wavelets, wave
% packets of finite length, at a range of scales (fringe wavelengths) and
% orientations. At every pixel it picks the scale and orientation with the
% strongest response, the "ridge", and takes the phase there. Because the
% analysis is local and adapts to the local fringe wavelength, it copes
% better with steep surfaces and sharp features. The cost is one
% transform per scale and orientation, so it is considerably slower.
%
% The wavelet scales are set by waveletWavelengths, the range of fringe
% wavelengths in RAW pixels, and waveletDivisions, the number of scales in
% that range. The range should bracket the fringe period of the reference
% (printed below) with room for the compression and stretching of the
% fringes over the surface. waveletAngles sets the orientations in degrees;
% one angle matched to the pattern is enough here, and every additional
% angle adds another set of transforms at every scale.
%
% To compare the two, both are run on a short stretch of frames and the
% phase is stored. setStoredOutputs chooses what is kept: here phase
% only, because the elevation model is introduced in section 10. The run
% time of each solve is stored in diagnostics.solveTime. By default, the
% solver uses the Fourier method.

fprintf("Reference fringe period: %.1f RAW pixels\n", vtxWaves.fringe.periodOrg)

testFrame = 501;
testRange = [501 520];      % frames used for the comparisons that follow
vtxWaves.setSolveRange(testRange(1), testRange(2));
vtxWaves.setStoredOutputs(elevation=false, phase=true);

vtxWaves.setDemodulation(method="fourier");
vtxWaves.solve();
phaseFourier = vtxWaves.outputSubarray("phase", testRange(1), testRange(2));
timeFourier = vtxWaves.diagnostics.solveTime;

vtxWaves.setDemodulation(method="wavelet", waveletWavelengths=[28 36], ...
                         waveletDivisions=10, waveletAngles=-90);
vtxWaves.solve();
phaseWavelet = vtxWaves.outputSubarray("phase", testRange(1), testRange(2));
timeWavelet = vtxWaves.diagnostics.solveTime;

fprintf("Fourier: %.2f s per frame, wavelet: %.2f s per frame\n", ...
    timeFourier/diff(testRange + [0 1]), timeWavelet/diff(testRange + [0 1]))

k = testFrame - testRange(1) + 1;
phaseDiff = phaseWavelet(:,:,k) - phaseFourier(:,:,k);

figure
tiledlayout(1, 3)
plotField(nexttile, xOut, yOut, phaseFourier(:,:,k), "Fourier", "\Delta\phi (rad)");
plotField(nexttile, xOut, yOut, phaseWavelet(:,:,k), "Wavelet", "\Delta\phi (rad)");
plotField(nexttile, xOut, yOut, phaseDiff, "Wavelet minus Fourier", "rad", [-0.05 0.05]);

fprintf("Median absolute phase difference between Fourier and Wavelet methods: %.3f rad\n", ...
                                median(abs(phaseDiff), 'all'))
fprintf("Mean absolute phase difference between Fourier and Wavelet methods: %.3f rad\n", ...
                                mean(abs(phaseDiff), 'all'))


% The difference between the two methods is very small here and is mostly
% seen at the edges. A slightly larger offset from the unwrap margin in the
% vertical direction for the output domain would prevent these edge effects
% leaking into the saved results. Though the wavelet method is slower than
% Fourier, it is particular useful when the surface slopes are steep or
% when the images very noisy.


%% 8. Unwrapping and fault detection
% Demodulation returns the phase wrapped into (-pi, pi]. Unwrapping removes
% the 2*pi jumps by assuming the true phase changes by less than pi between
% neighboring pixels, and adding or subtracting 2*pi wherever it jumps by
% more. The default ("1D") unwraps along one image direction and then the
% other. It is fast, but a single bad pixel, for example from debris, a
% bubble or a reflection, can be mistaken for a jump, and the error then
% runs on as a 2*pi step along the rest of the line.
%
% The solver checks every frame for such faults. They show up as narrow
% lines where the phase gradient is far larger than in their surroundings.
% Faults above a length threshold are recorded, frame by frame, in
% diagnostics.flaggedRegions as bounding boxes in pixel coordinates of the
% output arrays. Nothing is changed in the data. The record tells you where
% to look, or which regions to mask in later analysis.
%
% The sensitivity of the detector is set with setFaultDetection. For
% example, residualFactor sets how far the gradient must exceed its local
% median to count as a fault, and lengthToFlagPeriods how long a fault must
% be to be recorded. All lengths are specified as a multiple of the fringe
% period:
%   vtxWaves.setFaultDetection(residualFactor=500, lengthToFlagPeriods=0.5);
%   vtxWaves.setFaultDetection();     % restore the defaults
%
% With the 2D unwrapping method, a frame whose faults exceed a length
% threshold is unwrapped again with a slower two-dimensional algorithm that
% unwraps the most reliable pixels first and so routes around bad ones
% (Herraez et al., Appl. Opt. 41, 2002). This 2D unwrapping can be the most
% time consuming part of the entire pipeline for large images. You can
% switch to it with:
%   vtxWaves.setUnwrapMethod("2D")


%% 9. Correcting the phase
% Unwrapping only fixes the phase of each pixel relative to its
% neighbors. The phase map as a whole can still be off by 2*pi*n for an
% unknown integer n, and since the elevation follows from the phase, the
% whole surface would then be shifted up or down. Phase correction finds
% n for every frame. There are two strategies implemented in this solver.
%
% Temporal correction (the default) compares each frame with the previous
% one at a set of sample points and picks the n that minimizes the change
% between frames. It needs no extra information, but it assumes the surface
% moves by much less than half a fringe-equivalent height between frames at
% a majority of the sample points, and it inherits n from the first frame.
% In addition, a single failure shifts every later frame by some multiple
% of 2*pi.
%
% Spatial correction tracks one fringe peak along a few lines across the
% image. The shift of that peak from its position in the reference image,
% divided by the fringe period, gives the true phase at the peak with no
% 2*pi ambiguity. Given that unwrapping produces a continuous phase field
% (if successful), a single absolute phase value in the domain is enough to
% correct the entire phase field. With the spatial method, each frame is
% corrected on its own and a failure does not propagate. It requires that
% the tracked peak stays visible and is identified correctly in every
% frame. Peaks are counted from the edge of the full image
% (vtxWaves.phaseCorrOpts.startEdge), so the fringes between that edge and
% the tracked peak must be detected reliably too. By default the solver
% picks the first peak that lies at least k*lambda pixels inside the output
% domain, where k is 2.5 by default, and lambda is the fringe period. k can
% be set with
%   setPhaseCorrection(peakClearance=k)
% To set the tracked fringe peak by hand and change the number of sample
% points (9 by default), call
%   setPhaseCorrection(trackedPeak=n, numSamples=20)
%
% The three runs below test the phase correction methods. The mean phase
% over the output window is plotted for every frame and a failed correction
% shows up as a jump of 2*pi between neighboring frames. Observe how the
% temporal method fails because its first frame, which runs without any
% phase correction, becomes the reference for all that follow. In contrast,
% the spatial method performs correctly throughout.

testRange = [232 265]; 
vtxWaves.setSolveRange(testRange(1), testRange(2));
vtxWaves.setDemodulation(method="fourier");

vtxWaves.setPhaseCorrection(method="none");
vtxWaves.solve();
phaseNone = vtxWaves.outputSubarray("phase", testRange(1), testRange(2));

vtxWaves.setPhaseCorrection(method="temporal");
vtxWaves.solve();
phaseTemporal = vtxWaves.outputSubarray("phase", testRange(1), testRange(2));

vtxWaves.setPhaseCorrection(method="spatial");
vtxWaves.solve();
phaseSpatial = vtxWaves.outputSubarray("phase", testRange(1), testRange(2));
fprintf("Tracked peak: %d, counted from the %s edge\n", ...
    vtxWaves.phaseCorrOpts.trackedPeakInd, vtxWaves.phaseCorrOpts.startEdge)

testFrames = testRange(1):testRange(2);
figure
plot(testFrames, squeeze(mean(phaseNone, [1 2])), "o-", ...
     testFrames, squeeze(mean(phaseTemporal,  [1 2])), "x-", ...
     testFrames, squeeze(mean(phaseSpatial,     [1 2])), "k^--" , 'LineWidth', 1.5)
xlabel("Frame number"); ylabel("Mean \Delta\phi (rad)")
legend("none", "temporal", "spatial", 'Location', 'best')
ylim([-6*pi 4*pi])
grid on
set(gca, 'fontsize', 15)

% The showTrackedPeaks() function shows how the spatial peak detection
% performs for a range of selected frames and is useful when
% troubleshooting failed peak detections. The reference line with its
% detected peaks and the tracked one highlighted, then the same line in a
% few frames (here every 10th test frame). Check that the highlighted peak
% is the same physical fringe in every panel. By default the middle
% correction line is shown. A fourth argument to the function picks another
% line.
%       vtxWaves.showTrackedPeaks(testRange(1), testRange(2), 6);

% Frames where the correction could not be determined reliably are listed
% in diagnostics.phaseCorrectionFailures and left uncorrected. Note, the
% absence of a timestep in the failures list does not guarantee that the
% phase correction was successful. It just means that a majority of the
% sample points were in agreement on the phase correction.
failures = vtxWaves.diagnostics.phaseCorrectionFailures;
fprintf("Spatial phase correction failed in %d test frames\n", numel(failures))


%% 10. From phase to elevation
% For a guide on calibrating the profilometry system, see the calibration
% script in this folder. Here, we use the calibration file
% profCalibrationVtxWaves.mat in the sampleData folder. Note the projected
% fringe periods are different between this case and the calibration case.
% The calibration images were recorded with a projected fringe period of 8
% projector pixels, while this recording used 12, so the polynomial
% coefficients used here were first converted to account for this
% difference. The phase shift that a given height produces is proportional
% to the carrier frequency of the pattern, so the same height gives 8/12 as
% much phase shift with the 12-pixel pattern. Typically, one would
% calibrate and perform measurements with the same fringe period.
%
% The calibrated region is drawn in yellow below. This domain extends well
% beyond it, so the smoothed coefficients are used. The two runs that
% follow show the difference: the pixelwise model leaves the uncalibrated
% part of the domain empty.

testRange = [150, 150];
calib = load(polyCalibFile);
inCalib = ~isnan(calib.coeffMatOrg(:,:,end));
[rr, cc] = find(inCalib);
calibRegion = [min(cc), min(rr), max(cc) - min(cc), max(rr) - min(rr)];

vtxWaves.showDomains();
hold on
rectangle(Position=calibRegion, EdgeColor="y", LineStyle="--", LineWidth=2)
hold off
title("Calibrated region (yellow)")

vtxWaves.setElevModelPoly(polyCalibFile);
vtxWaves.setStoredOutputs(elevation=true, phase=false);
vtxWaves.setSolveRange(testRange(1), testRange(2));

vtxWaves.setPolyPixelwise(true);
vtxWaves.solve();
elevPixelwise = vtxWaves.outputSubarray("elevation", testRange(1), testRange(2));

vtxWaves.setPolyPixelwise(false);
vtxWaves.solve();
elevSmoothed = vtxWaves.outputSubarray("elevation", testRange(1), testRange(2));

figure
tiledlayout(1, 2)

firstClims = plotField(nexttile, xOut, yOut, elevPixelwise, ...
                "Pixelwise coefficients", "\eta (mm)");

plotField(nexttile, xOut, yOut, elevSmoothed,  ...
                "Smoothed coefficients (extrapolated)",  "\eta (mm)", firstClims);


%% 11. Parallax correction
% The camera sees the surface along slanted rays. Where the surface is
% raised or lowered, the ray through a given pixel meets it at a point
% shifted sideways from where it would meet the reference plane, so the
% height measured at that pixel belongs to a slightly different position.
% The shift grows with the elevation and with the angle of the ray, so it
% increases with distance from the camera's optical axis (approximately the
% middle of the image in our case).
%
% With a pinhole camera calibration, setParallaxCorrection(true) traces
% every pixel's ray to the measured surface, finds the true (x,y) position
% of each height sample, and interpolates the samples back onto the regular
% output grid. Due to the fact that the ratio between the surface elevation
% and camera height is small here (h/L = 0.01), the parallax effect is
% minimal.

testFrame = 500;
vtxWaves.setSolveRange(testFrame, testFrame);

vtxWaves.setParallaxCorrection(false);
vtxWaves.setDemodulation(method="wavelet");
vtxWaves.setStoredOutputs(elevation=true, phase=false);

vtxWaves.solve();
elevNotCorrected = vtxWaves.outputSubarray("elevation", testFrame, testFrame);

vtxWaves.setParallaxCorrection(true);
vtxWaves.solve();
elevCorrected = vtxWaves.outputSubarray("elevation", testFrame, testFrame);

figure
tiledlayout(1, 2)
plotField(nexttile, xOut, yOut, elevCorrected, "Parallax corrected", "\eta (mm)");
plotField(nexttile, xOut, yOut, elevCorrected - elevNotCorrected, ...
                            "Corrected minus uncorrected", "mm");


%% 12. Preprocessing options
% None of the options in this section are needed for this dataset, but
% they are worth knowing.
%
% Frames can be discarded automatically when their mean intensity in the
% computational domain exceeds a threshold, for example from a laser pulse
% (when combining with particle image velocimetry). Discarded frames are
% stored as NaN, so frame numbering is kept. To activate automatic
% discarding of frames:
%   vtxWaves.setPreprocessing(discardThreshold=2000);  % =[] switches it off
%
% Burnt (stuck or hot) camera pixels can be filled from their neighbors
% before anything else is done. The mask is a logical array the size of
% the raw, unrectified camera image, true at the burnt pixels:
%   vtxWaves.setPreprocessing(burntPixelMask=mask, burntPixelInterpMethod="fastLinear");
%
% Background removal (dcSubtraction, on by default) and smoothing of the
% reference image along the fringes (refImgSmoothing) are switched the
% same way:
%   vtxWaves.setPreprocessing(dcSubtraction=true, refImgSmoothing=false);


%% 13. Resolution, memory and streaming
% Two settings control resolution independently. resizeFactor (section 3)
% sets the resolution of the computation and so its cost.
% resizeFactorOutput sets how much the unwrapped phase is downsampled
% before it is stored, and so the size of everything that comes out of a
% run. Both are changed with setGeometry:
%   vtxWaves.setGeometry(resizeFactor=0.5, resizeFactorOutput=1);
%
% The results are held in memory in single precision, one array per
% stored quantity. For this case the full run needs about:

stackSize = size(vtxWaves.worldCoords.mesh.x);
fprintf("%.2f GB per stored quantity for %d frames\n", ...
    4*prod(stackSize)*nFrames/1e9, nFrames)

% Long recordings at high resolution can exceed available memory. To avoid
% this, stream results to disk with enableStreaming(bufferLength), where
% bufferLength is the maximum number of frames held in memory. Each full
% buffer is written to the result files in the save folder during the run.
% The files are opened only while being written, so an interrupted run
% does not corrupt the files unless it occurs during the write process.
%   vtxWaves.setSaveFolder(resultsFolder);
%   vtxWaves.enableStreaming(200);   % hold 200 frames in memory
%   vtxWaves.disableStreaming();     % back to holding everything in memory
%
% With streaming on, only the last buffer remains in memory after the run.


%% 14. The full run
% The settings chosen above are restated here, so this section can also be
% run directly after section 3. Toggle streaming to disk on and aff with
% the useStreaming switch. This run takes about 6 minutes with the fourier
% demodulation method and 13 minutes with the wavelet settings below on an
% Intel Core i9-12900 desktop computer.

solveRange = [1 600];
vtxWaves.setSolveRange(solveRange(1), solveRange(2));
vtxWaves.setDemodulation(method="fourier");
% vtxWaves.setDemodulation(method="wavelet", waveletWavelengths=[28 36], ...
%                          waveletDivisions=10, waveletAngles=-90);
vtxWaves.setUnwrapMethod("1D");
vtxWaves.setPhaseCorrection(method="spatial");
vtxWaves.setElevModelPoly(polyCalibFile);
vtxWaves.setPolyPixelwise(false);
vtxWaves.setParallaxCorrection(false);
vtxWaves.setStoredOutputs(elevation=true, phase=true);

vtxWaves.setSaveFolder(resultsFolder);

useStreaming = false;
if useStreaming
    vtxWaves.enableStreaming(50);
else
    vtxWaves.disableStreaming();
end

vtxWaves.solve();
fprintf("Solved %d frames in %.0f s\n", diff(solveRange) + 1, vtxWaves.diagnostics.solveTime)


%% 15. Saving the results
% writeCase writes everything into <saveFolder>/<caseID>:
%
%   vtxWaves.mat              the solver object itself, without the data
%   vtxWaves_elevation.h5     the elevation field in units of (mm)
%   vtxWaves_phase.h5         the phase field in units of (rad)
%   coordinateMesh.mat        x and y world coordinates (mm) of every pixel
%                             of the result arrays
%   vtxWaves_diagnostics.mat  flagged regions and phase-correction failures
%   vtxWaves_log.txt          every setting of the run, in plain text
%
% With streaming on, the phase and elevation result files already exist and
% writeCase() only checks that they are complete before writing the rest

vtxWaves.writeCase();
caseFolder = fullfile(resultsFolder, "vtxWaves");
type(fullfile(caseFolder, "vtxWaves_log.txt"))


%% 16. Inspecting the results
% animate() plays the elevation as a shaded surface. The options used here
% are camFPS, the acquisition frame rate, which turns the frame counter
% into a time stamp, DataAspectRatio, which here exaggerates the vertical
% scale so that millimeter-scale features are visible over a domain
% hundreds of millimeters wide, and smoothing, a light Gaussian filter for
% display only. subtractMean=true would remove the time-averaged surface
% and show only the fluctuations. If the run was streamed, use the data
% option of section 19 instead, since only the last buffer is in memory.
nFrames = diff(vtxWaves.procOpts.solveRange) + 1;
vtxWaves.animate(1, nFrames, 15, camFPS=frameRate, lightPosition=[1 1 5], ...
                 DataAspectRatio=[10 10 1], smoothing=true);

% probe extracts the time series at a point given in world coordinates.
% With patchSize it combines the values in a square or circle around the
% point, which reduces noise. Here the centre of the output window is
% probed, both at a single point and averaged over a 30 mm diameter circle.
[xOut, yOut] = vtxWaves.outputCoordinates();
probePoint = [mean(xOut, "all"), mean(yOut, "all")];
etaPoint = vtxWaves.probe(probePoint);
etaPatch = vtxWaves.probe(probePoint, patchSize=30, shape="circle", operation="mean");

t = frameTimes(1:nFrames);
figure
plot(t, etaPoint, t, etaPatch)
xlabel("t (s)"); ylabel("\eta (mm)")
legend("single point", "30 mm circle")
title(sprintf("Elevation at (x,y) = (%.0f, %.0f) mm", probePoint))
set(gca, 'fontsize', 15)

%% 17. Exporting a video
% exportVideo takes the same display options as animate and writes an MP4
% file named after the case. With fast=true each frame is captured straight
% from the screen. Without it, each frame is first exported as a
% high-resolution image and the video is assembled at the end, which looks
% better but takes considerably longer.

animationFPS = 15;
vtxWaves.exportVideo(1, nFrames, animationFPS, saveAddr=caseFolder, fast=true, ...
                     camFPS=frameRate, DataAspectRatio=[5 5 1], smoothing=true);


%% 18. Reloading a saved case
% The saved object holds every setting but no images or results. Loading
% it gives back a solver that can be inspected, modified and solved
% again. The images are read again from the stored paths, so the data must
% still be where it was.

clear vtxWaves
saved = load(fullfile(caseFolder, "vtxWaves.mat"));
vtxWaves = saved.vtxWaves;
fprintf("Reloaded case %s: frames %d to %d, %s demodulation\n", vtxWaves.caseID, ...
    vtxWaves.procOpts.solveRange, vtxWaves.demodOpts.method)

% For example, to solve a different stretch of frames at full output
% resolution and also keep the phase:
%   vtxWaves.setSolveRange(100, 200);
%   vtxWaves.setGeometry(resizeFactorOutput=1);
%   vtxWaves.setStoredOutputs(elevation=true, phase=true);
%   vtxWaves.solve();

%% 19. Reading the results without the solver
% The result files are meant to be used on their own. Each HDF5 file holds
% one dataset of size (rows, columns, frames) with attributes giving its
% units and the dataset frame number of its first frame. FtpSolver.readData
% reads all frames or a range of them. Note that the frame numbers provided
% to readData() count from 1 at the start of the file. For example, if
% solveRange is set to [501 600] and the results are saved to a file, then
% the firstFrame and lastFrame arguments given to readData() to read the
% full 100 frames are (1, 100), not (501, 600). In Python (h5py), the same
% dataset appears with its dimensions in reverse order, [frames, columns,
% rows].

elevFile = fullfile(caseFolder, "vtxWaves_elevation.h5");
fprintf("Units: %s, first frame of the file is dataset frame %d\n", ...
    h5readatt(elevFile, "/elevation", "units"), ...
    h5readatt(elevFile, "/elevation", "firstTimestep"))

[eta, fullSize] = FtpSolver.readData(elevFile, 1, diff(solveRange) + 1);
fprintf("The file holds %d x %d x %d values\n", fullSize)

mesh = load(fullfile(caseFolder, "coordinateMesh.mat"));    % xMesh, yMesh
diagnosticsFile = load(fullfile(caseFolder, "vtxWaves_diagnostics.mat"));

% animate also accepts data from outside the object in the form of a cell
% array containing the data and, optionally, its coordinates.
vtxWaves.animate(1, 150, 15, data={eta, mesh.xMesh, mesh.yMesh}, ...
    camFPS=frameRate, DataAspectRatio=[10 10 1]);
% without the coordinate arrays
%       vtxWaves.animate(1, 100, 15, data={eta}, camFPS=frameRate, DataAspectRatio=[10 10 1]);

%% Local helper functions

function plottedClims = plotField(ax, xMesh, yMesh, field, titleText, colorbarLabel, clims)
    if ~exist('clims', 'var')
        clims = [];
    end
    pcolor(ax, xMesh, yMesh, field, 'EdgeColor','none')
    axis(ax, "image")
    cb = colorbar(ax);
    cb.Label.String = colorbarLabel;
    if ~isempty(clims)
        clim(ax, clims)
    else
        plottedClims = get(ax, 'CLim');
    end
    xlabel(ax, "x (mm)"); ylabel(ax, "y (mm)")
    title(ax, titleText)
    set(ax, 'fontsize', 15)
end

function w = hannWindow(n)
    if n == 1
        w = 1;
    else
        w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n - 1)));
    end
end