classdef FtpSolver < handle
% FtpSolver  Fourier Transform Profilometry for time-resolved surface elevation.
%
%   Reconstructs surface elevation fields from sequences of fringe-pattern
%   images using the Fourier Transform Profilometry (FTP) method. Supports
%   both standard FT demodulation and continuous wavelet demodulation, two
%   phase-to-height calibration models, and two phase correction strategies
%   (spatial and temporal).
%
%   WORKFLOW
%     profCase = FtpSolver(caseID, ...
%                             refPath="path/to/refImage.tiff", ...
%                             dataPath="path/to/data.tiff", ...
%                             resizeFactor=0.7, ...
%                             compRect=[100 150 1500 1200], ...
%                             camCalibPath="path/to/calibration.mat" ...
%                             );
%     profCase.drawCompDomain();
%     profCase.setSaveFolder("path/to/results");
%     profCase.solve();
%     profCase.animate(1, 100, framesPerSecond);
%     profCase.writeCase();
%
%   KEY METHODS
%     solve             - Run the full processing pipeline
%     drawCompDomain    - Interactively define the computational domain
%     calibratePoly     - Calibrate phase-to-height model
%     probe             - Extract time series at a spatial location
%     animate           - Visualize the reconstructed surface
%     setSaveFolder     - Choose where writeCase (and streamed output) writes
%     writeCase         - Save results, settings log and object
%     readData          - Read frames back from a saved result file
%
%   DEPENDENCIES
%     xml2struct, rectifyImagePinhole, readDavisCalibration
%     unwrap2D (only for setUnwrapMethod("2D"))
%     Davis readimx library (for .set and .im7 file formats)
%     Image Processing Toolbox
%     Curve Fitting Toolbox
%     Signal Processing Toolbox
%     Computer Vision Toolbox
%     Wavelet Toolbox
%
%   REFERENCE
%     Semati, A., Shankaran, A., Smeltzer, B.K. et al. Simultaneous
%     free-surface profilometry and subsurface velocimetry with fringe
%     projection and PIV. Exp Fluids 67, 120 (2026).
%
% Author: Ali Semati
% April 2024;

%------------- BEGIN CODE --------------

properties (SetAccess = private)
    caseID              % object name
    source              % data attributes
    imgData             % reference and current images
    phaseData           % computed phase
    elevData            % computed surface elevation
    fringe              % fringe pattern parameters
    demodOpts           % demodulation settings
    phaseCorrOpts       % phase correction settings and state
    procOpts            % preprocessing and pipeline settings and flags
    loopState           % timestep bookkeeping
    cameraCalib         % camera calibration data
    worldCoords         % spatial coordinates
    elevModel           % phase to elevation conversion parameters
    outputConfig        % output/save control
    diagnostics         % diagnostics of the run
end

%% ----------- Setup and calibration -----------
methods (Access = public)
    function obj = FtpSolver(caseID, opts)
        arguments
            caseID              {mustBeTextScalar}              = "default"
            opts.refPath        {mustBeTextScalar}              = ""
            opts.dataPath       {mustBeTextScalar}              = ""
            opts.mmPerPixel     (1,1) double                    = NaN
            opts.period         (1,1) double                    = NaN
            opts.normalAxis     string {mustBeMember(opts.normalAxis, ["X", "Y"])} = "X"
            opts.resizeFactor   (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.resizeFactor, 1)} = 1
            opts.compRect       double                          = []
            opts.camCalibPath   {mustBeTextScalar}              = ""
            opts.camCalibNum    (1,1) double                    = 1
            opts.dataCamInd     (1,1) {mustBeInteger}           = 1
            opts.refCamInd      (1,1) {mustBeInteger}           = 1
            opts.imgRotAngle    (1,1) double                    = 0
            opts.burntPixelMask                                 = []
        end

        if ~isvarname(caseID)
            error("caseID must be a valid variable name.")
        end

        obj.caseID          = char(caseID);
        obj.source          = FtpSolver.defaultSource();
        obj.imgData         = FtpSolver.defaultImgData();
        obj.phaseData       = FtpSolver.defaultPhaseData();
        obj.elevData        = FtpSolver.defaultElevData();
        obj.fringe          = FtpSolver.defaultFringe();
        obj.phaseCorrOpts   = FtpSolver.defaultPhaseCorrOpts();
        obj.demodOpts       = FtpSolver.defaultDemodOpts();
        obj.procOpts        = FtpSolver.defaultProcOpts();
        obj.loopState       = FtpSolver.defaultLoopState();
        obj.cameraCalib     = FtpSolver.defaultCameraCalib();
        obj.worldCoords     = FtpSolver.defaultWorldCoords();
        obj.elevModel       = FtpSolver.defaultElevModel();
        obj.outputConfig    = FtpSolver.defaultOutputConfig();
        obj.diagnostics     = FtpSolver.defaultDiagnostics();


        obj.procOpts.resizeFactor       = opts.resizeFactor;
        obj.procOpts.compRectOrg        = opts.compRect;
        obj.procOpts.imgRotAngle        = opts.imgRotAngle;
        obj.procOpts.burntPixelMask     = logical(opts.burntPixelMask);
        obj.source.dataCamInd           = opts.dataCamInd;
        obj.source.refCamInd            = opts.refCamInd;

        if opts.refPath == ""
            return   % bare construction used by loadobj()
        end

        obj.readRefImage(opts.refPath)
        obj.initDataset(opts.dataPath)

        if ~isnan(opts.mmPerPixel)
            obj.worldCoords.scaling.X.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.Y.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.mmPerPixelOrg   = opts.mmPerPixel;
        end

        if strlength(opts.camCalibPath) > 0
            obj.readCameraCalibration(opts.camCalibPath, opts.camCalibNum);
        end

        obj.checkBurntPixelMask(obj.procOpts.burntPixelMask);
        if ~isempty(obj.procOpts.burntPixelMask)
            [obj.procOpts.burntPixelRows, ...
                obj.procOpts.burntPixelCols] = find(obj.procOpts.burntPixelMask);
        end

        obj.fringe.periodOrg = opts.period;
        obj.fringe.period = opts.period;

        obj.fringe.normalAxis = opts.normalAxis;
        if strcmpi(opts.normalAxis, 'X')
            obj.fringe.normalVec = [1, eps];
            obj.phaseCorrOpts.startEdge = "left";
        elseif strcmpi(opts.normalAxis, 'Y')
            obj.fringe.normalVec = [eps, 1];
            obj.phaseCorrOpts.startEdge = "top";
        end


        obj.rectifyRefImage();

        % auto-detect the fringe period and pattern axis when none given
        if isnan(opts.period)
            img = imcrop(obj.imgData.refRaw, obj.procOpts.compRectOrg);
            [period, normalAxis, tiltAngle, normalVec] = obj.analyzeFringe(img);
            obj.fringe.periodOrg = period;
            obj.fringe.normalAxis = normalAxis;
            obj.fringe.normalVec = normalVec;
            obj.fringe.tiltAngle = tiltAngle;
            if strcmpi(normalAxis, 'Y')
                obj.phaseCorrOpts.startEdge = "top";
            else
                obj.phaseCorrOpts.startEdge = "left";
            end
            obj.updateGeometry();
        end
        if strcmpi(obj.fringe.normalAxis, "Y")
            obj.setDemodulation(waveletAngles=-90);
        end
    end

    function setSolveRange(obj, r1, r2)
    % Set the frame range processed by solve().
    % Call without arguments to process all available images.
    %   Inputs:
    %           r1 - first frame
    %           r2 - last frame
        arguments
            obj
            r1 (1,1) double {mustBeInteger, mustBePositive} = obj.source.fullRange(1)
            r2 (1,1) double {mustBeInteger, mustBePositive} = obj.source.fullRange(2);
        end

        if r2 > obj.source.fullRange(2) || r1 > r2
            error("Invalid input. Available data range is [1 %g]", obj.source.fullRange(2));
        end
        obj.procOpts.solveRange = [r1 r2];
    end

    function setGeometry(obj, opts)
    % Set the computational domain and the resize factors.
    %   setGeometry(compRect=[100 150 1500 1200], resizeFactor=0.7, ...
    %                                             resizeFactorOutput=0.5)
    %   Options not named keep their current setting.
    %   compRect            - computational domain in the RAW frame,
    %                         [x y width height]; [] uses the full image.
    %                         Restores automatic selection of the tracked
    %                         fringe peak.
    %   resizeFactor        - factor (between 0 and 1) applied to the images
    %                         before processing to reduce computational cost
    %   resizeFactorOutput  - factor (between 0 and 1) applied after the
    %                         phase is calculated to reduce storage size
        arguments
            obj
            opts.compRect            {mustBeNumeric}
            opts.resizeFactor        (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.resizeFactor, 1)}
            opts.resizeFactorOutput (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.resizeFactorOutput, 1)}
        end
        geometryChanged = false;

        if isfield(opts, 'compRect')
            if ~isempty(opts.compRect) && numel(opts.compRect) ~= 4
                error("compRect must be [x y width height], or [] for the full image.")
            end
            obj.procOpts.compRectOrg = reshape(double(opts.compRect), 1, []);
            obj.procOpts.outputRectOrg = [];
            obj.phaseCorrOpts.trackedPeakInd = [];
            geometryChanged = true;
        end
        if isfield(opts, 'resizeFactor')
            if obj.procOpts.resizeFactor ~= opts.resizeFactor
                obj.procOpts.resizeFactor = opts.resizeFactor;
                geometryChanged = true;
            end
        end
        if isfield(opts, 'resizeFactorOutput')
            if obj.procOpts.resizeFactorOutput ~= opts.resizeFactorOutput
                obj.procOpts.resizeFactorOutput = opts.resizeFactorOutput;
                geometryChanged = true;
            end
        end
        if geometryChanged
            obj.updateGeometry();
        end
    end

    function drawCompDomain(obj)
    % Define the computational domain interactively.
       % Shows the raw reference image and waits for a rectangle to be
       % drawn. The result becomes compRectOrg, the output rectangle is
       % reset to one unwrap margin inside it, and the geometry is refreshed.

        if isempty(obj.imgData.refScaled)
            error("Missing reference image.")
        end

        figure()
        tempHndl = imshow(obj.imgData.refRaw, []);

        while true
            rect = drawrectangle(tempHndl.Parent);
            if isvalid(rect)
                rect.Position = round(rect.Position);
                if ~isempty(rect.Position)
                    obj.procOpts.compRectOrg = round(rect.Position);
                    obj.procOpts.outputRectOrg = [];   % refilled by updateGeometry
                    break
                end
            else
                break
            end
        end

        if isvalid(tempHndl)
            close(tempHndl.Parent.Parent)
        end

        obj.phaseCorrOpts.trackedPeakInd = [];
        obj.updateGeometry();
    end

    function setOutputRectFromMargin(obj, offsetX, offsetY)
    % Set the output domain as an offset from the unwrap margin in RAW pixels.
        arguments
            obj
            offsetX (1,1) double {mustBeNonnegative, mustBeInteger} = 0;
            offsetY (1,1) double {mustBeNonnegative, mustBeInteger} = 0;
        end
        obj.procOpts.outputRectOrg = obj.marginOutputRect([offsetX offsetY]);
        obj.phaseCorrOpts.trackedPeakInd = [];
        obj.updateGeometry();
    end

    function setUnwrapMargin(obj, marginX, marginY)
    % Set the margins, in RAW-frame pixels, that are ignored by the unwrapping algorithm.
    %   setUnwrapMargin() restores the default, which follows the fringe
    %   period: one period along the normal axis and 0.2 periods along the
    %   other axis.
        arguments
            obj
            marginX {mustBeNumeric, mustBeScalarOrEmpty, mustBeNonnegative, mustBeInteger} = []
            marginY {mustBeNumeric, mustBeScalarOrEmpty, mustBeNonnegative, mustBeInteger} = []
        end
        if nargin == 1
            obj.procOpts.unwrapMarginOrg = [];
        elseif nargin == 3 && ~isempty(marginX) && ~isempty(marginY)
            obj.procOpts.unwrapMarginOrg = double([marginX, marginY]);
        else
            error("Give both margins, or none to restore the default.")
        end
        obj.updateGeometry();
    end

    function setPreprocessing(obj, opts)
    % Set the image preprocessing options.
    %   setPreprocessing(dcSubtraction=true, discardThreshold=2000)
    %   Options not named keep their current setting.
    %   dcSubtraction          - remove the background (DC) intensity from
    %                            each image before demodulation
    %                            (default true)
    %   refImgSmoothing        - smooth the reference image along the
    %                            fringes (default false)
    %   discardThreshold       - frames whose mean intensity inside the
    %                            computational domain exceeds this value
    %                            are stored as NaN; [] disables (default)
    %   burntPixelMask         - logical mask of burnt pixels, the size of
    %                            the raw, unrectified image; [] for none
    %   burntPixelInterpMethod - how burnt pixels are filled: "fastLinear"
    %                            (4-neighbor average, default) or "smooth"
    %                            (regionfill)
        arguments
            obj
            opts.dcSubtraction          (1,1) logical
            opts.refImgSmoothing        (1,1) logical
            opts.discardThreshold       {mustBeNumeric, mustBeScalarOrEmpty}
            opts.burntPixelMask         {mustBeNumericOrLogical}
            opts.burntPixelInterpMethod (1,1) string {mustBeMember(opts.burntPixelInterpMethod, ["fastLinear", "smooth"])}
        end
        if isfield(opts, 'dcSubtraction')
            obj.procOpts.dcSubtraction = opts.dcSubtraction;
        end
        if isfield(opts, 'refImgSmoothing')
            obj.procOpts.refImgSmoothing = opts.refImgSmoothing;
        end
        if isfield(opts, 'discardThreshold')
            obj.procOpts.discardThreshold = double(opts.discardThreshold);
        end
        if isfield(opts, 'burntPixelMask')
            obj.checkBurntPixelMask(opts.burntPixelMask);
            obj.procOpts.burntPixelMask = logical(opts.burntPixelMask);
            [obj.procOpts.burntPixelRows, ...
             obj.procOpts.burntPixelCols] = find(obj.procOpts.burntPixelMask);
        end
        if isfield(opts, 'burntPixelInterpMethod')
            obj.procOpts.burntPixelInterpMethod = opts.burntPixelInterpMethod;
        end

        maskChanged   = isfield(opts, 'burntPixelMask');
        methodMatters = isfield(opts, 'burntPixelInterpMethod') ...
                                   && ~isempty(obj.procOpts.burntPixelMask);
        if (maskChanged || methodMatters) && ~isempty(obj.imgData.refRaw)
            obj.readRefImage(obj.source.refPath);
            obj.rectifyRefImage();
        end
    end

    function setParallaxCorrection(obj, state)
    % Switch the parallax error correction on or off (default off).

        arguments
            obj
            state (1,1) logical
        end
        obj.procOpts.parallaxCorrection = state;
        if state && ~isfield(obj.cameraCalib, 'intrinsicsMatlab')
            warning("FtpSolver:parallaxCorrectionWithoutPinhole", ...
                "No pinhole camera calibration is loaded. Parallax error " + ...
                "correction is skipped until one is read.")
        end
    end

    function setDemodulation(obj, opts)
    % Set the demodulation method and its parameters.
    %   setDemodulation(method="fourier", filterWidthFrac=0.6)
    %   setDemodulation(method="wavelet", waveletWavelengths=[10 16], ...
    %                   waveletDivisions=20, waveletAngles=0)
    %   Options not named keep their current setting.
    %   method             - "fourier" (default) or "wavelet"
    %   Fourier method:
    %   filterWidthFrac    - width of the bandpass filter as a fraction of
    %                        the carrier frequency, between 0 and 1
    %                        (default 0.6)
    %   Wavelet method:
    %   waveletWavelengths - [min max] fringe wavelength, in RAW-frame
    %                        pixels, spanned by the wavelet scales
    %   waveletDivisions   - number of scales between min and max
    %   waveletAngles      - wavelet orientations in degrees, mapped to
    %                        (-135, 45] for consistency with the Fourier
    %                        method (default 0)
    %   waveletSigma       - width of the Morlet wavelet (default 0.6)
    %   waveletGamma       - anisotropy of the Morlet wavelet (default 1)
        arguments
            obj
            opts.method             (1,1) string {mustBeMember(opts.method, ["fourier", "wavelet"])}
            opts.filterWidthFrac    (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.filterWidthFrac, 1)}
            opts.waveletWavelengths (1,2) double {mustBePositive}
            opts.waveletDivisions   (1,1) double {mustBeInteger, mustBePositive}
            opts.waveletAngles      double {mustBeVector}
            opts.waveletSigma       (1,1) double {mustBePositive}
            opts.waveletGamma       (1,1) double {mustBePositive}
        end
        if isfield(opts, 'method')
            obj.demodOpts.method = opts.method;
        end
        if isfield(opts, 'filterWidthFrac')
            obj.demodOpts.filterWidthFrac = opts.filterWidthFrac;
        end
        if isfield(opts, 'waveletWavelengths')
            if opts.waveletWavelengths(1) > opts.waveletWavelengths(2)
                error("waveletWavelengths must be [min max].")
            end
            obj.demodOpts.waveletWavelengths = opts.waveletWavelengths;
        end
        if isfield(opts, 'waveletDivisions')
            obj.demodOpts.waveletDivisions = opts.waveletDivisions;
        end
        if isfield(opts, 'waveletAngles')
            obj.demodOpts.waveletAngles = mod(opts.waveletAngles - 45, -180) + 45;
        end
        if isfield(opts, 'waveletSigma')
            obj.demodOpts.waveletSigma = opts.waveletSigma;
        end
        if isfield(opts, 'waveletGamma')
            obj.demodOpts.waveletGamma = opts.waveletGamma;
        end

        % the wavelet scales are derived from these in updateGeometry
        if isfield(opts, 'waveletWavelengths') || isfield(opts, 'waveletDivisions')
            obj.updateGeometry();
        end
    end

    function setUnwrapMethod(obj, method)
    % Set the unwrapping algorithm.
    %   "1D"   - unwrap row-wise, then column-wise (default)
    %   "2D"   - as "1D", but when the 1-D result contains phase faults
    %            longer than lengthFor2DPeriods fringe periods (see
    %            setFaultDetection) it is replaced by the result of unwrap2D
    %   "none" - leave the phase wrapped
        arguments
            obj
            method (1,1) string {mustBeMember(method, ["1D", "2D", "none"])}
        end
        obj.procOpts.unwrapMethod = method;
    end

    function setFaultDetection(obj, opts)
    % Set how unwrapAndFlag detects and records phase faults.
    %   setFaultDetection(residualFactor=500, lengthToFlagPeriods=0.5)
    %   Options not named keep their current setting.
    %   setFaultDetection() restores all defaults.
    %   A pixel is marked as part of a fault when the phase gradient
    %   magnitude there differs from its median-smoothed value by more
    %   than residualFactor times the median of that difference. Lengths
    %   are in fringe periods. A fault's length is set as half the number
    %   of pixels it marks. See defaultProcOpts() for default values.
    %
    %   medianWindowPeriods - width of the median filter. The
    %                         filter is made odd and at least 5 pixels wide
    %
    %   residualFactor      - threshold as a multiple of the median
    %                         difference
    %
    %   lengthFor2DPeriods  - with setUnwrapMethod("2D"), fault length above
    %                         which unwrap2D is used
    %
    %   lengthToFlagPeriods - fault length above which the faults of a frame
    %                         are recorded in diagnostics.flaggedRegions
    %
    %   dilationPeriods     - padding added around recorded faults
    %
    %   closePeriods        - radius of the closing that merges nearby
    %                         faults into one region
    %
        arguments
            obj
            opts.medianWindowPeriods (1,1) double {mustBePositive}
            opts.residualFactor      (1,1) double {mustBePositive}
            opts.lengthFor2DPeriods  (1,1) double {mustBeNonnegative}
            opts.lengthToFlagPeriods (1,1) double {mustBeNonnegative}
            opts.dilationPeriods     (1,1) double {mustBeNonnegative}
            opts.closePeriods        (1,1) double {mustBeNonnegative}
        end
        if isempty(fieldnames(opts))
            defaults = FtpSolver.defaultProcOpts();
            obj.procOpts.faultDetection = defaults.faultDetection;
            return
        end
        for name = string(fieldnames(opts))'
            obj.procOpts.faultDetection.(name) = opts.(name);
        end
    end

    function setPhaseCorrection(obj, opts)
    % Set how the 2*n*pi offset of each frame is fixed.
    %   setPhaseCorrection(method="spatial", trackedPeak=3)
    %   Options not named keep their current setting.
    %   method - "temporal" (default): match each frame's phase to the
    %               previous frame. Fails if the acquisition frequency is
    %               too low. A failure at a single timestep propagates to
    %               the end.
    %            "spatial": track a fringe peak on the correction lines and
    %               compare it to the same peak in the reference image.
    %               Requires the tracked peak to remain visible throughout
    %               the recording.
    %            "none": leave the unwrapped phase as is.
    %   numSamples         - number of points (temporal) or lines (spatial)
    %                        used for the phase correction algorithm
    %
    %   Spatial method:
    %   trackedPeak        - index of the tracked fringe peak, counted from
    %                        startEdge; [] selects it automatically
    %   peakClearance      - the automatic selection picks the first peak
    %                        at least this many fringe periods inside the
    %                        output window, measured from its startEdge
    %                        side
    %   startEdge          - edge of the full image (not the computational
    %                        domain) from which the peaks are counted:
    %                        "left" or "right" when fringe.normalAxis is
    %                        "X", "top" or "bottom" when it is "Y"
    %   peakMinFrac        - minimum peak height, as a fraction of the
    %                        median intensity on the line
    %   peakProminenceFrac - minimum peak prominence, as a fraction of the
    %                        5-95 percentile intensity range on the line
    %
    %   Naming startEdge, peakClearance, peakMinFrac, numSamples or
    %   peakProminenceFrac restores automatic peak selection, unless
    %   trackedPeak is named in the same call.
        arguments
            obj
            opts.method             (1,1) string {mustBeMember(opts.method, ["spatial", "temporal", "none"])}
            opts.numSamples         (1,1) double {mustBeInteger, mustBePositive}
            opts.trackedPeak        {mustBeNumeric, mustBeScalarOrEmpty, mustBeInteger, mustBePositive}
            opts.peakClearance      (1,1) double {mustBePositive}
            opts.startEdge          (1,1) string {mustBeMember(opts.startEdge, ["left", "right", "top", "bottom"])}
            opts.peakMinFrac        (1,1) double {mustBeNonnegative}
            opts.peakProminenceFrac (1,1) double {mustBeNonnegative}
        end
        if isfield(opts, 'method')
            obj.phaseCorrOpts.method = opts.method;
        end
        if isfield(opts, 'numSamples')
            obj.phaseCorrOpts.numSamples = opts.numSamples;
            obj.placeCorrectionLines();
        end

        if isfield(opts, 'startEdge')
            if strcmpi(obj.fringe.normalAxis, 'Y')
                validEdges = ["top", "bottom"];
            else
                validEdges = ["left", "right"];
            end
            if ~any(strcmpi(opts.startEdge, validEdges))
                error("startEdge must be ""%s"" or ""%s"" when fringe.normalAxis is ""%s"".", ...
                    validEdges(1), validEdges(2), obj.fringe.normalAxis)
            end
            obj.phaseCorrOpts.startEdge = opts.startEdge;
        end
        if isfield(opts, 'peakClearance')
            obj.phaseCorrOpts.peakClearance = opts.peakClearance;
        end
        if isfield(opts, 'peakMinFrac')
            obj.phaseCorrOpts.peakMinFrac = opts.peakMinFrac;
        end
        if isfield(opts, 'peakProminenceFrac')
            obj.phaseCorrOpts.peakPromFrac = opts.peakProminenceFrac;
        end

        if isfield(opts, 'trackedPeak')
            obj.phaseCorrOpts.trackedPeakInd = double(opts.trackedPeak);
        elseif any(isfield(opts, {'startEdge', 'numSamples', 'peakClearance', ...
                                  'peakMinFrac', 'peakProminenceFrac'}))
            obj.phaseCorrOpts.trackedPeakInd = [];
        end
        obj.refreshAutoPeak();
    end

    function setElevModelTakeda(obj, L, d)
        % Select Takeda's phase-to-height model.
        %   In:  L - camera height above the reference plane (mm)
        %        d - camera-projector separation (mm)

        obj.elevModel.type = 'takeda';
        obj.elevModel.L = L;
        obj.elevModel.d = d;
    end

    function setElevModelPoly(obj, path)
        % Select a polynomial phase-to-height model and load calibration data from a .mat file.
        %   Input:
        %   path - .mat file written by writePolyCalibration()

        path = FtpSolver.absolutePath(path);
        pixelwise = true;
        if isfield(obj.elevModel, 'poly') && isfield(obj.elevModel.poly, 'pixelwise')
            pixelwise = obj.elevModel.poly.pixelwise;
        end

        obj.elevModel.type = 'polynomial';
        obj.elevModel.path = path;
        obj.elevModel.poly = load(path);
        obj.elevModel.poly.pixelwise = pixelwise;
    end

    function setPolyPixelwise(obj, state)
    % Choose pixelwise or smoothed polynomial coefficients.
    %   Input: state - true for the raw pixelwise fit, false for the
    %                  second-order surfaces fitted to the coefficient planes

        if strcmpi(obj.elevModel.type, 'polynomial') && isfield(obj.elevModel, 'poly')
            obj.elevModel.poly.pixelwise = logical(state);
        else
            warning("Polynomial model not found. Call setElevModelPoly(path) first.")
        end
    end

    function setStoredOutputs(obj, opts)
    % Choose which quantities solve() keeps in memory and writeCase() writes.
    %   setStoredOutputs(elevation=true, phase=false, images=false)
    %   Quantities not named keep their current setting.
    %   Elevation and phase are each written to their own HDF5 file.
    %   Images are kept in memory (imgData.stack) for debugging only:
    %   they are never written, and cannot be stored when streaming.
        arguments
            obj
            opts.elevation  (1,1) logical = obj.outputConfig.elevStorage
            opts.phase      (1,1) logical = obj.outputConfig.phaseStorage
            opts.images     (1,1) logical = obj.outputConfig.imageStorage
        end
        obj.outputConfig.elevStorage   = opts.elevation;
        obj.outputConfig.phaseStorage  = opts.phase;
        obj.outputConfig.imageStorage = opts.images;
    end

    function enableStreaming(obj, bufferLength)
    % Write results to disk during solve() instead of holding the whole run in memory.
    %   Holds at most bufferLength timesteps in memory. solve() creates the
    %   result files in <saveFolder>/<caseID> at the start of the run and
    %   writes each buffer to them at its frame offset as soon as it is
    %   full. Requires setSaveFolder() before solve(). Images cannot be
    %   stored.
        arguments
            obj
            bufferLength (1,1) double {mustBeInteger, mustBePositive}
        end
        obj.outputConfig.bufferLength = bufferLength;
        obj.outputConfig.streaming = true;
    end

    function disableStreaming(obj)
    % Hold the entire output in memory.
        obj.outputConfig.streaming = false;
    end

    function setSaveFolder(obj, folder)
    % Set the folder that results are saved to.
    %   writeCase() writes into the subfolder <folder>/<caseID>, which is
    %   created if needed. With streaming on, solve() also writes the
    %   result files there while it runs.
        arguments
            obj
            folder {mustBeTextScalar}
        end
        
        if ~isfolder(folder)
            mkdir(folder)
        end
        % store the absolute path
        listing = dir(folder);
        obj.outputConfig.saveFolder = string(listing(1).folder);
    end

    function setProgressReport(obj, state)
    % Switch the per-frame progress report of solve() on or off (default on).
        arguments
            obj
            state (1,1) logical
        end
        obj.outputConfig.reportProgressEnabled = state;
    end

    function prepForCalibration(obj)
    % Set properties in preparation for calibration.
        obj.setGeometry(resizeFactor=1, resizeFactorOutput=1);
        obj.setStoredOutputs(phase=true, elevation=false);
    end

    function meanAbsErr = calibrateTakeda(obj, heightVec, opts)
    % Fit L and d in Takeda's phase-to-height formula.
    %
    %   Fits the two geometric parameters of
    %   h = L*dPhi ./ (dPhi - 2*pi*f0*d),    f0 = 1/(period*mmPerPixel)
    %   to the calibration stack by weighted nonlinear least squares over
    %   the calibration region, within optional user-specified bounds. L is
    %   the camera height above the reference plane and d is the signed
    %   camera-projector separation, both in units of mm. With image axes
    %   defined pointing right and down, d is positive if the projector
    %   sits towards the positive image axis with respect to the camera.
    %   Returns the mean absolute error of each plane.
    %
    %   Options:
    %       excludeHeights - logical mask over heightVec, true = drop plane
    %       weights        - per-plane fit weights
    %       boundsL        - [min max] bounds for L
    %       boundsD        - [min max] bounds for d
    %
    %   Stores obj.elevModel.L and obj.elevModel.d and sets
    %   obj.elevModel.type = 'takeda'.
        arguments
            obj
            heightVec double            = []            % in units of mm
            opts.excludeHeights         = false(size(heightVec))
            opts.weights double         = ones(size(heightVec))
            opts.boundsL (1,2) double   = [0 Inf]
            opts.boundsD                = []
            opts.plot (1,1) logical     = true;
        end

        [phaseCrop, heightVec, weights, geom] = obj.prepareCalibrationData( ...
            heightVec, opts.excludeHeights, opts.weights);

        % carrier frequency on the reference plane (1/mm)
        f0 = 1 / (obj.fringe.period * obj.worldCoords.scaling.mmPerPixel);

        % take the per-plane mean of the phase
        phaseVec = squeeze(mean(phaseCrop, [1, 2], 'omitmissing'));
        phaseVec = double(phaseVec(:));

        % initial guesses from the linear inverse model
        %   1/h = 1/L - (2*pi*f0*d/L) * (1/dPhi)
        Xlin = [ones(numel(phaseVec), 1), 1./phaseVec];
        ab   = (weights.*Xlin) \ (weights.*(1./heightVec));
        L0   = 1/ab(1);
        d0   = -ab(2)*L0 / (2*pi*f0);

        if isempty(opts.boundsD)
            opts.boundsD = sort([0 sign(d0)*Inf]);
        end

        % If the linear estimate is not finite or lies outside the bounds, start
        % instead from the midpoint of the bounds. For an unbounded parameter,
        % start L at ten times the largest calibration height and d at L/4,
        % clipped to the bounds.
        if all(isfinite(opts.boundsL))
            Lfall = mean(opts.boundsL);
        else
            Lfall = min(max(10*max(abs(heightVec)), opts.boundsL(1)), opts.boundsL(2));
        end
        if all(isfinite(opts.boundsD))
            dfall = mean(opts.boundsD);
        else
            dfall = min(max(Lfall/4, opts.boundsD(1)), opts.boundsD(2));
        end
        if ~isfinite(L0) || L0 < opts.boundsL(1) || L0 > opts.boundsL(2)
            L0 = Lfall;
        end
        if ~isfinite(d0) || d0 < opts.boundsD(1) || d0 > opts.boundsD(2)
            d0 = dfall;
        end

        % Coefficient order [L, d] follows the argument order of the handle,
        % and f0 is fixed at the value it has when the handle is created.
        % fit() minimizes sum(w.*r.^2), whereas calibratePoly minimizes
        % sum((w.*r).^2), so the weights are squared to make the two agree.
        ft = fittype(@(L, d, dPhi) L*dPhi ./ (dPhi - 2*pi*f0*d), 'independent', 'dPhi');
        fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Lower', [opts.boundsL(1), opts.boundsD(1)], ...
                        'Upper', [opts.boundsL(2), opts.boundsD(2)], ...
                        'StartPoint', [L0, d0], ...
                        'Weights', weights.^2);
        takedaFit = fit(phaseVec, heightVec, ft, fo);

        L = takedaFit.L;
        d = takedaFit.d;

        % Diagnostics
        elevCrop = FtpSolver.evalTakedaElev(double(phaseCrop), L, d, f0);

        midY = floor((size(phaseCrop, 1) + 1) / 2);
        midX = floor((size(phaseCrop, 2) + 1) / 2);
        phaseSample = squeeze(double(phaseCrop(midY, midX, :)));

        fitPhase  = linspace(min(phaseSample) - 1, max(phaseSample) + 1, 1000);
        fitHeight = FtpSolver.evalTakedaElev(fitPhase, L, d, f0);

        if opts.plot
            obj.plotCalibrationDiagnostics(phaseCrop, elevCrop, fitPhase, ...
                fitHeight, heightVec, geom.calibRect, "Takeda model", {});
        end

        fprintf("Takeda calibration: L = %.5g mm, d = %.5g mm, mean absolute error = %.3g mm\n", ...
            L, d, mean(abs(elevCrop - reshape(heightVec, 1, 1, [])), 'all', 'omitmissing'));

        errorMat = elevCrop - reshape(heightVec, 1, 1, []);
        meanAbsErr = squeeze(mean(abs(errorMat), [1 2]));
        meanAbsErr = meanAbsErr(:);

        % store the model
        obj.elevModel.L = L;
        obj.elevModel.d = d;
        obj.elevModel.type = 'takeda';
    end

    function meanAbsErr = calibratePoly(obj, polyOrder, heightVec, opts)
    % Fit a per-pixel polynomial phase-to-height model.
    %
    %   For every pixel in the calibration region, fits h = a0 + a1*dPhi +
    %   a2*dPhi^2 + ... + a_polyOrder*dPhi^polyOrder (a0 fitted only when
    %   constOffset = true, zero otherwise) by weighted least squares
    %   against the calibration stack in phaseData.stack. Returns the mean
    %   absolute error of each plane when pixelwise coefficients are used.
    %
    %   Stores in obj.elevModel.poly:
    %       coeffMatOrg     - pixelwise coefficients, RAW frame, NaN-padded
    %       surfaceCoeffs   - 6 x nPlanes array of second-order surface
    %                         coefficients for each plane of coeffMatOrg
        arguments
            obj
            polyOrder (1,1) double = 2
            heightVec double = []
            opts.excludeHeights = false(size(heightVec))
            opts.weights double = ones(size(heightVec))
            opts.constOffset (1,1) logical = false
            opts.plot (1,1) logical = true;
        end

        if isfield(obj.elevModel, 'path')
            obj.elevModel = rmfield(obj.elevModel, 'path');
        end

        [phaseCrop, heightVec, weights, geom] = obj.prepareCalibrationData( ...
            heightVec, opts.excludeHeights, opts.weights);

        % Pixelwise weighted least-squares fit
        [rows, cols, nPlanes] = size(phaseCrop);
        phaseVec  = reshape(double(phaseCrop), rows*cols, nPlanes);
        weightMat = diag(weights);
        hWeighted = weightMat * heightVec;

        nCoeff    = polyOrder + opts.constOffset;
        coeffFlat = zeros(rows*cols, nCoeff);

        for iPx = 1:rows*cols
            % Design matrix: [1 (optional), dPhi, dPhi^2, ...], one row per plane
            Amat = phaseVec(iPx, :)' .^ (1:polyOrder);        % nPlanes x polyOrder
            if opts.constOffset
                Amat = [ones(nPlanes, 1), Amat];
            end
            coeffFlat(iPx, :) = (weightMat * Amat) \ hWeighted;
        end

        coeffMat = reshape(coeffFlat, rows, cols, []);

        % Keep the coefficient stack indexed by power 0..polyOrder even when
        % no constant offset was fitted: prepend an all-zero 0th-order plane.
        if ~opts.constOffset
            coeffMat = cat(3, zeros(rows, cols), coeffMat);
        end

        % Embed pixelwise coefficients into the RAW frame (NaN-padded)
        coeffMatOrg = NaN(geom.rawSize(1), geom.rawSize(2), size(coeffMat, 3));
        coeffMatOrg(geom.embedRows, geom.embedCols, :) = coeffMat;

        % Smooth (poly22) coefficients
        surfaceCoeffs   = FtpSolver.fitCoeffSurfaces(coeffMatOrg);
        smoothCoeffCrop = FtpSolver.evalCoeffSurfaces(surfaceCoeffs, ...
            geom.pMeshX_crop, geom.pMeshY_crop);

        % Diagnostics
        elevCropSmooth = FtpSolver.evalPolyElev(phaseCrop, smoothCoeffCrop);
        midY = floor((size(phaseCrop, 1) + 1) / 2);
        midX = floor((size(phaseCrop, 2) + 1) / 2);
        phaseSample = squeeze(phaseCrop(midY, midX, :));

        fitPhase = linspace(min(phaseSample) - 1, max(phaseSample) + 1, 1000);
        fitHeightSmooth = FtpSolver.evalPolyElev(fitPhase, ...
                                        smoothCoeffCrop(midY, midX, :));

        elevCropPixelwise = FtpSolver.evalPolyElev(phaseCrop, coeffMat);
        fitHeightPixelwise = FtpSolver.evalPolyElev(fitPhase, ...
                                coeffMat(midY, midX, :));

        if opts.plot
            cLims = obj.plotCalibrationDiagnostics(phaseCrop, elevCropSmooth, ...
                fitPhase, fitHeightSmooth, heightVec, ...
                geom.calibRect, "Smoothed coefficients", {});

            obj.plotCalibrationDiagnostics(phaseCrop, elevCropPixelwise, ...
                        fitPhase, fitHeightPixelwise, heightVec, ...
                        geom.calibRect, "Pixelwise coefficients", cLims);
        end

        errorMat = elevCropPixelwise - reshape(heightVec, 1, 1, []);
        meanAbsErr = squeeze(mean(abs(errorMat), [1 2]));
        meanAbsErr = meanAbsErr(:);

        % Store the model
        obj.elevModel.type = 'polynomial';
        obj.elevModel.poly.coeffMatOrg  = coeffMatOrg;
        obj.elevModel.poly.surfaceCoeffs = surfaceCoeffs;
        obj.elevModel.poly.pixelwise = true;
    end

    function writePolyCalibration(obj, path)
    % Save the polynomial calibration to a .mat file.
    %   The resized working copy (coeffMat) is dropped; only the RAW-frame
    %   coefficients and the fitted surfaces are written.
    %   Input:  path - write path of .mat calibration file

        temp = struct;
        temp.coeffMatOrg = obj.elevModel.poly.coeffMatOrg;
        temp.surfaceCoeffs = obj.elevModel.poly.surfaceCoeffs;
        save(path, '-struct', 'temp')
        obj.elevModel.path = FtpSolver.absolutePath(path);
    end
end

%% ----------- Solve, results and saving -----------
methods (Access = public)
    function solve(obj)
    % Run the full processing pipeline over solveRange.
       % Prepares the reference image, then for every frame: load,
       % preprocess, demodulate, unwrap, phase-correct, convert to
       % elevation and store. Frames flagged as invalid are stored as NaN.
       % With streaming on, the result files are created in the save
       % folder at the start, and solve() writes each buffer to the result
       % files as soon as it is full.

        t0 = tic();

        fprintf("Case %s\n", obj.caseID)

        obj.prepareRun();

        for timestep = obj.procOpts.solveRange(1):obj.procOpts.solveRange(2)
            obj.loopState.timestep = timestep;
            obj.loopState.stackInd = obj.loopState.stackInd + 1;

            if obj.outputConfig.streaming && obj.loopState.stackInd > obj.outputConfig.bufferLength
                obj.flushBuffer();
                obj.loopState.stackInd = obj.loopState.stackInd - obj.outputConfig.bufferLength;
                obj.loopState.bufferIndex = obj.loopState.bufferIndex + 1;
            end
            obj.reportProgress();

            obj.imgData.curr = obj.readImage(timestep);
            obj.loopState.discardFrame = obj.preprocessCurrentImage();

            if obj.loopState.discardFrame
                obj.storeData();     % write NaN frame
                continue
            end

            obj.demodulate();
            obj.unwrapAndFlag();
            obj.correctPhase();
            obj.calculateElevation();
            obj.storeData();

            obj.loopState.isFirstTimestep = false;
        end
        fprintf('\nDone\n')

        if obj.outputConfig.streaming
            obj.flushBuffer();
        end
        obj.diagnostics.solveTime = toc(t0);
    end

    function writeCase(obj)
    % Save results, settings and the object to <saveFolder>/<caseID>.
    %   Files written:
    %     <caseID>_elevation.h5     - surface elevation (mm), if stored
    %     <caseID>_phase.h5         - phase (rad), if stored
    %     coordinateMesh.mat        - xMesh, yMesh (mm) of the result arrays
    %     <caseID>_diagnostics.mat  - possible unwrapping faults (flagged
    %                                 regions, per frame) and frames where
    %                                 phase correction failed
    %     <caseID>_log.txt          - solver settings and run summary
    %     <caseID>.mat              - the FtpSolver object, without data
    %   With streaming on, the .h5 files are written by solve(); here
    %   they are only checked for missing frames. Read result files with
    %   FtpSolver.readData().
        folder = obj.caseFolder();
        if ~isfolder(folder)
            mkdir(folder)
        end

        temp.(obj.caseID) = obj;
        save(fullfile(folder, obj.caseID + ".mat"), '-struct', 'temp')

        % with streaming on, the result files are written during solve()
        % and are only checked here
        % otherwise the whole run is in memory and is written now
        solveRange = obj.procOpts.solveRange;
        for quantity = obj.storedQuantities()
            if obj.outputConfig.streaming
                obj.checkResultFile(quantity)
            else
                frames = obj.outputSubarray(quantity, solveRange(1), solveRange(2));
                if isempty(frames)
                    warning("No %s data in memory. Run solve() first.", quantity)
                    continue
                end
                obj.createResultFile(quantity, size(frames, 1:3));
                obj.writeFrames(quantity, frames, 1);
            end
        end

        obj.writeCoordinateMesh()
        obj.writeDiagnostics()
        obj.writeLog()
    end

    function savedObj = saveobj(obj)
        savedObj.caseID             = obj.caseID;
        savedObj.source             = obj.source;
        savedObj.imgData            = obj.defaultImgData();
        savedObj.phaseData          = obj.defaultPhaseData();
        savedObj.elevData           = obj.defaultElevData();
        savedObj.fringe             = obj.fringe;
        savedObj.phaseCorrOpts      = obj.phaseCorrOpts;
        savedObj.demodOpts          = obj.demodOpts;
        savedObj.procOpts           = obj.procOpts;
        savedObj.loopState          = obj.loopState;
        savedObj.cameraCalib        = obj.cameraCalib;
        savedObj.worldCoords        = obj.worldCoords;
        savedObj.worldCoords.mesh   = [];
        savedObj.elevModel          = obj.elevModel;
        savedObj.outputConfig       = obj.outputConfig;
        savedObj.diagnostics        = obj.diagnostics;
        if isfield(savedObj.elevModel, 'poly')
            savedObj.elevModel.poly.coeffMatOrg = [];
            savedObj.elevModel.poly.coeffMat    = [];
        end
    end

    function arr = outputSubarray(obj, quantity, firstTimestep, lastTimestep)
        % OUTPUT-frame subarray of the elevation or phase stack; [] if not stored.
        arguments
            obj
            quantity      (1,1) string {mustBeMember(quantity, ["elevation", "phase"])}
            firstTimestep (1,1) double {mustBeInteger, mustBePositive}
            lastTimestep  (1,1) double {mustBeInteger, mustBePositive}
        end

        solveRange = obj.procOpts.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = firstTimestep - solveRange(1) + 1:lastTimestep - solveRange(1) + 1;

        arr = [];
        if quantity == "elevation" && ~isempty(obj.elevData.stack)
            arr = obj.elevData.stack(rowInds, colInds, tRange);
        elseif quantity == "phase" && ~isempty(obj.phaseData.stack)
            arr = obj.phaseData.stack(rowInds, colInds, tRange);
        end
    end

    function [xMesh, yMesh] = outputCoordinates(obj)
        % Return world coordinates of the OUTPUT domain.
        [rowInds, colInds] = obj.outputInds();
        xMesh = obj.worldCoords.mesh.x(rowInds, colInds);
        yMesh = obj.worldCoords.mesh.y(rowInds, colInds);
    end

    function tseries = probe(obj, location, opts)
    % Time series of elevation or phase at a location in world coordinates.
    %   probe([x y]) returns the elevation at the grid point nearest to
    %   (x, y). With patchSize=d it combines the grid points in a square
    %   of side d, or a circle of diameter d, centered on (x, y).
    %
    %   Examples:
    %     probe([120 45])
    %     probe([120 45], patchSize=5, shape="circle", operation="max")
    %     probe([120 45], quantity="phase")
    %
    %   Inputs:
    %        location - [x y], world units (mm)
    %        Name-value options:
    %           quantity  - "elevation" (default) or "phase"
    %           patchSize - side length or diameter of the patch (mm);
    %                      0 (default) for the nearest grid point
    %           shape     - "square" (default) or "circle"
    %           operation - "mean" (default), "min" or "max" over the
    %                      patch, ignoring NaN
    %   Outputs:
    %           tseries  - (N x 1) vector, one value per stored frame
        arguments
            obj
            location       (1,2) double {mustBeFinite}
            opts.quantity  (1,1) string {mustBeMember(opts.quantity, ["elevation", "phase"])} = "elevation"
            opts.patchSize (1,1) double {mustBeNonnegative, mustBeFinite} = 0
            opts.shape     (1,1) string {mustBeMember(opts.shape, ["square", "circle"])} = "square"
            opts.operation (1,1) string {mustBeMember(opts.operation, ["mean", "min", "max"])} = "mean"
        end

        switch opts.quantity
            case "elevation"
                stack = obj.elevData.stack;
            case "phase"
                stack = obj.phaseData.stack;
        end
        if isempty(stack)
            error("No stored %s data. Call setStoredOutputs(%s=true) and run solve().", ...
                opts.quantity, opts.quantity)
        end

        meshX = obj.worldCoords.mesh.x;
        meshY = obj.worldCoords.mesh.y;
        if ~isequal(size(stack, [1 2]), size(meshX))
            error("The %s stack does not match the world-coordinate mesh.", opts.quantity)
        end

        dx = meshX - location(1);
        dy = meshY - location(2);
        % spacing between neighboring grid points
        spacing = hypot(meshX(1,2) - meshX(1,1), meshY(1,2) - meshY(1,1));
        nFrames = size(stack, 3);

        % single point: nearest grid node
        if opts.patchSize == 0
            [dist, ind] = min(hypot(dx, dy), [], 'all');
            if dist > spacing
                error("Location (%g, %g) is outside the measurement area.", ...
                    location(1), location(2))
            end
            [row, col] = ind2sub(size(meshX), ind);
            tseries = reshape(stack(row, col, :), nFrames, 1);
            return
        end

        % patch: grid points inside the square or circle around (x, y)
        halfSize = opts.patchSize / 2;
        switch opts.shape
            case "square"
                inPatch = max(abs(dx), abs(dy)) <= halfSize;
            case "circle"
                inPatch = hypot(dx, dy) <= halfSize;
        end
        if ~any(inPatch, 'all')
            error("The patch contains no grid points: it is outside the " + ...
                "measurement area or smaller than the grid spacing (%g).", spacing)
        end
        if any(inPatch([1 end], :), 'all') || any(inPatch(:, [1 end]), 'all')
            warning("The patch reaches the edge of the measurement area " + ...
                "and may be truncated.")
        end

        % one row per grid point, one column per frame
        % inPatch(:) lists the grid points in the same order, so it selects the patch rows
        vals = reshape(stack, [], nFrames);
        vals = vals(inPatch(:), :);

        switch opts.operation
            case "mean"
                tseries = mean(vals, 1, "omitnan");
            case "min"
                tseries = min(vals, [], 1, "omitnan");
            case "max"
                tseries = max(vals, [], 1, "omitnan");
        end
        tseries = tseries(:);
    end

    function refImage = getRefImage(obj)
    % Return the rectified and preprocessed reference image.
        if strlength(obj.source.refPath) > 0 && isempty(obj.imgData.refRaw)
            obj.reloadInputs();
        end
        obj.preprocessRefImage();
        refImage = obj.imgData.ref;
    end

    function imageCell = getImages(obj, firstTimestep, lastTimestep)
    % Rectify, preprocess and return the image sequence in a cell array.
        arguments
            obj
            firstTimestep (1,1) double = obj.procOpts.solveRange(1)
            lastTimestep (1,1) double = obj.procOpts.solveRange(2)
        end

        imageCell = cell(1, lastTimestep - firstTimestep + 1);

        for timestep = firstTimestep:lastTimestep
            obj.imgData.curr = obj.readImage(timestep);
            obj.preprocessCurrentImage();
            imageCell{timestep - firstTimestep + 1} = obj.imgData.curr;
        end
    end

    function imageCell = getRectifiedImages(obj, firstTimestep, lastTimestep)
    % Rectify and return the image sequence in a cell array without preprocessing.
        arguments
            obj
            firstTimestep (1,1) double = obj.procOpts.solveRange(1)
            lastTimestep (1,1) double = obj.procOpts.solveRange(2)
        end

        imageCell = cell(1, lastTimestep - firstTimestep + 1);

        for timestep = firstTimestep:lastTimestep
            img = obj.readImage(timestep);
            img = obj.rectify(img);
            imageCell{timestep - firstTimestep + 1} = img;
        end
    end
end

methods (Static)
    function [data, fullSize] = readData(path, firstFrame, lastFrame)
    % Read elevation or phase frames from a result file written by writeCase.
    %   data = FtpSolver.readData(path) reads every frame.
    %   data = FtpSolver.readData(path, firstFrame, lastFrame) reads
    %   frames firstFrame to lastFrame, counted from 1 at the start of the
    %   file. The file's firstTimestep attribute gives the dataset frame
    %   index of frame 1. Called without a path, it asks for the file.
    %   Outputs:
    %       data     - [rows x cols x frames] single array
    %       fullSize - size of the whole array in the file
        arguments
            path        {mustBeTextScalar} = ""
            firstFrame  (1,1) double {mustBeInteger, mustBePositive} = 1
            lastFrame   double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
        end

        if strlength(path) == 0
            [file, folder] = uigetfile('*.h5');
            if isequal(file, 0)
                data = [];
                fullSize = [];
                return
            end
            path = fullfile(folder, file);
        end

        info = h5info(path);
        if numel(info.Datasets) ~= 1
            error("readData: %s is not an FtpSolver result file " + ...
                "(expected one dataset, found %d).", path, numel(info.Datasets))
        end
        dataset = "/" + info.Datasets.Name;
        fullSize = info.Datasets.Dataspace.Size;
        fullSize(end + 1:3) = 1;

        if isempty(lastFrame)
            lastFrame = fullSize(3);
        end
        if firstFrame > lastFrame || lastFrame > fullSize(3)
            error("readData: requested frames %d to %d, but the file " + ...
                "holds frames 1 to %d.", firstFrame, lastFrame, fullSize(3))
        end

        % frames after framesWritten were never filled (interrupted run)
        attrNames = string({info.Datasets.Attributes.Name});
        if any(attrNames == "framesWritten")
            framesWritten = h5readatt(path, dataset, "framesWritten");
            if lastFrame > framesWritten
                warning("readData: %s holds only %d written frames; " + ...
                    "later frames read as 0.", path, framesWritten)
            end
        end

        data = h5read(path, dataset, [1 1 firstFrame], ...
            [fullSize(1:2), lastFrame - firstFrame + 1]);
    end

    function newObj = loadobj(fileObj)
        if isstruct(fileObj)
            newObj = FtpSolver("temp");
            newObj.caseID           = fileObj.caseID;
            newObj.source           = fileObj.source;
            newObj.fringe           = fileObj.fringe;
            newObj.phaseCorrOpts    = fileObj.phaseCorrOpts;
            newObj.demodOpts        = fileObj.demodOpts;
            newObj.procOpts         = fileObj.procOpts;
            newObj.loopState        = fileObj.loopState;
            newObj.cameraCalib      = fileObj.cameraCalib;
            newObj.worldCoords      = fileObj.worldCoords;
            newObj.elevModel        = fileObj.elevModel;
            newObj.outputConfig     = fileObj.outputConfig;
            newObj.diagnostics      = fileObj.diagnostics;
        else
            newObj = fileObj;
        end
    end
end

%% ----------- Visualization -----------
methods (Access = public)
    function showDomains(obj)
    % Draw the domain rectangles on the reference image.
        figure("WindowStyle","normal", 'Units', 'normalized', ...
               'Position', [0.1 0.1 0.5 0.6]);
        imHndl = imshow(obj.imgData.refRaw, []);
        obj.overlayDomainRectangles(imHndl.Parent)
    end

    function showTrackedPeaks(obj, firstTimestep, lastTimestep, frameSkip, lineNum)
    % Diagnostic plot of phase-correction peak detection.
    %  Shows the reference line with its detected peaks and the tracked
    %  peak highlighted, then the same line in the requested images.
        arguments
            obj
            firstTimestep (1,1) double {mustBeInteger, mustBePositive}
            lastTimestep  (1,1) double {mustBeInteger, mustBePositive}
            frameSkip     (1,1) double {mustBeInteger, mustBePositive} = 1
            lineNum       (1,1) double {mustBeInteger, mustBePositive} ...
                                = ceil(length(obj.phaseCorrOpts.lineInds)/2)
        end
        
        if lineNum > length(obj.phaseCorrOpts.lineInds)
            error("Selected line number exceeds the number of lines used " ...
                   + "by the spatial phase correction algorithm. Call " ...
                   + "setPhaseCorrection(numSamples=N) to change the " ...
                   + "number of tracked lines.")
        end

        lineInd = obj.phaseCorrOpts.lineInds(lineNum);
        [optInd, ~, refPeakLocsCell] = obj.findOptimalPeakInd(obj.phaseCorrOpts.trackedPeakInd);
        peakInd = obj.phaseCorrOpts.trackedPeakInd;

        if isempty(peakInd)
            if isnan(optInd)
                error("Cannot determine a peak to highlight: phaseCorrOpts.trackedPeakInd is " + ...
                    "empty and no eligible peak was found. Check that the reference " + ...
                    "image and fringe period are set.")
            end
            peakInd = optInd;
            fprintf("phaseCorrOpts.trackedPeakInd is empty; showing the automatic choice (%g).\n", optInd)
        end

        refPeakLocs = refPeakLocsCell{lineNum};
        if isempty(refPeakLocs) || peakInd > length(refPeakLocs)
            error("Peak %g was not detected on phase-correction line %g of the reference image.", ...
                peakInd, lineInd)
        end

        if strcmpi(obj.fringe.normalAxis, 'X')
            refLine = double(obj.imgData.refRaw(lineInd, :));
        else
            refLine = double(obj.imgData.refRaw(:, lineInd));
        end
        refPeakVals = refLine(refPeakLocs);

        % detect the same line in the requested images
        targetImages = firstTimestep:frameSkip:lastTimestep;
        N = length(targetImages);

        currPeaksCell = cell(1, N);
        currLineCell = cell(1, N);

        for n = 1:N
            timestep = targetImages(n);
            img = obj.readImage(timestep);
            img = obj.rectify(img);

            [peakVals, peakLocs] = detectFringePeaks(obj, img, lineInd);

            if strcmpi(obj.phaseCorrOpts.startEdge, 'right') || strcmpi(obj.phaseCorrOpts.startEdge, 'bottom')
                peakLocs = flip(peakLocs);
                peakVals = flip(peakVals);
            end

            currPeaksCell{n} = [peakLocs(:), peakVals(:)];

            if strcmpi(obj.fringe.normalAxis, 'X')
                currLine = double(img(lineInd, :));
            else
                currLine = double(img(:, lineInd));
            end

            currLineCell{n} = currLine;
        end

        % plot
        figure;
        tiledlayout(3, ceil(N/3) + 1, 'TileSpacing','compact', 'Padding', 'compact')
        nexttile([3 1])
        plot(refLine)
        hold on
        scatter(refPeakLocs, 1.01*refPeakVals, 'v', 'filled')
        scatter(refPeakLocs(peakInd), 1.01*refPeakVals(peakInd), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
        hold off
        xlabel(obj.fringe.normalAxis + " (pixels)")
        ylabel("Counts")
        if strcmpi(obj.phaseCorrOpts.startEdge, 'left') || strcmpi(obj.phaseCorrOpts.startEdge, 'top')
            xLimit = [1 round(1.2*refPeakLocs(peakInd))];
        else
            xLimit = [round(0.9*refPeakLocs(peakInd)) length(refLine)];
        end
        yLimitRef = [0.8*min(refLine(xLimit(1):xLimit(2))) 1.2*max(refLine(xLimit(1):xLimit(2)))];

        xlim(xLimit)
        ylim(yLimitRef)
        title("Reference: x_p = " + num2str(refPeakLocs(peakInd)) )
        set(gca, 'fontsize', 12)

        for i = 1:N
            if i == 1                
                yLimitData = [0.8*min(refLine(xLimit(1):xLimit(2))) 1.2*max(currLineCell{1}(xLimit(1):xLimit(2)))];
            end
            nexttile
            plot(currLineCell{i})
            hold on
            scatter(currPeaksCell{i}(:,1), 1.05*currPeaksCell{i}(:,2), 'v', 'filled')
            if size(currPeaksCell{i}, 1) >= peakInd
                scatter(currPeaksCell{i}(peakInd, 1), 1.05*currPeaksCell{i}(peakInd,2), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
                title("Frame " + num2str(targetImages(i)) + ": x_p = " + num2str(currPeaksCell{i}(peakInd, 1)) )
            else
                title("Frame " + num2str(targetImages(i)) + ": peak " + num2str(peakInd) + " not detected")
            end
            hold off
            xlim(xLimit)
            ylim(yLimitData)
        end
    end

    function animate(obj, firstTimestep, lastTimestep, framerate, opts)
    % Animate the reconstructed surface or phase.
    %   Plots the surface elevation or phase (or an array passed in opts.data) frame
    %   by frame, with optional plane/mean/constant subtraction, NaN-frame
    %   interpolation, smoothing and a time or frame annotation.
    %
    %   Example: animate(1, 200, 24, quantity="phase")
    %   In:  firstTimestep, lastTimestep - frame range ([] = to the end)
    %        framerate          - playback rate (fps)
    %        opts               - name-value display options (data, ZLim,
    %                             view, colormap, camFPS, subtractPlane,
    %                             subtractMean, smoothing, ...)

        arguments
            obj
            firstTimestep                         = [];
            lastTimestep                          = [];
            framerate               (1,1) double  = 30;
            opts.quantity           (1,1) string ...
                {mustBeMember(opts.quantity, ["elevation", "phase"])} ...
                                                  = "elevation";
            opts.camFPS             (1,1) double  = 1;
            opts.data                     cell    = {};
            opts.showTime           (1,1) logical = true;
            opts.subtractPlane      (1,1) logical = false;
            opts.subtractMean       (1,1) logical = false;
            opts.subtractConst                    = [];
            opts.smoothing          (1,1) logical = false;
            opts.smoothingSigma     (1,1) double  = 2;
            opts.annotStartTime     (1,1) double  = 0;
            opts.fillNanFrames      (1,1) logical = false;
            opts.ZLim               (1,2) double  = [0 0];
            opts.lightPosition      (1,3) double  = [1 1 5];
            opts.FontSize           (1,1) double  = 16;
            opts.DataAspectRatio    (1,3) double  = [2 2 1];
            opts.view               (1,2) double  = [23.2 23.8]
            opts.frameStep          (1,1) double  = 1;
            opts.figPosition        (1,4) double  = [0.05 0.13 0.6 0.6];
            opts.timeFormat               string  = "";
            opts.colormap                         = [];
        end

        [dataArray, xMesh, yMesh, meshUnits] = ...
                    animationData(obj, firstTimestep, lastTimestep, opts);

        pauseTime = 1/framerate;

        if isequal(opts.ZLim, [0 0])
            zMin = prctile(dataArray(1:4:end,1:4:end,:), 0.01, 'all');
            zMax = prctile(dataArray(1:4:end,1:4:end,:), 99.99, 'all');
        else
            zMin = opts.ZLim(1);
            zMax = opts.ZLim(2);
        end

        midHeight = (zMin + zMax)/2;
        amp = (zMax - zMin)/2;

        figHndl = figure('Units', 'normalized', 'WindowStyle', 'normal', ...
                        'Position', opts.figPosition);
        axHndl = axes();

        surfHndl = surfl(axHndl, xMesh, yMesh, dataArray(:,:,1), 'light', 'EdgeColor', 'none');
        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % surfHndl(1).AmbientStrength = 0.4;
        surfHndl(2).Style = 'infinite';
        surfHndl(2).Position = opts.lightPosition;

        xlim([min(xMesh, [], 'all'), max(xMesh, [], 'all')])
        ylim([min(yMesh, [], 'all'), max(yMesh, [], 'all')])
        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        cbHndl = colorbar;

        if strcmpi(opts.quantity, "elevation")
            cbHndl.Label.String = 'Elevation (mm)';
        elseif strcmpi(opts.quantity, 'phase')
            cbHndl.Label.String = '\Delta\phi (rad)';
        end
        clim(axHndl, [midHeight - 0.95*amp, midHeight + 0.95*amp])

        xlabel(axHndl, "x " + meshUnits)
        ylabel(axHndl, "y " + meshUnits)
        if strcmpi(opts.quantity, "elevation")
            zlabel(axHndl, 'z (mm)')
        elseif strcmpi(opts.quantity, "phase")
            zlabel(axHndl, '\Delta\phi (rad)')
        end

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(figHndl, 'Position', opts.figPosition)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        cbHndl.Position = [0.917,0.11,0.02,0.815];

        if opts.showTime
            annotHndl = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.timeFormat) ~= 0
                timeFormat = opts.timeFormat;
            else
                if opts.camFPS == 1
                    timeFormat = "Frame = %05g";
                else
                    timeFormat = "t = %.3f s";
                end
            end
            timeString = sprintf(timeFormat, opts.annotStartTime);
            set(annotHndl, 'String', timeString)
        end

        set(gca, 'DataAspectRatio', opts.DataAspectRatio);

        pause(pauseTime)

        for i = 2:size(dataArray, 3)
            if ~ishghandle(figHndl)
                break
            end
            surfHndl(1).ZData = dataArray(:,:,i);

            if opts.showTime
                if opts.camFPS == 1
                    timeString = sprintf(timeFormat, opts.annotStartTime + (i - 1)*opts.frameStep);
                else
                    timeString = sprintf(timeFormat, opts.annotStartTime + (i - 1)*opts.frameStep/opts.camFPS);
                end
                set(annotHndl, 'String', timeString)
            end

            pause(pauseTime)
        end
    end

    function exportVideo(obj, firstTimestep, lastTimestep, framerate, opts)
    % Animate and export a video of the reconstructed surface or phase.
    % Plots the surface elevation or phase (or an array passed in opts.data) frame
    % by frame, with optional plane/mean/constant subtraction, NaN-frame
    % interpolation, smoothing and a time or frame annotation. Exports to
    % .mp4 format frame by frame if using the opts.fast option, or exports
    % each frame to a .png image in a temporary folder for higher quality
    % if opts.fast is off (default) and animates the images into a video at
    % the end.
    %
    % Example: exportVideo(1, 200, 24, quantity="phase")
    %   Inputs:
    %        firstTimestep, lastTimestep - frame range ([] = to the end)
    %        framerate          - video framerate (fps)
    %        opts               - name-value display options (data, ZLim,
    %                             view, colormap, camFPS, subtractPlane,
    %                             subtractMean, smoothing, ...)
        arguments
            obj
            firstTimestep                         = [];
            lastTimestep                          = [];
            framerate               (1,1) double  = 30;
            opts.quantity           (1,1) string ...
                {mustBeMember(opts.quantity, ["elevation", "phase"])} ...
                                                  = "elevation";
            opts.camFPS             (1,1) double  = 1;
            opts.data                     cell    = {};
            opts.saveAddr                 string  = [];
            opts.overwriteFile      (1,1) logical = false;
            opts.fast               (1,1) logical = false;
            opts.showTime           (1,1) logical = true;
            opts.subtractPlane      (1,1) logical = false;
            opts.subtractMean       (1,1) logical = false;
            opts.subtractConst                    = [];
            opts.smoothing          (1,1) logical = false;
            opts.smoothingSigma     (1,1) double  = 2;
            opts.annotStartTime     (1,1) double  = 0;
            opts.fillNanFrames      (1,1) logical = false;
            opts.ZLim               (1,2) double  = [0 0];
            opts.lightPosition      (1,3) double  = [1 1 5];
            opts.FontSize           (1,1) double  = 16;
            opts.DataAspectRatio    (1,3) double  = [2 2 1];
            opts.view               (1,2) double  = [23.2 23.8]
            opts.frameStep          (1,1) double  = 1;
            opts.imageResolution    (1,1) double  = 150;
            opts.figPosition        (1,4) double  = [0.05 0.13 0.6 0.6];
            opts.timeFormat               string  = "";
            opts.colormap                         = [];
            opts.cleanup            (1,1) logical = false;
        end    

        if isempty(opts.saveAddr)
            opts.saveAddr = uigetdir();
            if opts.saveAddr == 0
                return
            end
        end

        if ~isfolder(opts.saveAddr)
            mkdir(opts.saveAddr)
        end

        [dataArray, xMesh, yMesh, meshUnits] = animationData(obj, firstTimestep, ...
                                        lastTimestep, opts);

        videoFilename = fullfile(opts.saveAddr, (obj.caseID + ".mp4"));

        if ~opts.overwriteFile
            videoSuffix = 2;
            while isfile(videoFilename)
                videoFilename = sprintf('%s_%g.mp4', fullfile(opts.saveAddr, obj.caseID), videoSuffix);
                videoSuffix = videoSuffix + 1;
            end
        end

        v = VideoWriter(videoFilename, 'MPEG-4');
        v.FrameRate = framerate;
        v.Quality = 100;
        open(v)

        figHndl = figure('Units', 'normalized', 'WindowStyle', 'normal', ...
                      'Position', opts.figPosition);
        axHndl = axes;

        if ~opts.fast
            mkdir(fullfile(opts.saveAddr, "tempImages_" + obj.caseID));
        end

        if isequal(opts.ZLim, [0 0])
            zMin = prctile(dataArray(1:4:end,1:4:end,:), 0.01, 'all');
            zMax = prctile(dataArray(1:4:end,1:4:end,:), 99.99, 'all');
        else
            zMin = opts.ZLim(1);
            zMax = opts.ZLim(2);
        end

        midHeight = (zMin + zMax)/2;
        amp = (zMax - zMin)/2;

        surfHndl = surfl(axHndl, xMesh, yMesh, dataArray(:,:,1), 'light', ...
                    'EdgeColor', 'none');

        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % surfHndl(1).AmbientStrength = 0.4;
        surfHndl(2).Style = 'infinite';
        surfHndl(2).Position = opts.lightPosition;

        xlim([min(xMesh, [], 'all'), max(xMesh, [], 'all')])
        ylim([min(yMesh, [], 'all'), max(yMesh, [], 'all')])

        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        cbHndl = colorbar;
        if strcmpi(opts.quantity, "elevation")
            cbHndl.Label.String = 'Elevation (mm)';
        elseif strcmpi(opts.quantity, "phase")
            cbHndl.Label.String = '\Delta\phi (rad)';
        end
        clim(axHndl, [midHeight - 0.95*amp, midHeight + 0.95*amp])

        xlabel(axHndl, "x " + meshUnits)
        ylabel(axHndl, "y " + meshUnits)
        if strcmpi(opts.quantity, "elevation")
            zlabel(axHndl, 'z (mm)')
        elseif strcmpi(opts.quantity, "phase")
            zlabel(axHndl, '\Delta\phi (rad)')
        end

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        cbHndl.Position = [0.917,0.11,0.02,0.815];

        if opts.showTime
            annotHndl = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.timeFormat) ~= 0
                timeFormat = opts.timeFormat;
            else
                if opts.camFPS == 1
                    timeFormat = "Frame = %05g";
                else
                    timeFormat = "t = %.3f s";
                end
            end
            timeString = sprintf(timeFormat, opts.annotStartTime);
            set(annotHndl, 'String', timeString)
        end

        set(gca, 'DataAspectRatio', opts.DataAspectRatio);

        if opts.fast
            frame = getframe(figHndl);
            writeVideo(v, frame)
        else
            exportgraphics(figHndl, ...
                fullfile(opts.saveAddr, "tempImages_" + obj.caseID, ...
                           sprintf("tempImage_%05g.png", 1)), ...
                           'Resolution', opts.imageResolution)
        end

        for i = 2:size(dataArray, 3)
            if ~ishghandle(figHndl)
                break
            end
            surfHndl(1).ZData = dataArray(:,:,i);

            if opts.showTime
                if opts.camFPS == 1
                    timeString = sprintf(timeFormat, opts.annotStartTime + (i - 1)*opts.frameStep);
                else
                    timeString = sprintf(timeFormat, opts.annotStartTime + (i - 1)*opts.frameStep/opts.camFPS);
                end
                set(annotHndl, 'String', timeString)
            end

            if opts.fast
                frame = getframe(figHndl);
                writeVideo(v, frame)
            else
                exportgraphics(figHndl, ...
                    fullfile(opts.saveAddr, "tempImages_" + obj.caseID, ...
                               sprintf("tempImage_%05g.png", i)), ...
                               'Resolution', opts.imageResolution)
            end
        end

        % gather images and create a video
        if ~opts.fast
            filenames = dir(fullfile(opts.saveAddr, "tempImages_" + obj.caseID, "*.png"));

            for i = 1:length(filenames)
                frame = imread(fullfile(filenames(i).folder, filenames(i).name));
                if size(frame, 2) > 1920
                    frame = imresize(frame, 1920/size(frame,2));
                end
                writeVideo(v, frame)
            end
            close(v)

            if opts.cleanup
                rmdir(fullfile(opts.saveAddr, "tempImages_" + obj.caseID), 's')
            end
        else
            close(v)
        end
    end

    function animateImages(obj, firstTimestep, lastTimestep, framerate)
    % Preprocess then animate the image sequence.
        arguments
            obj
            firstTimestep (1,1) double = obj.procOpts.solveRange(1)
            lastTimestep (1,1) double = obj.procOpts.solveRange(2)
            framerate (1,1) double = 5
        end
        figHndl = figure();
        axHndl = axes(figHndl);

        pauseTime = 1/framerate;

        obj.imgData.curr = obj.readImage(firstTimestep);
        obj.preprocessCurrentImage();

        [values, edges] = histcounts(obj.imgData.curr, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));
        minThresh = edges(find(values > 0.001, 1, 'first'));

        imHndl = imshow(obj.imgData.curr, [minThresh maxThresh], 'Parent', axHndl);
        title(num2str(firstTimestep))
        pause(pauseTime)

        for timestep = firstTimestep + 1:lastTimestep
            if ~ishghandle(figHndl)
                break
            end

            obj.imgData.curr = obj.readImage(timestep);
            obj.preprocessCurrentImage();

            set(imHndl, 'CData', obj.imgData.curr)
            title(axHndl, num2str(timestep))
            pause(pauseTime)
        end
    end

    function animateRectifiedImages(obj, firstTimestep, lastTimestep, framerate, opts)
    % Rectify (if possible) the image sequence and then animate without additional preprocessing.
        arguments
            obj
            firstTimestep (1,1) double = obj.procOpts.solveRange(1)
            lastTimestep (1,1) double = obj.procOpts.solveRange(2)
            framerate (1,1) double = 5
            opts.showDomains (1,1) logical = true
        end
        figHndl = figure();
        axHndl = axes(figHndl);

        pauseTime = 1/framerate;

        img = obj.readImage(firstTimestep);

        [values, edges] = histcounts(img, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));

        for timestep = firstTimestep:lastTimestep
            if ~ishghandle(figHndl)
                break
            end

            img = obj.readImage(timestep);
            img = obj.rectify(img);

            if timestep == firstTimestep
                imHndl = imshow(img, [0 maxThresh], 'Parent', axHndl);
                if opts.showDomains
                    obj.overlayDomainRectangles(imHndl.Parent)
                end
            else
                set(imHndl, 'CData', img)
            end

            title(axHndl, num2str(timestep))
            pause(pauseTime)
        end
    end
end

%% ----------- Internal: solve() pipeline -----------
methods (Access = private)
    function prepareRun(obj)
    % Reset loop and and prepare everything solve() needs before its loop.
    %  Validates the settings, resets the loop state and diagnostics, and
    %  refreshes the geometry (reloading the images and camera mesh of a
    %  case loaded from a .mat file). Then prepares, peak-tracks and
    %  demodulates the reference image, resamples the polynomial
    %  calibration onto the computational grid, sets up the parallax
    %  geometry, and allocates the result stacks. If streaming is turned
    %  on, the results files are created.

        if obj.outputConfig.elevStorage && strcmpi(obj.elevModel.type, 'none')
            error("Elevation output is on but no elevation model is set. " + ...
                   "Call setElevModelTakeda or setElevModelPoly, or " + ...
                   "setStoredOutputs(elevation=false).")
        end

        if isempty(obj.fringe.periodOrg) || isnan(obj.fringe.periodOrg)
            error("Missing data: fringe.period")
        end

        if strcmpi(obj.demodOpts.method, 'wavelet') ...
                && (isempty(obj.demodOpts.waveletWavelengths) ...
                    || isempty(obj.demodOpts.waveletDivisions))
            error("The wavelet method needs a range of wavelet scales. Call " + ...
                  "setDemodulation(waveletWavelengths=[min max], waveletDivisions=n), " + ...
                  "with the wavelengths in RAW-frame pixels.")
        end

        if ~obj.outputConfig.imageStorage ...
                && ~obj.outputConfig.phaseStorage ...
                && ~obj.outputConfig.elevStorage ...
                && (diff(obj.procOpts.solveRange) + 1 > 1)
            error("All storage switches are off. Call setStoredOutputs().")
        end

        if obj.outputConfig.streaming
            if strlength(obj.outputConfig.saveFolder) == 0
                error("Streaming writes the results during solve() and " + ...
                    "needs a save folder. Call setSaveFolder(folder) first.")
            end
            if obj.outputConfig.imageStorage
                error("Images cannot be stored when streaming. Call " + ...
                    "setStoredOutputs(images=false) or disableStreaming().")
            end
        end

        obj.loopState.stackInd = 0;
        obj.loopState.bufferIndex = 1;
        obj.loopState.isFirstTimestep = true;
        obj.elevData.curr = [];
        obj.phaseData.curr = [];
        obj.diagnostics.flaggedRegions = struct('timestep', {}, 'rects', {});
        obj.diagnostics.phaseCorrectionFailures = struct('timestep', {}, 'candidates', {});
        obj.diagnostics.solveTime = [];

        rf = obj.procOpts.resizeFactor;
        rf_o = obj.procOpts.resizeFactorOutput;
        Nt = diff(obj.procOpts.solveRange) + 1;
        bufferLength = obj.outputConfig.bufferLength;

        % refresh coordinate mesh if running a case loaded from a .mat file
        if ~isfield(obj.worldCoords.mesh, 'xOrg')
            if strlength(obj.cameraCalib.path) > 0 && strcmpi(obj.cameraCalib.type, 'pinhole')
                obj.readCameraCalibration(obj.cameraCalib.path, obj.cameraCalib.camNum);
            end
        end
        % refresh images if running a case loaded from a .mat file
        if strlength(obj.source.refPath) > 0 && isempty(obj.imgData.refRaw)
            obj.reloadInputs();
        end
        obj.updateGeometry();

        obj.preprocessRefImage();
        obj.locateRefPeaks();
        obj.demodulateRef();

        % initialize polynomial calibration model
        [Ny, Nx] = size(obj.imgData.refRaw);
        if obj.outputConfig.elevStorage && strcmpi(obj.elevModel.type, 'polynomial')
            if obj.elevModel.poly.pixelwise
                coeffSrc = obj.elevModel.poly.coeffMatOrg;
            else
                [xMesh, yMesh] = meshgrid(1:Nx, 1:Ny);
                coeffSrc = FtpSolver.evalCoeffSurfaces( ...
                    obj.elevModel.poly.surfaceCoeffs, xMesh, yMesh);
            end
            coeffMat = FtpSolver.shrinkArr(coeffSrc, rf, 'bilinear', Antialiasing=false);


            compRect3 = [obj.procOpts.compRect(1:2), 1, ...
                obj.procOpts.compRect(3:4), ...
                size(obj.elevModel.poly.coeffMatOrg, 3) - 1 ...
                ];
            coeffMat = imcrop3(coeffMat, compRect3);
            coeffMat = FtpSolver.shrinkArr(coeffMat, rf_o, 'bilinear', Antialiasing=false);
            obj.elevModel.poly.coeffMat = coeffMat;
        end

        % ray geometry of the 2-D parallax error correction
        obj.worldCoords.mesh.parallax = [];
        if obj.procOpts.parallaxCorrection && obj.outputConfig.elevStorage ...
                && isfield(obj.cameraCalib, 'intrinsicsMatlab') ...
                && min(size(obj.imgData.ref)) >= 2
            obj.prepareParallaxGeometry();
        end

        % preallocate stacks
        if obj.outputConfig.streaming && Nt > bufferLength
            obj.allocateStacks(bufferLength)
        else
            obj.allocateStacks(Nt)
        end

        % when streaming, the result files are created now and filled one
        % buffer at a time during the run
        if obj.outputConfig.streaming
            [rowInds, colInds] = obj.outputInds();
            for quantity = obj.storedQuantities()
                obj.createResultFile(quantity, [numel(rowInds), numel(colInds), Nt]);
            end
        end
    end

    function allocateStacks(obj, nTimesteps)
    % Preallocate the result stacks for nTimesteps.
    %   Allocates a NaN-filled single-precision stack for every quantity
    %   whose output switch is on and clears the others.

        [NyImg, NxImg] = size(obj.imgData.ref);
        rf_o = obj.procOpts.resizeFactorOutput;
        stackSize = FtpSolver.shrunkArrSize([NyImg, NxImg], rf_o);
        NyStack = stackSize(1);
        NxStack = stackSize(2);

        obj.phaseData.stack = [];
        obj.elevData.stack  = [];
        obj.imgData.stack   = [];

        if obj.outputConfig.phaseStorage
            obj.phaseData.stack = nan(NyStack, NxStack, nTimesteps, 'single');
        end
        if obj.outputConfig.elevStorage
            obj.elevData.stack = nan(NyStack, NxStack, nTimesteps, 'single');
        end
        if obj.outputConfig.imageStorage
            obj.imgData.stack = nan(NyImg, NxImg, nTimesteps, 'single');
        end
    end

    function preprocessRefImage(obj)
    % Prepare the reference image.
    %  Crops the scaled reference to the computational domain, removes the
    %  DC component and (optionally) applies a smoothing filter along the
    %  fringe.

        obj.imgData.ref = imcrop(obj.imgData.refScaled, obj.procOpts.compRect);

        if obj.procOpts.dcSubtraction
            obj.imgData.ref = obj.subtractDC(obj.imgData.ref, ...
                obj.fringe.period);
        end

        if obj.procOpts.refImgSmoothing
            if strcmpi(obj.fringe.normalAxis, 'Y')
                obj.imgData.ref = imgaussfilt(obj.imgData.ref, [1e-06 5]);
            else
                obj.imgData.ref = imgaussfilt(obj.imgData.ref, [5 1e-06]);
            end
        end
    end

    function locateRefPeaks(obj)
    % Prepare the spatial phase correction from the reference image.
    %  Selects or validates phaseCorrOpts.trackedPeakInd and measures the
    %  local fringe wavelength and reference peak position on every
    %  correction line.
        nLines = length(obj.phaseCorrOpts.lineInds);
        obj.phaseCorrOpts.refPeakLocs = nan(1, nLines);
        obj.phaseCorrOpts.localPeriod = nan(1, nLines);

        if strcmpi(obj.phaseCorrOpts.method, 'spatial')
            [optInd, isSafeLine, peakLocsCell] = ...
                obj.findOptimalPeakInd(obj.phaseCorrOpts.trackedPeakInd);

            if isempty(obj.phaseCorrOpts.trackedPeakInd)
                if isnan(optInd)
                    error("Spatial phase correction: no peak lies inside the output " + ...
                        "domain with a peak clearance of %.3g fringe periods from the " + ...
                        "'%s' edge on any phase-correction line. Enlarge the output " + ...
                        "domain, or lower peakClearance or adjust peakMinFrac and " + ...
                        "peakProminenceFrac with setPhaseCorrection().", ...
                        obj.phaseCorrOpts.peakClearance, obj.phaseCorrOpts.startEdge)
                end
                obj.phaseCorrOpts.trackedPeakInd = optInd;

            elseif any(~isSafeLine)
                lineStr = strjoin(string(obj.phaseCorrOpts.lineInds(~isSafeLine)), ', ');
                if isnan(optInd)
                    warning("Peak %g is missing or does not lie safely inside the output domain " + ...
                        "on phase-correction line(s) %s, and no peak satisfies the peak clearance " + ...
                        "of %.3g fringe periods from the '%s' edge. Enlarge the output domain, " + ...
                        "or lower peakClearance or adjust peakMinFrac and peakProminenceFrac " + ...
                        "with setPhaseCorrection().", ...
                        obj.phaseCorrOpts.trackedPeakInd, lineStr, obj.phaseCorrOpts.peakClearance, obj.phaseCorrOpts.startEdge)
                else
                    warning("Peak %g is missing or does not lie safely inside the output domain " + ...
                        "on phase-correction line(s) %s. Optimal peak number: %g (closest usable " + ...
                        "peak to the '%s' edge with a peak clearance of %.3g fringe periods). " + ...
                        "Call setPhaseCorrection(trackedPeak=%g), or setPhaseCorrection(trackedPeak=[]) " + ...
                        "to select the peak automatically.", ...
                        obj.phaseCorrOpts.trackedPeakInd, lineStr, optInd, obj.phaseCorrOpts.startEdge, ...
                        obj.phaseCorrOpts.peakClearance, optInd)
                end
            end

            % local wavelength and reference peak position on each line
            peakInd = obj.phaseCorrOpts.trackedPeakInd;
            for i = 1:nLines
                peaksLine = peakLocsCell{i};

                % the wavelength stencil needs a peak on each side of
                % peakInd; leave NaN (skipped by correctPhase) when the
                % required peaks do not exist
                if peakInd < 2 || peakInd > length(peaksLine) - 1
                    continue
                end

                obj.phaseCorrOpts.localPeriod(i) = abs( peaksLine(peakInd + 1) ...
                    - peaksLine(peakInd - 1) ) / 2;
                obj.phaseCorrOpts.refPeakLocs(i) = peaksLine(peakInd);

                if abs(obj.phaseCorrOpts.localPeriod(i) - obj.fringe.periodOrg)/obj.fringe.periodOrg > 0.5
                    warning("Relative difference between local period" + ...
                        " calculated for spatial phase correction and fringe period is greater than 50%.")
                end
            end
        end
    end

    function demodulateRef(obj)
        if strcmpi(obj.demodOpts.method, 'fourier')
            obj.demodulateRefFourier();
        elseif strcmpi(obj.demodOpts.method, 'wavelet')
            obj.demodulateRefWavelet();
        else
            error("Invalid demodulation method")
        end
    end

    function demodulateRefFourier(obj)
    % Precompute all reference-dependent quantities for FT demodulation.
    %
    % Sets:
    %   obj.phaseData.bandpassFilter   - super-Gaussian band-pass centered on carrier
    %   obj.phaseData.refComplexCoeffs - band-passed complex reference signal,
    %                                    hI0 = ifft2(filter .* fft2(I0))
    %   obj.phaseData.fourierPadding   - padding around the image in
    %                                    [pixelsY, pixelsX] format
        
        [nyRef, nxRef] = size(obj.imgData.ref);

        alongFringePad = ceil(2/obj.demodOpts.filterWidthFrac * obj.fringe.period);
        acrossFringePad = ceil(1/obj.demodOpts.filterWidthFrac * obj.fringe.period);

        if strcmpi(obj.fringe.normalAxis, 'Y')
            obj.phaseData.fourierPadding = ...
                [acrossFringePad * (nyRef > 1), alongFringePad * (nxRef > 1)];
        else
            obj.phaseData.fourierPadding = ...
                [alongFringePad * (nyRef > 1), acrossFringePad * (nxRef > 1)];
        end
    
        refPadded = obj.padImage(obj.imgData.ref);
        [ny, nx] = size(refPadded);

        kxv = 2*pi/nx * [0:floor(nx/2), -ceil(nx/2)+1:-1];
        kyv = 2*pi/ny * [0:floor(ny/2), -ceil(ny/2)+1:-1];
        [kx, ky] = meshgrid(kxv, kyv);

        % carrier from the detected fringe properties
        nVec = obj.fringe.normalVec;
        omega = 2*pi/obj.fringe.period;
        kxG = omega * nVec(1);
        kyG = omega * nVec(2);

        kr = sqrt((kx - kxG).^2 + (ky - kyG).^2);
        w  = omega * obj.demodOpts.filterWidthFrac;
        obj.phaseData.bandpassFilter = exp(-(kr/w).^8);

        obj.phaseData.refComplexCoeffs = obj.cropPadding( ...
            ifft2(obj.phaseData.bandpassFilter .* fft2(refPadded)));
    end

    function demodulateRefWavelet(obj)
    % Wavelet ridge of the reference image.
    %   Runs the Morlet CWT over demodOpts.waveletScales/waveletAngles,
    %   picks the maximum-amplitude (scale, angle) per pixel, smooths that
    %   ridge and interpolates the complex coefficients between the four
    %   neighboring (scale, angle) planes. Result is stored in
    %   phaseData.refComplexCoeffs.

        [Ny, Nx] = size(obj.imgData.ref);
        % Demodulate reference image
        W4D = obj.morletCWT(obj.imgData.ref, ...
            obj.demodOpts.waveletScales, obj.demodOpts.waveletAngles, ...
            obj.demodOpts.waveletSigma, obj.demodOpts.waveletGamma, ...
            obj.fringe.normalAxis);
        W4D_abs = abs(W4D);

        % For each pixel, find the (scale, angle) at which the amplitude of
        % the wavelet transform is maximal, then smooth the ridge before
        % extracting the phase
        nScales = size(W4D, 3);
        nAngles = size(W4D, 4);   % 1 for 1-D input regardless of waveletAngles
        sigma_s = 2;              % Smoothing for scale map
        sigma_a = 2;              % Smoothing for angle map

        % --- Decompose and smooth the index maps ---
        [~, maxInd] = max(W4D_abs, [], [3 4], 'linear');
        [~, ~, scale_idx_raw, angle_idx_raw] = ind2sub(size(W4D), maxInd);


        % Apply Gaussian smoothing to get fractional scale and angle maps
        s_smooth = imgaussfilt(double(scale_idx_raw), sigma_s);
        a_smooth = imgaussfilt(double(angle_idx_raw), sigma_a);

        % Clamp to stay within valid array bounds [1, N-0.001]
        s_smooth = max(1, min(nScales - 0.001, s_smooth));
        a_smooth = max(1, min(nAngles - 0.001, a_smooth));

        % Bilinear interpolation
        [I, J] = ndgrid(1:Ny, 1:Nx);

        % Find the 4 integer neighbors for every pixel
        s_low = floor(s_smooth); s_high = min(nScales, s_low + 1);
        a_low = floor(a_smooth); a_high = min(nAngles, a_low + 1);

        % Calculate weights
        u = s_smooth - s_low; % scale weight
        v = a_smooth - a_low; % angle weight

        % Get linear indices for the 4 neighbors in the 4D array
        idx00 = sub2ind(size(W4D), I, J, s_low,  a_low);  % Low Scale, Low Angle
        idx10 = sub2ind(size(W4D), I, J, s_high, a_low);  % High Scale, Low Angle
        idx01 = sub2ind(size(W4D), I, J, s_low,  a_high); % Low Scale, High Angle
        idx11 = sub2ind(size(W4D), I, J, s_high, a_high); % High Scale, High Angle

        % Bilinear blend of complex coefficients
        W_interp = (1-u).*(1-v).*W4D(idx00) + ... % Point (0,0)
            u .*(1-v).*W4D(idx10) + ... % Point (1,0)
            (1-u).* v  .*W4D(idx01) + ... % Point (0,1)
            u .* v  .*W4D(idx11);       % Point (1,1)

        obj.phaseData.refComplexCoeffs = W_interp;
    end

    function imgPadded = padImage(obj, img)
        % Pad an image before Fourier demodulation.
        %   Applies symmetric padding in the direction along the fringe
        %   pattern and zero padding across it.

        if obj.procOpts.dcSubtraction
            padValue = 0;
        else
            padValue = mean(img, 'all');
        end

        if strcmpi(obj.fringe.normalAxis, 'Y')
            imgPadded = padarray(img, [obj.phaseData.fourierPadding(1), 0], padValue);
            imgPadded = padarray(imgPadded, [0, obj.phaseData.fourierPadding(2)], 'symmetric');
        else
            imgPadded = padarray(img, [0, obj.phaseData.fourierPadding(2)], padValue);
            imgPadded = padarray(imgPadded, [obj.phaseData.fourierPadding(1), 0], 'symmetric');
        end
    end

    function arr = cropPadding(obj, arr)
        % Remove the padding added by padImage().
        pad = obj.phaseData.fourierPadding;
        arr = arr(pad(1) + 1:end - pad(1), pad(2) + 1:end - pad(2));
    end

    function prepareParallaxGeometry(obj)
    % Ray geometry of the 2-D parallax error correction.
    %   The camera ray through STACK node p meets the reference plane at
    %   node p itself, and the plane z = h at the world position
    %   (x + h*a, y + h*b), where (x, y) is the node's position and (a, b)
    %   is the sideways shift of the ray per mm of elevation. The mesh is a
    %   regular grid, so world positions convert to grid positions
    %   (fractional row and column) by an affine map, and in grid units the
    %   sample stored at node p lies at p + h*slope(p), where slope is in
    %   grid pixels per mm of elevation. Its row and column components are
    %   stored in worldCoords.mesh.parallax.
        mesh = obj.worldCoords.mesh;
        [ny, nx] = size(mesh.x);

        K_inv = inv(obj.cameraCalib.intrinsicsMatlab.K);
        R     = obj.cameraCalib.extrinsics.R;
        t     = obj.cameraCalib.extrinsics.Translation';
        pix = [mesh.xPixel(:), mesh.yPixel(:), ones(nx*ny, 1)]';  % 3 x N homogeneous pixels
        RKuv  = R' * (K_inv * pix);     % ray directions
        rinvT = R' * t;                 % camera-position term
        zSign = -sign(rinvT(3));        % see correctParallax1D

        a  = reshape(zSign * RKuv(1, :) ./ RKuv(3, :), ny, nx);
        b  = reshape(zSign * RKuv(2, :) ./ RKuv(3, :), ny, nx);

        % affine map [x; y] = A*[col - 1; row - 1] + [x(1,1); y(1,1)]
        A = [mesh.x(1, 2) - mesh.x(1, 1), mesh.x(2, 1) - mesh.x(1, 1)
             mesh.y(1, 2) - mesh.y(1, 1), mesh.y(2, 1) - mesh.y(1, 1)];
        Ainv = inv(A);

        obj.worldCoords.mesh.parallax = struct( ...
            'slopeRow', Ainv(2,1)*a + Ainv(2,2)*b, ...
            'slopeCol', Ainv(1,1)*a + Ainv(1,2)*b);
    end

    function discardFrame = preprocessCurrentImage(obj)
    % Bring imgData.curr into the computational frame.
    %  Interpolates bad pixels, applies the camera calibration
    %  (polynomial dewarp or pinhole undistort + rectify), flags the
    %  frame for discarding, locates the spatial phase-correction peaks
    %  in the RAW frame, then resizes, crops and removes the DC
    %  component.

        obj.imgData.curr = obj.rectify(obj.imgData.curr);
        compRectOrg = obj.procOpts.compRectOrg;

        discardFrame = false;
        % check if image should be discarded
        if ~isempty(obj.procOpts.discardThreshold)
            imageMean = mean(imcrop(obj.imgData.curr, compRectOrg), 'all', 'omitmissing');
            if imageMean > obj.procOpts.discardThreshold
                discardFrame = true;
            end
        end

        if strcmpi(obj.phaseCorrOpts.method, 'spatial') && ~discardFrame
            obj.phaseCorrOpts.currPeakLocs = zeros(1, length(obj.phaseCorrOpts.lineInds));
            for i = 1:length(obj.phaseCorrOpts.lineInds)
                [~, peakLocs] = detectFringePeaks(obj, obj.imgData.curr, obj.phaseCorrOpts.lineInds(i));

                if strcmpi(obj.phaseCorrOpts.startEdge, 'right') || strcmpi(obj.phaseCorrOpts.startEdge, 'bottom')
                    peakLocs = flip(peakLocs);
                end
                if length(peakLocs) >= obj.phaseCorrOpts.trackedPeakInd
                    obj.phaseCorrOpts.currPeakLocs(i) = peakLocs(obj.phaseCorrOpts.trackedPeakInd);
                else
                    obj.phaseCorrOpts.currPeakLocs(i) = nan;
                end
            end
        end
        if obj.procOpts.resizeFactor ~= 1
            obj.imgData.curr = FtpSolver.shrinkArr(obj.imgData.curr, obj.procOpts.resizeFactor);
        end

        if ~isempty(obj.procOpts.compRect)
            obj.imgData.curr = imcrop(obj.imgData.curr, obj.procOpts.compRect);
        end

        % store the image without background subtraction if it is invalid
        if obj.procOpts.dcSubtraction && ~discardFrame
            obj.imgData.curr = obj.subtractDC(obj.imgData.curr, ...
                obj.fringe.period);
        end
    end

    function demodulate(obj)
        if strcmpi(obj.demodOpts.method, 'fourier')
            obj.demodulateFourier();
        elseif strcmpi(obj.demodOpts.method, 'wavelet')
            obj.demodulateWavelet();
        else
            error("Invalid demodulation method")
        end
    end

    function demodulateFourier(obj)
    % Find the phase of the current image relative to the reference.
    % The band-pass filter and the complex reference coefficients are
    % precomputed by demodulateRefFourier. Only the current image is
    % processed here.

        % Band-pass and inverse transform
        currPadded = obj.padImage(obj.imgData.curr);
        hI = obj.cropPadding( ...
            ifft2(obj.phaseData.bandpassFilter .* fft2(currPadded)));

        % Phase difference relative to the reference
        delPhi = angle(hI .* conj(obj.phaseData.refComplexCoeffs));

        if strcmpi(obj.phaseCorrOpts.method, 'temporal')
            obj.phaseData.prev = obj.phaseData.curr;
        end

        obj.phaseData.curr = delPhi;
    end

    function demodulateWavelet(obj)
    % Wavelet phase of the current frame.
    %  Same ridge extraction and bilinear (scale, angle) interpolation as
    %  demodulateRefWavelet, applied to imgData.curr; phaseData.curr is
    %  set to the phase relative to the reference coefficients.

        [Ny, Nx] = size(obj.imgData.curr);

         W4D = obj.morletCWT(obj.imgData.curr, ...
             obj.demodOpts.waveletScales, obj.demodOpts.waveletAngles, ...
             obj.demodOpts.waveletSigma, obj.demodOpts.waveletGamma, ...
             obj.fringe.normalAxis);
         W4D_abs = abs(W4D);

         % For each pixel, find the (scale, angle) at which the amplitude of
         % the wavelet transform is maximal, then smooth the ridge before
         % extracting the phase.
         nScales = size(W4D, 3);
         nAngles = size(W4D, 4);   % 1 for 1-D input regardless of waveletAngles
         sigma_s = 2;              % Smoothing for scale map
         sigma_a = 2;              % Smoothing for angle map

         % --- Decompose and smooth the index maps ---
         [~, maxInd] = max(W4D_abs, [], [3 4], 'linear');
         [~, ~, scale_idx_raw, angle_idx_raw] = ind2sub(size(W4D), maxInd);

         % Apply Gaussian smoothing to get fractional scale and angle maps
         s_smooth = imgaussfilt(double(scale_idx_raw), sigma_s);
         a_smooth = imgaussfilt(double(angle_idx_raw), sigma_a);

         % Clamp to stay within valid array bounds [1, N-0.001]
         s_smooth = max(1, min(nScales - 0.001, s_smooth));
         a_smooth = max(1, min(nAngles - 0.001, a_smooth));

         % Bilinear interpolation
         [I, J] = ndgrid(1:Ny, 1:Nx);

         % Find the 4 integer neighbors for every pixel
         s_low = floor(s_smooth); s_high = min(nScales, s_low + 1);
         a_low = floor(a_smooth); a_high = min(nAngles, a_low + 1);

         % Calculate weights
         u = s_smooth - s_low; % scale weight
         v = a_smooth - a_low; % angle weight

         % Get linear indices for the 4 neighbors in the 4D array
         idx00 = sub2ind(size(W4D), I, J, s_low,  a_low);  % Low Scale, Low Angle
         idx10 = sub2ind(size(W4D), I, J, s_high, a_low);  % High Scale, Low Angle
         idx01 = sub2ind(size(W4D), I, J, s_low,  a_high); % Low Scale, High Angle
         idx11 = sub2ind(size(W4D), I, J, s_high, a_high); % High Scale, High Angle

         % Bilinear blend of complex coefficients
         W_interp = (1-u).*(1-v).*W4D(idx00) + ... % Point (0,0)
             u .*(1-v).*W4D(idx10) + ... % Point (1,0)
             (1-u).* v  .*W4D(idx01) + ... % Point (0,1)
             u .* v  .*W4D(idx11);       % Point (1,1)

         if strcmpi(obj.phaseCorrOpts.method, 'temporal')
             obj.phaseData.prev = obj.phaseData.curr;
         end

         obj.phaseData.curr = angle(W_interp .* conj(obj.phaseData.refComplexCoeffs));
    end

    function unwrapAndFlag(obj)
    % Unwrap the current phase map, flag possible faults, resize to STACK.
    %   Unwraps the part of phaseData.curr inside the unwrap margin, first
    %   along dimension 1 and then along dimension 2. With unwrapMethod
    %   "2D", the result is replaced by that of unwrap2D when its faults
    %   are longer than procOpts.faultDetection.lengthFor2DPeriods; unwrap2D
    %   unwraps from scratch and uses its input only to rank pixel
    %   reliability. Faults that remain are recorded in
    %   diagnostics.flaggedRegions. With unwrapMethod "none" the phase
    %   stays wrapped and nothing is flagged. In every case the phase map
    %   is then resized to the STACK frame with the margin set to NaN.
        [rows, cols, mRow, mCol] = obj.interiorInds();
        interior = obj.phaseData.curr(rows, cols);

        if obj.procOpts.unwrapMethod ~= "none"
            interior = unwrap(interior, [], 1);
            interior = unwrap(interior, [], 2);
            faults = obj.detectPhaseFaults(interior);

            if obj.procOpts.unwrapMethod == "2D" && FtpSolver.faultLength(faults) ...
                    > obj.procOpts.faultDetection.lengthFor2DPeriods * obj.fringe.period
                interior = unwrap2D(interior);
                faults = obj.detectPhaseFaults(interior);
            end

            obj.recordFlaggedRegions(faults, mRow, mCol);
        end

        obj.phaseData.curr = obj.resizeToStack(interior, mRow, mCol);
    end

    function [rows, cols, mRow, mCol] = interiorInds(obj)
    % Rows and columns of the COMP-frame phase map inside the unwrap margin.
    %   mRow and mCol are the margins applied along dimensions 1 and 2. A
    %   dimension too short to hold the margin on both sides (1-D data)
    %   gets none.
        [ny, nx] = size(obj.phaseData.curr);
        margin = obj.procOpts.unwrapMargin;      % [mX mY], COMP
        mRow = margin(2) * (ny > 2*margin(2));
        mCol = margin(1) * (nx > 2*margin(1));
        rows = mRow + 1:ny - mRow;
        cols = mCol + 1:nx - mCol;
    end

    function faults = detectPhaseFaults(obj, phase)
    % Mask of pixels where the phase gradient jumps as at a 2*pi fault.
    %   Compares the gradient magnitude with a median-filtered copy of
    %   itself and marks the pixels whose difference exceeds
    %   faultDetection.residualFactor times the median difference. The
    %   median window spans faultDetection.medianWindowPeriods fringe
    %   periods. It is odd, so the filter is centered, and at least 5 pixels
    %   wide, so that a fault line two pixels wide cannot fill half of it.
        settings = obj.procOpts.faultDetection;
        window = obj.periodsToPixels(settings.medianWindowPeriods, 5);
        window = 2*floor(window/2) + 1;
        gradMag = imgradient(phase);
        residual = abs(gradMag - medfilt2(gradMag, [window window], 'symmetric'));
        faults = residual > settings.residualFactor * median(residual, 'all');
    end

    function recordFlaggedRegions(obj, faults, mRow, mCol)
    % Record bounding boxes of phase faults in diagnostics.flaggedRegions.
    %   Nothing is recorded when the faults are no longer than
    %   faultDetection.lengthToFlagPeriods. Otherwise the faults are
    %   padded by faultDetection.dilationPeriods, faults close to
    %   each other are merged (closing with radius faultDetection.closePeriods),
    %   and the bounding box of each merged region is converted from the
    %   unwrap interior (COMP frame) to pixel coordinates of the result
    %   arrays (OUTPUT frame) and clipped to them. The current frame is
    %   recorded only if a box survives the clipping.
        settings = obj.procOpts.faultDetection;
        if FtpSolver.faultLength(faults) <= settings.lengthToFlagPeriods * obj.fringe.period
            return
        end
        faults = imdilate(faults, strel("disk", obj.periodsToPixels(settings.dilationPeriods, 0)));
        faults = imclose(faults, strel("disk", obj.periodsToPixels(settings.closePeriods, 0)));
        regions = regionprops(faults, 'BoundingBox');
        rects = cat(1, regions.BoundingBox);
        if isempty(rects)
            return
        end

        % unwrap interior (COMP) to OUTPUT frame
        rf_o = obj.procOpts.resizeFactorOutput;
        outputRect = round(obj.procOpts.outputRect);
        rects(:, 1:2) = FtpSolver.rescaleCoord(rects(:, 1:2) + [mCol, mRow], rf_o) ...
                        - outputRect(1:2) + 1;
        rects(:, 3:4) = rects(:, 3:4) * rf_o;

        % clip to the result arrays, whose far edge lies half a pixel beyond
        % the last pixel center
        bottomRight = min(rects(:, 1:2) + rects(:, 3:4), outputRect(3:4) + 1.5);
        rects(:, 1:2) = max(0.5, rects(:, 1:2));
        rects(:, 3:4) = bottomRight - rects(:, 1:2);
        rects(any(rects(:, 3:4) <= 0, 2), :) = [];

        if ~isempty(rects)
            obj.diagnostics.flaggedRegions(end + 1) = struct( ...
                'timestep', obj.loopState.timestep, 'rects', rects);
        end
    end

    function phase = resizeToStack(obj, interior, mRow, mCol)
    % Resize the unwrap interior to the STACK frame, with the margin NaN.
    %   The margin is first filled with the nearest interior values, so the
    %   resize blends only interior phase into the pixels at the interior
    %   edge. Resized pixels whose centers lie in the margin are then set
    %   to NaN.
        rf_o = obj.procOpts.resizeFactorOutput;
        phase = padarray(interior, [mRow, mCol], 'replicate');
        [ny, nx] = size(phase);
        phase = FtpSolver.shrinkArr(phase, rf_o, 'bilinear', Antialiasing=false);

        % resized pixel k is centered at (k - 0.5)/rf_o + 0.5 in the input
        [nyS, nxS] = size(phase);
        rowCenters = FtpSolver.rescaleCoord(1:nyS, 1/rf_o);
        colCenters = FtpSolver.rescaleCoord(1:nxS, 1/rf_o);
        phase(rowCenters < mRow + 0.5 | rowCenters > ny - mRow + 0.5, :) = NaN;
        phase(:, colCenters < mCol + 0.5 | colCenters > nx - mCol + 0.5) = NaN;
    end

    function n = periodsToPixels(obj, nPeriods, minPixels)
    % Whole number of COMP-frame pixels spanning nPeriods fringe periods.
        if nargin < 3
            minPixels = 1;
        end
        n = max(minPixels, round(nPeriods * obj.fringe.period));
    end

    function correctPhase(obj)
    % Remove 2*pi ambiguities from the current phase map.
    %   Frames where the estimate is unreliable are listed in
    %   diagnostics.phaseCorrectionFailures and left uncorrected.

        switch obj.phaseCorrOpts.method
            case "temporal"
                nJumps = obj.temporalJumpCount();
            case "spatial"
                nJumps = obj.spatialJumpCount();
            otherwise
                return
        end
        obj.phaseData.curr = obj.phaseData.curr + 2*pi*nJumps;
    end

    function nJumps = temporalJumpCount(obj)
    % Calculates the 2*pi jumps necessary to correct this frame's phase using the previous frame.
    %   Counts the whole 2*pi jumps between this frame and the previous one
    %   at phaseCorrOpts.numSamples interior points. Returns zero and records
    %   the frame in diagnostics.phaseCorrectionFailures if no jump count
    %   is shared by more than half the points. The first frame has no
    %   predecessor and sets the reference for the rest of the run.

        nJumps = 0;
        if obj.loopState.isFirstTimestep
            return
        end
        numSamples = obj.phaseCorrOpts.numSamples;
        [Ny, Nx] = size(obj.phaseData.curr);

        margin = ceil(obj.procOpts.unwrapMargin * obj.procOpts.resizeFactorOutput);
        % no margin along a dimension too short to have one (1-D data)
        margin(2) = margin(2) * (Ny > 2*margin(2));
        margin(1) = margin(1) * (Nx > 2*margin(1));

        rows = floor(linspace(1 + margin(2), Ny - margin(2), numSamples + 2));
        cols = floor(linspace(1 + margin(1), Nx - margin(1), numSamples + 2));
        idx = sub2ind([Ny, Nx], rows(2:end - 1), cols(2:end - 1));
        candidates = round((obj.phaseData.prev(idx) - obj.phaseData.curr(idx)) / (2*pi));
        [nJumps, frequency] = mode(candidates);
        if frequency <= numSamples/2
            obj.recordPhaseCorrectionFailure(candidates);
            nJumps = 0;
        end
    end

    function nJumps = spatialJumpCount(obj)
    % Calculates the 2*pi jumps necessary to correct this frame's phase using the tracked fringe peak.
    %   On each correction line, the shift of the tracked peak from its
    %   reference position gives the true phase difference at the peak's
    %   current position. Comparing it with the demodulated phase there
    %   gives a whole number of 2*pi jumps for that line. This function
    %   returns the count that more than half the lines agree on. If no
    %   line gives a usable count, or fewer than a majority of the lines
    %   agree on the jump, the frame is recorded in
    %   diagnostics.phaseCorrectionFailures and the function returns zero.
    %   The recorded candidates hold one count per correction line and show
    %   NaN where the line could not be used.

        nJumps = 0;
        lineInds = obj.phaseCorrOpts.lineInds(:);
        candidates = nan(1, numel(lineInds));

        currPeakLocs   = obj.phaseCorrOpts.currPeakLocs;
        refPeakLocs    = obj.phaseCorrOpts.refPeakLocs;
        localPeriod = obj.phaseCorrOpts.localPeriod;

        % lines on which the tracked peak was found in both images
        lines = find(~isnan(currPeakLocs) & ~isnan(refPeakLocs) & ~isnan(localPeriod));
        if isempty(lines)
            obj.recordPhaseCorrectionFailure(candidates);
            return
        end

        % position of each peak in the phase map (STACK frame)
        switch upper(obj.fringe.normalAxis)
            case "X"
                rowCol = obj.RAW2STACK([lineInds(lines), currPeakLocs(lines).']);
            case "Y"
                rowCol = obj.RAW2STACK([currPeakLocs(lines).', lineInds(lines)]);
            otherwise
                error("Pattern normal axis not defined.")
        end

        % keep the peaks that fall inside the phase map
        [Ny, Nx] = size(obj.phaseData.curr);
        inside = rowCol(:, 1) >= 1 & rowCol(:, 1) <= Ny & ...
                 rowCol(:, 2) >= 1 & rowCol(:, 2) <= Nx;
        lines  = lines(inside);
        rowCol = rowCol(inside, :);

        % jumps per line: true phase from the peak shift minus measured phase
        truePhase = 2*pi * (refPeakLocs(lines) - currPeakLocs(lines)) ./ localPeriod(lines);
        measuredPhase = obj.phaseData.curr(sub2ind([Ny, Nx], rowCol(:, 1), rowCol(:, 2)));
        candidates(lines) = round((truePhase(:) - measuredPhase(:)) / (2*pi));

        counts = candidates(~isnan(candidates));
        if isempty(counts)
            obj.recordPhaseCorrectionFailure(candidates);
            return
        end
        [nJumps, frequency] = mode(counts);
        if frequency <= numel(lineInds)/2
            obj.recordPhaseCorrectionFailure(candidates);
            nJumps = 0;
        end
    end

    function recordPhaseCorrectionFailure(obj, candidates)
    % Add the current frame to diagnostics.phaseCorrectionFailures.
        obj.diagnostics.phaseCorrectionFailures(end + 1) = struct( ...
            'timestep', obj.loopState.timestep, 'candidates', candidates);
    end

    function calculateElevation(obj)
    % Calculates the surface elevation based on the selected model.
        if ~obj.outputConfig.elevStorage || strcmpi(obj.elevModel.type, 'none')
            return
        end

        if strcmpi(obj.elevModel.type, 'takeda')
            period = obj.fringe.period;
            pixelPitch = obj.worldCoords.scaling.mmPerPixel;

            obj.elevData.curr = ...
                FtpSolver.evalTakedaElev(obj.phaseData.curr, obj.elevModel.L, ...
                                        obj.elevModel.d, (period*pixelPitch)^-1);
        elseif strcmpi(obj.elevModel.type, 'polynomial')
            obj.elevData.curr = ...
                FtpSolver.evalPolyElev(obj.phaseData.curr, ...
                                       obj.elevModel.poly.coeffMat);
        end

        if obj.procOpts.parallaxCorrection
            obj.correctParallax()
        end
    end

    function correctParallax(obj)
    % Resample elevData.curr onto the nominal world grid.
    %   A pinhole camera sees each surface point along a slanted ray, so
    %   the pixel nominally associated with world position (x, y) actually
    %   samples the surface at a laterally shifted position whenever the
    %   surface is displaced from the reference plane z = 0. This method
    %   interpolates the elevation at the nominal mesh nodes from the true
    %   sample positions using fixed point iteration. Requires a pinhole
    %   calibration.

        if ~isfield(obj.cameraCalib, 'intrinsicsMatlab')
            return
        end

        % The unwrap margin is NaN after unwrapAndFlag and is left out of
        % the resampling together with any other invalid pixels. Without
        % unwrapping it is not masked and enters the resampling.
        if strcmpi(obj.procOpts.unwrapMethod, 'none') && obj.loopState.isFirstTimestep
            warning("FtpSolver:parallaxCorrectionWithoutUnwrap", ...
                "Parallax error correction is on but unwrapping is disabled. " + ...
                "The unwrap margin is not masked and is included in the resampling.")
        end

        if min(size(obj.imgData.ref)) < 2
            obj.correctParallax1D();
            return
        end

        tol        = 1e-3;  % convergence tolerance, grid pixels
        maxIter    = 30;
        activeFrac = 0.02;  % share of moving nodes below which only those are iterated
        edgeTol    = 0.05;  % grid pixels, applied to pixels on the edges

        elev = obj.elevData.curr;
        geom = obj.worldCoords.mesh.parallax;
        [ny, nx] = size(elev);
        valid = isfinite(elev);
        if ~any(valid, 'all')
            return
        end

        % NaN nodes take the nearest valid elevation, only for dRow and dCol
        elevFilled = elev;
        if ~all(valid, 'all')
            [~, nearest] = bwdist(valid);
            elevFilled = elev(nearest);
        end
        % row and column distance (in grid steps) from each node to the
        % point where its elevation was measured
        dRow = elevFilled .* geom.slopeRow;
        dCol = elevFilled .* geom.slopeCol;
        dInterp = griddedInterpolant({1:ny, 1:nx}, cat(3, dRow, dCol), 'linear', 'none');

        % the first guess is the node itself, so the first update uses dRow
        % and dCol directly without interpolation
        [qCol, qRow] = meshgrid(1:nx, 1:ny);
        pRow = qRow - dRow;
        pCol = qCol - dCol;

        % all nodes together while many are still moving
        active = [];
        k = 0;
        while k < maxIter
            k = k + 1;
            % dRow and dCol interpolated at the current guesses
            D = dInterp(pRow, pCol);
            newRow = qRow - D(:, :, 1);
            newCol = qCol - D(:, :, 2);
            moving = abs(newRow - pRow) >= tol | abs(newCol - pCol) >= tol;
            pRow = newRow;
            pCol = newCol;
            if nnz(moving) <= activeFrac * ny * nx
                active = find(moving);
                break
            end
        end

        % then the few that are still moving
        while ~isempty(active) && k < maxIter
            k = k + 1;
            D = dInterp(pRow(active), pCol(active));
            newRow = qRow(active) - D(:, 1);
            newCol = qCol(active) - D(:, 2);
            moving = abs(newRow - pRow(active)) >= tol | abs(newCol - pCol(active)) >= tol;
            pRow(active) = newRow;
            pCol(active) = newCol;
            active = active(moving);
        end

        elevInterp = griddedInterpolant({1:ny, 1:nx}, elev, 'linear', 'none');
        elevOut = elevInterp(pRow, pCol);

        % samples on the edge of the data
        retry = find(isnan(elevOut) & isfinite(pRow) & isfinite(pCol));
        r = pRow(retry);
        c = pCol(retry);
        r = r + (round(r) - r) .* (abs(round(r) - r) < edgeTol);
        c = c + (round(c) - c) .* (abs(round(c) - c) < edgeTol);
        elevOut(retry) = elevInterp(r, c);

        obj.elevData.curr = elevOut;
    end

    function correctParallax1D(obj)
    % Parallax error correction of 1-D data (see correctParallax).
    %   Computes the true position of every sample from its camera ray and
    %   interpolates the elevations back onto the nominal mesh.

        elev = obj.elevData.curr;
        u    = obj.worldCoords.mesh.xPixel;
        v    = obj.worldCoords.mesh.yPixel;

        % Back-project every pixel to a world-frame camera ray
        K_inv = inv(obj.cameraCalib.intrinsicsMatlab.K);
        R     = obj.cameraCalib.extrinsics.R;
        t     = obj.cameraCalib.extrinsics.Translation';

        pix   = [u(:), v(:), ones(numel(u), 1)]';   % 3 x N homogeneous pixels
        RKuv  = R' * (K_inv * pix);                 % 3 x N ray directions
        rinvT = R' * t;                             % camera-position term

        % elevData is positive when the surface is raised above the
        % reference plane. If the world z-axis points away from the camera,
        % a raised surface has negative world z. Multiplying by
        % -sign(rinvT(3)), which is +1 when z points towards the camera,
        % gives the correct sign either way.
        zSign = -sign(rinvT(3));

        kField = (zSign * elev(:)' + rinvT(3)) ./ RKuv(3, :);
        res    = RKuv .* kField - rinvT;    % 3 x N true world positions

        X_corr = res(1, :);
        Y_corr = res(2, :);

        % handle NaN values
        ok = isfinite(X_corr) & isfinite(Y_corr) & isfinite(elev(:)');

        elev = elev(:);
        if strcmpi(obj.fringe.normalAxis, 'X')
            obj.elevData.curr = interp1(X_corr(ok), elev(ok), ...
                                        obj.worldCoords.mesh.x, 'linear');
        else
            obj.elevData.curr = interp1(Y_corr(ok), elev(ok), ...
                                        obj.worldCoords.mesh.y, 'linear');
        end
    end

    function storeData(obj)
    % Store current elevation, phase or image to corresponding stack.
        if obj.outputConfig.phaseStorage
            if obj.loopState.discardFrame
                obj.phaseData.stack(:,:,obj.loopState.stackInd) = nan;
            else
                obj.phaseData.stack(:,:,obj.loopState.stackInd) = obj.phaseData.curr;
            end
        end

        if obj.outputConfig.elevStorage
            if obj.loopState.discardFrame
                obj.elevData.stack(:,:,obj.loopState.stackInd) = nan;
            else
                obj.elevData.stack(:,:,obj.loopState.stackInd) = obj.elevData.curr;
            end
        end

        if obj.outputConfig.imageStorage
            obj.imgData.stack(:,:,obj.loopState.stackInd) = obj.imgData.curr;
        end
    end

    function reportProgress(obj)
        if ~obj.outputConfig.reportProgressEnabled
            return
        end
        textstring = sprintf('Processing image %5d...', obj.loopState.timestep);
        if obj.loopState.timestep > obj.procOpts.solveRange(1)
            textstringOld = sprintf('Processing image %5d...', obj.loopState.timestep - 1);
            fprintf(repmat('\b', 1, numel(textstringOld)))
            fprintf(textstring)
        else
            % first timestep
            fprintf(textstring)
        end
    end

    function flushBuffer(obj)
    % Write the buffer held in memory to the result files.
    %   Writes the output subarray of each stored stack into its result
    %   file at the buffer's frame offset, then reallocates the stacks for
    %   the next buffer unless this was the last one.

        solveRange = obj.procOpts.solveRange;
        bufferIndex = obj.loopState.bufferIndex;
        bufferLength = obj.outputConfig.bufferLength;
        nTimesteps = diff(solveRange) + 1;
        nBuffers = ceil(nTimesteps / bufferLength);
        firstFrame = (bufferIndex - 1)*bufferLength + 1;     % frame index in the file

        if obj.outputConfig.elevStorage
            nInBuffer = size(obj.elevData.stack, 3);
        else
            nInBuffer = size(obj.phaseData.stack, 3);
        end

        for quantity = obj.storedQuantities()
            frames = obj.outputSubarray(quantity, solveRange(1), ...
                solveRange(1) + nInBuffer - 1);
            obj.writeFrames(quantity, frames, firstFrame);
        end

        if bufferIndex == nBuffers
            return
        end

        % the last buffer holds whatever timesteps remain
        if bufferIndex == nBuffers - 1
            nextBufferLength = nTimesteps - bufferIndex*bufferLength;
        else
            nextBufferLength = bufferLength;
        end
        obj.allocateStacks(nextBufferLength)
    end

    function refreshAutoPeak(obj)
    % Fill an empty phaseCorrOpts.trackedPeakInd with the automatic choice, if one exists.
        if isempty(obj.phaseCorrOpts.trackedPeakInd)
            optInd = obj.findOptimalPeakInd();
            if ~isnan(optInd)
                obj.phaseCorrOpts.trackedPeakInd = optInd;
            end
        end
    end

    function [optInd, isSafeLine, peakLocsCell] = findOptimalPeakInd(obj, peakInd)
    % Find which phase-correction peaks can be tracked inside the output domain.
    %
    %   [optInd, isSafeLine, peakLocsCell] = findOptimalPeakInd(obj, peakInd)
    %
    %   Detects the reference peak train on every phase-correction line and
    %   determines which peaks lie inside the output domain while keeping a
    %   clearance of phaseCorrOpts.peakClearance fringe periods from the
    %   starting edge, so that motion of the surface cannot move the peak
    %   out of view. The first and last detected peaks are never eligible
    %   because the local-wavelength stencil in preprocessRefImage needs a
    %   neighbor on each side.
    %
    %   optInd        - recommended trackedPeakInd: the eligible peak
    %                   closest to the starting edge, maximized over the
    %                   lines so the choice is valid on all of them. NaN
    %                   when no line has an eligible peak, or when the
    %                   inputs listed under "prerequisites" below are not
    %                   set.
    %   isSafeLine    - per-line logical: true where the candidate peakInd
    %                   is eligible. All false when peakInd is omitted or
    %                   empty.
    %   peakLocsCell  - detected peak positions per line, RAW-frame
    %                   coordinates, flipped for startEdge 'right'/'bottom'
    %

        if nargin < 2
            peakInd = [];
        end

        nLines = length(obj.phaseCorrOpts.lineInds);
        optInd = nan;
        isSafeLine = false(1, nLines);
        peakLocsCell = cell(1, nLines);

        % prerequisites: reference image, correction lines, fringe period,
        % pattern normal axis, and output-domain geometry must all exist
        if isempty(obj.imgData.refRaw) || nLines == 0 ...
                || ~isfield(obj.fringe, 'periodOrg') ...
                || isempty(obj.fringe.periodOrg) ...
                || isnan(obj.fringe.periodOrg) ...
                || strlength(string(obj.fringe.normalAxis)) == 0 ...
                || isempty(obj.procOpts.compRectOrg) ...
                || isempty(obj.procOpts.outputRectOrg)
            return
        end

        % safe window along the pattern normal axis, RAW frame
        compRectOrg = obj.procOpts.compRectOrg;
        outputRectOrg = obj.procOpts.outputRectOrg;
        if strcmpi(obj.fringe.normalAxis, 'X')
            winStart = compRectOrg(1) + outputRectOrg(1) - 1;
            winEnd   = winStart + outputRectOrg(3);
        else
            winStart = compRectOrg(2) + outputRectOrg(2) - 1;
            winEnd   = winStart + outputRectOrg(4);
        end

        safetyDist = obj.phaseCorrOpts.peakClearance * obj.fringe.periodOrg;

        % the starting edge is winEnd for 'right'/'bottom' (where the peak
        % trains are flipped) and winStart otherwise
        if strcmpi(obj.phaseCorrOpts.startEdge, 'right') || strcmpi(obj.phaseCorrOpts.startEdge, 'bottom')
            safeLo = winStart;
            safeHi = winEnd - safetyDist;
        else
            safeLo = winStart + safetyDist;
            safeHi = winEnd;
        end

        % detect and evaluate every line
        optIndsPerLine = nan(1, nLines);
        for i = 1:nLines
            [~, peakLocs] = detectFringePeaks(obj, obj.imgData.refRaw, obj.phaseCorrOpts.lineInds(i));

            if strcmpi(obj.phaseCorrOpts.startEdge, 'right') || strcmpi(obj.phaseCorrOpts.startEdge, 'bottom')
                peakLocs = flip(peakLocs);
            end

            peakLocsCell{i} = peakLocs;

            % inside the safe window, with a neighbor on each side
            eligible = peakLocs >= safeLo & peakLocs <= safeHi;
            if ~isempty(eligible)
                eligible(1) = false;
                eligible(end) = false;
            end

            firstEligible = find(eligible, 1, 'first');
            if ~isempty(firstEligible)
                optIndsPerLine(i) = firstEligible;
            end

            if ~isempty(peakInd) && peakInd >= 2 && peakInd <= length(peakLocs) - 1
                isSafeLine(i) = eligible(peakInd);
            end
        end

        if ~all(isnan(optIndsPerLine))
            optInd = max(optIndsPerLine);   % NaN lines are ignored by max
        end
    end

    function [peakVals, peakLocs] = detectFringePeaks(obj, img, lineInd)
        compRectOrg = obj.procOpts.compRectOrg;
        if strcmpi(obj.fringe.normalAxis, 'X')
            cropStart = compRectOrg(1);
            cropEnd = compRectOrg(1) + compRectOrg(3);
            phaseCorrLine = double(img(lineInd, :));
        else
            cropStart = compRectOrg(2);
            cropEnd = compRectOrg(2) + compRectOrg(4);
            phaseCorrLine = double(img(:, lineInd));
        end

        segment = phaseCorrLine(cropStart:cropEnd);
        [minHeight, minProm] = FtpSolver.peakThresholds(segment, ...
            obj.phaseCorrOpts.peakMinFrac, obj.phaseCorrOpts.peakPromFrac);

        [peakVals, peakLocs] = findpeaks(phaseCorrLine, ...
            'MinPeakDistance', 0.6*obj.fringe.periodOrg, ...
            'MinPeakProminence', minProm, ...
            'MinPeakHeight', minHeight);
    end
end

%% ----------- Internal: input, geometry, output and plotting -----------
methods (Access = private)
    function readRefImage(obj, refPath)
    % Read the reference image into imgData.refRaw.
    %   Averages over the whole set if given a Davis .set file.
        refPath = FtpSolver.absolutePath(refPath);
        obj.source.refPath = refPath;
        [~, ~, ext] = fileparts(refPath);

        if strcmpi(ext, ".set")
            img = obj.readDavisFrame(refPath, obj.source.refCamInd, 'avg');
        elseif strcmpi(ext, ".im7")
            img = obj.readDavisFrame(refPath, obj.source.refCamInd);
        elseif ismember(lower(extractAfter(ext, ".")), [imformats().ext])
            img = imread(refPath);
        else
            error("Invalid file type for reference image.")
        end
        obj.imgData.refRaw = double(img);
    end

    function initDataset(obj, dataPath)
    % Record where the frames of the dataset are.
    %   dataPath - a Davis .set file, or any one file of an image sequence.
    %              Every other file with the same extension in that folder
    %              is treated as a frame, in file-name order.
        dataPath = FtpSolver.absolutePath(dataPath);
        [folder, ~, ext] = fileparts(dataPath);
        obj.source.dataPath = dataPath;

        if strcmpi(ext, ".set")
            obj.source.dataFormat = "set";
            obj.source.frameFiles = strings(0, 1);
            nFrames = lvsetsize(char(dataPath));
        else
            if strcmpi(ext, ".im7")
                obj.source.dataFormat = "im7";
            elseif ismember(lower(extractAfter(ext, ".")), [imformats().ext])
                obj.source.dataFormat = "image";
            else
                error("Invalid file type for dataset.")
            end
            obj.source.frameFiles = obj.listFrameFiles(folder, ext);
            nFrames = numel(obj.source.frameFiles);
        end

        obj.source.fullRange = [1, nFrames];
        if isempty(obj.procOpts.solveRange)
            obj.procOpts.solveRange = obj.source.fullRange;
        end
    end

    function files = listFrameFiles(obj, folder, ext)
    % Full paths of the frame files in folder, in file-name order.
        % keep the extension's original case so the pattern matches on
        % case-sensitive file systems
        listing = dir(fullfile(folder, "*" + ext));
        listing = listing(~[listing.isdir]);
        files = fullfile(string({listing.folder}), string({listing.name}))';

        % warn if the reference image is in the data sequence
        isRef = strcmpi(files, obj.source.refPath);
        if any(isRef)
            warning("The reference image is included as frame %d in " + ...
                     "the image sequence.", find(isRef))
        end

        if isempty(files)
            error("No %s files found in %s.", ext, folder)
        end

        % names of different lengths usually mean unpadded frame numbers,
        % which text order puts out of sequence (1, 10, 2, ...)
        names = string({listing(~isRef).name});
        if numel(unique(strlength(names))) > 1
            warning("Frame file names differ in length, so file-name order " + ...
                "may not be frame order. Zero-pad the frame numbers.")
        end
    end

    function img = readImage(obj, timestep)
    % Read frame timestep of the dataset as a double image.
        switch obj.source.dataFormat
            case "set"
                img = obj.readDavisFrame(obj.source.dataPath, obj.source.dataCamInd, timestep);
            case "im7"
                img = obj.readDavisFrame(obj.source.frameFiles(timestep), obj.source.dataCamInd);
            case "image"
                img = imread(obj.source.frameFiles(timestep));
            otherwise
                error("No dataset loaded.")
        end
        img = double(img);
    end

    function reloadInputs(obj)
    % Restore the images of a case loaded from a .mat.
    %  Image data and polynomial coefficients are not saved with the object,
    %  so this re-reads the reference image and the dataset listing from
    %  the stored paths, rebuilds the geometry, and reloads the polynomial
    %  calibration if the case used one.

        obj.readRefImage(obj.source.refPath);
        obj.initDataset(obj.source.dataPath);

        obj.rectifyRefImage();
        if strcmpi(obj.elevModel.type, 'polynomial') && isfield(obj.elevModel, 'path')
            obj.setElevModelPoly(obj.elevModel.path)
        end
    end

    function readCameraCalibration(obj, calibPath, camNumber)
        % Read camera calibration from .mat or .xml file.
        [~, ~, ext] = fileparts(calibPath);

        if strcmpi(ext, '.mat')
            obj.cameraCalib = load(calibPath);
            obj.cameraCalib.type = "pinhole";

            imSize = obj.cameraCalib.intrinsics.ImageSize;
            [~, pTransform, scaling, xMesh, yMesh] = rectifyImagePinhole(zeros(imSize), obj.cameraCalib, 0, obj.procOpts.imgRotAngle);
            obj.cameraCalib.pTransform = pTransform;
            obj.cameraCalib.scaling = scaling;

            obj.worldCoords.mesh.xOrg = xMesh;
            obj.worldCoords.mesh.yOrg = yMesh;

            obj.storePinholeScaling(scaling);

        elseif strcmpi(ext, '.xml')
            calibCell = readDavisCalibration(calibPath);
            calibData = calibCell{camNumber,1};
            calibType = calibCell{camNumber,2};

            if strcmpi(calibType, "pinhole")
                obj.cameraCalib = calibData;

                imSize = obj.cameraCalib.intrinsics.ImageSize;
                [~, pTransform, scaling, xMesh, yMesh] = rectifyImagePinhole(zeros(imSize), obj.cameraCalib, 0, obj.procOpts.imgRotAngle);
                obj.cameraCalib.pTransform = pTransform;
                obj.cameraCalib.scaling = scaling;

                obj.worldCoords.mesh.xOrg = xMesh;
                obj.worldCoords.mesh.yOrg = yMesh;

                obj.storePinholeScaling(scaling);

            elseif strcmpi(calibType, "polynomial")
                polyCalib = calibData{1};
                orgSize = polyCalib.orgSize;
                [Gx, Gy] = meshgrid(0:orgSize(1) - 1,0:orgSize(2) - 1);

                obj.cameraCalib.A = polyCalib.A;
                obj.cameraCalib.B = polyCalib.B;
                obj.cameraCalib.origin = polyCalib.origin;
                obj.cameraCalib.orgSize = orgSize;
                obj.cameraCalib.dwSize = polyCalib.dwSize;
                obj.cameraCalib.offset = polyCalib.offset;
                obj.cameraCalib.pixelsPerMM = polyCalib.pixelsPerMM;
                obj.cameraCalib.Gx = Gx;
                obj.cameraCalib.Gy = Gy;

                obj.cameraCalib.GSx = polyCalib.GSx;
                obj.cameraCalib.GSy = polyCalib.GSy;

                mmPerPixelX = polyCalib.scaling.X.Slope;
                mmPerPixelY = polyCalib.scaling.Y.Slope;
                offSetX = polyCalib.scaling.X.Offset;
                offSetY = polyCalib.scaling.Y.Offset;

                obj.worldCoords.scaling.mmPerPixelOrg = abs(mmPerPixelX);
                obj.worldCoords.scaling.X.Slope = mmPerPixelX;
                obj.worldCoords.scaling.X.SlopeOrg = mmPerPixelX;
                obj.worldCoords.scaling.X.Unit = 'mm';
                obj.worldCoords.scaling.Y.Slope = mmPerPixelY;
                obj.worldCoords.scaling.Y.SlopeOrg = mmPerPixelY;
                obj.worldCoords.scaling.Y.Unit = 'mm';
                obj.worldCoords.scaling.X.Offset = offSetX;
                obj.worldCoords.scaling.Y.Offset = offSetY;
            end
            obj.cameraCalib.type = calibType;
        end

        obj.cameraCalib.path = calibPath;
        obj.cameraCalib.camNum = camNumber;
    end

    function storePinholeScaling(obj, scalingStruct)
        if obj.procOpts.imgRotAngle == 0
            obj.worldCoords.scaling.X.Offset = scalingStruct.X.Offset;
            obj.worldCoords.scaling.X.Slope = scalingStruct.X.Slope;
            obj.worldCoords.scaling.X.SlopeOrg = scalingStruct.X.Slope;
            obj.worldCoords.scaling.Y.Offset = scalingStruct.Y.Offset;
            obj.worldCoords.scaling.Y.Slope = scalingStruct.Y.Slope;
            obj.worldCoords.scaling.Y.SlopeOrg = scalingStruct.Y.Slope;
        else
            obj.worldCoords.scaling.X.Offset = NaN;
            obj.worldCoords.scaling.X.Slope = NaN;
            obj.worldCoords.scaling.X.SlopeOrg = NaN;
            obj.worldCoords.scaling.Y.Offset = NaN;
            obj.worldCoords.scaling.Y.Slope = NaN;
            obj.worldCoords.scaling.Y.SlopeOrg = NaN;
        end
        obj.worldCoords.scaling.mmPerPixelOrg = scalingStruct.mmPerPixel;
    end

    function rectifyRefImage(obj)
        obj.imgData.refRaw = obj.rectify(obj.imgData.refRaw);
        obj.updateGeometry();
    end

    function imgRectified = rectify(obj, img)
        if ~isempty(obj.procOpts.burntPixelMask)
            img = obj.interpBurntPixels(img);
        end

        if strcmpi(obj.cameraCalib.type, 'polynomial')
            imgRectified = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, double(img), obj.cameraCalib.GSx, obj.cameraCalib.GSy, 'linear');
            imgRectified(isnan(imgRectified)) = 0;
        elseif strcmpi(obj.cameraCalib.type, 'pinhole')
            imgRectified = undistortImage(img, obj.cameraCalib.intrinsicsMatlab);
            imgRectified = imwarp(imgRectified, obj.cameraCalib.pTransform);
        else            % no camera calibration
            imgRectified = img;
        end
    end

    function checkBurntPixelMask(obj, mask)
    % Throw an error if the size of the burnt pixel mask does not match the image.
        % size [rows cols] of the camera image, from the calibration.
        switch lower(obj.cameraCalib.type)
            case "pinhole"
                sz = obj.cameraCalib.intrinsics.ImageSize;
            case "polynomial"
                sz = flip(obj.cameraCalib.orgSize);     % stored as [width height]
            otherwise
                if ~isempty(obj.imgData.refRaw)
                    sz = size(obj.imgData.refRaw);
                else
                    sz = [];
                end
        end

        if ~isempty(mask) && ~isempty(sz) && ~isequal(size(mask), sz)
            error("FtpSolver:burntPixelMaskSize", ...
                "burntPixelMask is %dx%d but the camera image is %dx%d.", ...
                size(mask), sz)
        end
    end

    function imgOut = interpBurntPixels(obj, img)
        if strcmpi(obj.procOpts.burntPixelInterpMethod, 'fastLinear')
            imgOut = obj.fastFill(img, obj.procOpts.burntPixelRows, ...
                obj.procOpts.burntPixelCols);
        elseif strcmpi(obj.procOpts.burntPixelInterpMethod, 'smooth')
            imgOut = regionfill(img, obj.procOpts.burntPixelMask);
        end
    end

    function updateGeometry(obj)
    % Recompute all resizeFactor-derived quantities
        rf = obj.procOpts.resizeFactor;
        rf_o = obj.procOpts.resizeFactorOutput;

        obj.fringe.period = obj.fringe.periodOrg * rf;

        % Morlet wavelet scales (Omega0 = 6) in the COMP frame, from the
        % wavelength range given in the RAW frame
        if ~isempty(obj.demodOpts.waveletWavelengths) ...
                && ~isempty(obj.demodOpts.waveletDivisions)
            scaleRange = obj.demodOpts.waveletWavelengths * 6/(2*pi) * rf;
            obj.demodOpts.waveletScales = linspace(scaleRange(1), scaleRange(2), ...
                                                obj.demodOpts.waveletDivisions);
        end

        if ~isempty(obj.worldCoords.scaling)
            obj.worldCoords.scaling.X.Slope = obj.worldCoords.scaling.X.SlopeOrg / rf;
            obj.worldCoords.scaling.Y.Slope = obj.worldCoords.scaling.Y.SlopeOrg / rf;
            obj.worldCoords.scaling.mmPerPixel = obj.worldCoords.scaling.mmPerPixelOrg / rf;
        end

        % an empty margin is filled using the fringe period
        % one period across the fringes and 0.2 periods along them
        p = obj.fringe.periodOrg;
        if isempty(obj.procOpts.unwrapMarginOrg) && ~isempty(p) && ~isnan(p)
            mrg = ceil([1, 0.2] * p);       % [across fringe, along fringe]
            if strcmpi(obj.fringe.normalAxis, 'Y')
                mrg = flip(mrg);
            end
            obj.procOpts.unwrapMarginOrg = mrg;         % [mX mY], RAW
        end
        obj.procOpts.unwrapMargin = floor(obj.procOpts.unwrapMarginOrg * rf);   % COMP

        if isempty(obj.procOpts.compRectOrg)
            [Ny, Nx] = size(obj.imgData.refRaw);
            obj.procOpts.compRectOrg = [1, 1, Nx - 1, Ny - 1];
        end

        if ~isempty(obj.imgData.refRaw)
            obj.imgData.refScaled = FtpSolver.shrinkArr(obj.imgData.refRaw, rf);
        end

        compRectOrg = round(obj.procOpts.compRectOrg);
        if isfield(obj.imgData, 'refRaw') && ~isempty(obj.imgData.refRaw)
            [NyRaw, NxRaw] = size(obj.imgData.refRaw);
            compRectOrg(1) = max(1, compRectOrg(1));
            compRectOrg(2) = max(1, compRectOrg(2));
            if (compRectOrg(1) + compRectOrg(3)) > NxRaw
                compRectOrg(3) = NxRaw - compRectOrg(1);
            end
            if (compRectOrg(2) + compRectOrg(4)) > NyRaw
                compRectOrg(4) = NyRaw - compRectOrg(2);
            end

            obj.procOpts.compRectOrg = compRectOrg;
        end

        [NyScaled, NxScaled] = size(obj.imgData.refScaled);
        compRect = round(compRectOrg * rf);
        compRect(1:2) = max(1, compRect(1:2));
        compRect(3) = min(NxScaled - compRect(1), compRect(3));
        compRect(4) = min(NyScaled - compRect(2), compRect(4));
        obj.procOpts.compRect = compRect;

        % the dimensions of the resized and cropped reference image is
        % required to initialize arrays
        % obj.imgData.ref is reinitialized in preprocessRefImage()
        obj.imgData.ref = imcrop(obj.imgData.refScaled, obj.procOpts.compRect);
        stackSize = FtpSolver.shrunkArrSize(size(obj.imgData.ref), rf_o);

        % An empty output rectangle is filled one unwrap margin inside the
        % computational domain; a non-empty one is left as it is. While the
        % margin is unknown (fringe period not yet auto-detected) it stays
        % empty and the full domain stands in for the derived rectangle.
        if isempty(obj.procOpts.outputRectOrg) && ~isempty(obj.procOpts.unwrapMarginOrg)
            obj.procOpts.outputRectOrg = obj.marginOutputRect(0);
        end
        outputRectOrg = obj.procOpts.outputRectOrg;
        if isempty(outputRectOrg)
            outputRectOrg = [1 1 compRectOrg(3) compRectOrg(4)];
        end
        
        s = rf*rf_o;
        firstInd = ceil(FtpSolver.rescaleCoord(outputRectOrg(1:2), s));
        lastInd  = floor(FtpSolver.rescaleCoord(outputRectOrg(1:2) + outputRectOrg(3:4), s));
        outputRect = [firstInd, lastInd - firstInd];

        outputRect(1:2) = max(1, outputRect(1:2));
        outputRect(3) = min(stackSize(2) - outputRect(1), outputRect(3));
        outputRect(4) = min(stackSize(1) - outputRect(2), outputRect(4));
        obj.procOpts.outputRect = outputRect;

        obj.placeCorrectionLines();

        [Ny, Nx] = size(obj.imgData.refRaw);
        % initialize coordinate mesh
        if isfield(obj.worldCoords.mesh, 'xOrg')
            xOrg = obj.worldCoords.mesh.xOrg;
            yOrg = obj.worldCoords.mesh.yOrg;
        else
            x0 = obj.worldCoords.scaling.X.Offset;
            y0 = obj.worldCoords.scaling.Y.Offset;
            xSlopeOrg = obj.worldCoords.scaling.X.SlopeOrg;
            ySlopeOrg = obj.worldCoords.scaling.Y.SlopeOrg;

            [xOrg, yOrg] = meshgrid(x0 + xSlopeOrg*(0:Nx-1), ...
                                    y0 + ySlopeOrg*(0:Ny-1));
        end
        
        xComp = FtpSolver.shrinkArr(xOrg, rf, 'bilinear', Antialiasing=false);
        yComp = FtpSolver.shrinkArr(yOrg, rf, 'bilinear', Antialiasing=false);
        xComp = imcrop(xComp, obj.procOpts.compRect);
        yComp = imcrop(yComp, obj.procOpts.compRect);

        obj.worldCoords.mesh.x = FtpSolver.shrinkArr(xComp, rf_o, 'bilinear', ...
            Antialiasing=false);
        obj.worldCoords.mesh.y = FtpSolver.shrinkArr(yComp, rf_o, 'bilinear', ...
            Antialiasing=false);

        % check if pinhole calibration model is available
        if isfield(obj.cameraCalib, 'pTransform')
            [xMeshPixel, yMeshPixel] = meshgrid(1:Nx, 1:Ny);
            pTransform = obj.cameraCalib.pTransform;
            RB = obj.cameraCalib.scaling.imageRB;
            u = xMeshPixel + RB.XWorldLimits(1) - 0.5;
            v = yMeshPixel + RB.YWorldLimits(1) - 0.5;
            pixelsOrg = pTransform.transformPointsInverse([u(: ) v(:)]);

            xPixelOrg = reshape(pixelsOrg(:,1), size(u));
            yPixelOrg = reshape(pixelsOrg(:,2), size(v));

            xPixelComp = FtpSolver.shrinkArr(xPixelOrg, rf, 'bilinear', Antialiasing=false);
            yPixelComp = FtpSolver.shrinkArr(yPixelOrg, rf, 'bilinear', Antialiasing=false);

            xPixelComp = imcrop(xPixelComp, obj.procOpts.compRect);
            yPixelComp = imcrop(yPixelComp, obj.procOpts.compRect);

            obj.worldCoords.mesh.xPixel = FtpSolver.shrinkArr(xPixelComp, rf_o, ...
                'bilinear', Antialiasing=false);
            obj.worldCoords.mesh.yPixel = FtpSolver.shrinkArr(yPixelComp, rf_o, ...
                'bilinear', Antialiasing=false);
        end

        % refresh the auto-selected phase-correction peak
        obj.refreshAutoPeak();
    end

    function rect = marginOutputRect(obj, offset)
    % Output rectangle placed unwrap margin + offset inside the
    % computational domain, in RAW pixels relative to compRectOrg.
    % Along an axis too short to hold the inset (1-D data) no inset is
    % applied.
        inset  = obj.procOpts.unwrapMarginOrg + offset;     % [x y]
        extent = obj.procOpts.compRectOrg(3:4);
        inset  = inset .* (extent + 1 > 2*inset);
        rect   = [inset + 1, extent - 2*inset];
    end

    function placeCorrectionLines(obj)
        compRectOrg = obj.procOpts.compRectOrg;
        numSamples = obj.phaseCorrOpts.numSamples;
        if strcmpi(obj.fringe.normalAxis, 'Y')
            targetInds = round(linspace(compRectOrg(1), compRectOrg(1) + compRectOrg(3), numSamples + 2));
            obj.phaseCorrOpts.lineInds = unique(targetInds(2:end - 1))';
        else                                                    % normal vector to pattern is parallel to X axis
            targetInds = round(linspace(compRectOrg(2), compRectOrg(2) + compRectOrg(4), numSamples + 2));
            obj.phaseCorrOpts.lineInds = unique(targetInds(2:end - 1))';
        end
    end

    function [rowInds, colInds] = outputInds(obj)
    % Row and column indices of the OUTPUT domain specified in STACK-frame pixel coordinates.
        r = round(obj.procOpts.outputRect);
        [Ny, Nx] = size(obj.worldCoords.mesh.x);
        rowInds = r(2) : min(Ny, r(2) + r(4));
        colInds = r(1) : min(Nx, r(1) + r(3));
    end

    function output = RAW2STACK(obj, input)
    % Convert (row, col) values from RAW to STACK frame.
        rf = obj.procOpts.resizeFactor;
        rf_o = obj.procOpts.resizeFactorOutput;
        compRect = obj.procOpts.compRect;

        rows = input(:,1);
        cols = input(:,2);

        % RAW to SCALED
        rowsScaled = FtpSolver.rescaleCoord(rows, rf);
        colsScaled = FtpSolver.rescaleCoord(cols, rf);

        % SCALED to COMP
        rowsCrop = rowsScaled - compRect(2) + 1;
        colsCrop = colsScaled - compRect(1) + 1;

        % COMP to STACK
        rowsStack = FtpSolver.rescaleCoord(rowsCrop, rf_o);
        colsStack = FtpSolver.rescaleCoord(colsCrop, rf_o);

        output = round([rowsStack, colsStack]);
    end

    function createResultFile(obj, quantity, dims)
    % Create an empty result file for one quantity, replacing any old one.
    %   The dataset is single precision and is not compressed. Attributes
    %   record the units, the dataset frame index of the file's first
    %   frame, the number of frames written so far and the dimension order.
        folder = obj.caseFolder();
        if ~isfolder(folder)
            mkdir(folder)
        end
        path = obj.resultFilePath(quantity);
        if isfile(path)
            delete(path)
        end

        dataset = "/" + quantity;
        h5create(path, dataset, dims, 'Datatype', 'single');

        if quantity == "elevation"
            units = "mm";
        else
            units = "rad";
        end
        h5writeatt(path, dataset, "units", units);
        h5writeatt(path, dataset, "firstTimestep", obj.procOpts.solveRange(1));
        h5writeatt(path, dataset, "framesWritten", 0);
        h5writeatt(path, dataset, "dimensionOrder", "MATLAB: [row, column, frame]. " + ...
            "Row-major readers such as h5py: [frame, column, row].");
    end

    function writeFrames(obj, quantity, frameBlock, firstFrame)
    % Write a block of frames into a result file and record the progress.
    %  Writes the N frames in frameBlock to frames firstFrame to
    %  firstFrame + N - 1 of the result file, which must already exist
    %  (createResultFile).
        path = obj.resultFilePath(quantity);
        dataset = "/" + quantity;
        h5write(path, dataset, frameBlock, [1 1 firstFrame], size(frameBlock, 1:3));
        h5writeatt(path, dataset, "framesWritten", firstFrame + size(frameBlock, 3) - 1);
    end

    function checkResultFile(obj, quantity)
    % Warn if a result file is missing or holds fewer frames than solveRange.
        path = obj.resultFilePath(quantity);
        if ~isfile(path)
            warning("Result file %s is missing. Run solve() first.", path)
            return
        end
        nFrames = diff(obj.procOpts.solveRange) + 1;
        framesWritten = h5readatt(path, "/" + quantity, "framesWritten");
        if framesWritten < nFrames
            warning("%s holds %d of %d frames; the run did not finish.", ...
                path, framesWritten, nFrames)
        end
    end

    function quantities = storedQuantities(obj)
    % Names of the result quantities written to file: "elevation", "phase".
        quantities = strings(1, 0);
        if obj.outputConfig.elevStorage
            quantities(end + 1) = "elevation";
        end
        if obj.outputConfig.phaseStorage
            quantities(end + 1) = "phase";
        end
    end

    function folder = caseFolder(obj)
    % Folder that this case's files are written to: <saveFolder>/<caseID>.
        if strlength(obj.outputConfig.saveFolder) == 0
            error("No save folder set. Call setSaveFolder(folder) first.")
        end
        folder = fullfile(obj.outputConfig.saveFolder, obj.caseID);
    end

    function path = resultFilePath(obj, quantity)
    % Path of the result file of one quantity.
        path = fullfile(obj.caseFolder(), obj.caseID + "_" + quantity + ".h5");
    end

    function writeCoordinateMesh(obj)
    % Save the world coordinates (mm) of the result arrays.
        [xMesh, yMesh] = obj.outputCoordinates();
        save(fullfile(obj.caseFolder(), "coordinateMesh.mat"), 'xMesh', 'yMesh')
    end

    function writeDiagnostics(obj)
    % Save diagnostics to a .mat file.
        flaggedRegions = obj.diagnostics.flaggedRegions;
        phaseCorrectionFailures = obj.diagnostics.phaseCorrectionFailures;
        rectFormat = ['Each row of rects is [x y width height] in pixel ' ...
            'coordinates of the result arrays (x = column, y = row), in the ' ...
            'regionprops BoundingBox convention: a box covers columns ' ...
            'x+0.5 to x+width-0.5 and rows y+0.5 to y+height-0.5. timestep ' ...
            'is the frame index in the dataset; frame k of the result files ' ...
            'is timestep firstTimestep + k - 1.'];
        candidateFormat = ['candidates in phaseCorrectionFailures holds ' ...
                           'whole-2*pi jump counts, one per sample '...
                           'point (temporal) or per correction line ' ...
                           '(spatial, in lineInds order), with NaN ' ...
                           'where unusable.'];
        save(fullfile(obj.caseFolder(), obj.caseID + "_diagnostics.mat"), ...
                            'flaggedRegions', 'rectFormat', ...
                            'phaseCorrectionFailures', 'candidateFormat')
    end

    function writeLog(obj)
    % Write the solver settings and a run summary to <caseID>_log.txt.
        sections = cell(0, 2);

        % ---- run ----
        solveRange = obj.procOpts.solveRange;
        nFrames = diff(solveRange) + 1;
        solveTime = obj.diagnostics.solveTime;
        nFlaggedFrames = numel(obj.diagnostics.flaggedRegions);
        entries = {
            "First and last frame solved",   solveRange
            "Number of frames",              nFrames
            "Processing time (s)",           round(solveTime, 1)
            };
        if ~isempty(solveTime)
            entries(end + 1, :) = {"Time per frame (s)", ...
                            round(solveTime/nFrames, 3, "significant")};
        end
        entries(end + 1, :) = {"Frames with phase correction failures", ...
            numel(obj.diagnostics.phaseCorrectionFailures)};
        entries(end + 1, :) = {"Frames with flagged regions", nFlaggedFrames};
        sections(end + 1, :) = {"Run", entries};

        % ---- input ----
        entries = {
            "Reference image",          obj.source.refPath
            "Data",                     obj.source.dataPath
            "Data format",              obj.source.dataFormat
            "Frames available",         obj.source.fullRange
            "Data camera index",        obj.source.dataCamInd
            "Reference camera index",   obj.source.refCamInd
            };
        sections(end + 1, :) = {"Input", entries};

        % ---- camera calibration ----
        calib = obj.cameraCalib;
        entries = {"Type", calib.type};
        if ~strcmpi(calib.type, "none")
            entries(end + 1, :) = {"File", calib.path};
            if endsWith(calib.path, ".xml", 'IgnoreCase', true)
                entries(end + 1, :) = {"Camera number", calib.camNum};
            end
        end
        entries(end + 1, :) = {"Image rotation (deg)", obj.procOpts.imgRotAngle};
        entries(end + 1, :) = {"mm per pixel (RAW)", obj.worldCoords.scaling.mmPerPixelOrg};
        sections(end + 1, :) = {"Camera calibration", entries};

        % ---- fringe pattern ----
        entries = {
            "Period (RAW px)",      obj.fringe.periodOrg
            "Normal axis",          obj.fringe.normalAxis
            "Tilt angle (deg)",     obj.fringe.tiltAngle
            };
        sections(end + 1, :) = {"Fringe pattern", entries};

        % ---- geometry ----
        entries = {
            "Computational domain [x y w h] (RAW px)",  obj.procOpts.compRectOrg
            "Output domain [x y w h] (RAW px, relative to computational domain)", ...
                obj.procOpts.outputRectOrg
            "Resize factor",                             obj.procOpts.resizeFactor
            "Output resize factor",                      obj.procOpts.resizeFactorOutput
            };
        if isstruct(obj.worldCoords.mesh) && isfield(obj.worldCoords.mesh, 'x')
            [rowInds, colInds] = obj.outputInds();
            entries(end + 1, :) = {"Result array size [rows cols]", ...
                [numel(rowInds), numel(colInds)]};
        end
        sections(end + 1, :) = {"Geometry", entries};

        % ---- preprocessing ----
        nBurnt = nnz(obj.procOpts.burntPixelMask);
        entries = {
            "DC subtraction",                       obj.procOpts.dcSubtraction
            "Reference image smoothing",            obj.procOpts.refImgSmoothing
            "Discard threshold (mean intensity)",   obj.procOpts.discardThreshold
            "Burnt pixels",                         nBurnt
            };
        if nBurnt > 0
            entries(end + 1, :) = {"Burnt pixel fill method", obj.procOpts.burntPixelInterpMethod};
        end
        entries(end + 1, :) = {"Parallax correction", obj.procOpts.parallaxCorrection};
        sections(end + 1, :) = {"Preprocessing", entries};

        % ---- demodulation ----
        demod = obj.demodOpts;
        entries = {"Method", demod.method};
        if strcmpi(demod.method, "fourier")
            entries(end + 1, :) = {"Filter width (fraction of carrier)", demod.filterWidthFrac};
        else
            entries = [entries; {
                "Wavelengths [min max] (RAW px)",   demod.waveletWavelengths
                "Divisions",                        demod.waveletDivisions
                "Scales (COMP px)",                 demod.waveletScales
                "Angles (deg)",                     demod.waveletAngles
                "Morlet sigma",                     demod.waveletSigma
                "Morlet gamma",                     demod.waveletGamma
                }];
        end
        sections(end + 1, :) = {"Demodulation", entries};

        % ---- unwrapping ----
        entries = {
            "Method",                   obj.procOpts.unwrapMethod
            "Margin [x y] (RAW px)",    obj.procOpts.unwrapMarginOrg
            };
        % fault detection runs only when the phase is unwrapped
        if obj.procOpts.unwrapMethod ~= "none"
            faultSettings = obj.procOpts.faultDetection;
            entries = [entries; {
                "Fault median window (periods)",        faultSettings.medianWindowPeriods
                "Fault residual threshold (x median)",  faultSettings.residualFactor
                "Fault length to flag (periods)",       faultSettings.lengthToFlagPeriods
                "Flag padding (periods)",               faultSettings.dilationPeriods
                "Flag merge radius (periods)",          faultSettings.closePeriods
                }];
            if obj.procOpts.unwrapMethod == "2D"
                entries(end + 1, :) = {"Fault length for unwrap2D (periods)", ...
                    faultSettings.lengthFor2DPeriods};
            end
        end
        sections(end + 1, :) = {"Unwrapping", entries};

        % ---- phase correction ----
        pCorr = obj.phaseCorrOpts;
        entries = {"Method", pCorr.method};
        if strcmpi(pCorr.method, "temporal")
            entries(end + 1, :) = {"Sample points", pCorr.numSamples};
        elseif strcmpi(pCorr.method, "spatial")
            entries = [entries; {
                "Correction lines requested",                   pCorr.numSamples
                "Tracked peak",                                 pCorr.trackedPeakInd
                "Start edge",                                   pCorr.startEdge
                "Tracked peak clearance (periods)",             pCorr.peakClearance
                "Minimum peak height (fraction of median)",     pCorr.peakMinFrac
                "Minimum peak prominence (fraction of range)",  pCorr.peakPromFrac
                "Correction lines (RAW px)",                    pCorr.lineInds
                }];
        end
        sections(end + 1, :) = {"Phase correction", entries};

        % ---- elevation model ----
        model = obj.elevModel;
        entries = {"Type", model.type};
        if strcmpi(model.type, "takeda")
            entries = [entries; {"L (mm)", model.L; "d (mm)", model.d}];
        elseif strcmpi(model.type, "polynomial")
            if isfield(model, 'path')
                entries(end + 1, :) = {"Calibration file", model.path};
            else
                entries(end + 1, :) = {"Calibration file", "none (fitted in this session, not saved)"};
            end
            if model.poly.pixelwise
                entries(end + 1, :) = {"Coefficients", "pixelwise"};
            else
                entries(end + 1, :) = {"Coefficients", "smoothed (poly22 fits)"};
            end
            if ~isempty(model.poly.coeffMatOrg)
                entries(end + 1, :) = {"Polynomial order", size(model.poly.coeffMatOrg, 3) - 1};
            end
        end
        sections(end + 1, :) = {"Elevation model", entries};

        % ---- output ----
        entries = {
            "Elevation stored",     obj.outputConfig.elevStorage
            "Phase stored",         obj.outputConfig.phaseStorage
            "Streamed to disk",     obj.outputConfig.streaming
            };
        if obj.outputConfig.streaming
            entries(end + 1, :) = {"Buffer length (frames)", obj.outputConfig.bufferLength};
        end
        entries(end + 1, :) = {"Save folder", obj.outputConfig.saveFolder};
        sections(end + 1, :) = {"Output", entries};

        % ---- format, names aligned within each section ----
        solverFile = string(which("FtpSolver"));
        solverInfo = dir(solverFile);
        txt = [
            "FtpSolver log, case " + obj.caseID
            "Written:     " + string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"))
            "MATLAB:      " + string(version)
            "Solver file: " + solverFile + " (modified " ...
                + string(datetime(solverInfo.datenum, "ConvertFrom", "datenum", ...
                "Format", "yyyy-MM-dd HH:mm")) + ")"
            ];
        for s = 1:size(sections, 1)
            txt(end + 1) = "";
            txt(end + 1) = "[" + sections{s, 1} + "]";
            entries = sections{s, 2};
            nameWidth = max(strlength(string(entries(:, 1)))) + 2;
            for e = 1:size(entries, 1)
                txt(end + 1) = "  " + pad(entries{e, 1} + ":", nameWidth) ...
                    + FtpSolver.logValue(entries{e, 2});
            end
        end

        writelines(txt, fullfile(obj.caseFolder(), obj.caseID + "_log.txt"));
    end

    function [phaseCrop, heightVec, weights, geom] = ...
            prepareCalibrationData(obj, heightVec, excludeHeights, weights)
    % Check and trim the calibration data and compute the calibration-region geometry.
    %
    %   Outputs:
    %     phaseCrop - phase stack restricted to the calibration region and
    %                 to the retained heights (range applied, excluded and
    %                 zero heights removed)
    %     heightVec - matching column vector of heights
    %     weights   - matching column vector of fit weights
    %     geom      - geometry struct: calibRect (crop frame), rawRect (RAW
    %                 frame), embedRows/embedCols, pixel meshes, compRect,
    %                 rawSize

        if obj.procOpts.resizeFactor ~= 1 || obj.procOpts.resizeFactorOutput ~= 1
            error("Run solve() with resizeFactor and resizeFactorOutput set to 1")
        end
        if size(obj.phaseData.stack, 3) ~= length(heightVec)
            error("Size of phase array incompatible with length of heightVec.")
        end

        heightVec = heightVec(~excludeHeights);
        heightVec = heightVec(:);
        weights   = weights(~excludeHeights);
        weights   = weights(:);

        % crop interior inside the unwrap margin
        mg       = obj.procOpts.unwrapMargin;      % [mX mY], COMP (= RAW when rf = 1)
        compRect = obj.procOpts.compRect;

        geom.compRect  = compRect;
        geom.calibRect = [mg + 1, compRect(3:4) - 2*mg];              % crop frame
        geom.rawRect   = [compRect(1:2) + mg, compRect(3:4) - 2*mg];  % RAW frame
        geom.embedRows = geom.rawRect(2) : geom.rawRect(2) + geom.rawRect(4);
        geom.embedCols = geom.rawRect(1) : geom.rawRect(1) + geom.rawRect(3);

        [Ny, Nx] = size(obj.imgData.refRaw);
        geom.rawSize = [Ny, Nx];
        [geom.pMeshX, geom.pMeshY] = meshgrid(1:Nx, 1:Ny);
        geom.pMeshX_crop = imcrop(geom.pMeshX, geom.rawRect);
        geom.pMeshY_crop = imcrop(geom.pMeshY, geom.rawRect);

        % crop the phase stack to region
        calibRect3 = [geom.calibRect(1:2), 1, geom.calibRect(3:4), ...
                      size(obj.phaseData.stack, 3) - 1];
        phaseCrop  = imcrop3(obj.phaseData.stack, calibRect3);
        phaseCrop(:,:,excludeHeights) = [];

        % drop zero height (reference) if present
        zeroInd = find(heightVec == 0);
        heightVec(zeroInd) = [];
        weights(zeroInd)   = [];
        phaseCrop(:,:,zeroInd) = [];
    end

    function cLims = plotCalibrationDiagnostics(obj, phaseCrop, elevCrop, ...
            fitPhase, fitHeight, heightVec, calibRect, figTitle, cLims)
        % Plot calibration diagnostics for one set of coefficients.
        %
        %   Layout (2x4, column-major): mean absolute error vs height with std
        %   bars; center-pixel phase-height data and fitted curve; four error
        %   maps at evenly spaced heights. Pass cLims = {} on the first call to
        %   capture color limits; pass the returned cLims on the second call so
        %   both figures share them.

        % errors over the calibration region
        errorMat = elevCrop - reshape(heightVec, 1, 1, []);
        errVec   = squeeze(mean(abs(errorMat), [1 2]));
        stdVec   = squeeze(std(errorMat, 0, [1 2]));

        % center-pixel fit curve
        midY = floor((size(phaseCrop, 1) + 1) / 2);
        midX = floor((size(phaseCrop, 2) + 1) / 2);
        phaseSample = squeeze(phaseCrop(midY, midX, :));

        % world coordinates of the calibration region
        X = imcrop(obj.worldCoords.mesh.x, calibRect);
        Y = imcrop(obj.worldCoords.mesh.y, calibRect);

        % plot
        figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
        tHndl = tiledlayout(2, 4, 'TileIndexing', 'columnmajor');
        txt = title(tHndl, figTitle);
        txt.FontSize = 18;

        % Mean absolute error vs height
        nexttile([1 2])
        errorbar(heightVec, errVec, stdVec, '-o', 'linewidth', 1.5);
        xlim([min(heightVec) - 2, max(heightVec) + 2])
        xlabel("Height (mm)")
        ylabel("Mean absolute error (mm)")
        set(gca, 'fontsize', 15)

        % Center-pixel calibration curve
        nexttile([1 2])
        plot(phaseSample, heightVec, 'x', 'linewidth', 1.5)
        hold on
        plot(fitPhase, fitHeight, 'linewidth', 1)
        hold off
        xlim([min(fitPhase) - 1, max(fitPhase) + 1])
        xlabel("\Delta\phi (rad)")
        ylabel("Height (mm)")
        set(gca, 'fontsize', 15)
        legend('Measurements', 'Model', 'Location', "best")

        % Error maps at four heights, color limits shared across figures
        nMaps = min(4, numel(heightVec));
        showErrInd = unique(round(linspace(1, numel(heightVec), nMaps)));
        applyLims = ~isempty(cLims);
        for k = 1:numel(showErrInd)
            nexttile
            imagesc(X(1,:), Y(:,1), errorMat(:,:, showErrInd(k)))
            set(gca, 'fontsize', 15)
            cbHndl = colorbar;
            cbHndl.Label.String = 'Error (mm)';
            if applyLims
                clim(cLims{k});
            else
                cLims{k} = cbHndl.Limits;
            end
            xlabel('X (mm)')
            ylabel('Y (mm)')
            title(sprintf('H = %.2f mm', heightVec(showErrInd(k))))
        end
    end

    function [dataArray, xMesh, yMesh, meshUnits] = ...
                    animationData(obj, firstTimestep, lastTimestep, opts)
        % Frames and coordinates shown by animate and exportVideo.
        %   data     - opts.data of the caller: {array} or {array, xMesh, yMesh};
        %              {} uses the stored result instead
        %   quantity - "elevation" or "phase", the stored result to use when data is {}
        %   An empty firstTimestep or lastTimestep means the first or last
        %   available frame.
        
        data = opts.data;
        quantity = opts.quantity;
        meshUnits = "(mm)";
        if ~isempty(data)
            dataArray = data{1};
            if isempty(firstTimestep)
                firstTimestep = 1;
            end
            if isempty(lastTimestep)
                lastTimestep = size(dataArray, 3);
            end
            dataArray = dataArray(:, :, firstTimestep:lastTimestep);
            if numel(data) > 1
                xMesh = data{2};
                yMesh = data{3};
            else
                [xMesh, yMesh] = meshgrid(1:size(dataArray, 2), 1:size(dataArray, 1));
                meshUnits = "(pixels)";
            end
        end

        if isempty(data)
            if strcmpi(quantity, "elevation")
                stack = obj.elevData.stack;
            else
                stack = obj.phaseData.stack;
            end
            if isempty(firstTimestep)
                firstTimestep = obj.procOpts.solveRange(1);
            end
            if isempty(lastTimestep)
                lastTimestep = obj.procOpts.solveRange(1) + size(stack, 3) - 1;
            end
    
            if obj.outputConfig.streaming
                % only the last buffer is in memory, and it ends at solveRange(2)
                solveRange = obj.procOpts.solveRange;
                bufferStart = solveRange(2) - size(stack, 3) + 1;
                shown = bufferStart + [firstTimestep, lastTimestep] - solveRange(1);
                warning("FtpSolver:streamedFrames", ...
                    "Streaming is on, so only the last buffer (timesteps %d to %d) " + ...
                    "is in memory. Timesteps are counted from the start of that buffer: " + ...
                    "the frames shown are timesteps %d to %d. Use FtpSolver.readData() " + ...
                    "to read other frames from the result file.", ...
                    bufferStart, solveRange(2), shown(1), shown(2))
            end

            dataArray = obj.outputSubarray(quantity, firstTimestep, lastTimestep);
            [xMesh, yMesh] = obj.outputCoordinates();

            if isempty(dataArray)
                error("No %s data to show. Call setStoredOutputs(%s=true) and run " + ...
                    "solve(), and check that firstTimestep <= lastTimestep.", ...
                    quantity, quantity)
            end
        end

        % check for nan frames and interpolate
        if opts.fillNanFrames
            probeVec = dataArray(floor(end/2), floor(end/2), :);
            nanTimes = find(isnan(probeVec));
            nanTimes(nanTimes == 1) = [];
            nanTimes(nanTimes == size(dataArray, 3)) = [];

            for i = 1:length(nanTimes)
                dataArray(:,:,nanTimes(i)) = 0.5*(dataArray(:,:,nanTimes(i) - 1) + dataArray(:,:,nanTimes(i) + 1));
            end
        end

        if ~isempty(opts.subtractConst)
            dataArray = dataArray - opts.subtractConst;
        end

        if opts.subtractPlane
            dataArray = obj.subtractPlane(dataArray);
        end

        if opts.subtractMean
            dataArray = dataArray - mean(dataArray, 3, 'omitmissing');
        end

        if opts.smoothing
            dataArray = imgaussfilt(dataArray, opts.smoothingSigma);
        end
    end

    function overlayDomainRectangles(obj, ax)
        % Draw the domains on axes showing a rectified full-resolution image.
        %  Shades the unwrap margin, whose pixels are NaN in the results,
        %  in blue, and outlines the computational domain (red) and the
        %  output domain (green). Does nothing if no computational domain
        %  is set.
        compRectOrg = obj.procOpts.compRectOrg;
        if isempty(compRectOrg)
            return
        end
    
        textFontSize = 14;
    
        % pixel-edge boundary [left top right bottom] of the computational
        % domain, whose pixel centers run from x0 to x0 + w
        outer = [compRectOrg(1:2) - 0.5, compRectOrg(1:2) + compRectOrg(3:4) + 0.5];
    
        % margin in RAW pixels; none along an axis too short to hold it
        % (1-D data), as in marginOutputRect
        mg = obj.procOpts.unwrapMarginOrg;     % [mX mY]
        if isempty(mg)
            mg = [0 0];
        end
        extent = compRectOrg(3:4);
        mg = mg .* (extent + 1 > 2*mg);
        inner = outer + [mg, -mg];
    
        % margin band: four strips between the outer and inner boundaries
        vertices = [outer(1) outer(2); outer(3) outer(2); outer(3) outer(4); outer(1) outer(4)
                    inner(1) inner(2); inner(3) inner(2); inner(3) inner(4); inner(1) inner(4)];
        faces = [1 2 6 5      % top
                 2 3 7 6      % right
                 3 4 8 7      % bottom
                 4 1 5 8];    % left
        patch(ax, 'Faces', faces, 'Vertices', vertices, 'FaceColor', 'b', ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none');
    
        rectangle(ax, 'Position', [outer(1:2), outer(3:4) - outer(1:2)], ...
            'EdgeColor', 'r', 'LineWidth', 2);
    
        hasOutput = ~isempty(obj.procOpts.outputRectOrg);
        if hasOutput
            outputRectOrg = obj.procOpts.outputRectOrg;
            firstPixel = compRectOrg(1:2) + outputRectOrg(1:2) - 1;   % absolute
            outputEdges = [firstPixel - 0.5, outputRectOrg(3:4) + 1];
            rectangle(ax, 'Position', outputEdges, 'EdgeColor', 'g', 'LineWidth', 1.5);
        end
    
        text(ax, outer(1), outer(2), 'Computational domain ', ...
            'Color', 'r', 'FontSize', textFontSize, 'VerticalAlignment', 'bottom');
        text(ax, inner(3), inner(2), 'Unwrap margin ', ...
            'Color', 'b', 'FontSize', textFontSize, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        if hasOutput
            text(ax, outputEdges(1), outputEdges(2) + outputEdges(4), ...
                ' Output domain ', 'Color', 'g', 'FontSize', textFontSize, ...
                'VerticalAlignment', 'bottom');
        end
end
end

%% ----------- Property initialization and static helpers -----------
methods (Access = private, Static)
    function output = defaultSource()
        output.dataPath   = "";
        output.dataFormat = "";
        output.frameFiles = strings(0, 1);
        output.refPath    = "";
        output.fullRange  = [];
        output.dataCamInd = 1;
        output.refCamInd  = 1;
    end

    function output = defaultImgData()
        output.refRaw = [];
        output.refScaled = [];
        output.ref = [];
        output.curr = [];
    end

    function output = defaultPhaseData()
        output.curr = [];
        output.prev = [];
        output.stack = [];
        output.refComplexCoeffs = [];
        output.bandpassFilter = [];
        output.fourierPadding = [0, 0];
    end

    function output = defaultElevData()
        output.stack = [];
        output.curr = [];
    end

    function output = defaultFringe()
        output.normalVec = [];
        output.normalAxis = "";
        output.periodOrg = [];
        output.period = [];
        output.tiltAngle = [];
    end

    function output = defaultDemodOpts()
        output.method = "fourier";
        output.waveletScales = [];
        output.waveletAngles = 0;
        output.waveletWavelengths = [];   % [min max], RAW px
        output.waveletDivisions = [];
        output.waveletSigma = 0.6;
        output.waveletGamma = 1;
        output.filterWidthFrac = 0.6;
    end

    function output = defaultPhaseCorrOpts()
        output.method = "temporal";
        output.lineInds = [];
        output.startEdge = "left";
        output.numSamples = 9;
        output.peakClearance = 2.5;
        output.trackedPeakInd = [];
        output.peakPromFrac = 1/4;
        output.peakMinFrac = 1/3;
    end

    function output = defaultProcOpts()
        output.solveRange = [];
        output.unwrapMethod = "1D";
        output.unwrapMarginOrg = [];        % [marginX marginY], RAW frame
        output.unwrapMargin = [];           % [marginX marginY], COMP frame
        output.imgRotAngle = 0;             % in degrees and positive clockwise
        output.parallaxCorrection = false;
        output.discardThreshold = [];
        output.dcSubtraction = true;
        output.refImgSmoothing = false;
        output.burntPixelInterpMethod = 'fastLinear';
        output.burntPixelMask = [];
        output.burntPixelRows = [];
        output.burntPixelCols = [];
        output.compRectOrg = [];
        output.compRect = [];
        output.outputRectOrg = [];
        output.outputRect = [];
        output.resizeFactor = 1;
        output.resizeFactorOutput = 0.5;

        % fault detection in unwrapAndFlag
        % lengths are in fringe periods
        output.faultDetection.medianWindowPeriods = 1;    % median filter on the phase gradient magnitude
        output.faultDetection.residualFactor      = 1000; % residual threshold, as a multiple of its median
        output.faultDetection.lengthFor2DPeriods  = 0.6;  % fault length above which unwrap2D runs (method "2D")
        output.faultDetection.lengthToFlagPeriods = 0.2;  % fault length above which faults are recorded
        output.faultDetection.dilationPeriods     = 1;    % padding around recorded faults
        output.faultDetection.closePeriods        = 1.2;  % closing radius that merges nearby faults
    end

    function output = defaultLoopState()
        output.stackInd = 0;
        output.timestep = [];
        output.isFirstTimestep = false;
        output.bufferIndex = 1;
        output.discardFrame = false;
    end

    function output = defaultCameraCalib()
        output.type = "none";
        output.path = "";
        output.X.Unit = 'pixel';
        output.X.SlopeOrg = 1;
        output.X.Offset = 0;
        output.Y.Unit = 'pixel';
        output.Y.SlopeOrg = 1;
        output.Y.Offset = 0;
    end

    function output = defaultWorldCoords()
        output.mesh = [];
        output.scaling.X.SlopeOrg = 1;
        output.scaling.X.Offset = 0;
        output.scaling.Y.SlopeOrg = 1;
        output.scaling.Y.Offset = 0;
        output.scaling.mmPerPixelOrg = 1;
    end

    function output = defaultElevModel()
        output.type = 'none';
    end

    function output = defaultOutputConfig()
        output.phaseStorage = false;
        output.imageStorage = false;
        output.elevStorage = true;
        output.reportProgressEnabled = true;
        output.streaming = false;
        output.bufferLength = 1000;
        output.saveFolder = "";
    end

    function output = defaultDiagnostics()
    output.phaseCorrectionFailures = struct('timestep', {}, 'candidates', {});
        output.flaggedRegions = struct('timestep', {}, 'rects', {});
        output.solveTime = [];
    end

    function out = rescaleCoord(x, scale)
    % Map a coordinate between frames related by imresize.
        % scale about the image edge
        out = (x - 0.5)*scale + 0.5;
    end

    function out = shrinkArr(A, scale, varargin)
    % Shrink an array with imresize, dropping output pixels that lie outside it.
        out = imresize(A, scale, varargin{:});
        n = FtpSolver.shrunkArrSize(size(A, [1 2]), scale);
        out = out(1:n(1), 1:n(2), :);
    end

    function n = shrunkArrSize(sz, scale)
    % Rows and columns returned by shrinkArr for an array of size sz.
        n = max(1, floor(scale*sz + (1 - scale)/2 + 1e-9));
    end

    function imgOut = fastFill(img, interpRows, interpCols)
    % Replace pixels with values interpolated from neighbors.
        imgOut = img;

        for i = 1:length(interpRows)
            r = interpRows(i);
            c = interpCols(i);
            % simple 4-neighbor average
            % assumes bad pixels aren't on the image edge
            imgOut(r, c, :) = (img(r-1, c, :) + img(r+1, c, :) + ...
                img(r, c-1, :) + img(r, c+1, :)) / 4;
        end
    end

    function len = faultLength(faults)
    % Approximate total length, in pixels, of the faults in a mask.
        % a 2*pi fault marks a line two pixels wide.
        len = nnz(faults) / 2;
    end

    function delta = quadPeakOffset(fm, f0, fp)
    % Sub-bin offset of a spectral peak from its neighbors
    %   Fits a parabola through log-magnitudes at bins [-1, 0, +1] and
    %   returns the vertex position in [-0.5, 0.5]
        a = log(max(fm, realmin));
        b = log(max(f0, realmin));
        c = log(max(fp, realmin));
        denom = a - 2*b + c;
        if denom >= 0        % not a local max in log domain
            delta = 0;
        else
            delta = min(0.5, max(-0.5, 0.5*(a - c)/denom));
        end
    end

    function [minHeight, minProm] = peakThresholds(lineVals, peakMinFrac, peakPromFrac)
    % Fringe mean level and peak-to-trough amplitude from the intensity distribution.
        medianLevel = median(lineVals, 'omitmissing');
        amplitude = prctile(lineVals, 95) - prctile(lineVals, 5);
        minHeight = peakMinFrac  * medianLevel;
        minProm   = peakPromFrac * amplitude;
    end

    function [image, scaling] = readDavisFrame(addr, camInd, timestep)
    % Read one frame (or the set average) from a Davis file.
    %   Undoes the zero padding and area-of-interest offset that Davis
    %   applies, using the RealFrameSize/AOIused/CameraMaxNx attributes.
    %   Inputs:
    %       addr     - path to the .set or .im7 file
    %       camInd   - camera index in the buffer
    %       timestep - set index, 0/omitted for a single image, or
    %                   'avg' to average the whole set
    %   Outputs:
    %       image   - image array
    %       scaling - Davis scaling struct of the last frame read

        addr = convertStringsToChars(addr);

        if ~exist('timestep', 'var')
            timestep = 0;
        elseif strcmpi(timestep, 'avg')
            timestep = 1:lvsetsize(addr);
        end

        for i = 1:length(timestep)
            if timestep(i) ~= 0
                temp = readimx(addr, timestep(i));
            else
                temp = readimx(addr);
            end

            image = temp.Frames{camInd}.Components{1}.Planes{1}';

            % Crop image if Davis has padded it with zeros
            rfsIdx = [];
            % reposition image if AOI used
            aoiIdx = [];
            maxNxIdx = [];
            for k = 1:length(temp.Frames{camInd}.Attributes)
                if strcmpi(temp.Frames{camInd}.Attributes{k}.Name, 'RealFrameSize')
                    rfsIdx = k;
                end
                if strcmpi(temp.Frames{camInd}.Attributes{k}.Name, 'AOIused')
                    aoiIdx = k;
                end
                if strcmpi(temp.Frames{camInd}.Attributes{k}.Name, 'CameraMaxNx')
                    maxNxIdx = k;
                end
            end
            if ~isempty(rfsIdx)
                realFrameX = temp.Frames{camInd}.Attributes{rfsIdx}.Value(1);
                realFrameY = temp.Frames{camInd}.Attributes{rfsIdx}.Value(2);
                if ~isequal(size(image), [realFrameY, realFrameX])
                    image = image(1:realFrameY, 1:realFrameX);
                end
            end

            if ~isempty(aoiIdx)
                roiX = temp.Frames{camInd}.Attributes{aoiIdx}.Value(1) + 1;
                roiY = temp.Frames{camInd}.Attributes{aoiIdx}.Value(2) + 1;
                binX = temp.Frames{camInd}.Attributes{aoiIdx}.Value(3);
                binY = temp.Frames{camInd}.Attributes{aoiIdx}.Value(4);
                Nx = str2double(temp.Frames{camInd}.Attributes{maxNxIdx}.Value);
                Ny = str2double(temp.Frames{camInd}.Attributes{maxNxIdx + 1}.Value);

                if binX ~=  1 || binY ~= 1
                    if roiX ~= 1 || roiY ~= 1
                        error("readDavisFrame() cannot handle images " +  ...
                            "that have an area of interest and are also binned.")
                    end
                    imageFullSize = uint16(zeros(realFrameY, realFrameX));
                else
                    imageFullSize = uint16(zeros(Ny, Nx));
                end
                imageFullSize(roiY : roiY +  realFrameY - 1, roiX : roiX + realFrameX - 1) = image;
                image = imageFullSize;
            end

            if i == 1
                imageSum = zeros(size(image)); % initialize on first iteration
            end
            imageSum = imageSum + double(image);
        end

        if length(timestep) > 1
            image = imageSum/length(timestep);
        end

        scaling = temp.Frames{camInd}.Scales;
        scaling.X.SlopeOrg = scaling.X.Slope;
        scaling.Y.SlopeOrg = scaling.Y.Slope;
    end

    function p = absolutePath(p)
    % Full path of an existing file, resolving relative paths and "..".
        info = dir(p);
        if numel(info) ~= 1 || info.isdir
            error("File not found: %s", p)
        end
        p = string(fullfile(info.folder, info.name));
    end

    function s = logValue(v)
    % Format one setting for the log file. Large arrays are summarized.
        if isempty(v)
            s = "none";
        elseif isstring(v) || ischar(v)
            s = strjoin(string(v), ", ");
        elseif islogical(v) && isscalar(v)
            if v
                s = "yes";
            else
                s = "no";
            end
        elseif isnumeric(v) || islogical(v)
            if numel(v) > 100
                s = "array of size " + mat2str(size(v));
            else
                if isvector(v)
                    v = reshape(v, 1, []);
                end
                s = string(mat2str(double(v), 6));
            end
        else
            s = "(" + class(v) + ")";
        end
    end

    function surfaceCoeffs = fitCoeffSurfaces(coeffMatOrg)
    % Fit a second-order surface to each pixelwise coefficient plane.
        % Each plane is fitted by least squares over its finite
        % (calibrated) pixels with 
        %       c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2, 
        % where x and y are the RAW pixel column and row indices. 
        % Input:
        %   coeffMatOrg   - pixelwise coefficients, RAW frame, NaN outside
        %                   the calibrated region
        % Output:
        % surfaceCoeffs - 6 x nPlanes; column i holds
        %                 [p00; p10; p01; p20; p11; p02] of plane i
        [Ny, Nx, nPlanes] = size(coeffMatOrg);
        [xMesh, yMesh] = meshgrid(1:Nx, 1:Ny);
        surfaceCoeffs = zeros(6, nPlanes);
        for i = 1:nPlanes
            c = coeffMatOrg(:,:,i);
            valid = isfinite(c);
            A = FtpSolver.quadBasis(xMesh(valid), yMesh(valid));
            % scale each column to unit length
            colNorms = vecnorm(A);
            surfaceCoeffs(:, i) = ((A ./ colNorms) \ c(valid)) ./ colNorms';
        end
    end

    function coeffPlanes = evalCoeffSurfaces(surfaceCoeffs, xMesh, yMesh)
    % Evaluate the fitted coefficient surfaces at RAW pixel coordinates.
        %   Inputs:
        %       surfaceCoeffs - 6 x nPlanes, from fitCoeffSurfaces
        %       xMesh, yMesh  - RAW pixel column and row indices
        %   Output: 
        %       coeffPlanes   - size(xMesh) x nPlanes
        coeffPlanes = reshape(FtpSolver.quadBasis(xMesh(:), yMesh(:)) * surfaceCoeffs, ...
                              [size(xMesh), size(surfaceCoeffs, 2)]);
    end

    function A = quadBasis(x, y)
    % Terms of the second-order surface, one column each: 1, x, y, x^2, x*y, y^2.
        A = [ones(size(x)), x, y, x.^2, x.*y, y.^2];
    end
end

%% ----------- Standalone functions -----------
methods (Static)
    function img = subtractDC(img, period)
    % Remove the low-frequency (DC) component from an image.
    %   Applies a double box filter of width ~4 fringe periods and
    %   subtracts the result.
    %
        filterSize = 2*floor(2*period) + 1;
        dc = imboxfilt(img, filterSize, 'padding', 'symmetric');
        dc = imboxfilt(dc, filterSize, 'padding', 'symmetric');
        img = double(img) - double(dc);
    end

    function data = subtractPlane(data, downscaleFactor)
    % Fit and subtract a plane from each timestep slice of data array.
    % Downscales data array in spatial dimensions by downscaleFactor.

        if ~exist('downscaleFactor', 'var')
            downscaleFactor = 0.2;
        end
        [Ny, Nx, Nt] = size(data);
        [xGrid, yGrid] = meshgrid(1:Nx, 1:Ny);
        xGrid_coarse = imresize(xGrid, downscaleFactor, 'bilinear', Antialiasing=false);
        yGrid_coarse = imresize(yGrid, downscaleFactor, 'bilinear', Antialiasing=false);

        for i = 1:Nt
            % fit plane to data and subtract
            temp = double(imresize(data(:,:,i), downscaleFactor, ...
                          'bilinear', Antialiasing=false));
            [xOut, yOut, zOut] = prepareSurfaceData(xGrid_coarse, yGrid_coarse, temp);
            if numel(zOut) < 3
                continue
            end
            sf = fit([xOut, yOut], zOut, 'poly11');
            fittedPlane = sf.p00 + xGrid*sf.p10 + yGrid*sf.p01;
            data(:,:,i) = data(:,:,i) - fittedPlane;
        end
    end

    function [period, normalAxis, tiltAngle, normalVec] = analyzeFringe(img)
    % Estimate fringe period and orientation via FFT
    %
    %   Applies a Hann window to the image and locates the dominant
    %   spectral peak outside the low-frequency core.
    %
    %   period      - fringe period in pixels
    %   normalAxis  - 'X' if the pattern normal is closer to the x axis,
    %                 'Y' otherwise
    %   tiltAngle   - signed angle in degrees, in [-45, 45], between the
    %                 pattern normal (carrier wave vector) and the selected
    %                 axis; 0 when the fringes are perfectly aligned.
    %                 Equivalently, the tilt of the fringe lines away from
    %                 the perpendicular axis. Positive when the normal is
    %                 rotated from the selected positive axis toward the
    %                 other positive axis.
    %   normalVec   - fringe normal unit vector [vx, vy]
    %
        img = double(img);
        [ny, nx] = size(img);

        img = FtpSolver.subtractDC(img, max(ny,nx)/8);

        % suppress spectral leakage
        win = hann(ny) * hann(nx)';
        F = abs(fft2(img .* win));

        % integer frequency grids: bin value = cycles across the crop
        px = [0:floor(nx/2), -ceil(nx/2)+1:-1];
        py = [0:floor(ny/2), -ceil(ny/2)+1:-1];
        [PX, PY] = meshgrid(px, py);

        % search one half-plane and exclude the low-frequency core
        % the carrier must complete at least minCycles fringes
        minCycles = 5;
        mask = (PX > 0 | (PX == 0 & PY > 0)) & (PX.^2 + PY.^2 >= minCycles^2);

        Fm = F;
        Fm(~mask) = 0;
        [peakVal, ind] = max(Fm(:));
        if peakVal == 0
            error("analyzeFringe: no spectral peak found.")
        end
        [iPk, jPk] = ind2sub(size(Fm), ind);

        % quality check
        if peakVal < 10*median(F(mask))
            warning("analyzeFringe: dominant spectral peak is weak " + ...
                "(SNR ~%.3g). The detected period may be unreliable; " + ...
                "verify against the image.", peakVal/median(F(mask)))
        end

        % sub-bin refinement: quadratic interpolation of the log-magnitude
        % through the peak and its two neighbors, per frequency direction
        dx = FtpSolver.quadPeakOffset( ...
                F(iPk, mod(jPk - 2, nx) + 1), peakVal, F(iPk, mod(jPk, nx) + 1));
        dy = FtpSolver.quadPeakOffset( ...
                F(mod(iPk - 2, ny) + 1, jPk), peakVal, F(mod(iPk, ny) + 1, jPk));

        fx = (px(jPk) + dx)/nx;      % cycles per pixel
        fy = (py(iPk) + dy)/ny;

        period = 1/hypot(fx, fy);

        if abs(fx) >= abs(fy)
            normalAxis = 'X';       % fx > 0 by construction of the half-plane
            tiltAngle = atand(fy/fx);
        else
            normalAxis = 'Y';
            if fy < 0                % report the conjugate with fy > 0
                fx = -fx;
                fy = -fy;
            end
            tiltAngle = atand(fx/fy);
        end
        normalVec = [fx, fy]/norm([fx, fy]);
    end

    function W4D = morletCWT(img, scales, anglesDeg, sigma, gamma, fringeNormalAxis)
    % Continuous Morlet wavelet transform of a 1-D or 2-D input.
    %   Returns the complex CWT coefficients of IMG as an array of size
    %   [size(IMG,1), size(IMG,2), numel(SCALES), nAngles], L1-normalized.
    %
    %   2-D input: Wavelet Toolbox CWTFT2 with the anisotropic Morlet
    %   wavelet {Omega0 = 6, sigma, gamma}
    %
    %   1-D input (row or column vector): the analytic Morlet wavelet is
    %   built directly in the Fourier domain (k0 = 6, width set by sigma);
    %   anglesdeg and gamma are ignored and nAngles = 1.
    %   If fringeNormalAxis ('X' or 'Y') is provided, pads the array in
    %   that direction only and in both directions otherwise.

        if ~exist('fringeNormalAxis', 'var')
            fringeNormalAxis = '';
        end
        padSize = ceil(3 * sigma * max(scales));

        if min(size(img)) == 1
            % 1-D signal ----------------
            isCol = iscolumn(img);
            sig = padarray(img(:).', [0 padSize], 0);      % work on a padded row
            N = numel(sig);
            fSig = fft(sig);

            % Wavenumber grid
            k = (2*pi/N) * [0:floor(N/2), -ceil(N/2)+1:-1];

            k0 = 6;   % Morlet central wavenumber
            W4D = complex(zeros(1, N, numel(scales)));
            for iS = 1:numel(scales)
                a = scales(iS);
                % Analytic Morlet in the frequency domain, L1-normalized
                psiHat = exp(-sigma^2 * (a*k - k0).^2 / 2);
                W4D(1, :, iS) = ifft(fSig .* psiHat);
            end
            W4D = W4D(1, padSize+1:end-padSize, :);

            if isCol
                W4D = permute(W4D, [2 1 3 4]);
            end
        else
            % 2-D image ----------------
            if strcmpi(fringeNormalAxis, 'X')
                paddedImg = padarray(img, [0 padSize], 0);
                paddedImg = padarray(paddedImg, [padSize 0], 'symmetric');
            elseif strcmpi(fringeNormalAxis, 'Y')
                paddedImg = padarray(img, [padSize 0], 0);
                paddedImg = padarray(paddedImg, [0 padSize], 'symmetric');
            else
                paddedImg = padarray(img, [padSize padSize], 0);
            end
            cwtOut = cwtft2(paddedImg, ...
                wavelet = {"morlet", {6, sigma, gamma}}, ...
                scales  = scales, ...
                angles  = deg2rad(anglesDeg), ...
                norm    = "L1");

            % cfs is [rows x cols x 1 x nScales x nAngles]
            % drop the singleton plane dimension and crop the padding.
            W4D = permute(cwtOut.cfs, [1 2 4 5 3]);
            W4D = W4D(padSize+1:end-padSize, padSize+1:end-padSize, :, :);
        end
    end

    function elev = evalPolyElev(phase, coeffs)
    % Evaluate the polynomial phase-to-height model.
    % h = sum_i coeffs(:,:,i) .* phase.^(i-1)
    % coeffs planes are broadcast over the time dimension of phase.
        elev = zeros(size(phase));
        for i = 1:size(coeffs, 3)
            elev = elev + phase.^(i - 1) .* coeffs(:,:,i);
        end
    end

    function elev = evalTakedaElev(phaseDiff, L, d, f0)
    % Evaluate the Takeda phase-to-height model.
        elev = L*phaseDiff ./ (phaseDiff - 2*pi*f0*d);
    end

    function makeVideo(filename, framerate)
    % Read all images in directory and convert into animation.
        addr = uigetdir();
        if addr == 0
            return
        end

        filenames = dir(fullfile(addr, "*.png"));

        videoFilename = fullfile(addr, (filename + ".mp4"));
        videoSuffix = 2;
        while isfile(videoFilename)
            videoFilename = sprintf('%s_%g.mp4', fullfile(addr, filename), videoSuffix);
            videoSuffix = videoSuffix + 1;
        end

        v = VideoWriter(videoFilename,  'MPEG-4');
        v.FrameRate = framerate;
        open(v)

        for i = 1:length(filenames)
            frame = imread(fullfile(filenames(i).folder, filenames(i).name));
            if size(frame, 2) > 1920
                frame = imresize(frame, 1920/size(frame,2));
            end

            writeVideo(v, frame)
        end
        close(v)
    end
end

end