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
%                             cropRect=[100 150 1500 1200], ...
%                             camCalibPath="path/to/calibration.mat" ...
%                             );
%     profCase.setROI();
%     profCase.solve();
%     profCase.animate(1, 100, framesPerSecond);
%     profCase.writeCase(outputDir);
%
%   KEY METHODS
%     solve             - Run the full processing pipeline
%     setROI            - Interactively define the computational domain
%     calibratePoly     - Calibrate phase-to-height model
%     probe             - Extract time series at a spatial location
%     animate           - Visualize the reconstructed surface
%     writeCase         - Save results
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
%     Statistics and Machine Learning Toolbox
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
    imgData             % image buffers
    phaseData           % computed phase
    elevData            % computed surface elevation
    fringe              % fringe pattern parameters
    demodOpts           % demodulation settings        
    pCorrOpts           % phase correction settings and state
    prcOpts             % preprocessing and pipeline settings and flags
    loopState           % timestep bookkeeping
    cameraCalib         % camera calibration data
    worldCoords         % spatial coordinates
    elevModel           % phase to elevation conversion parameters
    outputConfig        % output/save control
    postData            % outlier coordinates
end

%% ----------- Setup and solve() -----------
methods (Access = public)
    function obj = FtpSolver(caseID, opts)
        arguments
            caseID              {mustBeTextScalar}              = "default"
            opts.refPath        {mustBeTextScalar}              = ""
            opts.dataPath       {mustBeTextScalar}              = ""
            opts.mmPerPixel     (1,1) double                    = NaN
            opts.period         (1,1) double                    = NaN
            opts.normAxis       string {mustBeMember(opts.normAxis, ["X", "Y"])} = "X"
            opts.resizeFactor   (1,1) double {mustBePositive}   = 1
            opts.cropRect       double                          = []
            opts.camCalibPath   {mustBeTextScalar}              = ""
            opts.camCalibNum    (1,1) double                    = 1    
            opts.dataCamInd     (1,1) {mustBeInteger}           = 1
            opts.refCamInd      (1,1) {mustBeInteger}           = 1
            opts.imgRotAngle    (1,1) double                    = 0
            opts.burntPixelMask                                = []
        end

        if ~isvarname(caseID)
            error("caseID must be a valid variable name.")
        end

        obj.caseID      = char(caseID);
        obj.source      = FtpSolver.defaultSource();
        obj.imgData     = FtpSolver.defaultImgData();
        obj.phaseData   = FtpSolver.defaultPhaseData();
        obj.elevData    = FtpSolver.defaultElevData();
        obj.fringe      = FtpSolver.defaultFringe();
        obj.pCorrOpts   = FtpSolver.defaultPCorrOpts();
        obj.demodOpts   = FtpSolver.defaultDemodOpts();
        obj.prcOpts     = FtpSolver.defaultPrcOpts();
        obj.loopState   = FtpSolver.defaultLoopState();
        obj.cameraCalib = FtpSolver.defaultCameraCalib();
        obj.worldCoords = FtpSolver.defaultWorldCoords();
        obj.elevModel   = FtpSolver.defaultElevModel();
        obj.outputConfig = FtpSolver.defaultOutputConfig();            


        obj.prcOpts.resizeFactor        = opts.resizeFactor;
        obj.prcOpts.cropRectOrg         = opts.cropRect;
        obj.prcOpts.imgRotAngle         = opts.imgRotAngle;
        obj.prcOpts.burntPixelMask     = opts.burntPixelMask;
        obj.source.dataCamInd           = opts.dataCamInd;
        obj.source.refCamInd            = opts.refCamInd;

        if opts.refPath == ""
            return   % bare construction used by loadobj()
        end
        

        % load reference image
        obj.loadRefImage(opts.refPath)
        
        % load data
        obj.loadDataset(opts.dataPath)

        if ~isnan(opts.mmPerPixel)
            obj.worldCoords.scaling.X.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.Y.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.mmPerPixelOrg   = opts.mmPerPixel;
        end

        if strlength(opts.camCalibPath) > 0
            obj.readCameraCalibration(opts.camCalibPath, opts.camCalibNum);
        end   

        obj.fringe.periodOrg = opts.period;
        obj.fringe.period = opts.period;

        obj.fringe.normAxis = opts.normAxis;
        if strcmpi(opts.normAxis, 'X')
            obj.fringe.normVec = [1, eps];
        elseif strcmpi(opts.normAxis, 'Y')
            obj.fringe.normVec = [eps, 1];
        end

        if ~isempty(obj.prcOpts.burntPixelMask)
            [obj.prcOpts.burntPixelRows, ...
             obj.prcOpts.burntPixelCols] = find(obj.prcOpts.burntPixelMask);
        end
                    
        obj.scaleAndTransformRef();

        % auto-detect the fringe period and pattern axis when none given
        if isnan(opts.period)
            img = imcrop(obj.imgData.refRaw, obj.prcOpts.cropRectOrg);
            [period, normAxis, tiltAngle, normVec] = obj.analyzeFringe(img);
            obj.fringe.periodOrg = period;
            obj.fringe.normAxis = normAxis;
            obj.fringe.normVec = normVec;
            obj.fringe.tiltAngle = tiltAngle;
            if strcmpi(normAxis, 'Y')
                obj.pCorrOpts.startEdge = 'top';
            end
            obj.updateGeometry();
        end
    end

    function solve(obj)
        % Run the full processing pipeline over solveRange.
        % Prepares the reference image, then for every frame: load,
        % preprocess, demodulate, unwrap, phase-correct, convert to
        % elevation and store. Frames flagged as invalid are stored as NaN,
        % and with chunked output on, each full chunk is written to a
        % temporary file. These files are merged when the user calls writeCase().

        t0 = tic();
        
        fprintf("Case %s\n", obj.caseID)

        obj.initParams();       
        obj.preprocessRefImage();
        obj.demodulateRef();

        for timestep = obj.prcOpts.solveRange(1):obj.prcOpts.solveRange(2)
            obj.loopState.timestep = timestep;
            obj.loopState.stackInd = obj.loopState.stackInd + 1;

            if obj.outputConfig.chunkedOutput && obj.loopState.stackInd > obj.outputConfig.chunkLength
                obj.writeChunk()
                obj.loopState.stackInd = obj.loopState.stackInd - obj.outputConfig.chunkLength;
                obj.loopState.chunkIndex = obj.loopState.chunkIndex + 1;
            end
            obj.reportProgress();
            
            obj.imgData.curr = obj.loadImage(timestep);
            obj.loopState.invalidFrame = obj.preprocessCurrentImage();

            if obj.loopState.invalidFrame
                obj.storeData()     % write NaN frame
                continue
            end
            
            obj.demodulate()
            obj.unwrapAndFlag();
            obj.correctPhase();
            obj.calculateElevation();
            obj.storeData();

            obj.loopState.isFirstTimestep = false;
        end
        fprintf('\nDone\n')

        if obj.outputConfig.chunkedOutput
            obj.writeChunk()
        end
        obj.loopState.processTime = toc(t0);
    end
    
    function calibrateTakeda(obj, heightVec, opts)
    % Fit L and d in Takeda's phase-to-height formula.
    %
    %   Fits the two geometric parameters of
    %       h = L*dPhi ./ (2*pi*f0*d + dPhi),    f0 = 1/(period*mmPerPixel)
    %   to the calibration stack by weighted nonlinear least squares over
    %   the calibration region, within user-specified bounds. L is the
    %   camera height above the reference plane and d the camera-projector
    %   separation, both in the units of heightVec.
    %
    %   Options:
    %       excludeHeights - logical mask over heightVec, true = drop plane
    %       weights        - per-plane fit weights
    %       boundsL        - [min max] bounds for L, default [0 Inf]
    %       boundsD        - [min max] bounds for d, default [0 Inf]
    %
    %   Stores obj.elevModel.L and obj.elevModel.d and sets
    %   obj.elevModel.type = 'takeda'.
        arguments
            obj
            heightVec double = []
            opts.excludeHeights = false(size(heightVec))
            opts.weights double = ones(size(heightVec))
            opts.boundsL (1,2) double = [0 Inf]
            opts.boundsD (1,2) double = [0 Inf]
        end
    
        [phaseCrop, heightVec, weights, geom] = obj.prepareCalibrationData( ...
            heightVec, opts.excludeHeights, opts.weights);
    
        % Carrier frequency on the reference plane (1/mm)
        f0 = 1 / (obj.fringe.period * obj.worldCoords.scaling.mmPerPixel);
    
        % Assemble fit data
        [rows, cols, nPlanes] = size(phaseCrop);
        phaseVec = reshape(double(phaseCrop), rows*cols, nPlanes);
    
        % Only two parameters are fitted, so a subsample of pixels is enough
        % Take every stride-th pixel to keep at most maxFitPoints samples
        % All planes are kept
        maxFitPoints = 2e5;
        stride  = max(1, ceil(rows*cols*nPlanes / maxFitPoints));
        phaseVec = phaseVec(1:stride:end, :);
        nPix = size(phaseVec, 1);
    
        xData = phaseVec(:);                    % plane-major stacking
        yData = repelem(heightVec, nPix);       % matching height per sample
        wData = repelem(weights, nPix);
    
        valid = isfinite(xData);                % drop NaN pixels (discarded
        xData = xData(valid);                   % frames, masked regions)
        yData = yData(valid);
        wData = wData(valid);
    
        % Start point from the linearized model
        %   1/h = 1/L + (2*pi*f0*d/L) * (1/dPhi)
        Xlin = [ones(numel(xData), 1), 1./xData];
        ab   = (wData.*Xlin) \ (wData.*(1./yData));
        L0   = 1/ab(1);
        d0   = ab(2)*L0 / (2*pi*f0);
    
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
        ft = fittype(@(L, d, dPhi) L*dPhi ./ (2*pi*f0*d + dPhi), 'independent', 'dPhi');
        fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Lower', [opts.boundsL(1), opts.boundsD(1)], ...
                        'Upper', [opts.boundsL(2), opts.boundsD(2)], ...
                        'StartPoint', [L0, d0], ...
                        'Weights', wData.^2);
        takedaFit = fit(xData, yData, ft, fo);
    
        L = takedaFit.L;
        d = takedaFit.d;
    
        % Diagnostics
        elevCrop = FtpSolver.evalTakedaElev(double(phaseCrop), L, d, f0);
    
        midY = floor((size(phaseCrop, 1) + 1) / 2);
        midX = floor((size(phaseCrop, 2) + 1) / 2);
        phaseSample = squeeze(double(phaseCrop(midY, midX, :)));
    
        fitPhase  = linspace(min(phaseSample) - 1, max(phaseSample) + 1, 1000);
        fitHeight = FtpSolver.evalTakedaElev(fitPhase, L, d, f0);
    
        obj.plotCalibrationDiagnostics(phaseCrop, elevCrop, fitPhase, fitHeight, ...
            heightVec, geom.calibRect, "Takeda model", {});
    
        fprintf("Takeda calibration: L = %.5g mm, d = %.5g mm, mean absolute error = %.3g mm\n", ...
            L, d, mean(abs(elevCrop - reshape(heightVec, 1, 1, [])), 'all', 'omitmissing'));
    
        % Store the model
        obj.elevModel.L = L;
        obj.elevModel.d = d;
        obj.elevModel.type = 'takeda';
    end

    function calibratePoly(obj, polyOrder, heightVec, opts)
    % Fit a per-pixel polynomial phase-to-height model.
    %
    %   For every pixel in the calibration region, fits
    %       h = a0 + a1*dPhi + a2*dPhi^2 + ... + a_polyOrder*dPhi^polyOrder
    %   (a0 fitted only when constOffset = true, zero otherwise) by weighted
    %   least squares against the calibration stack in phaseData.stack.
    %
    %   Stores in obj.elevModel.poly:
    %       coeffMatOrg  - pixelwise coefficients, RAW frame, NaN-padded
    %       fittedPlanes - poly22 fits of each coefficient over the region
     
        arguments
            obj
            polyOrder (1,1) double = 2
            heightVec double = []
            opts.excludeHeights = false(size(heightVec))
            opts.weights double = ones(size(heightVec))
            opts.constOffset (1,1) logical = false
        end
     
        [phaseCrop, heightVec, weights, geom] = obj.prepareCalibrationData( ...
            heightVec, opts.excludeHeights, opts.weights);
        % ================================================================
        % Pixelwise weighted least-squares fit
        % ================================================================
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
        % ================================================================
        % Smooth (poly22) coefficients
        % ================================================================
        fittedPlanes = cell(1, size(coeffMat, 3));
        smoothCoeffCrop = zeros(size(coeffMat));
        for i = 1:size(coeffMat, 3)
            [xData, yData, zData] = prepareSurfaceData( ...
                geom.pMeshX_crop, geom.pMeshY_crop, coeffMat(:,:,i));
            fittedPlanes{i} = fit([xData, yData], zData, fittype('poly22'));
            smoothCoeffCrop(:,:,i) = fittedPlanes{i}(geom.pMeshX_crop, ...
                                                        geom.pMeshY_crop);
        end
        % ================================================================
        % Diagnostics
        % ================================================================
        elevCropSmooth = FtpSolver.evalPolyElev(phaseCrop, smoothCoeffCrop);
        midY = floor((size(phaseCrop, 1) + 1) / 2);
        midX = floor((size(phaseCrop, 2) + 1) / 2);
        phaseSample = squeeze(phaseCrop(midY, midX, :));

        fitPhase = linspace(min(phaseSample) - 1, max(phaseSample) + 1, 1000);
        fitHeightSmooth = FtpSolver.evalPolyElev(fitPhase, ...
                                        smoothCoeffCrop(midY, midX, :));

        cLims = obj.plotCalibrationDiagnostics(phaseCrop, elevCropSmooth, ...
                        fitPhase, fitHeightSmooth, heightVec, ...
                        geom.calibRect, "Smoothed coefficients", {});
        
        elevCropPixelwise = FtpSolver.evalPolyElev(phaseCrop, coeffMat);
        fitHeightPixelwise = FtpSolver.evalPolyElev(fitPhase, ...
            coeffMat(midY, midX, :));

        obj.plotCalibrationDiagnostics(phaseCrop, elevCropPixelwise, ...
                        fitPhase, fitHeightPixelwise, heightVec, ...
                        geom.calibRect, "Pixelwise coefficients", cLims);
        % ================================================================
        % Store the model
        % ================================================================
        obj.elevModel.type = 'poly';
        obj.elevModel.poly.coeffMatOrg  = coeffMatOrg;
        obj.elevModel.poly.fittedPlanes = fittedPlanes;
        obj.elevModel.pixelwise = true;
    end

    function writePolyCalibration(obj, addr)
    % Save the polynomial calibration to a .mat file.
    %   The resized working copy (coeffMat) is dropped; only the RAW-frame
    %   coefficients and the fitted planes are written.
    %   In:  addr - output .mat file path

        temp = obj.elevModel.poly;
        if isfield(temp, 'coeffMat')
            temp = rmfield(temp, 'coeffMat');
        end
        save(addr, '-struct', 'temp')
    end

    function setStoredOutputs(obj, opts)
    % Choose which quantities solve() keeps in memory and writeCase() writes.
    %   setStoredOutputs(elev=true, phase=false, images=false)
    %   Quantities not named keep their current setting.
        arguments
            obj
            opts.elev   (1,1) logical = obj.outputConfig.elevOutput
            opts.phase  (1,1) logical = obj.outputConfig.phaseOutput
            opts.images (1,1) logical = obj.outputConfig.imageOutput
        end
        obj.outputConfig.elevOutput   = opts.elev;
        obj.outputConfig.phaseOutput  = opts.phase;
        obj.outputConfig.imageOutput = opts.images;
    end

    function setElevModelTakeda(obj, L, d)
        % Select Takeda's phase-to-height model.
        %   In:  L - camera height above the reference plane (mm)
        %        d - camera-projector separation (mm)

        obj.elevModel.type = 'takeda';
        obj.elevModel.L = L;
        obj.elevModel.d = d;
    end

    function setElevModelPoly(obj, addr)
        % Select a polynomial phase-to-height model and load calibration data from a .mat file.
        %   Input:
        %   addr - .mat file written by writePolyCalibration()

        obj.elevModel.type = 'poly';
        obj.elevModel.poly = load(addr);
        obj.elevModel.polyAddr = addr;
        if ~isfield(obj.elevModel, 'pixelwise')
            obj.elevModel.pixelwise = true;
        end
    end

    function setPolyPixelwise(obj, state)
    % Choose pixelwise or smoothed polynomial coefficients.
    %   In:  state - true for the raw pixelwise fit, false for the poly22
    %                smoothed coefficient planes

        if strcmpi(obj.elevModel.type, 'poly') && isfield(obj.elevModel, 'poly')
            obj.elevModel.pixelwise = logical(state);
        else
            warning("Polynomial model not found. Run setElevModelPoly(addr) first.")
        end
    end

    function setUnwrapMargin(obj, val)
    % Set the margin in RAW-frame pixels that is ignored by the unwrapping algorithm.
        obj.prcOpts.unwrapMarginOrg = val;
        obj.updateGeometry();
    end

    function setUnwrapMethod(obj, method)
    % Set the unwrapping algorithm.
        arguments 
            obj
            method (1,1) string {mustBeMember(method, ["1D", "2D"])}
        end
        obj.prcOpts.unwrapMethod = method;
    end

    function setSolveRange(obj, r1, r2)
    % Set the frame range processed by solve().
    %   In:  r1 - first frame, or 'full' to use the whole dataset
    %        r2 - last frame (ignored when r1 is 'full')

        if strcmpi(r1, 'full')
            obj.prcOpts.solveRange = obj.source.fullRange;
            return
        end
        if r2 > obj.source.fullRange(2) || r1 < 1
            error("Invalid input. Available data range is [1 %g]", obj.source.fullRange(2));
        end
        obj.prcOpts.solveRange = [r1 r2];

    end
    
    function setROI(obj)
    % Define the computational domain interactively.
    % Shows the raw reference image and waits for a rectangle to be
    % drawn. The result becomes cropRectOrg, the display rectangle is
    % set one unwrap margin inside it, and the geometry is refreshed.

        if isempty(obj.imgData.refScaled)
            error("Missing reference image.")
        end
        mg = obj.prcOpts.unwrapMarginOrg;

        figure()
        tempHndl = imshow(obj.imgData.refRaw, []);

        while true
            rect = drawrectangle(tempHndl.Parent);
            rect.Position = round(rect.Position);
            if isvalid(rect)
                if ~isempty(rect.Position)
                    obj.prcOpts.cropRectOrg = round(rect.Position);
                    obj.prcOpts.cropRectDisplayOrg = [mg + 1, mg + 1, ... 
                                            obj.prcOpts.cropRectOrg(3) - 2*mg, ...
                                            obj.prcOpts.cropRectOrg(4) - 2*mg];
                    break
                end
            else
                break
            end
        end

        if isvalid(tempHndl)
            close(tempHndl.Parent.Parent)
        end

        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry();
    end
    
    function setCropRect(obj, cropRect)
    % Specify the computational domain at full image resolution in [x y width height] format.
        obj.prcOpts.cropRectOrg = cropRect;
        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry(); 
    end

    function setDemodMethod(obj, method)
    % Set the demodulation algorithm.
        arguments 
            obj
            method (1,1) string {mustBeMember(method, ["fourier", "wavelet"])}
        end
        obj.demodOpts.method = method;
    end

    function setFourierFiltWidthFrac(obj, frac)
    % Set the width of the bandpass filter as a fraction of the carrier frequency (between 0 and 1).
        arguments 
            obj
            frac (1,1) double
        end
        obj.demodOpts.filtWidthFrac = frac;
    end
    
    function setPhaseCorrMethod(obj, method)
    % Choose how the 2*n*pi offset of each frame is fixed.
    %   setPhaseCorrMethod("temporal") - (default) Match each frame's 
    %       phase to the previous frame. Fails if the acquisition 
    %       frequency is too low. Failure at a single timestep will
    %       propagate to the end. 
    %   setPhaseCorrMethod("spatial")  - Track a fringe peak on the
    %       correction lines and align it with the reference. Requires
    %       the first fringe peak to remain visible throughout the
    %       recording.
    %   setPhaseCorrMethod("none")     - Leave the unwrapped phase as is.

        arguments 
            obj
            method (1,1) string {mustBeMember(method, ["spatial", "temporal", "none"])}
        end
        obj.pCorrOpts.method = method;
    end

    function setTrackedPeakInd(obj, peakInd)
    % Set the fringe peak that is tracked by the spatial phase correction algorithm.
        arguments
            obj
            peakInd {mustBeNumeric, mustBeScalarOrEmpty, mustBeInteger, mustBePositive}
        end
        obj.pCorrOpts.peakInd = double(peakInd);
    end

    function setTrackedPeakSafetyFactor(obj, factor)
    % Set the margin that is used to automatically determine the tracked fringe peak.
    % If no fringe peak is set by the user, the spatial phase correction
    % algorithm finds the closest peak within factor*fringePeriod pixels to
    % the edge and tracks that. 
        arguments
            obj
            factor (1,1) double {mustBePositive}
        end
        obj.pCorrOpts.edgeSafetyFactor = factor;
    end

    function setDisplayFromMargin(obj, offset)
    % Set the display crop rectangle as an offset from the unwrap margin.
        mrg = obj.prcOpts.unwrapMarginOrg;
        obj.prcOpts.cropRectDisplayOrg = [mrg + 1 + offset, mrg + 1 + offset, ... 
                                obj.prcOpts.cropRectOrg(3) - 2*mrg - 2*offset, ...
                                obj.prcOpts.cropRectOrg(4) - 2*mrg - 2*offset];
        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry();
    end

    function setResizeFactor(obj, factor)
    % Set the image resize factor (between 0 and 1) to reduce the computational cost of processing.
        arguments
            obj
            factor (1,1) double {mustBePositive}
        end
        obj.prcOpts.resizeFactor = factor;
        obj.updateGeometry();
    end

    function setResizeFactorDisplay(obj, factor)
    % Set the data resize factor (between 0 and 1) that is applied after the phase is calculated. Reduces storage size. 
        arguments
            obj
            factor (1,1) double {mustBePositive}
        end
        obj.prcOpts.resizeFactorDisplay = factor;
        obj.updateGeometry();
    end

    function setWaveletScales(obj, minWavelength, maxWavelength, divisions)
    % set wavelet scales based on pattern wavelength in pixels
    % wavelength should be for the raw, unscaled image
        scaleMin = minWavelength*6/(2*pi) * obj.prcOpts.resizeFactor;
        scaleMax = maxWavelength*6/(2*pi) * obj.prcOpts.resizeFactor;
        obj.demodOpts.scaleList = linspace(scaleMin, scaleMax, divisions);
    end

    function setWaveletAngles(obj, angleList)
    % Restrict angles to the range (-135, 45] for consistency with FT method.
        obj.demodOpts.angleList = mod(angleList - 45, -180) + 45;
    end
    
    function setBurntPixelMask(obj, mask)
        % Store a logical mask which indicates burnt pixels. 
        % Mask has to be the same size as the raw, unrectified image.
        obj.prcOpts.burntPixelMask = mask;
        [obj.prcOpts.burntPixelRows, ...
            obj.prcOpts.burntPixelCols] = find(mask);
    end
    
    function calibModeOn(obj)
    % Set properties in preparation for calibration.  
        obj.setResizeFactor(1);
        obj.setResizeFactorDisplay(1);
        obj.outputConfig.phaseOutput = true;
        obj.outputConfig.elevOutput = false;
    end

    function calibModeOff(obj)
    % Switch off storage of phase and switch on storage of elevation. 
        obj.outputConfig.phaseOutput = false;
        obj.outputConfig.elevOutput = true;
    end

    function readCameraCalibration(obj, calibPath, camNumber)
    % Read camera calibration from .mat or .xml file. 
        [~, ~, ext] = fileparts(calibPath);

        if strcmpi(ext, '.mat')
            obj.cameraCalib = load(calibPath);
            obj.cameraCalib.type = "Pinhole";

            imSize = obj.cameraCalib.intrinsics.ImageSize;
            [~, pTransform, scaling, xMesh, yMesh] = rectifyImagePinhole(zeros(imSize), obj.cameraCalib, 0, obj.prcOpts.imgRotAngle);
            obj.cameraCalib.pTransform = pTransform;
            obj.cameraCalib.scaling = scaling;

            obj.worldCoords.mesh.xOrg = xMesh;
            obj.worldCoords.mesh.yOrg = yMesh;

            obj.storePinholeScaling(scaling);

        elseif strcmpi(ext, '.xml')
            calibCell = readDavisCalibration(calibPath);
            calibData = calibCell{camNumber,1};
            calibType = calibCell{camNumber,2};

            if strcmpi(calibType, "Pinhole")
                obj.cameraCalib = calibData;

                imSize = obj.cameraCalib.intrinsics.ImageSize;
                [~, pTransform, scaling, xMesh, yMesh] = rectifyImagePinhole(zeros(imSize), obj.cameraCalib, 0, obj.prcOpts.imgRotAngle);
                obj.cameraCalib.pTransform = pTransform;
                obj.cameraCalib.scaling = scaling;
                
                obj.worldCoords.mesh.xOrg = xMesh;
                obj.worldCoords.mesh.yOrg = yMesh;

                obj.storePinholeScaling(scaling);
    
            elseif strcmpi(calibType, "Polynomial")
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

    function enableChunkedOutput(obj, chunkFolder, chunkLength)
    % Hold results in memory in chunks and flush each to disk.
    %   Store only chunks of chunkLength timesteps in memory. Temporary 
    %   chunk files are saved to path chunkFolder/TEMP_chunks_<caseID>
    %   and are merged into a single file when the user calls writeCase().
    %   Temporary files are deleted after they are merged.
        arguments
            obj
            chunkFolder string {mustBeFolder}
            chunkLength (1,1) double {mustBeInteger, mustBePositive}
        end
        obj.outputConfig.chunkLength = chunkLength;
        obj.outputConfig.chunkFolder = chunkFolder;
        obj.outputConfig.chunkedOutput = true;
    end

    function disableChunkedOutput(obj)
    % Store entire output in memory.
        obj.outputConfig.chunkedOutput = false;
    end

    function writeCase(obj, addr)
    % Save FtpSolver object, data binaries and coordinate mesh to a folder named <caseID>.
        if ~isfolder(fullfile(addr, obj.caseID))
            mkdir(fullfile(addr, obj.caseID))
        end

        temp.(obj.caseID) = obj;
        save(fullfile(addr, obj.caseID, obj.caseID + ".mat"), '-struct', 'temp')
        obj.writeResultFiles(fullfile(addr, obj.caseID))
    end

    function [optInd, isSafeLine, peakLocsCell] = findOptimalPeakInd(obj, peakInd)
    % Find which phase-correction peaks can be tracked inside the output frame.
    %
    %   [optInd, isSafeLine, peaksLineCell] = findOptimalPeakInd(obj, peakInd)
    %
    %   Detects the reference peak train on every phase-correction line and
    %   determines which peaks lie inside the output (display) frame while
    %   keeping a safety distance of pCorrOpts.edgeSafetyFactor pattern
    %   periods from the starting edge, so that motion of the surface
    %   cannot move the peak out of view. The first and last detected
    %   peaks are never eligible because the local-wavelength stencil in
    %   preprocessRefImage needs a neighbor on each side.
    %
    %   optInd        - recommended pCorrOpts.peakInd: the eligible peak
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

        nLines = length(obj.pCorrOpts.lineInds);
        optInd = nan;
        isSafeLine = false(1, nLines);
        peakLocsCell = cell(1, nLines);

        % prerequisites: reference image, correction lines, fringe period,
        % pattern normal axis, and output-frame geometry must all exist
        if isempty(obj.imgData.refRaw) || nLines == 0 ...
                || ~isfield(obj.fringe, 'periodOrg') ...
                || isempty(obj.fringe.periodOrg) ...
                || isnan(obj.fringe.periodOrg) ...
                || strlength(string(obj.fringe.normAxis)) == 0 ...
                || isempty(obj.prcOpts.cropRectOrg) ...
                || isempty(obj.prcOpts.cropRectDisplayOrg)
            return
        end

        % safe window along the pattern normal axis, RAW frame
        cropRectOrg = obj.prcOpts.cropRectOrg;
        dispRectOrg = obj.prcOpts.cropRectDisplayOrg;
        if strcmpi(obj.fringe.normAxis, 'X')
            winStart = cropRectOrg(1) + dispRectOrg(1) - 1;
            winEnd   = winStart + dispRectOrg(3);
        else
            winStart = cropRectOrg(2) + dispRectOrg(2) - 1;
            winEnd   = winStart + dispRectOrg(4);
        end

        safetyDist = obj.pCorrOpts.edgeSafetyFactor * obj.fringe.periodOrg;

        % the starting edge is winEnd for 'right'/'bottom' (where the peak
        % trains are flipped) and winStart otherwise
        if strcmpi(obj.pCorrOpts.startEdge, 'right') || strcmpi(obj.pCorrOpts.startEdge, 'bottom')
            safeLo = winStart;
            safeHi = winEnd - safetyDist;
        else
            safeLo = winStart + safetyDist;
            safeHi = winEnd;
        end

        % detect and evaluate every line
        optIndsPerLine = nan(1, nLines);
        for i = 1:nLines
            [~, peakLocs] = detectFringePeaks(obj, obj.imgData.refRaw, obj.pCorrOpts.lineInds(i));

            if strcmpi(obj.pCorrOpts.startEdge, 'right') || strcmpi(obj.pCorrOpts.startEdge, 'bottom')
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
    
    function updateGeometry(obj)
    % Recompute all resizeFactor-derived quantities
        rf = obj.prcOpts.resizeFactor;
        rf_o = obj.prcOpts.resizeFactorDisplay;

        obj.fringe.period = obj.fringe.periodOrg * rf;

        if ~isempty(obj.worldCoords.scaling)
            obj.worldCoords.scaling.X.Slope = obj.worldCoords.scaling.X.SlopeOrg / rf;
            obj.worldCoords.scaling.Y.Slope = obj.worldCoords.scaling.Y.SlopeOrg / rf;
            obj.worldCoords.scaling.mmPerPixel = obj.worldCoords.scaling.mmPerPixelOrg / rf;
        end

        obj.prcOpts.unwrapMargin = round(obj.prcOpts.unwrapMarginOrg * rf * rf_o);

        if isempty(obj.prcOpts.cropRectOrg)
            [Ny, Nx] = size(obj.imgData.refRaw);
            obj.prcOpts.cropRectOrg = [1, 1, Nx - 1, Ny - 1];
        end

        if ~isempty(obj.imgData.refRaw)
            obj.imgData.refScaled = imresize(obj.imgData.refRaw, rf);
        end

        cropRectOrg = round(obj.prcOpts.cropRectOrg);
        if isfield(obj.imgData, 'refRaw') && ~isempty(obj.imgData.refRaw)
            [NyRaw, NxRaw] = size(obj.imgData.refRaw);
            cropRectOrg(1) = max(1, cropRectOrg(1));
            cropRectOrg(2) = max(1, cropRectOrg(2));
            if (cropRectOrg(1) + cropRectOrg(3)) > NxRaw
                cropRectOrg(3) = NxRaw - cropRectOrg(1);
            end
            if (cropRectOrg(2) + cropRectOrg(4)) > NyRaw
                cropRectOrg(4) = NyRaw - cropRectOrg(2);
            end

            obj.prcOpts.cropRectOrg = cropRectOrg;
        end

        [NyScaled, NxScaled] = size(obj.imgData.refScaled);
        cropRect = round(cropRectOrg * rf);
        cropRect(1:2) = max(1, cropRect(1:2));
        cropRect(3) = min(NxScaled - cropRect(1), cropRect(3));
        cropRect(4) = min(NyScaled - cropRect(2), cropRect(4));
        obj.prcOpts.cropRect = cropRect;

        % the dimensions of the resized and cropped reference image is
        % required to initialize arrays
        % obj.imgData.ref is reinitialized in preprocessRefImage()
        obj.imgData.ref = imcrop(obj.imgData.refScaled, obj.prcOpts.cropRect);
        stackSize = ceil(rf_o*size(obj.imgData.ref));
        
        if isempty(obj.prcOpts.cropRectDisplayOrg)
            obj.prcOpts.cropRectDisplayOrg = [1 1 cropRectOrg(3) cropRectOrg(4)];
        end
        cropRectDisplay(1:2) = ceil(obj.prcOpts.cropRectDisplayOrg(1:2) * rf * rf_o);
        cropRectDisplay(3:4) = floor(obj.prcOpts.cropRectDisplayOrg(3:4) * rf * rf_o);

        cropRectDisplay(1:2) = max(1, cropRectDisplay(1:2));
        cropRectDisplay(3) = min(stackSize(2) - cropRectDisplay(1), cropRectDisplay(3));
        cropRectDisplay(4) = min(stackSize(1) - cropRectDisplay(2), cropRectDisplay(4));
        obj.prcOpts.cropRectDisplay = cropRectDisplay;

        obj.initPhaseCorrection();

        % initialize coordinate mesh
        [ny, nx, ~] = size(obj.imgData.refScaled);
        [Ny, Nx] = size(obj.imgData.refRaw);

        if isfield(obj.worldCoords.mesh, 'xOrg')
            xComp = imresize(obj.worldCoords.mesh.xOrg, ...
                                    rf, 'bilinear', Antialiasing=false);
            yComp = imresize(obj.worldCoords.mesh.yOrg, ...
                                    rf, 'bilinear', Antialiasing=false);

        else
            x0 = obj.worldCoords.scaling.X.Offset;
            xSlope = obj.worldCoords.scaling.X.Slope;
            y0 = obj.worldCoords.scaling.Y.Offset;
            ySlope = obj.worldCoords.scaling.Y.Slope;

            [xComp, yComp] = meshgrid(x0:xSlope:x0 + xSlope*(nx - 1), ...
                                        y0:ySlope:y0 + ySlope*(ny - 1));
        end

        xComp = imcrop(xComp, obj.prcOpts.cropRect);
        yComp = imcrop(yComp, obj.prcOpts.cropRect);

        obj.worldCoords.mesh.x = imresize(xComp, rf_o, 'bilinear', ...
            Antialiasing=false);
        obj.worldCoords.mesh.y = imresize(yComp, rf_o, 'bilinear', ...
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

            xPixelComp = imresize(xPixelOrg, rf, 'bilinear', Antialiasing=false);
            yPixelComp = imresize(yPixelOrg, rf, 'bilinear', Antialiasing=false);

            xPixelComp = imcrop(xPixelComp, obj.prcOpts.cropRect);
            yPixelComp = imcrop(yPixelComp, obj.prcOpts.cropRect);

            obj.worldCoords.mesh.xPixel = imresize(xPixelComp, rf_o, ...
                'bilinear', Antialiasing=false);
            obj.worldCoords.mesh.yPixel = imresize(yPixelComp, rf_o, ...
                'bilinear', Antialiasing=false);
        end

        % refresh the auto-selected phase-correction peak 
        if isempty(obj.pCorrOpts.peakInd)
            optInd = obj.findOptimalPeakInd();
            if ~isnan(optInd)
                obj.pCorrOpts.peakInd = optInd;
            end
        end
    end
end

%% ----------- Inspection and Visualization -----------
methods (Access = public)
    function showROI(obj)
        % Draw the domain rectangles on the reference image.
        % Overlays the computational domain (red), the unwrap margin (blue)
        % and the display domain (green) for visual checking.

        figure;
        image = obj.imgData.refRaw;
        image = rescale(image);

        if ~isempty(obj.prcOpts.cropRectOrg)
            image = insertShape(image, 'rectangle', obj.prcOpts.cropRectOrg, ...
                'LineWidth', 5, 'ShapeColor', 'red');
            image = insertText(image, obj.prcOpts.cropRectOrg(1:2), ...
                'Computational domain', 'FontColor', 'red', 'FontSize', 20);
            % draw unwrap margin
            mg = obj.prcOpts.unwrapMarginOrg;
            unwrapRect = obj.prcOpts.cropRectOrg + [mg, mg, -2*mg, -2*mg];
            image = insertShape(image, 'rectangle', unwrapRect, ...
                'LineWidth', 5, 'ShapeColor', 'blue');
            image = insertText(image, unwrapRect(1:2) + [unwrapRect(3) 0], ...
                'Unwrap margin', 'FontColor', 'blue', ...
                'FontSize', 20, 'AnchorPoint', 'RightTop');
        end

        if ~isempty(obj.prcOpts.cropRectDisplayOrg)
            absoluteCropDisplay = obj.prcOpts.cropRectDisplayOrg;
            absoluteCropDisplay(1:2) = absoluteCropDisplay(1:2) + obj.prcOpts.cropRectOrg(1:2);
            image = insertShape(image, 'rectangle', absoluteCropDisplay, ...
                'LineWidth', 2, 'ShapeColor', 'green');
            image = insertText(image, [absoluteCropDisplay(1), ...
                absoluteCropDisplay(2) ...
                + absoluteCropDisplay(4)], ...
                'Display domain', 'FontColor', ...
                'green', 'FontSize', 20, ...
                'AnchorPoint', 'LeftBottom');
        end
        imshow(image, [])
    end

    function [subarray, xCrop, yCrop] = getElevDisplaySubarray(obj, firstTimestep, lastTimestep)
    % Crop the elevation stack to display frame and return.
        solveRange = obj.prcOpts.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = firstTimestep - solveRange(1) + 1:lastTimestep - solveRange(1) + 1; 

        if ~isempty(obj.elevData.stack)
            subarray = obj.elevData.stack(rowInds, colInds, tRange);
        else
            subarray = 0;
        end
        xCrop = obj.worldCoords.mesh.x(rowInds, colInds);
        yCrop = obj.worldCoords.mesh.y(rowInds, colInds);
    end

    function [subarray, xCrop, yCrop] = getPhaseDisplaySubarray(obj, firstTimestep, lastTimestep)
    % Crop the phase stack to display frame and return.
        solveRange = obj.prcOpts.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = firstTimestep - solveRange(1) + 1:lastTimestep - solveRange(1) + 1; 

        if ~isempty(obj.phaseData.stack)
            subarray = obj.phaseData.stack(rowInds, colInds, tRange);
        else
            subarray = 0;
        end 
        xCrop = obj.worldCoords.mesh.x(rowInds, colInds);
        yCrop = obj.worldCoords.mesh.y(rowInds, colInds);
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
    %           quantity  - "elev" (default) or "phase"
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
            opts.quantity  (1,1) string {mustBeMember(opts.quantity, ["elev", "phase"])} = "elev"
            opts.patchSize (1,1) double {mustBeNonnegative, mustBeFinite} = 0
            opts.shape     (1,1) string {mustBeMember(opts.shape, ["square", "circle"])} = "square"
            opts.operation (1,1) string {mustBeMember(opts.operation, ["mean", "min", "max"])} = "mean"
        end

        switch opts.quantity
            case "elev"
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

    function animate(obj, firstTimestep, lastTimestep, framerate, opts)
    % Animate the reconstructed surface or phase.
    %   Plots the surface elevation or phase (or an array passed in opts.data) frame
    %   by frame, with optional plane/mean/constant subtraction, NaN-frame
    %   interpolation, smoothing and a time or frame annotation.
    %    
    %   Example usage: animate(1, 200, 24, quantity="phase")
    %   In:  firstTimestep, lastTimestep - frame range ([] = to the end)
    %        framerate          - playback rate (fps)
    %        opts               - name-value display options (data, ZLim,
    %                             view, colormap, camFPS, subtractPlane,
    %                             subtractMean, smoothing, ...)

        arguments
            obj
            firstTimestep = [];
            lastTimestep = [];
            framerate (1,1) double = 30;
            opts.quantity (1,1) string ...
                    {mustBeMember(opts.quantity, ["elev", "phase"])}= "elev";
            opts.camFPS (1,1) double = 1;
            opts.data cell = {};
            opts.dispTime (1,1) logical = true;
            opts.subtractPlane (1,1) logical = false;
            opts.subtractMean (1,1) logical = false;
            opts.subtractConst = [];
            opts.smoothing (1,1) logical = false;
            opts.smoothingSigma (1,1) double = 2;
            opts.annotTStart (1,1) double = 0;
            opts.nanInterp (1,1) logical = false;
            opts.ZLim (1,2) double = [0 0];
            opts.lightPosition (1,3) double = [1 1 5];
            opts.FontSize (1,1) double = 16;
            opts.DataAspectRatio = [2 2 1];
            opts.view = [23.2 23.8]
            opts.nominalSkip (1,1) double = 1;
            opts.figPosition (1,4) double = [0.05 0.13 0.6 0.6];
            opts.tFormatString string = "";
            opts.colormap = [];
        end  

        if ~isempty(opts.data)
            dataArray = opts.data{1};

            if isempty(firstTimestep)
                firstTimestep = 1;
            end

            if isempty(lastTimestep)
                lastTimestep = size(dataArray, 3);
            end
            dataArray = dataArray(:,:,firstTimestep:lastTimestep);
            if length(opts.data) > 1
                xCrop = opts.data{2};
                yCrop = opts.data{3};
            else
                [xCrop, yCrop] = meshgrid(1:size(dataArray, 2), 1:size(dataArray, 1));
            end
        else
            if isempty(firstTimestep)
                firstTimestep = obj.prcOpts.solveRange(1);
            end

            if isempty(lastTimestep)
                if strcmpi(opts.quantity, 'elev')
                    stackSize = size(obj.elevData.stack, 3);
                else
                    stackSize = size(obj.phaseData.stack, 3);
                end
                lastTimestep = obj.prcOpts.solveRange(1) + stackSize - 1;
            end
            if strcmpi(opts.quantity, "elev")
                [dataArray, xCrop, yCrop] = getElevDisplaySubarray(obj, firstTimestep, lastTimestep);
            elseif strcmpi(opts.quantity, "phase")
                [dataArray, xCrop, yCrop] = getPhaseDisplaySubarray(obj, firstTimestep, lastTimestep);
            end

        end

        pauseTime = 1/framerate;
    
        % check for nan frames and interpolate
        if opts.nanInterp
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

        surfHndl = surfl(axHndl, xCrop, yCrop, dataArray(:,:,1), 'light', 'EdgeColor', 'none');
        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % surfHndl(1).AmbientStrength = 0.4;
        surfHndl(2).Style = 'infinite';
        surfHndl(2).Position = opts.lightPosition;

        xlim([min(xCrop, [], 'all'), max(xCrop, [], 'all')])
        ylim([min(yCrop, [], 'all'), max(yCrop, [], 'all')])
        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        cbHndl = colorbar;

        if strcmpi(opts.quantity, 'elev')
            cbHndl.Label.String = 'Elevation (mm)';
        elseif strcmpi(opts.quantity, 'phase')
            cbHndl.Label.String = '\Delta\phi (rad)';
        end
        clim(axHndl, [midHeight - 0.95*amp, midHeight + 0.95*amp])

        xlabel(axHndl, 'x (mm)')
        ylabel(axHndl, 'y (mm)')
        if strcmpi(opts.quantity, 'elev')
            zlabel(axHndl, 'z (mm)')
        elseif strcmpi(opts.quantity, 'phase')
            zlabel(axHndl, '\Delta\phi (rad)')
        end

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(figHndl, 'Position', opts.figPosition)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        cbHndl.Position = [0.917,0.11,0.02,0.815];

        if opts.dispTime
            annotHndl = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.tFormatString) ~= 0
                tFormatString = opts.tFormatString;
            else
                if opts.camFPS == 1
                    tFormatString = "Frame = %05g";
                else
                    tFormatString = "t = %.3f s";
                end
            end
            timeString = sprintf(tFormatString, opts.annotTStart);
            set(annotHndl, 'String', timeString)
        end

        set(gca, 'DataAspectRatio', opts.DataAspectRatio);

        pause(pauseTime)
    
        for i = 2:size(dataArray, 3)
            if ~ishghandle(figHndl)
                break
            end
            surfHndl(1).ZData = dataArray(:,:,i);

            if opts.dispTime
                if opts.camFPS == 1                    
                    timeString = sprintf(tFormatString, opts.annotTStart + (i - 1)*opts.nominalSkip);
                else
                    timeString = sprintf(tFormatString, opts.annotTStart + (i - 1)*opts.nominalSkip/opts.camFPS);
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
    % Example usage: exportVideo(1, 200, 24, quantity="phase")
    %   Inputs:
    %        firstTimestep, lastTimestep - frame range ([] = to the end)
    %        framerate          - video framerate (fps)
    %        opts               - name-value display options (data, ZLim,
    %                             view, colormap, camFPS, subtractPlane,
    %                             subtractMean, smoothing, ...)
        arguments
            obj
            firstTimestep = [];
            lastTimestep = [];
            framerate (1,1) double = 30;
            opts.quantity (1,1) string ...
                    {mustBeMember(opts.quantity, ["elev", "phase"])}= "elev";
            opts.camFPS (1,1) double = 1;
            opts.data cell = {};
            opts.saveAddr string = [];
            opts.fast (1,1) logical = false;
            opts.dispTime (1,1) logical = true;
            opts.subtractPlane (1,1) logical = false;
            opts.subtractMean (1,1) logical = false;
            opts.subtractConst = [];
            opts.smoothing (1,1) logical = false;
            opts.smoothingSigma (1,1) double = 2;
            opts.annotTStart (1,1) double = 0;
            opts.nanInterp (1,1) logical = false;
            opts.cleanup (1,1) logical = false;
            opts.ZLim (1,2) double = [0 0];
            opts.lightPosition (1,3) double = [1 1 5];
            opts.FontSize (1,1) double = 16;
            opts.DataAspectRatio = [2 2 1];
            opts.view = [23.2 23.8]
            opts.overwriteFile (1,1) logical = false;
            opts.nominalSkip (1,1) double = 1;
            opts.imageResolution (1,1) double = 150;
            opts.figPosition (1,4) double = [0.05 0.13 0.6 0.6];
            opts.tFormatString string = "";
            opts.colormap = [];
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

        if ~isempty(opts.data)
            dataArray = opts.data{1};
            if isempty(firstTimestep)
                firstTimestep = 1;
            end
            if isempty(lastTimestep)
                lastTimestep = size(dataArray, 3);
            end
            dataArray = dataArray(:,:,firstTimestep:lastTimestep);
            if length(opts.data) > 1
                xCrop = opts.data{2};
                yCrop = opts.data{3};
            else
                [xCrop, yCrop] = meshgrid(1:size(dataArray, 2), 1:size(dataArray, 1));
            end
        else
            if isempty(firstTimestep)
                firstTimestep = obj.prcOpts.solveRange(1);
            end
            if isempty(lastTimestep)
                if strcmpi(opts.quantity, 'elev')
                    stackSize = size(obj.elevData.stack, 3);
                else
                    stackSize = size(obj.phaseData.stack, 3);
                end
                lastTimestep = obj.prcOpts.solveRange(1) + stackSize - 1;
            end

            if strcmpi(opts.quantity, "elev")
                [dataArray, xCrop, yCrop] = getElevDisplaySubarray(obj, firstTimestep, lastTimestep);
            elseif strcmpi(opts.quantity, "phase")
                [dataArray, xCrop, yCrop] = getPhaseDisplaySubarray(obj, firstTimestep, lastTimestep);
            end
        end
               
        % check for nan frames and interpolate
        if opts.nanInterp
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

        surfHndl = surfl(axHndl, xCrop, yCrop, dataArray(:,:,1), 'light', ...
                    'EdgeColor', 'none');

        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % surfHndl(1).AmbientStrength = 0.4;
        surfHndl(2).Style = 'infinite';
        surfHndl(2).Position = opts.lightPosition;

        xlim([min(xCrop, [], 'all'), max(xCrop, [], 'all')])
        ylim([min(yCrop, [], 'all'), max(yCrop, [], 'all')])

        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        cbHndl = colorbar;
        if strcmpi(opts.quantity, 'elev')
            cbHndl.Label.String = 'Elevation (mm)';
        elseif strcmpi(opts.quantity, 'phase')
            cbHndl.Label.String = '\Delta\phi (rad)';
        end
        clim(axHndl, [midHeight - 0.95*amp, midHeight + 0.95*amp])

        xlabel(axHndl, 'x (mm)')
        ylabel(axHndl, 'y (mm)')
        if strcmpi(opts.quantity, 'elev')
            zlabel(axHndl, 'z (mm)')
        elseif strcmpi(opts.quantity, 'phase')
            zlabel(axHndl, '\Delta\phi (rad)')
        end

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        cbHndl.Position = [0.917,0.11,0.02,0.815];

        if opts.dispTime
            annotHndl = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.tFormatString) ~= 0
                tFormatString = opts.tFormatString;
            else
                if opts.camFPS == 1
                    tFormatString = "Frame = %05g";
                else
                    tFormatString = "t = %.3f s";
                end
            end
            timeString = sprintf(tFormatString, opts.annotTStart);
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
            
            if opts.dispTime
                if opts.camFPS == 1                    
                    timeString = sprintf(tFormatString, opts.annotTStart + (i - 1)*opts.nominalSkip);
                else
                    timeString = sprintf(tFormatString, opts.annotTStart + (i - 1)*opts.nominalSkip/opts.camFPS);
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
            firstTimestep (1,1) double = obj.prcOpts.solveRange(1)
            lastTimestep (1,1) double = obj.prcOpts.solveRange(2)
            framerate (1,1) double = 5
        end
        figHndl = figure(120);
        
        pauseTime = 1/framerate;

        obj.imgData.curr = obj.loadImage(firstTimestep);
        obj.preprocessCurrentImage();

        [values, edges] = histcounts(obj.imgData.curr, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));
        minThresh = edges(find(values > 0.001, 1, 'first'));
        
        imHndl = imshow(obj.imgData.curr, [minThresh maxThresh]);
        title(num2str(firstTimestep))                
        pause(pauseTime)

        for timestep = firstTimestep + 1:lastTimestep
            if ~ishghandle(figHndl)
                break
            end

            obj.imgData.curr = obj.loadImage(timestep);
            obj.preprocessCurrentImage();

            set(imHndl, 'CData', obj.imgData.curr)
            title(num2str(timestep))                
            pause(pauseTime)
        end
    end

    function animateRawImages(obj, firstTimestep, lastTimestep, framerate)
    % Rectify (if possible) the image sequence and then animate without additional preprocessing.
        arguments
            obj
            firstTimestep (1,1) double = obj.prcOpts.solveRange(1)
            lastTimestep (1,1) double = obj.prcOpts.solveRange(2)
            framerate (1,1) double = 5
        end
        figHndl = figure();
        
        pauseTime = 1/framerate;

        img = obj.loadImage(firstTimestep);

        [values, edges] = histcounts(img, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));

        for timestep = firstTimestep:lastTimestep
            if ~ishghandle(figHndl)
                break
            end

            img = obj.loadImage(timestep);
            img = obj.rectify(img);
   
            if timestep == firstTimestep
                imHndl = imshow(img, [0 maxThresh]);
            else
                set(imHndl, 'CData', img)
            end
            title(num2str(timestep))                
            pause(pauseTime)
        end
    end

    function imageCell = getRawImages(obj, firstTimestep, lastTimestep)
    % Rectify and return the image sequence in a cell array without preprocessing.
        arguments
            obj
            firstTimestep (1,1) double = obj.prcOpts.solveRange(1)
            lastTimestep (1,1) double = obj.prcOpts.solveRange(2)
        end

        imageCell = cell(1, lastTimestep - firstTimestep + 1);

        for timestep = firstTimestep:lastTimestep
            img = obj.loadImage(timestep);
            img = obj.rectify(img);
            imageCell{timestep - firstTimestep + 1} = img;
        end            
    end

    function imageCell = getImages(obj, firstTimestep, lastTimestep)
    % Preprocess and return the image sequence in a cell array.
        arguments
            obj
            firstTimestep (1,1) double = obj.prcOpts.solveRange(1)
            lastTimestep (1,1) double = obj.prcOpts.solveRange(2)
        end

        imageCell = cell(1, lastTimestep - firstTimestep + 1);

        for timestep = firstTimestep:lastTimestep
            obj.imgData.curr = obj.loadImage(timestep);
            obj.preprocessCurrentImage();
            imageCell{timestep - firstTimestep + 1} = obj.imgData.curr;
        end            
    end

    function drawPeaks(obj, firstTimestep, lastTimestep, skip)
    % Diagnostic plot of phase-correction peak detection.
    % Shows the reference line with its detected peaks and the tracked
    % peak highlighted, then the same line in the requested images.
        if ~exist('skip', 'var')
            skip = 1;
        end

        lineNum = min(2, length(obj.pCorrOpts.lineInds));
        lineInd = obj.pCorrOpts.lineInds(lineNum);

        [optInd, ~, refPeakLocsCell] = obj.findOptimalPeakInd(obj.pCorrOpts.peakInd);

        peakInd = obj.pCorrOpts.peakInd;
        if isempty(peakInd)
            if isnan(optInd)
                error("Cannot determine a peak to highlight: pCorrOpts.peakInd is " + ...
                    "empty and no eligible peak was found. Check that the reference " + ...
                    "image, fringe period, and crop rectangles are set.")
            end
            peakInd = optInd;
            fprintf("pCorrOpts.peakInd is empty; showing the automatic choice (%g).\n", optInd)
        end

        refPeakLocs = refPeakLocsCell{lineNum};
        if isempty(refPeakLocs) || peakInd > length(refPeakLocs)
            error("Peak %g was not detected on phase-correction line %g of the reference image.", ...
                peakInd, lineInd)
        end

        if strcmpi(obj.fringe.normAxis, 'X')
            refLine = double(obj.imgData.refRaw(lineInd, :));
        else
            refLine = double(obj.imgData.refRaw(:, lineInd));
        end
        refPeakVals = refLine(refPeakLocs);

        % detect the same line in the requested images
        targetImages = firstTimestep:skip:lastTimestep;
        N = length(targetImages);

        currPeaksCell = cell(1, N);
        currLineCell = cell(1, N);

        for n = 1:N
            timestep = targetImages(n);
            img = obj.loadImage(timestep);
            img = obj.rectify(img);

            [peakVals, peakLocs] = detectFringePeaks(obj, img, lineInd);

            if strcmpi(obj.pCorrOpts.startEdge, 'right') || strcmpi(obj.pCorrOpts.startEdge, 'bottom')
                peakLocs = flip(peakLocs);
                peakVals = flip(peakVals);
            end

            currPeaksCell{n} = [peakLocs(:), peakVals(:)];

            if strcmpi(obj.fringe.normAxis, 'X')
                currLine = double(img(lineInd, :));
            else
                currLine = double(img(:, lineInd));
            end

            currLineCell{n} = currLine;
        end

        % plot
        figure;
        tiledlayout(3, ceil(N/3) + 1)
        nexttile([3 1])
        plot(refLine)
        hold on
        scatter(refPeakLocs, 1.01*refPeakVals, 'v', 'filled')
        scatter(refPeakLocs(peakInd), 1.01*refPeakVals(peakInd), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
        hold off
        if strcmpi(obj.pCorrOpts.startEdge, 'left') || strcmpi(obj.pCorrOpts.startEdge, 'top')
            xLimit = [1 round(1.5*refPeakLocs(peakInd))];
        else
            xLimit = [round(0.9*refPeakLocs(peakInd)) length(refLine)];
        end
        yLimit = [0.8*min(refLine(xLimit(1):xLimit(2))) 1.2*max(refLine(xLimit(1):xLimit(2)))];

        xlim(xLimit)
        ylim(yLimit)
        title("Reference: x_p = " + num2str(refPeakLocs(peakInd)) )
        for i = 1:N
            nexttile
            plot(currLineCell{i})
            hold on
            scatter(currPeaksCell{i}(:,1), 1.05*currPeaksCell{i}(:,2), 'v', 'filled')
            if size(currPeaksCell{i}, 1) >= peakInd
                scatter(currPeaksCell{i}(peakInd, 1), 1.05*currPeaksCell{i}(peakInd,2), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
                title(num2str(targetImages(i)) + ": x_p = " + num2str(currPeaksCell{i}(peakInd, 1)) )
            else
                title(num2str(targetImages(i)) + ": peak " + num2str(peakInd) + " not detected")
            end
            hold off
            xlim(xLimit)
            ylim(yLimit)
        end
    end
end

%% ----------- Internal methods -----------
methods (Access = private)
    function initParams(obj)
    % Reset loop state and preallocate arrays used by solve().
    %   Refreshes the geometry (and the camera mesh of a reloaded case),
    %   preallocates the elevation/phase stacks for the whole run or
    %   for one chunk and resamples the polynomial calibration onto the
    %   computational grid.

        if obj.outputConfig.elevOutput && strcmpi(obj.elevModel.type, 'none')
            error("Elevation output is on but no elevation model is set." + ...
                   "Call setElevModelTakeda or setElevModelPoly, or " + ...
                   "setStoredOutputs(elev=false).")
        end

        if isempty(obj.fringe.periodOrg) || isnan(obj.fringe.periodOrg)
            error("Missing data: fringe.period")
        end

        if ~obj.outputConfig.imageOutput ...
                && ~obj.outputConfig.phaseOutput ...
                && ~obj.outputConfig.elevOutput ...
                && (diff(obj.prcOpts.solveRange) + 1 > 1)
            error("All storage switches are off. Check outputConfig.")
        end

        obj.loopState.stackInd = 0;
        obj.loopState.chunkIndex = 1;
        obj.loopState.isFirstTimestep = true;
        obj.elevData.curr = [];
        obj.phaseData.curr = [];
        obj.pCorrOpts.warningTimesteps = [];

        rf = obj.prcOpts.resizeFactor;
        rf_o = obj.prcOpts.resizeFactorDisplay;
        Nt = diff(obj.prcOpts.solveRange) + 1;
        chunkLength = obj.outputConfig.chunkLength;

        % refresh coordinate mesh if running a case loaded from a .mat file
        if isfield(obj.cameraCalib, 'addr') && ~isfield(obj.worldCoords.mesh, 'xOrg')
            if strlength(obj.cameraCalib.path) > 0 && strcmpi(obj.cameraCalib.type, 'Pinhole')
                obj.readCameraCalibration(obj.cameraCalib.path, obj.cameraCalib.camNum);
            end
        end
        % refresh images if running a case loaded from a .mat file
        if strlength(obj.source.refPath) > 0 && isempty(obj.imgData.refRaw)
            obj.initImages();
        end
        obj.updateGeometry();

        % preallocate stacks
        if obj.outputConfig.chunkedOutput && Nt > chunkLength
            obj.allocateStacks(chunkLength)
        else
            obj.allocateStacks(Nt)
        end
        
        % initialize polynomial calibration model
        [Ny, Nx] = size(obj.imgData.refRaw);
        if obj.outputConfig.elevOutput && strcmpi(obj.elevModel.type, 'poly')
            if obj.elevModel.pixelwise
                coeffMat = imresize(obj.elevModel.poly.coeffMatOrg, ...
                                    rf, 'bilinear', Antialiasing=false);
            else
                [xMesh, yMesh] = meshgrid(1:Nx, 1:Ny);
                fPlanes = obj.elevModel.poly.fittedPlanes;
                coeffMat = imresize(zeros(Ny, Nx, length(fPlanes)), ...
                                    rf, 'bilinear', Antialiasing=false);
                for i = 1:length(fPlanes)
                    coeffVals = fPlanes{i}(xMesh, yMesh);
                    coeffMat(:,:,i) = imresize(coeffVals, ...
                                    rf, 'bilinear', Antialiasing=false);
                end
            end

            cropRect3 = [obj.prcOpts.cropRect(1:2), 1, ...
                obj.prcOpts.cropRect(3:4), ...
                size(obj.elevModel.poly.coeffMatOrg, 3) - 1 ...
                ];
            coeffMat = imcrop3(coeffMat, cropRect3);
            coeffMat = imresize(coeffMat, rf_o, 'bilinear', Antialiasing=false);
            obj.elevModel.poly.coeffMat = coeffMat;
        end

        obj.postData.phaseAnomalies = struct('timestep', {}, 'rects', {});
    end

    function allocateStacks(obj, nTimesteps)
    % Preallocate the result stacks for nTimesteps.
    %   Allocates a NaN-filled single-precision stack for every quantity
    %   whose output switch is on and clears the others.
    
        [NyImg, NxImg] = size(obj.imgData.ref);
        rf_o = obj.prcOpts.resizeFactorDisplay;
        NyStack = ceil(rf_o * NyImg);
        NxStack = ceil(rf_o * NxImg);
    
        obj.phaseData.stack = [];
        obj.elevData.stack  = [];
        obj.imgData.stack   = [];
    
        if obj.outputConfig.phaseOutput
            obj.phaseData.stack = nan(NyStack, NxStack, nTimesteps, 'single');
        end
        if obj.outputConfig.elevOutput
            obj.elevData.stack = nan(NyStack, NxStack, nTimesteps, 'single');
        end
        if obj.outputConfig.imageOutput
            obj.imgData.stack = nan(NyImg, NxImg, nTimesteps, 'single');
        end
    end

    function preprocessRefImage(obj)
    % Prepare the reference image and the spatial phase-correction reference.
    %   Crops the scaled reference to the computational domain, removes
    %   the DC component, and for spatial phase correction selects or
    %   validates pCorrOpts.peakInd and measures the local fringe
    %   wavelength and reference peak position on every correction line.

        obj.imgData.ref = imcrop(obj.imgData.refScaled, obj.prcOpts.cropRect);

        if obj.prcOpts.dcSubtraction
            obj.imgData.ref = obj.subtractDC(obj.imgData.ref, ...
                obj.fringe.period);
        end

        nLines = length(obj.pCorrOpts.lineInds);
        obj.pCorrOpts.peaksRef = nan(1, nLines);
        obj.pCorrOpts.localPeriod = nan(1, nLines);

        if strcmpi(obj.pCorrOpts.method, 'spatial')
            [optInd, isSafeLine, peaksLineCell] = ...
                obj.findOptimalPeakInd(obj.pCorrOpts.peakInd);

            if isempty(obj.pCorrOpts.peakInd)
                if isnan(optInd)
                    error("Spatial phase correction: no peak lies inside the output " + ...
                        "frame with a safety margin of %.3g fringe periods from the " + ...
                        "'%s' edge on any phase-correction line. Enlarge the display " + ...
                        "window, lower pCorrOpts.edgeSafetyFactor, or adjust the peak " + ...
                        "detection settings.", ...
                        obj.pCorrOpts.edgeSafetyFactor, obj.pCorrOpts.startEdge)
                end
                obj.pCorrOpts.peakInd = optInd;

            elseif any(~isSafeLine)
                lineStr = strjoin(string(obj.pCorrOpts.lineInds(~isSafeLine)), ', ');
                if isnan(optInd)
                    warning("Peak %g is missing or does not lie safely inside the output frame " + ...
                        "on phase-correction line(s) %s, and no peak satisfies the safety margin " + ...
                        "of %.3g fringe periods from the '%s' edge. Enlarge the display window, " + ...
                        "lower pCorrOpts.edgeSafetyFactor, or adjust the peak detection settings.", ...
                        obj.pCorrOpts.peakInd, lineStr, obj.pCorrOpts.edgeSafetyFactor, obj.pCorrOpts.startEdge)
                else
                    warning("Peak %g is missing or does not lie safely inside the output frame " + ...
                        "on phase-correction line(s) %s. Optimal peak number: %g (closest usable " + ...
                        "peak to the '%s' edge with a safety margin of %.3g fringe periods). " + ...
                        "Call setTrackedPeakInd(%g), or setTrackedPeakInd([]) to select the peak automatically.", ...
                        obj.pCorrOpts.peakInd, lineStr, optInd, obj.pCorrOpts.startEdge, ...
                        obj.pCorrOpts.edgeSafetyFactor, optInd)
                end
            end

            % local wavelength and reference peak position on each line
            peakInd = obj.pCorrOpts.peakInd;
            for i = 1:nLines
                peaksLine = peaksLineCell{i};

                % the wavelength stencil needs a peak on each side of
                % peakInd; leave NaN (skipped by correctPhase) when the
                % required peaks do not exist
                if peakInd < 2 || peakInd > length(peaksLine) - 1
                    continue
                end

                obj.pCorrOpts.localPeriod(i) = abs( peaksLine(peakInd + 1) ...
                    - peaksLine(peakInd - 1) ) / 2;
                obj.pCorrOpts.peaksRef(i) = peaksLine(peakInd);

                if abs(obj.pCorrOpts.localPeriod(i) - obj.fringe.periodOrg)/obj.fringe.periodOrg > 0.5
                    warning("Relative difference between local period" + ...
                        " calculated for spatial phase correction and fringe period is greater than 50%.")
                end
            end
        end

        if obj.prcOpts.smoothRefImg
            if strcmpi(obj.fringe.normAxis, 'Y')
                obj.imgData.ref = imgaussfilt(obj.imgData.ref, [1e-06 5]);
            else
                obj.imgData.ref = imgaussfilt(obj.imgData.ref, [5 1e-06]);
            end
        end

    end

    function invalidFrame = preprocessCurrentImage(obj)
    % Bring imgData.curr into the computational frame.
    %   Interpolates bad pixels, applies the camera calibration
    %   (polynomial dewarp or pinhole undistort + rectify), flags the
    %   frame for discarding, locates the spatial phase-correction peaks
    %   in the RAW frame, then resizes, crops and removes the DC
    %   component.

        obj.imgData.curr = obj.rectify(obj.imgData.curr);
        cropRectOrg = obj.prcOpts.cropRectOrg;

        invalidFrame = false;
        % check if image should be discarded
        if ~isempty(obj.prcOpts.discardThreshold)
            imageMean = mean(imcrop(obj.imgData.curr, cropRectOrg), 'all', 'omitmissing');
            if imageMean > obj.prcOpts.discardThreshold
                invalidFrame = true;
            end
        end

        if strcmpi(obj.pCorrOpts.method, 'spatial') && ~invalidFrame
            obj.pCorrOpts.peaksCurr = zeros(1, length(obj.pCorrOpts.lineInds));
            for i = 1:length(obj.pCorrOpts.lineInds)
                [~, peakLocs] = detectFringePeaks(obj, obj.imgData.curr, obj.pCorrOpts.lineInds(i));

                if strcmpi(obj.pCorrOpts.startEdge, 'right') || strcmpi(obj.pCorrOpts.startEdge, 'bottom')
                    peakLocs = flip(peakLocs);
                end
                if length(peakLocs) >= obj.pCorrOpts.peakInd
                    obj.pCorrOpts.peaksCurr(i) = peakLocs(obj.pCorrOpts.peakInd);
                else
                    obj.pCorrOpts.peaksCurr(i) = nan;
                end
            end
        end
        if obj.prcOpts.resizeFactor ~= 1
            obj.imgData.curr = imresize(obj.imgData.curr, obj.prcOpts.resizeFactor);
        end

        % if ~isempty(obj.rotationAngle)
        %     obj.imgData.curr = imrotate(obj.imgData.curr, obj.rotationAngle, "bicubic");
        % end

        if ~isempty(obj.prcOpts.cropRect)
            obj.imgData.curr = imcrop(obj.imgData.curr, obj.prcOpts.cropRect);
        end

        % store the image without background subtraction if it is invalid
        if obj.prcOpts.dcSubtraction && ~invalidFrame
            obj.imgData.curr = obj.subtractDC(obj.imgData.curr, ...
                obj.fringe.period);
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

    function demodulate(obj)
        if strcmpi(obj.demodOpts.method, 'fourier')
            obj.demodulateFourier();
        elseif strcmpi(obj.demodOpts.method, 'wavelet')
            obj.demodulateWavelet();
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

        [ny, nx] = size(obj.imgData.ref);

        kxv = 2*pi/nx * [0:floor(nx/2), -ceil(nx/2)+1:-1];
        kyv = 2*pi/ny * [0:floor(ny/2), -ceil(ny/2)+1:-1];
        [kx, ky] = meshgrid(kxv, kyv);
    
        % carrier from the detected fringe properties
        nVec = obj.fringe.normVec;
        omega = 2*pi/obj.fringe.period;
        kxG = omega * nVec(1);
        kyG = omega * nVec(2);
    
        kr = sqrt((kx - kxG).^2 + (ky - kyG).^2);
        w  = omega * obj.demodOpts.filtWidthFrac;
        obj.phaseData.bandpassFilter = exp(-(kr/w).^8);
    
        obj.phaseData.refComplexCoeffs = ...
            ifft2(obj.phaseData.bandpassFilter .* fft2(obj.imgData.ref));
    end

    function demodulateFourier(obj)
    % Find the phase of the current image relative to the reference.
    % The band-pass filter and the complex reference coefficients are 
    % precomputed by demodulateRefFourier. Only the current image is
    % processed here.

        % Band-pass and inverse transform
        hI = ifft2(obj.phaseData.bandpassFilter .* fft2(obj.imgData.curr));

        % Phase difference relative to the reference
        delPhi = angle(hI .* conj(obj.phaseData.refComplexCoeffs));

        if strcmpi(obj.pCorrOpts.method, 'temporal')
            obj.phaseData.prev = obj.phaseData.curr;
        end

        obj.phaseData.curr = delPhi;
    end

    function demodulateRefWavelet(obj)
    % Wavelet ridge of the reference image.
    %   Runs the Morlet CWT over demodOpts.scaleList/angleList, picks the
    %   maximum-amplitude (scale, angle) per pixel, smooths that ridge and
    %   interpolates the complex coefficients between the four
    %   neighboring (scale, angle) planes. Result is stored in
    %   phaseData.refComplexCoeffs.

        [Ny, Nx] = size(obj.imgData.ref);
        % Demodulate reference image
        W4D = obj.morletCWT(obj.imgData.ref, ...
            obj.demodOpts.scaleList, obj.demodOpts.angleList, ...
            obj.demodOpts.sigma, obj.demodOpts.gamma, obj.fringe.normAxis);
        W4D_abs = abs(W4D);

        % For each pixel, find the (scale, angle) at which the amplitude of
        % the wavelet transform is maximal, then smooth the ridge before
        % extracting the phase
        nScales = size(W4D, 3);
        nAngles = size(W4D, 4);   % 1 for 1-D input regardless of angleList
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

    function demodulateWavelet(obj)
    % Wavelet phase of the current frame.
    %   Same ridge extraction and bilinear (scale, angle) interpolation as
    %   demodulateRefWavelet, applied to imgData.curr; phaseData.curr is
    %   set to the phase relative to the reference coefficients.

        [Ny, Nx] = size(obj.imgData.curr);
    
         W4D = obj.morletCWT(obj.imgData.curr, ...
             obj.demodOpts.scaleList, obj.demodOpts.angleList, ...
             obj.demodOpts.sigma, obj.demodOpts.gamma, obj.fringe.normAxis);
         W4D_abs = abs(W4D);
    
         % For each pixel, find the (scale, angle) at which the amplitude of
         % the wavelet transform is maximal, then smooth the ridge before
         % extracting the phase.
         nScales = size(W4D, 3);
         nAngles = size(W4D, 4);   % 1 for 1-D input regardless of angleList
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
    
         if strcmpi(obj.pCorrOpts.method, 'temporal')
             obj.phaseData.prev = obj.phaseData.curr;
         end

         obj.phaseData.curr = angle(W_interp .* conj(obj.phaseData.refComplexCoeffs));
    end

    function scaleAndTransformRef(obj)
        obj.imgData.refRaw = obj.rectify(obj.imgData.refRaw);
        obj.updateGeometry();
    end

    function imgRectified = rectify(obj, img)
        if ~isempty(obj.prcOpts.burntPixelMask)
            img = obj.interpBurntPixels(img);
        end

        if strcmpi(obj.cameraCalib.type, 'Polynomial')
            imgRectified = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, double(img), obj.cameraCalib.GSx, obj.cameraCalib.GSy, 'linear');
            imgRectified(isnan(imgRectified)) = 0;
        elseif strcmpi(obj.cameraCalib.type, 'Pinhole')
            imgRectified = undistortImage(img, obj.cameraCalib.intrinsicsMatlab);
            imgRectified = imwarp(imgRectified, obj.cameraCalib.pTransform);
        else            % no camera calibration
            imgRectified = img;
        end
    end

    function [phaseCrop, heightVec, weights, geom] = prepareCalibrationData( ...
            obj, heightVec, excludeHeights, weights)
    % Check and trim the calibration data and compute
    %   the calibration-region geometry.
    %
    %   Returns:
    %     phaseCrop - phase stack restricted to the calibration region and
    %                 to the retained heights (range applied, excluded and
    %                 zero heights removed)
    %     heightVec - matching column vector of heights
    %     weights   - matching column vector of fit weights
    %     geom      - geometry struct: calibRect (crop frame), rawRect (RAW
    %                 frame), embedRows/embedCols, pixel meshes, cropRect,
    %                 rawSize
    
        if obj.prcOpts.resizeFactor ~= 1 || obj.prcOpts.resizeFactorDisplay ~= 1
            error("Run solve() with resizeFactor and resizeFactorDisplay set to 1")
        end
        if size(obj.phaseData.stack, 3) ~= length(heightVec)
            error("Size of phase array incompatible with length of heightVec.")
        end
   
        heightVec = heightVec(~excludeHeights);
        heightVec = heightVec(:);
        weights   = weights(~excludeHeights);
        weights   = weights(:);
        
        % crop interior inside the unwrap margin
        mg       = obj.prcOpts.unwrapMarginOrg;
        cropRect = obj.prcOpts.cropRect;
    
        geom.cropRect  = cropRect;
        geom.calibRect = [mg + 1, mg + 1, ...
            cropRect(3) - 2*mg, cropRect(4) - 2*mg];     % crop frame
        geom.rawRect   = [cropRect(1) + mg, cropRect(2) + mg, ...
            cropRect(3) - 2*mg, cropRect(4) - 2*mg];     % RAW frame
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

    function unwrapAndFlag(obj)
    % Unwrap the current phase map and flag anomalies.
    %   Unwraps phaseData.curr row-wise then column-wise inside the
    %   unwrap margin (optionally refining with unwrap2D), detects
    %   residual gradient outliers, records their bounding boxes in
    %   postData.phaseAnomalies for the current frame, and finally
    %   resizes the phase map to the output frame.

        if ~obj.prcOpts.unwrapEnabled
            obj.phaseData.curr = imresize(obj.phaseData.curr, ...
                obj.prcOpts.resizeFactorDisplay, ...
                'bilinear', Antialiasing=false);
            return
        end
        [ny,nx] = size(obj.phaseData.curr);
        sumOutliers = inf;
        rf_o = obj.prcOpts.resizeFactorDisplay;
        margin = round(obj.prcOpts.unwrapMargin / rf_o); 

        % no margin along a dimension too short to have one (1-D data)
        mRow = margin * (ny > 2*margin);
        mCol = margin * (nx > 2*margin);
        rows = mRow + 1:ny - mRow;
        cols = mCol + 1:nx - mCol;

        tempDeltaPhi = obj.phaseData.curr;

        tempDeltaPhi(rows, cols) = ...
            unwrap(tempDeltaPhi(rows, cols), [], 1);
        tempDeltaPhi(rows, cols) = ...
            unwrap(tempDeltaPhi(rows, cols), [], 2);

        if strcmpi(obj.prcOpts.unwrapMethod, '2D')
            gMag = imgradient(tempDeltaPhi(rows, cols));
            gMagSmooth = medfilt2(gMag, [8 8], 'symmetric');
            residual = abs(gMag - gMagSmooth);
            resMed = median(residual, 'all');
            outlierMask = residual > 1000*resMed;
            sumOutliers = sum(outlierMask, 'all');

            if sumOutliers > 10
                obj.phaseData.curr(rows, cols) = ...
                    unwrap2D(tempDeltaPhi(rows, cols));
            else
                obj.phaseData.curr = tempDeltaPhi;
            end
        else
            obj.phaseData.curr = tempDeltaPhi;
        end
        if sumOutliers > 0
            gMag = imgradient(obj.phaseData.curr(rows, cols));
            gMagSmooth = medfilt2(gMag, [8 8], 'symmetric');
            residual = abs(gMag - gMagSmooth);
            resMed = median(residual, 'all');
            outlierMask = residual > 1000*resMed;
            sumOutliers = sum(outlierMask, 'all');
         
            if sumOutliers > 3
                outlierMask = imdilate(outlierMask, strel("disk", 8));
                outlierMask = imclose(outlierMask, strel("disk", 10));
                regions = regionprops(outlierMask);


                rectArr = cat(1, regions([regions.Area] > 3).BoundingBox);

                if ~isempty(rectArr)
                    dispRect = round(obj.prcOpts.cropRectDisplay); 
                    rectArr(:, 1:2) = (rectArr(:, 1:2) + [mCol, mRow] - 0.5)*rf_o + 0.5;
                    rectArr(:, 3:4) = rectArr(:, 3:4)*rf_o;
                    rectArr(:, 1:2) = rectArr(:, 1:2) - dispRect(1:2) + 1;

                    % Clip the bottom-right to limits of the OUTPUT frame
                    bottomRightLims = dispRect(3:4) + 1.5;
                    bottomRight = min(rectArr(:, 1:2) + rectArr(:, 3:4), ...
                                      bottomRightLims);

                    % Clip the top-left coordinates to a minimum of 1
                    rectArr(:, 1:2) = max(0.5, rectArr(:, 1:2));


                    % Recalculate width and height based on the clipped top-left
                    rectArr(:, 3:4) = bottomRight - rectArr(:, 1:2);

                    rectArr(any(rectArr(:, 3:4) <= 0, 2),:) = [];
                    
                    if ~isempty(rectArr)
                        obj.postData.phaseAnomalies(end + 1) = struct( ...
                                'timestep', obj.loopState.timestep, ...
                                'rects', rectArr);
                    end
                end
            end
        end
        
        % The margin is still wrapped, so its values jump at the edge of the
        % unwrapped region. Replace it with the nearest unwrapped values so the
        % resize blends only unwrapped phase into the pixels at that edge.
        obj.phaseData.curr = padarray(obj.phaseData.curr(rows, cols), ...
            [mRow, mCol], 'replicate');

        obj.phaseData.curr = imresize(obj.phaseData.curr, rf_o, ...
                                        'bilinear', Antialiasing=false);

        % Mark as invalid every resized pixel whose centre lies in the margin.
        % Output pixel k is centred at (k - 0.5)/rf_o + 0.5 in the input.
        [nyS, nxS] = size(obj.phaseData.curr);
        rowCentres = FtpSolver.rescaleCoord(1:nyS, 1/rf_o);
        colCentres = FtpSolver.rescaleCoord(1:nxS, 1/rf_o);
        obj.phaseData.curr(rowCentres < mRow + 0.5 | rowCentres > ny - mRow + 0.5, :) = NaN;
        obj.phaseData.curr(:, colCentres < mCol + 0.5 | colCentres > nx - mCol + 0.5) = NaN;
    end

    function correctPhase(obj)
    % Remove 2*pi ambiguities from the current phase map.
    %   'temporal': adds the modal 2*pi jump between this frame and the
    %   previous one, measured at three interior points. 'spatial': uses
    %   the tracked fringe peak on each correction line to predict the
    %   true phase there and adds the modal 2*pi offset. Frames where the
    %   estimate is unreliable are listed in pCorrOpts.warningTimesteps
    %   and left uncorrected.

        if strcmpi(obj.pCorrOpts.method, 'none')
            return
        end
 
        lineInds = obj.pCorrOpts.lineInds;
 
        if ~obj.loopState.isFirstTimestep && strcmpi(obj.pCorrOpts.method, 'temporal')
            [Ny, Nx] = size(obj.phaseData.curr);
            i_ind = floor(linspace(1, Ny, 5));
            i_ind = i_ind(2:4);
            j_ind = floor(linspace(1, Nx, 5));
            j_ind = j_ind(2:4);
            linInd = sub2ind([Ny, Nx], i_ind, j_ind);
 
            temporalUnwrap = unwrap([obj.phaseData.prev(linInd); obj.phaseData.curr(linInd)]);
            phaseDiff = temporalUnwrap(2,:) - obj.phaseData.curr(linInd);
            tol = 1e-08;
            phaseDiff = round(phaseDiff/tol)*tol;
            [temporalJump, F] = mode(phaseDiff);
            if F == 1
                obj.pCorrOpts.warningTimesteps(end + 1, :) = [obj.loopState.timestep, temporalJump];
                temporalJump = 0;
            end
        end
 
        if strcmpi(obj.pCorrOpts.method, 'temporal')
            if obj.loopState.isFirstTimestep
                return
            end

            if temporalJump ~= 0
                obj.phaseData.curr = obj.phaseData.curr + temporalJump;
            end
        else                        % spatial or manual phase correction
            if ~isempty(obj.pCorrOpts.manualVals)
                obj.phaseData.curr = obj.phaseData.curr + obj.pCorrOpts.manualVals(obj.loopState.timestep);
                return
            end

            peaksCurr = obj.pCorrOpts.peaksCurr;
            peaksRef = obj.pCorrOpts.peaksRef;
            localPeriod = obj.pCorrOpts.localPeriod;
 
            % keep only the lines on which both the reference and the
            % current peak were detected (failures are marked with NaN)
            validLine = ~isnan(peaksCurr) & ~isnan(peaksRef) & ~isnan(localPeriod);
 
            if ~any(validLine)
                obj.pCorrOpts.warningTimesteps(end + 1, :) = ...
                    [obj.loopState.timestep, nan(1, length(lineInds))];
                return
            end

            peaksCurr = peaksCurr(validLine);
            peaksRef = peaksRef(validLine);
            localPeriod = localPeriod(validLine);
            lineIndsUsed = lineInds(validLine);
 
            if strcmpi(obj.fringe.normAxis, 'X')
                targetPoint = obj.RAW2STACK([lineIndsUsed, peaksCurr']);
            elseif  strcmpi(obj.fringe.normAxis, 'Y')
                targetPoint = obj.RAW2STACK([peaksCurr', lineIndsUsed]);
            else
                error("Pattern normal axis not defined.")
            end
 
            % drop target points that fall outside the computational domain
            [NyPhase, NxPhase] = size(obj.phaseData.curr);
            inBounds = ( targetPoint(:, 1) >= 1 & targetPoint(:, 1) <= NyPhase & ...
                         targetPoint(:, 2) >= 1 & targetPoint(:, 2) <= NxPhase ).';
 
            if ~any(inBounds)
                warning("All phase-correction peaks fall outside the computational " + ...
                    "domain in image %g; phase correction skipped for this frame.", ...
                    obj.loopState.timestep)
                obj.pCorrOpts.warningTimesteps(end + 1, :) = ...
                    [obj.loopState.timestep, nan(1, length(lineInds))];
                return
            end
 
            targetPoint = targetPoint(inBounds, :);
            peaksCurr = peaksCurr(inBounds);
            peaksRef = peaksRef(inBounds);
            localPeriod = localPeriod(inBounds);
 
            realDeltaPhi = ( peaksRef - peaksCurr ) ...
                                ./ localPeriod ...
                                *2*pi;
            targetInds = sub2ind(size(obj.phaseData.curr), ...
                                    targetPoint(:, 1), targetPoint(:, 2));
            calcDeltaPhi = obj.phaseData.curr(targetInds);
            difference = realDeltaPhi - calcDeltaPhi';
            correction = 2*pi*round(difference/2/pi);
            tol = 1e-08;
            correction = round(correction/tol)*tol;
 
            [M, F] = mode(correction);
 
            if F == 1 && length(correction) > 2
                corrRow = nan(1, length(lineInds));
                usedLines = find(validLine);
                corrRow(usedLines(inBounds)) = correction;
                obj.pCorrOpts.warningTimesteps(end + 1, :) = [obj.loopState.timestep, corrRow];
                correction = 0;
            else
                correction = M;
            end
 
            obj.phaseData.curr = obj.phaseData.curr + correction;
        end
    end

    function calculateElevation(obj)
    % Calculates the surface elevation based on the selected model.
        if ~obj.outputConfig.elevOutput || strcmpi(obj.elevModel.type, 'none')
            return
        end

        if strcmpi(obj.elevModel.type, 'takeda')
            period = obj.fringe.period;
            pixelPitch = obj.worldCoords.scaling.mmPerPixel;

            obj.elevData.curr = ...
                FtpSolver.evalTakedaElev(obj.phaseData.curr, obj.elevModel.L, ...
                                        obj.elevModel.d, (period*pixelPitch)^-1);
        elseif strcmpi(obj.elevModel.type, 'poly')
            obj.elevData.curr = ...
                FtpSolver.evalPolyElev(obj.phaseData.curr, ...
                                       obj.elevModel.poly.coeffMat);
        end
        
        if obj.prcOpts.lateralShiftCorrection
            obj.resampleAtTrueCoords()
        end
    end

    function resampleAtTrueCoords(obj)
    % Resample elevData.curr onto the nominal world grid.
    %
    %   A pinhole camera sees each surface point along a slanted ray, so the
    %   pixel nominally associated with world position (x, y) actually samples
    %   the surface at a laterally shifted position whenever the surface is
    %   displaced from the reference plane z = 0. This method computes the
    %   true sample positions from the camera geometry and resamples the
    %   height field back onto the nominal mesh. Requires a pinhole calibration.
        
        if ~isfield(obj.cameraCalib, 'intrinsicsMatlab')
            return
        end

        % Crop away the unwrap margin
        mrg = obj.prcOpts.unwrapMargin;

        elevCrop = FtpSolver.cropOutMargin(obj.elevData.curr, mrg);
        u        = FtpSolver.cropOutMargin(obj.worldCoords.mesh.xPixel, mrg);
        v        = FtpSolver.cropOutMargin(obj.worldCoords.mesh.yPixel, mrg);
    
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

        kField = (zSign * elevCrop(:)' + rinvT(3)) ./ RKuv(3, :);
        res    = RKuv .* kField - rinvT;    % 3 x N true world positions
        
        X_corr = res(1, :);
        Y_corr = res(2, :);

        % handle NaN values
        ok = isfinite(X_corr) & isfinite(Y_corr) & isfinite(elevCrop(:)');

        % 1-D data
        if min(size(obj.imgData.ref)) < 2
            elevCrop = elevCrop(:);
            if strcmpi(obj.fringe.normAxis, 'X')
                obj.elevData.curr = interp1(X_corr(ok), elevCrop(ok), ...
                                            obj.worldCoords.mesh.x, 'linear');
            else
                obj.elevData.curr = interp1(Y_corr(ok), elevCrop(ok), ...
                                            obj.worldCoords.mesh.y, 'linear');
            end
            return
        end
    
        % 2-D data
        F = scatteredInterpolant(X_corr(ok)', Y_corr(ok)', elevCrop(ok)', ...
                                 'linear', 'none');
        obj.elevData.curr = F(obj.worldCoords.mesh.x, obj.worldCoords.mesh.y);
    end

    function storeData(obj)
        if obj.outputConfig.phaseOutput
            if obj.loopState.invalidFrame
                obj.phaseData.stack(:,:,obj.loopState.stackInd) = nan;
            else
                obj.phaseData.stack(:,:,obj.loopState.stackInd) = obj.phaseData.curr;
            end
        end

        if obj.outputConfig.elevOutput
            if obj.loopState.invalidFrame
                obj.elevData.stack(:,:,obj.loopState.stackInd) = nan;
            else
                obj.elevData.stack(:,:,obj.loopState.stackInd) = obj.elevData.curr;
            end
        end          

        if obj.outputConfig.imageOutput
            obj.imgData.stack(:,:,obj.loopState.stackInd) = obj.imgData.curr;
        end
    end

    function reportProgress(obj)
        if ~obj.outputConfig.reportProgressEnabled
            return
        end
        textstring = sprintf('Processing image %5d...', obj.loopState.timestep);
        if obj.loopState.timestep > obj.prcOpts.solveRange(1)
            textstringOld = sprintf('Processing image %5d...', obj.loopState.timestep - 1);
            fprintf(repmat('\b', 1, numel(textstringOld)))
            fprintf(textstring)
        else
            % first timestep
            fprintf(textstring)
        end
    end

    function writeResultFiles(obj, addr)
        solveRange = obj.prcOpts.solveRange;
        if obj.outputConfig.chunkedOutput
            if obj.outputConfig.elevOutput
                obj.mergeChunkFiles("elevData", ...
                    fullfile(addr, obj.caseID + "_elevData.bin"));
            end
            if obj.outputConfig.phaseOutput
                obj.mergeChunkFiles("phaseData", ...
                    fullfile(addr, obj.caseID + "_phaseData.bin"));
            end
        else
            if obj.outputConfig.elevOutput
                filename = fullfile(addr, obj.caseID + "_elevData.bin");

                displayArr = obj.getElevDisplaySubarray(solveRange(1), solveRange(2));

                FtpSolver.writeTimeSeriesBinary(filename, displayArr);
            end

            if obj.outputConfig.phaseOutput
                filename = fullfile(addr, obj.caseID + "_phaseData.bin");

                [displayArr, ~, ~] = obj.getPhaseDisplaySubarray(solveRange(1), solveRange(2));

                FtpSolver.writeTimeSeriesBinary(filename, displayArr);
            end
        end
    
        [rowInds, colInds] = obj.outputInds();
        xMesh = obj.worldCoords.mesh.x(rowInds,colInds);
        yMesh = obj.worldCoords.mesh.y(rowInds,colInds);
        
        temp.xMesh = xMesh;
        temp.yMesh = yMesh;
        save(fullfile(addr, "coordinateMesh.mat"), '-struct', 'temp')
        
        if obj.outputConfig.chunkedOutput
            rmdir(obj.caseChunkFolder(), 's')
        end
    end
    
    function writeChunk(obj)
    % Write the chunk held in memory to the chunk folder.
    %   Writes the output subarray of the elevation and phase stacks for
    %   the current chunk to one file each, then reallocates the stacks
    %   for the next chunk unless this was the last one.

        folder = obj.caseChunkFolder();
        if ~isfolder(folder)
            mkdir(folder);
        end

        solveRange = obj.prcOpts.solveRange;
        chunkIndex = obj.loopState.chunkIndex;
        chunkLength = obj.outputConfig.chunkLength;
        nTimesteps = diff(solveRange) + 1;
        nChunks = ceil(nTimesteps / chunkLength);

        if obj.outputConfig.elevOutput
            chunkFilename = obj.chunkFilePath(chunkIndex, "elevData");

            chunkArr = obj.getElevDisplaySubarray(solveRange(1), ...
                solveRange(1) + size(obj.elevData.stack, 3) - 1);

            FtpSolver.writeTimeSeriesBinary(chunkFilename, chunkArr);
        end

        if obj.outputConfig.phaseOutput
            chunkFilename = obj.chunkFilePath(chunkIndex, "phaseData");

            chunkArr = obj.getPhaseDisplaySubarray(solveRange(1), ...
                solveRange(1) + size(obj.phaseData.stack, 3) - 1);

            FtpSolver.writeTimeSeriesBinary(chunkFilename, chunkArr);
        end

        if chunkIndex == nChunks
            return
        end

        % the last chunk holds whatever timesteps remain
        if chunkIndex == nChunks - 1
            nextChunkLength = nTimesteps - chunkIndex*chunkLength;
        else
            nextChunkLength = chunkLength;
        end
        obj.allocateStacks(nextChunkLength)
    end

    function folder = caseChunkFolder(obj)
        % Subfolder of chunkFolder holding this case's chunks.
        folder = fullfile(obj.outputConfig.chunkFolder, "TEMP_chunks_" + obj.caseID);
    end

    function path = chunkFilePath(obj, k, quantity)
        % Path of chunk k of one quantity.
        path = fullfile(obj.caseChunkFolder(), ...
            sprintf("chunk_%05d_%s.bin", k, quantity));
    end

    function mergeChunkFiles(obj, quantity, outFile)
        % Concatenate the chunk files of one quantity to a single file.
        %   Checks every chunk file of a chunked run, then writes them in
        %   chunk order to outFile as one binary in the readBinary format,
        %   with dimensions [Ny Nx nTimesteps]. Nothing is written if a 
        %   chunk is missing, has a different frame size, or the chunks do
        %   not add up to the timesteps in solveRange.
        %   Inputs:  
        %       quantity - "elevData" or "phaseData"
        %       outFile  - path of the merged binary

        nTimesteps = diff(obj.prcOpts.solveRange) + 1;
        nChunks = ceil(nTimesteps / obj.outputConfig.chunkLength);

        chunkFiles  = strings(nChunks, 1);
        chunkTimesteps = zeros(nChunks, 1);
        for k = 1:nChunks
            chunkFiles(k) = obj.chunkFilePath(k, quantity);
            if ~isfile(chunkFiles(k))
                error("mergeChunkFiles: chunk file %s is missing.", chunkFiles(k))
            end

            dims = FtpSolver.readBinaryHeader(chunkFiles(k));
            if k == 1
                frameSize = dims(1:end - 1);
            elseif ~isequal(dims(1:end - 1), frameSize)
                error("mergeChunkFiles: chunk %d has frame size %s, expected %s.", ...
                    k, mat2str(dims(1:end - 1)), mat2str(frameSize))
            end
            chunkTimesteps(k) = dims(end);
        end

        if sum(chunkTimesteps) ~= nTimesteps
            error("mergeChunkFiles: the chunks hold %d frames but solveRange has %d.", ...
                sum(chunkTimesteps), nTimesteps)
        end

        outID = fopen(outFile, 'w');
        if outID < 0
            error("mergeChunkFiles: cannot open %s for writing.", outFile)
        end
        closeOut = onCleanup(@() fclose(outID));

        FtpSolver.writeBinaryHeader(outID, [frameSize, nTimesteps])
        for k = 1:nChunks
            fwrite(outID, FtpSolver.readBinary(chunkFiles(k)), 'single');
        end
    end

    function initImages(obj)
        % Restore the image buffers of a case loaded from a .mat.
        %   Image data and polynomial coefficients are not saved with the
        %   object, so this re-reads the reference image and the dataset 
        %   listing from the stored paths, rebuilds the geometry, and
        %   reloads the polynomial calibration when the case used one.

        obj.loadRefImage(obj.source.refPath);
        obj.loadDataset(obj.source.dataPath);

        obj.scaleAndTransformRef();
        if strcmpi(obj.elevModel.type, 'poly') && isfield(obj.elevModel, 'polyAddr')
            obj.setElevModelPoly(obj.elevModel.polyAddr)
        end
    end

    function initPhaseCorrection(obj)      
        cropRectOrg = obj.prcOpts.cropRectOrg;
        if strcmpi(obj.fringe.normAxis, 'Y')
            targetInds = round(linspace(cropRectOrg(1), cropRectOrg(1) + cropRectOrg(3), 5));
            obj.pCorrOpts.lineInds = unique(targetInds(2:end - 1))';
        else                                                    % normal vector to pattern is parallel to X axis
            targetInds = round(linspace(cropRectOrg(2), cropRectOrg(2) + cropRectOrg(4), 5));
            obj.pCorrOpts.lineInds = unique(targetInds(2:end - 1))';
        end
    end

    function storePinholeScaling(obj, scalingStruct)
        if obj.prcOpts.imgRotAngle == 0
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

    function imgOut = interpBurntPixels(obj, img)
        if strcmpi(obj.prcOpts.interpBurntPixelMethod, 'fastLinear')
            imgOut = obj.fastFill(img, obj.prcOpts.burntPixelRows, ...
                obj.prcOpts.burntPixelCols);
        elseif strcmpi(obj.prcOpts.interpBurntPixelMethod, 'smooth')
            imgOut = regionfill(img, obj.prcOpts.burntPixelMask);
        end
    end

    function [rowInds, colInds] = outputInds(obj)
        % Row and column indices of the OUTPUT frame
        r = round(obj.prcOpts.cropRectDisplay);
        [Ny, Nx] = size(obj.worldCoords.mesh.x);
        rowInds = r(2) : min(Ny, r(2) + r(4));
        colInds = r(1) : min(Nx, r(1) + r(3));
    end

    function output = RAW2STACK(obj, input)
        % Convert (row, col) values from RAW to STACK frame
        rf = obj.prcOpts.resizeFactor;
        rf_o = obj.prcOpts.resizeFactorDisplay;
        cropRect = obj.prcOpts.cropRect;

        rows = input(:,1);
        cols = input(:,2);

        % RAW to SCALED
        rowsScaled = FtpSolver.rescaleCoord(rows, rf);
        colsScaled = FtpSolver.rescaleCoord(cols, rf);

        % SCALED to COMP
        rowsCrop = rowsScaled - cropRect(2) + 1;
        colsCrop = colsScaled - cropRect(1) + 1;

        % COMP to STACK
        rowsStack = FtpSolver.rescaleCoord(rowsCrop, rf_o);
        colsStack = FtpSolver.rescaleCoord(colsCrop, rf_o);

        output = round([rowsStack, colsStack]);
    end

    function [peakVals, peakLocs] = detectFringePeaks(obj, img, lineInd)
        cropRectOrg = obj.prcOpts.cropRectOrg;
        if strcmpi(obj.fringe.normAxis, 'X')
            cropStart = cropRectOrg(1);
            cropEnd = cropRectOrg(1) + cropRectOrg(3);
            phaseCorrLine = double(img(lineInd, :));
        else
            cropStart = cropRectOrg(2);
            cropEnd = cropRectOrg(2) + cropRectOrg(4);
            phaseCorrLine = double(img(:, lineInd));
        end

        segment = phaseCorrLine(cropStart:cropEnd);
        [minHeight, minProm] = FtpSolver.peakThresholds(segment, ...
            obj.pCorrOpts.peakMinFrac, obj.pCorrOpts.peakPromFrac);

        [peakVals, peakLocs] = findpeaks(phaseCorrLine, ...
            'MinPeakDistance', 0.6*obj.fringe.periodOrg, ...
            'MinPeakProminence', minProm, ...
            'MinPeakHeight', minHeight);
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
        legend('Measurements', 'Model')

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

    function img = loadImage(obj, timestep)
    % Read frame timestep of the dataset as a double image.
        switch obj.source.dataFormat
            case "set"
                img = obj.getDavisFrame(obj.source.dataPath, obj.source.dataCamInd, timestep);
            case "im7"
                img = obj.getDavisFrame(obj.source.frameFiles(timestep), obj.source.dataCamInd);
            case "image"
                img = imread(obj.source.frameFiles(timestep));
            otherwise
                error("No dataset loaded.")
        end
        img = double(img);
    end

    function loadRefImage(obj, refPath)
    % Read the reference image into imgData.refRaw.
    %   Averages over the whole set if given a Davis .set file.
        refPath = FtpSolver.absolutePath(refPath);
        obj.source.refPath = refPath;
        [~, ~, ext] = fileparts(refPath);
    
        if strcmpi(ext, ".set")
            img = obj.getDavisFrame(refPath, obj.source.refCamInd, 'avg');
        elseif strcmpi(ext, ".im7")
            img = obj.getDavisFrame(refPath, obj.source.refCamInd);
        elseif ismember(lower(extractAfter(ext, ".")), [imformats().ext])
            img = imread(refPath);
        else
            error("Invalid file type for reference image.")
        end
        obj.imgData.refRaw = double(img);
    end
    
    function loadDataset(obj, dataPath)
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
        if isempty(obj.prcOpts.solveRange)
            obj.prcOpts.solveRange = obj.source.fullRange;
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
            warning("The reference image is in the data folder and is " + ...
                "included as frame %d.", find(isRef))
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

    function savedObj = saveobj(obj)
        savedObj.caseID = obj.caseID;
        savedObj.source = obj.source;
        savedObj.imgData = obj.defaultImgData();         
        savedObj.phaseData = obj.defaultPhaseData();    
        savedObj.elevData = obj.defaultElevData();
        savedObj.fringe = obj.fringe;
        savedObj.pCorrOpts = obj.pCorrOpts;           
        savedObj.demodOpts = obj.demodOpts;
        savedObj.prcOpts = obj.prcOpts;
        savedObj.loopState = obj.loopState;
        savedObj.cameraCalib = obj.cameraCalib;         
        savedObj.worldCoords = obj.worldCoords;
        savedObj.worldCoords.mesh = [];
        savedObj.elevModel = obj.elevModel;
        savedObj.outputConfig = obj.outputConfig;
        savedObj.postData = obj.postData;
        if isfield(savedObj.elevModel, 'poly')
            savedObj.elevModel = rmfield(savedObj.elevModel, 'poly');
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
    end

    function output = defaultPCorrOpts()
        output.method = "temporal";
        output.lineInds = [];
        output.arr = [];
        output.startEdge = "left";
        output.edgeSafetyFactor = 2.5;
        output.peakInd = [];
        output.manualVals = [];
        output.peakPromFrac = 1/4;
        output.peakMinFrac = 1/3;
        output.warningTimesteps = [];
    end

    function output = defaultDemodOpts()
        output.method = "fourier";
        output.scaleList = [];
        output.angleList = 0;
        output.sigma = 0.6;
        output.gamma = 1;
        output.filtWidthFrac = 0.6;
    end

    function output = defaultElevData()
        output.stack = [];
        output.curr = [];
    end

    function output = defaultFringe()
        output.normVec = [];
        output.normAxis = "";
        output.periodOrg = [];
        output.period = [];
        output.tiltAngle = [];
    end

    function output = defaultPrcOpts()
        output.solveRange = [];
        output.unwrapEnabled = true;
        output.unwrapMethod = "1D";
        output.unwrapMarginOrg = 20;
        output.unwrapMargin = [];
        output.imgRotAngle = 0; % in degrees and positive clockwise
        output.lateralShiftCorrection = false;
        output.discardThreshold = [];
        output.dcSubtraction = true;
        output.smoothRefImg = false;
        output.interpBurntPixelMethod = 'fastLinear';
        output.burntPixelMask = [];
        output.burntPixelRows = [];
        output.burntPixelCols = [];
        output.cropRectOrg = [];
        output.cropRect = [];
        output.cropRectDisplayOrg = [];
        output.cropRectDisplay = [];
        output.resizeFactor = 1;
        output.resizeFactorDisplay = 0.5;
    end

    function output = defaultLoopState()
        output.stackInd = 0;
        output.timestep = [];
        output.isFirstTimestep = false;
        output.chunkIndex = 1;
        output.processTime = [];
        output.invalidFrame = false;
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
        output.phaseOutput = false;
        output.imageOutput = false;
        output.elevOutput = true;
        output.reportProgressEnabled = true;
        output.chunkedOutput = false;
        output.chunkLength = 1000;
        output.chunkFolder = "";
    end
     
    function writeTimeSeriesBinary(path, arr, nSpaceDims)
        % Write an array whose last dimension is time.
        %   Creates (or overwrites) the file at path. The header always holds
        %   nSpaceDims + 1 dimensions.
        arguments
            path {mustBeTextScalar}
            arr
            nSpaceDims (1,1) double {mustBeInteger, mustBePositive} = 2
        end

        if ndims(arr) > nSpaceDims + 1
            error("writeTimeSeriesBinary: array has %d dimensions, but nSpaceDims = %d allows at most %d.", ...
                ndims(arr), nSpaceDims, nSpaceDims + 1)
        end

        fID = fopen(path, 'w');
        if fID < 0
            error("writeTimeSeriesBinary: cannot open %s for writing.", path)
        end
        closeFile = onCleanup(@() fclose(fID));

        FtpSolver.writeBinaryHeader(fID, size(arr, 1:nSpaceDims + 1));
        fwrite(fID, arr, 'single');
    end

    function writeBinaryHeader(fID, dims)
        % Write the dimension header of a result binary.
        fwrite(fID, numel(dims), 'uint32');
        fwrite(fID, dims, 'uint32');
    end

    function dims = readBinaryHeader(path)
        % Dimensions stored in the header of a result binary.
        fID = fopen(path, 'r');
        if fID < 0
            error("readBinaryHeader: cannot open %s.", path)
        end
        closeIn = onCleanup(@() fclose(fID));
        nDims = fread(fID, 1, 'uint32');
        dims  = fread(fID, nDims, 'uint32')';
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

    function croppedArr = cropOutMargin(arr, margin)
        [Ny, Nx] = size(arr);

        if Ny > 2*margin
            y = margin + 1 : Ny - margin;
        else
            y = 1:Ny;
        end

        if Nx > 2*margin
            x = margin + 1 : Nx - margin;
        else
            x = 1:Nx;
        end

        croppedArr = arr(y, x);
    end

    function out = rescaleCoord(x, scale)
        % Map a coordinate between frames related by imresize.
        % scale about the image edge
        out = (x - 0.5)*scale + 0.5;
    end

    function [minHeight, minProm] = peakThresholds(lineVals, peakMinFrac, peakPromFrac)
        % Fringe mean level and peak-to-trough amplitude from the intensity
        % distribution
        medianLevel = median(lineVals, 'omitmissing');
        amplitude = prctile(lineVals, 95) - prctile(lineVals, 5);
        minHeight = peakMinFrac  * medianLevel;
        minProm   = peakPromFrac * amplitude;
    end

    function [image, scaling] = getDavisFrame(addr, camInd, timestep)
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
                        error("getDavisFrame() cannot handle images " +  ...
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
end

%% ----------- Public static helpers -----------
methods(Static)
    function newObj = loadobj(fileObj)
        if isstruct(fileObj)
            newObj = FtpSolver("temp");
            newObj.caseID = fileObj.caseID;
            newObj.source = fileObj.source;
            newObj.imgData = fileObj.imgData;         
            newObj.phaseData = fileObj.phaseData;    
            newObj.elevData = fileObj.elevData;
            newObj.fringe = fileObj.fringe;
            newObj.pCorrOpts = fileObj.pCorrOpts;           
            newObj.demodOpts = fileObj.demodOpts;
            newObj.prcOpts = fileObj.prcOpts;
            newObj.loopState = fileObj.loopState;
            newObj.cameraCalib = fileObj.cameraCalib;         
            newObj.worldCoords = fileObj.worldCoords;             
            newObj.elevModel = fileObj.elevModel;
            newObj.outputConfig = fileObj.outputConfig;
            newObj.postData = fileObj.postData;
        else
            newObj = fileObj;
        end
    end

    function [data, dataSizeAll] = readBinary(addr, seekTime)
    % Read back binary data file created by the class.
    % Inputs:
    %       addr -      Path to the binary file. 
    %       seekTime -  (1x2) vector of timesteps to return in the format
    %                   [firstTimestep lastTimestep]. Reads the full binary
    %                   if not specified.
    % Outputs:
    %       data -          The data array.
    %       dataSizeAll -   The dimensions of the the full binary array.
    %
        if ~exist('addr', 'var')
            [file, path] = uigetfile('*.bin');
            addr = fullfile(path, file);
        end
        fID = fopen(addr, 'rb');
        closeIn = onCleanup(@() fclose(fID));

        numDim = fread(fID, 1, 'uint32');
        dataSize = fread(fID, numDim, 'uint32');
        dataSize = dataSize';
        dataSizeAll = dataSize;

        if exist('seekTime', 'var')
            if dataSize(end) < seekTime(2)
                error("Requested timestep %g exceeds number of timesteps in array (%g).", seekTime(2), dataSize(end))
            end
            dataSize(end) = seekTime(2) - seekTime(1) + 1;
            % assumes single precision (4 bytes)
            pageElements = prod(dataSize(1:end - 1));
            skipLength = 4*pageElements*(seekTime(1) - 1);
            readElements = pageElements*(seekTime(2) - seekTime(1) + 1);
            % skip to seekTime(1)
            fseek(fID, skipLength, 'cof');
            data = fread(fID, readElements, '*single');
        else
            data = fread(fID, '*single');
        end

        data = reshape(data, dataSize);            
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

    function [period, normAxis, tiltAngle, normVec] = analyzeFringe(img)
    % Estimate fringe period and orientation via FFT
    %
    %   Applies a Hann window to the image and locates the dominant
    %   spectral peak outside the low-frequency core.
    %
    %   period      - fringe period in pixels
    %   normAxis - 'X' if the pattern normal is closer to the x axis,
    %                 'Y' otherwise
    %   tiltAngle   - signed angle in degrees, in [-45, 45], between the
    %                 pattern normal (carrier wave vector) and the selected
    %                 axis; 0 when the fringes are perfectly aligned.
    %                 Equivalently, the tilt of the fringe lines away from
    %                 the perpendicular axis. Positive when the normal is
    %                 rotated from the selected positive axis toward the
    %                 other positive axis.
    %   normVec     - fringe normal unit vector [vx, vy]
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
            normAxis = 'X';       % fx > 0 by construction of the half-plane
            tiltAngle = atand(fy/fx);
        else
            normAxis = 'Y';
            if fy < 0                % report the conjugate with fy > 0
                fx = -fx;
                fy = -fy;
            end
            tiltAngle = atand(fx/fy);
        end
        normVec = [fx, fy]/norm([fx, fy]);
    end

    function W4D = morletCWT(img, scales, anglesDeg, sigma, gamma, fringeNormAxis)
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
    %   If fringeNormAxis ('X' or 'Y') is provided, pads the array in that
    %   direction only and in both directions otherwise.

        if ~exist('fringeNormAxis', 'var')
            fringeNormAxis = '';
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
            if strcmpi(fringeNormAxis, 'X')
                paddedImg = padarray(img, [0 padSize], 0);
                paddedImg = padarray(paddedImg, [padSize 0], 'symmetric');
            elseif strcmpi(fringeNormAxis, 'Y')
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

            % cfs is [rows x cols x 1 x nScales x nAngles]; drop the
            % singleton plane dimension and crop the padding.
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
        elev = L*phaseDiff ./ (2*pi*f0*d + phaseDiff);
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
