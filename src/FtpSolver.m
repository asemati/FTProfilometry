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
%                             refAddr="path/to/refImage.tiff", ...
%                             dataAddr="path/to/data.tiff", ...
%                             resizeFactor=0.7, ...
%                             cropRect=[100 150 1500 1200], ...
%                             camCalibAddr="path/to/calibration.mat" ...
%                             );
%     obj.setROI();
%     obj.solve();
%     obj.animate(1, 100, framesPerSecond);
%     obj.writeCase(outputDir);
%
%   KEY METHODS
%     solve             - Run the full processing pipeline
%     setROI            - Interactively define the computational domain
%     calibratePoly     - Calibrate phase-to-height model
%     probeElev         - Extract time series at a spatial location
%     animate           - Visualize the reconstructed surface
%     writeCase         - Save results
%
%   DEPENDENCIES
%     xml2struct, rectifyImagePinhole, readDavisCalibration
%     unwrap2D (only for prcOpts.unwrapMethod = '2D')
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
% April 2024; Last revision: 17-August-2026

%------------- BEGIN CODE --------------
    
    properties (SetAccess = public)
        caseID          
        inputData           %
        imgState            % image buffers
        phaseData           % computed phase
        surfData            % computed surface elevation
        demodOpts           % demodulation settings        
        pCorrOpts           % phase correction settings and state
        prcOpts             % preprocessing and pipeline settings and flags
        loopState           % timestep bookkeeping
        cameraCalib         % camera calibration data
        worldCoords         % spatial coordinates
        surfParams          % phase to elevation conversion parameters
        outputConfig        % output/save control
        postData            % outlier info and results of user-defined function
    end
    
methods
    function obj = FtpSolver(caseID, opts)
        arguments
            caseID              {mustBeTextScalar}              = "default"
            opts.refAddr        {mustBeTextScalar}              = ""
            opts.dataAddr       {mustBeTextScalar}              = ""
            opts.mmPerPixel     (1,1) double                    = NaN
            opts.period         (1,1) double                    = NaN
            opts.patNormAxis    string {mustBeMember(opts.patNormAxis, ["X", "Y"])} = "X"
            opts.resizeFactor   (1,1) double {mustBePositive}   = 1
            opts.cropRect       double                          = []
            opts.camCalibAddr   {mustBeTextScalar}              = ""
            opts.camCalibNum    (1,1) double                    = 1    
            opts.imgFrameNum    (1,1) {mustBeInteger}           = 1
            opts.refImgFrameNum (1,1) {mustBeInteger}           = 1
            opts.imgRotAngle    (1,1) double                    = 0
            opts.interpPixelMask                                = []
        end

        if ~isvarname(caseID)
            error("caseID must be a valid variable name.")
        end

        obj.caseID      = char(caseID);
        obj.inputData   = FtpSolver.defaultInput();
        obj.imgState    = FtpSolver.defaultImgState();
        obj.phaseData   = FtpSolver.defaultPhaseData();
        obj.surfData    = FtpSolver.defaultSurfData();
        obj.pCorrOpts   = FtpSolver.defaultPCorrOpts();
        obj.demodOpts   = FtpSolver.defaultDemodOpts();
        obj.prcOpts     = FtpSolver.defaultPrcOpts();
        obj.loopState   = FtpSolver.defaultLoopState();
        obj.cameraCalib = FtpSolver.defaultCameraCalib();
        obj.worldCoords = FtpSolver.defaultWorldCoords();
        obj.outputConfig = FtpSolver.defaultOutputConfig();            


        obj.prcOpts.resizeFactor        = opts.resizeFactor;
        obj.prcOpts.cropRectOrg         = opts.cropRect;
        obj.prcOpts.imgRotAngle         = opts.imgRotAngle;
        obj.prcOpts.interpPixelMask     = opts.interpPixelMask;
        obj.inputData.imgFrameNum       = opts.imgFrameNum;
        obj.inputData.refImgFrameNum    = opts.refImgFrameNum;

        if opts.refAddr == ""
            return   % bare construction - used by loadobj and init()
        end
        

        % load reference image
        obj.loadRefImage(opts.refAddr)
        
        % load data
        obj.loadDataset(opts.dataAddr)

        if ~isnan(opts.mmPerPixel)
            obj.worldCoords.scaling.X.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.Y.SlopeOrg      = opts.mmPerPixel;
            obj.worldCoords.scaling.mmPerPixelOrg   = opts.mmPerPixel;
        end

        if strlength(opts.camCalibAddr) > 0
            obj.readCameraCalibration(opts.camCalibAddr, opts.camCalibNum);
        end   

        obj.surfParams.periodOrg = opts.period;
        obj.surfParams.period = opts.period;

        obj.prcOpts.patNormAxis = opts.patNormAxis;
        if strcmpi(opts.patNormAxis, 'X')
            obj.prcOpts.patNormVec = [1, eps];
        elseif strcmpi(opts.patNormAxis, 'Y')
            obj.prcOpts.patNormVec = [eps, 1];
        end

        if ~isempty(obj.prcOpts.interpPixelMask)
            [obj.prcOpts.interpPixelRows, ...
             obj.prcOpts.interpPixelCols] = find(obj.prcOpts.interpPixelMask);
        end
                    
        obj.scaleAndTransformRef();

        % auto-detect the fringe period and pattern axis when none given
        if isnan(opts.period)
            img = imcrop(obj.imgState.refRaw, obj.prcOpts.cropRectOrg);
            [period, patNormAxis, ~, normVec] = obj.analyzeFringe(img);
            obj.surfParams.periodOrg = period;
            obj.prcOpts.patNormAxis = patNormAxis;
            obj.prcOpts.patNormVec = normVec;
            obj.updateGeometry();
        end
    end

%%
    function loadRefImage(obj, refAddr)
        % Read the reference image into imgState.refRaw
        % Averages over the whole set if given a Davis .set file
        obj.inputData.refAddr = refAddr;
        [~, ~, refExt] = fileparts(refAddr);
        refExt = lower(extractAfter(refExt, '.'));
        obj.inputData.refFiletype = refExt;

        if strcmpi(refExt, 'set') 
            obj.imgState.refRaw = obj.getDavisFrame(obj.inputData.refAddr, obj.inputData.refImgFrameNum, 'avg');
        elseif strcmpi(refExt, 'im7')
            obj.imgState.refRaw = obj.getDavisFrame(obj.inputData.refAddr, obj.inputData.refImgFrameNum);
        elseif ismember(refExt, [imformats().ext])
            obj.imgState.refRaw = imread(obj.inputData.refAddr);
        else
            error("Invalid file type for reference image.")
        end

        obj.imgState.refRaw = double(obj.imgState.refRaw);

    end
%%
    function loadDataset(obj, dataAddr)
        %LOADDATASET Register the image sequence to be processed.
        %   Stores either the .set file path or a dir() listing of all
        %   .im7/image files in the same folder. Sets inputData.fullRange
        %   from the number of frames found and initializes solveRange when 
        %   empty.
        %   In:  dataAddr - path to a .set file or to one file of the sequence

        [filepath, ~, dataExt] = fileparts(dataAddr);
        dataExt = lower(extractAfter(dataExt, '.'));
        obj.inputData.dataFiletype = dataExt;

        if strcmpi(dataExt, 'set') 
            obj.inputData.dataAddr = convertStringsToChars(dataAddr);
            obj.inputData.fullRange = [1, lvsetsize(obj.inputData.dataAddr)];
        elseif strcmpi(dataExt, 'im7')
            obj.inputData.dataAddr = dir(fullfile(filepath, '*.im7'));
            obj.inputData.fullRange = [1, length(obj.inputData.dataAddr)];
        elseif ismember(dataExt, [imformats().ext])
            obj.inputData.dataAddr = dir(fullfile(filepath, "*." + dataExt));
            obj.inputData.fullRange = [1, length(obj.inputData.dataAddr)];
        else
            error("Invalid file type for dataset.")
        end

        if isempty(obj.inputData.solveRange)
            obj.inputData.solveRange = obj.inputData.fullRange;
        end

    end
%%
    function initImages(obj)
        %INITIMAGES Restore the image buffers of a case loaded from a .mat.
        %   Image data and polynomial coefficients are not not saved with 
        %   the object, so this re-reads the reference image and the 
        %   dataset listing from the stored paths, rebuilds the geometry, 
        %   and reloads the polynomial calibration when the case used one.

        obj.loadRefImage(obj.inputData.refAddr);

        if isstruct(obj.inputData.dataAddr)
            field = obj.inputData.dataAddr(1);
            addr = fullfile(field.folder, field.name);
            obj.loadDataset(addr)
        else
            obj.loadDataset(obj.inputData.dataAddr)
        end

        obj.scaleAndTransformRef();
        if isfield(obj.surfParams, 'profPolynomialAddr')
            obj.setProfModelPoly(obj.surfParams.profPolynomialAddr)
        end
    end
%%        
    function setUnwrapMargin(obj, val)
        % margin is in RAW-frame pixels
        obj.prcOpts.unwrapMarginOrg = val;
        obj.updateGeometry();
    end
%%        
    function solve(obj)
        %SOLVE Run the full processing pipeline over inputData.solveRange.
        %   Prepares the reference image, then for every frame: load,
        %   preprocess, demodulate, unwrap, phase-correct, convert to
        %   elevation and store. Handles frame discarding and block-wise
        %   writing to disk. Results end up in surfData/phaseData.

        t0 = tic();
        
        fprintf("Case %s\n", obj.caseID)

        obj.initParams();       
        obj.preprocessRefImage();
        obj.demodulateRef();

        for loopIndex = obj.inputData.solveRange(1):obj.inputData.solveRange(2)
            obj.loopState.surfLoopInd = loopIndex;
            obj.loopState.surfDataInd = obj.loopState.surfDataInd + 1;

            if obj.prcOpts.useSegmentation && obj.loopState.surfDataInd > obj.prcOpts.blockSize
                obj.writeTempBlock()
                obj.loopState.surfDataInd = obj.loopState.surfDataInd - obj.prcOpts.blockSize;
                obj.loopState.blockNumber = obj.loopState.blockNumber + 1;
            end
            obj.reportProgress();
            
            obj.loadNextImage();
            obj.preprocessCurrentImage();
            
            if obj.loopState.discardCurr
                obj.storeData()
                obj.loopState.discardCurr = false;
                continue
            end
            
            obj.demodulate()
            obj.unwrapAndFlag();
            obj.correctPhase();
            obj.calculateSurface();

            obj.userFunction();
            obj.storeData();

            obj.loopState.firstTimestep = false;
        end
        fprintf('\nDone\n')

        obj.writeTempBlock()
        obj.loopState.prcTime = toc(t0);
    end

%%
    function initParams(obj)
        %INITPARAMS Reset loop state and allocate everything solve() needs.
        %   Refreshes the geometry (and the camera mesh of a reloaded case),
        %   preallocates the surface/phase stacks for the whole run or for
        %   one block, resamples the polynomial calibration onto the
        %   computational grid, and optionally equalizes the reference
        %   exposure against the first data frame.

        obj.loopState.surfDataInd = 0;
        obj.loopState.blockNumber = 1;
        obj.loopState.firstTimestep = true;
        obj.surfData.curr = [];
        obj.phaseData.curr = [];
        obj.pCorrOpts.warningFrames = [];

        % refresh coordinate mesh if running a case loaded from a .mat file
        if isfield(obj.cameraCalib, 'addr') && ~isfield(obj.worldCoords.mesh, 'xOrg')
            if strlength(obj.cameraCalib.addr) > 0 && strcmpi(obj.cameraCalib.type, 'Pinhole')
                obj.readCameraCalibration(obj.cameraCalib.addr, obj.cameraCalib.camNum);
            end
        end
        % refresh images if running a case loaded from a .mat file
        if strlength(obj.inputData.refAddr) > 0 && isempty(obj.imgState.refRaw)
            obj.initImages();
        end
        obj.updateGeometry();

        rf = obj.prcOpts.resizeFactor;
        rf_o = obj.prcOpts.resizeFactorDisplay;
        [NyRef, NxRef] = size(obj.imgState.ref);
        NxStack = ceil(rf_o * NxRef);
        NyStack = ceil(rf_o * NyRef);
        Nt = diff(obj.inputData.solveRange) + 1;
        blockSize = obj.prcOpts.blockSize;

        if obj.prcOpts.useSegmentation && Nt > blockSize
            if obj.outputConfig.savePhase
                obj.phaseData.stack = single(nan(NyStack, NxStack, blockSize));
            end

            if obj.outputConfig.saveSurf
                obj.surfData.stack = single(nan(NyStack, NxStack, blockSize));
            end
        else
            if obj.outputConfig.savePhase
                obj.phaseData.stack = single(nan(NyStack, NxStack, Nt));
            end

            if obj.outputConfig.saveSurf
                obj.surfData.stack = single(nan(NyStack, NxStack, Nt));
            end
        end

        if ~isfield(obj.surfParams,'period')
            error("Missing data: surfParams")
        end
        
        % initialize polynomial calibration model
        [Ny, Nx] = size(obj.imgState.refRaw);
        if isfield(obj.surfParams, 'profPolynomial')
            if obj.surfParams.profPolynomial.pixelwise
                coeffMat = imresize(obj.surfParams.profPolynomial.coeffMatOrg, ...
                                    rf, 'bilinear', Antialiasing=false);
            else
                [xMesh, yMesh] = meshgrid(1:Nx, 1:Ny);
                fPlanes = obj.surfParams.profPolynomial.fittedPlanes;
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
                size(obj.surfParams.profPolynomial.coeffMatOrg, 3) - 1 ...
                ];
            coeffMat = imcrop3(coeffMat, cropRect3);
            coeffMat = imresize(coeffMat, rf_o, 'bilinear', Antialiasing=false);
            obj.surfParams.profPolynomial.coeffMat = coeffMat;
        end

        % normalize refImgRaw in case of a large difference in camera
        % exposure between ref and data frames
        if obj.prcOpts.equalizeExposure
            obj.loopState.surfLoopInd = obj.inputData.solveRange(1);
            obj.loadNextImage();
            cropRectOrg = obj.prcOpts.cropRectOrg;

            if isfield(obj.cameraCalib, 'GSx')
                obj.imgState.curr = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, ...
                                        double(obj.imgState.curr), ...
                                        obj.cameraCalib.GSx, ...
                                        obj.cameraCalib.GSy, ...
                                        'linear');
                obj.imgState.curr(isnan(obj.imgState.curr)) = 0;
            elseif isfield(obj.cameraCalib, 'intrinsicsMatlab')
                obj.imgState.curr = undistortImage(obj.imgState.curr, obj.cameraCalib.intrinsicsMatlab);
                obj.imgState.curr = imwarp(obj.imgState.curr, obj.cameraCalib.pTransform);
            end

            meanCurr = mean(imcrop(obj.imgState.curr, cropRectOrg), 'all', 'omitmissing');
            meanRef = mean(imcrop(obj.imgState.refRaw, cropRectOrg), 'all', 'omitmissing');
            ratio = meanCurr/meanRef;

            if ratio > 1.2 || ratio < 0.8
                warning("Rescaling reference image by k = %g due to a difference in exposure compared to dataset.", ratio)
                obj.imgState.refRaw = obj.imgState.refRaw*ratio;
                obj.imgState.refScaled = obj.imgState.refScaled*ratio;
            end
        end

        obj.postData.phaseAnomalies = struct('frame', {}, 'rects', {});
    end

%%
    function preprocessRefImage(obj)
        %PREPROCESSREFIMAGE Prepare the reference image and the spatial
        %   phase-correction reference.
        %   Crops the scaled reference to the computational domain, removes
        %   the DC component, and for spatial phase correction selects or
        %   validates pCorrOpts.peakInd and measures the local fringe
        %   wavelength and reference peak position on every correction line.

        obj.imgState.ref = imcrop(obj.imgState.refScaled, obj.prcOpts.cropRect);

        if ~isequal(size(obj.pCorrOpts.arr), size(obj.imgState.refRaw))
            error("pCorrOpts.arr and imgState.refRaw arrays are not the same size.")
        end

        if obj.prcOpts.dcSubtraction
            obj.imgState.ref = obj.subtractDC(obj.imgState.ref, ...
                                              obj.surfParams.period);
        end

        nLines = length(obj.pCorrOpts.lineInds);
        obj.pCorrOpts.peaksRef = nan(1, nLines);
        obj.pCorrOpts.wvLength = nan(1, nLines);

        if strcmp(obj.pCorrOpts.method, 'spatial')
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
                        "Set obj.pCorrOpts.peakInd = %g and rerun, or set it to [] to select " + ...
                        "the optimal peak automatically.", ...
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

                obj.pCorrOpts.wvLength(i) = abs( peaksLine(peakInd + 1) ...
                                - peaksLine(peakInd - 1) ) / 2;
                obj.pCorrOpts.peaksRef(i) = peaksLine(peakInd);

                if abs(obj.pCorrOpts.wvLength(i) - obj.surfParams.periodOrg)/obj.surfParams.periodOrg > 0.5
                    warning("Relative difference between local period" + ...
                    " calculated for spatial phase correction and fringe period is greater than 50%.")
                end
            end
        end

        if obj.prcOpts.smoothRefImg
            if strcmpi(obj.prcOpts.patNormAxis, 'Y')
                obj.imgState.ref = imgaussfilt(obj.imgState.ref, [1e-06 5]);
            else
                obj.imgState.ref = imgaussfilt(obj.imgState.ref, [5 1e-06]);
            end
        end

    end
%%
    function [optInd, isSafeLine, peaksLineCell] = findOptimalPeakInd(obj, peakInd)
        %FINDOPTIMALPEAKIND Evaluate phase-correction peaks against the output frame.
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
        %                   when it cannot be determined (required inputs not
        %                   in place yet, or no eligible peak on any line).
        %   isSafeLine    - per-line logical: true where the candidate peakInd
        %                   is eligible. All false when peakInd is omitted or
        %                   empty.
        %   peaksLineCell - detected peak positions per line, RAW-frame
        %                   coordinates, flipped for startEdge 'right'/'bottom'
        %
        %   Query only: does not modify pCorrOpts.peakInd. Note: updates
        %   pCorrOpts.peakMin as a side effect of peak detection.

        if nargin < 2
            peakInd = [];
        end

        nLines = length(obj.pCorrOpts.lineInds);
        optInd = nan;
        isSafeLine = false(1, nLines);
        peaksLineCell = cell(1, nLines);

        % prerequisites: reference image, correction lines, fringe period,
        % pattern normal axis, and output-frame geometry must all exist
        if isempty(obj.imgState.refRaw) || nLines == 0 ...
                || ~isfield(obj.surfParams, 'periodOrg') ...
                || isempty(obj.surfParams.periodOrg) ...
                || isnan(obj.surfParams.periodOrg) ...
                || strlength(string(obj.prcOpts.patNormAxis)) == 0 ...
                || isempty(obj.prcOpts.cropRectOrg) ...
                || isempty(obj.prcOpts.cropRectDisplayOrg)
            return
        end

        if isempty(obj.pCorrOpts.peakMin)
            obj.pCorrOpts.peakMin = ones(1, nLines)*mean(obj.imgState.refRaw, 'all')/3;
            obj.pCorrOpts.peakMin = round(obj.pCorrOpts.peakMin);
        elseif isscalar(obj.pCorrOpts.peakMin)
            obj.pCorrOpts.peakMin = ones(1, nLines)*obj.pCorrOpts.peakMin;
        end

        if length(obj.pCorrOpts.peakMin) ~= nLines
            obj.pCorrOpts.peakMin = ones(1, nLines)*mean(obj.imgState.refRaw, 'all')/3;
            obj.pCorrOpts.peakMin = round(obj.pCorrOpts.peakMin);
        end

        % safe window along the pattern normal axis, RAW frame
        cropRectOrg = obj.prcOpts.cropRectOrg;
        dispRectOrg = obj.prcOpts.cropRectDisplayOrg;
        if strcmp(obj.prcOpts.patNormAxis, 'X')
            winStart = cropRectOrg(1) + dispRectOrg(1) - 1;
            winEnd   = winStart + dispRectOrg(3);
        else
            winStart = cropRectOrg(2) + dispRectOrg(2) - 1;
            winEnd   = winStart + dispRectOrg(4);
        end

        safetyDist = obj.pCorrOpts.edgeSafetyFactor * obj.surfParams.periodOrg;

        % the starting edge is winEnd for 'right'/'bottom' (where the peak
        % trains are flipped) and winStart otherwise
        if strcmp(obj.pCorrOpts.startEdge, 'right') || strcmp(obj.pCorrOpts.startEdge, 'bottom')
            safeLo = winStart;
            safeHi = winEnd - safetyDist;
        else
            safeLo = winStart + safetyDist;
            safeHi = winEnd;
        end

        obj.pCorrOpts.peakMinHint = nan(1, nLines);
        obj.pCorrOpts.peakPromHint = nan(1, nLines);

        % detect and evaluate every line
        optIndsPerLine = nan(1, nLines);
        for i = 1:nLines
                if strcmp(obj.prcOpts.patNormAxis, 'X')
                    phaseCorrLine = double(obj.imgState.refRaw(obj.pCorrOpts.lineInds(i), :));
                else
                    phaseCorrLine = double(obj.imgState.refRaw(:, obj.pCorrOpts.lineInds(i)));
                end
            [~, initPeaks] = findpeaks(phaseCorrLine, ...
                'MinPeakDistance', 0.6*obj.surfParams.periodOrg, ...
                'MinPeakProminence', obj.pCorrOpts.peakProm);
            [~, troughX] = findpeaks(-phaseCorrLine, ...
                'MinPeakDistance', 0.6*obj.surfParams.periodOrg, ...
                'MinPeakProminence', obj.pCorrOpts.peakProm);
            peakVals = phaseCorrLine(initPeaks);
            troughVals = phaseCorrLine(troughX);
            vecLength = min(length(peakVals), length(troughVals));
            meanVal = (peakVals(1:vecLength) + troughVals(1:vecLength))/2;
            amplitude = (peakVals(1:vecLength) - troughVals(1:vecLength));

            obj.pCorrOpts.peakMinHint(i) = round(median(meanVal)/3);
            obj.pCorrOpts.peakPromHint(i) = round(median(amplitude)/4);

            [~, peaksLine] = findpeaks(phaseCorrLine, ...
                'MinPeakDistance', 0.6*obj.surfParams.periodOrg, ...
                'MinPeakProminence', obj.pCorrOpts.peakProm, ...
                'MinPeakHeight', obj.pCorrOpts.peakMin(i));

            if strcmp(obj.pCorrOpts.startEdge, 'right') || strcmp(obj.pCorrOpts.startEdge, 'bottom')
                peaksLine = flip(peaksLine);
            end

            peaksLineCell{i} = peaksLine;

            % inside the safe window, with a neighbor on each side
            eligible = peaksLine >= safeLo & peaksLine <= safeHi;
            if ~isempty(eligible)
                eligible(1) = false;
                eligible(end) = false;
            end

            firstEligible = find(eligible, 1, 'first');
            if ~isempty(firstEligible)
                optIndsPerLine(i) = firstEligible;
            end

            if ~isempty(peakInd) && peakInd >= 2 && peakInd <= length(peaksLine) - 1
                isSafeLine(i) = eligible(peakInd);
            end
        end

        if ~all(isnan(optIndsPerLine))
            optInd = max(optIndsPerLine);   % NaN lines are ignored by max
        end
    end

%%
    function preprocessCurrentImage(obj)
        %PREPROCESSCURRENTIMAGE Bring imgState.curr into the computational frame.
        %   Interpolates bad pixels, applies the camera calibration
        %   (polynomial dewarp or pinhole undistort + rectify), flags the
        %   frame for discarding, locates the spatial phase-correction peaks
        %   in the RAW frame, then resizes, crops and removes the DC
        %   component.

        if ~isempty(obj.prcOpts.interpPixelMask)
            obj.imgState.curr = obj.interpBadPixels(obj.imgState.curr);
        end
        
        if isfield(obj.cameraCalib, 'GSx')
            obj.imgState.curr = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, ...
                                    double(obj.imgState.curr), ...
                                    obj.cameraCalib.GSx, ...
                                    obj.cameraCalib.GSy, ...
                                    'linear');
            obj.imgState.curr(isnan(obj.imgState.curr)) = 0;
        elseif isfield(obj.cameraCalib, 'intrinsicsMatlab')
            obj.imgState.curr = undistortImage(obj.imgState.curr, obj.cameraCalib.intrinsicsMatlab);
            obj.imgState.curr = imwarp(obj.imgState.curr, obj.cameraCalib.pTransform);
        end

        % check if image should be discarded
        if ~isempty(obj.prcOpts.discardThreshold)
            cropRectOrg = obj.prcOpts.cropRectOrg;
            imageMean = mean(imcrop(obj.imgState.curr, cropRectOrg), 'all');
            if imageMean > obj.prcOpts.discardThreshold
                obj.loopState.discardCurr = true;
            end
        end
        
        if strcmp(obj.pCorrOpts.method, 'spatial')
            obj.pCorrOpts.peaksCurr = zeros(1, length(obj.pCorrOpts.lineInds));
            for i = 1:length(obj.pCorrOpts.lineInds)
                if strcmp(obj.prcOpts.patNormAxis, 'X')
                    phaseCorrLine = double(obj.imgState.curr(obj.pCorrOpts.lineInds(i), :));
                else
                    phaseCorrLine = double(obj.imgState.curr(:, obj.pCorrOpts.lineInds(i)));
                end
                [~, peaksLine] = findpeaks(phaseCorrLine, ...
                        'MinPeakDistance', 0.6*obj.surfParams.periodOrg, ...
                        'MinPeakProminence', obj.pCorrOpts.peakProm, ...
                        'MinPeakHeight', obj.pCorrOpts.peakMin(i));

                if strcmp(obj.pCorrOpts.startEdge, 'right') || strcmp(obj.pCorrOpts.startEdge, 'bottom')
                    peaksLine = flip(peaksLine);
                end
                if length(peaksLine) >= obj.pCorrOpts.peakInd
                    obj.pCorrOpts.peaksCurr(i) = peaksLine(obj.pCorrOpts.peakInd);
                else
                    obj.pCorrOpts.peaksCurr(i) = nan;
                end
            end
        end
        if obj.prcOpts.resizeFactor ~= 1
            obj.imgState.curr = imresize(obj.imgState.curr, obj.prcOpts.resizeFactor);
        end

        % if ~isempty(obj.rotationAngle)
        %     obj.imgState.curr = imrotate(obj.imgState.curr, obj.rotationAngle, "bicubic");
        % end

        if ~isempty(obj.prcOpts.cropRect)
            obj.imgState.curr = imcrop(obj.imgState.curr, obj.prcOpts.cropRect);
        end
        
        if obj.prcOpts.dcSubtraction
            obj.imgState.curr = obj.subtractDC(obj.imgState.curr, ...
                                               obj.surfParams.period);
        end
    end
%%
    function [image, scaling] = getDavisFrame(obj, addr, frameNum, time)
        %GETDAVISFRAME Read one frame (or the set average) from a Davis file.
        %   Undoes the zero padding and area-of-interest offset that Davis
        %   applies, using the RealFrameSize/AOIused/CameraMaxNx attributes.
        %   In:  addr     - path to the .set or .im7 file
        %        frameNum - frame number inside the buffer
        %        time     - set index, 0/omitted for a single image, or
        %                   'avg' to average the whole set
        %   Out: image    - image array
        %        scaling  - Davis scaling struct of the last frame read

            addr = convertStringsToChars(addr);

            if ~exist('time', 'var')
                time = 0;
            elseif strcmp(time, 'avg')
                time = 1:lvsetsize(addr);
            end

            for i = 1:length(time)
                if time(i) ~= 0
                    temp = readimx(addr, time(i));
                else
                    temp = readimx(addr);
                end
    
                image = temp.Frames{frameNum}.Components{1}.Planes{1}';

                % Crop image if Davis has padded it with zeros
                rfsIdx = [];
                % reposition image if AOI used
                aoiIdx = [];
                maxNxIdx = [];
                for k = 1:length(temp.Frames{frameNum}.Attributes)
                    if strcmp(temp.Frames{frameNum}.Attributes{k}.Name, 'RealFrameSize')
                        rfsIdx = k;
                    end
                    if strcmp(temp.Frames{frameNum}.Attributes{k}.Name, 'AOIused')
                        aoiIdx = k;
                    end
                    if strcmp(temp.Frames{frameNum}.Attributes{k}.Name, 'CameraMaxNx')
                        maxNxIdx = k;
                    end
                end
                if ~isempty(rfsIdx)
                    realFrameX = temp.Frames{frameNum}.Attributes{rfsIdx}.Value(1);
                    realFrameY = temp.Frames{frameNum}.Attributes{rfsIdx}.Value(2);
                    if ~isequal(size(image), [realFrameY, realFrameX])
                        image = image(1:realFrameY, 1:realFrameX);
                    end
                elseif obj.loopState.surfDataInd == 0
                    warning("getDavisFrame: Unable to find 'RealFrameSize' attribute in Davis image");
                end
    
                if ~isempty(aoiIdx)
                    roiX = temp.Frames{frameNum}.Attributes{aoiIdx}.Value(1) + 1;
                    roiY = temp.Frames{frameNum}.Attributes{aoiIdx}.Value(2) + 1;
                    binX = temp.Frames{frameNum}.Attributes{aoiIdx}.Value(3);
                    binY = temp.Frames{frameNum}.Attributes{aoiIdx}.Value(4);
                    Nx = str2double(temp.Frames{frameNum}.Attributes{maxNxIdx}.Value);
                    Ny = str2double(temp.Frames{frameNum}.Attributes{maxNxIdx + 1}.Value);

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

            if length(time) > 1
                image = imageSum/length(time);
            end

            scaling = temp.Frames{frameNum}.Scales;
            scaling.X.SlopeOrg = scaling.X.Slope;
            scaling.Y.SlopeOrg = scaling.Y.Slope;
        end
%%
    function setRange(obj, r1, r2)
        %SETRANGE Set the frame range processed by solve().
        %   In:  r1 - first frame, or 'full' to use the whole dataset
        %        r2 - last frame (ignored when r1 is 'full')

        if strcmp(r1, 'full')
            obj.inputData.solveRange = obj.inputData.fullRange;
            return
        end
        if r2 > obj.inputData.fullRange(2) || r1 < 1
            error("Invalid input. Available data range is [1 %g]", obj.inputData.fullRange(2));
        end
        obj.inputData.solveRange = [r1 r2];

    end

%%
    function demodulateRef(obj)
        if strcmpi(obj.demodOpts.method, 'FT')
            obj.demodulateRefFT();
        elseif strcmpi(obj.demodOpts.method, 'wavelet')
            obj.demodulateRefWavelet();
        else
            error("Invalid demodulation method")
        end
    end

%%
    function demodulate(obj)
        if strcmpi(obj.demodOpts.method, 'FT')
            obj.demodulateFT();
        elseif strcmpi(obj.demodOpts.method, 'wavelet')
            obj.demodulateWavelet();
        else
            error("Invalid demodulation method")
        end
    end
%%
    function demodulateRefFT(obj)
        % Precompute all reference-dependent quantities for FT demodulation.
        %
        % Sets:
        %   obj.phaseData.bandpassFilter   - super-Gaussian band-pass centered on carrier
        %   obj.phaseData.refComplexCoeffs - band-passed complex reference signal,
        %                                    hI0 = ifft2(filter .* fft2(I0))

        [ny, nx] = size(obj.imgState.ref);

        kxv = 2*pi/nx * [0:floor(nx/2), -ceil(nx/2)+1:-1];
        kyv = 2*pi/ny * [0:floor(ny/2), -ceil(ny/2)+1:-1];
        [kx, ky] = meshgrid(kxv, kyv);
    
        % carrier from the detected fringe properties
        nVec = obj.prcOpts.patNormVec;
        omega = 2*pi/obj.surfParams.period;
        kxG = omega * nVec(1);
        kyG = omega * nVec(2);
    
        kr = sqrt((kx - kxG).^2 + (ky - kyG).^2);
        w  = omega * obj.demodOpts.filtWidthFrac;
        obj.phaseData.bandpassFilter = exp(-(kr/w).^8);
    
        obj.phaseData.refComplexCoeffs = ...
            ifft2(obj.phaseData.bandpassFilter .* fft2(obj.imgState.ref));
    end

%%
    function demodulateFT(obj)
        % Find the phase of the current image relative to the reference.
        % All reference-dependent quantities (carrier location, band-pass
        % filter, complex reference coefficients) are precomputed by
        % demodulateRefFT; only the current image is transformed here.

        % Band-pass and inverse transform
        hI = ifft2(obj.phaseData.bandpassFilter .* fft2(obj.imgState.curr));

        % Phase difference relative to the reference
        delPhi = angle(hI .* conj(obj.phaseData.refComplexCoeffs));

        if strcmpi(obj.pCorrOpts.method, 'temporal')
            obj.phaseData.old = obj.phaseData.curr;
        end

        obj.phaseData.curr = delPhi;
    end

%%
    function demodulateRefWavelet(obj)
        %DEMODULATEREFWAVELET Wavelet ridge of the reference image.
        %   Runs the Morlet CWT over demodOpts.scaleList/angleList, picks the
        %   maximum-amplitude (scale, angle) per pixel, smooths that ridge and
        %   interpolates the complex coefficients between the four
        %   neighboring (scale, angle) planes. Result is stored in
        %   phaseData.refComplexCoeffs.

        [Ny, Nx] = size(obj.imgState.ref);
        % Demodulate reference image
        W4D = obj.morletCWT(obj.imgState.ref, ...
            obj.demodOpts.scaleList, obj.demodOpts.angleList, ...
            obj.demodOpts.sigma, obj.demodOpts.gamma, obj.prcOpts.patNormAxis);
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
     
        % freqMode = mode(scale_idx_raw, 'all');
        % obj.surfParams.W = 6/obj.demodOpts.scaleList(freqMode);  % 6 = Morlet Omega0    
    end
 %%
    function demodulateWavelet(obj)
        %DEMODULATEWAVELET Wavelet phase of the current frame.
        %   Same ridge extraction and bilinear (scale, angle) interpolation as
        %   demodulateRefWavelet, applied to imgState.curr; phaseData.curr is
        %   set to the phase relative to the reference coefficients.

        [Ny, Nx] = size(obj.imgState.curr);
    
         W4D = obj.morletCWT(obj.imgState.curr, ...
             obj.demodOpts.scaleList, obj.demodOpts.angleList, ...
             obj.demodOpts.sigma, obj.demodOpts.gamma, obj.prcOpts.patNormAxis);
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
             obj.phaseData.old = obj.phaseData.curr;
         end

         obj.phaseData.curr = angle(W_interp .* conj(obj.phaseData.refComplexCoeffs));
    end
%%
    function setWaveletScales(obj, minWavelength, maxWavelength, divisions)
        % set wavelet scales based on pattern wavelength in pixels
        % wavelength should be for the raw, unscaled image
        scaleMin = minWavelength*6/(2*pi) * obj.prcOpts.resizeFactor;
        scaleMax = maxWavelength*6/(2*pi) * obj.prcOpts.resizeFactor;
        obj.demodOpts.scaleList = linspace(scaleMin, scaleMax, divisions);
    end

%%
    function setWaveletAngles(obj, angleList)
        % restrict angles to the range (-135, 45] for consistency with FT method
        obj.demodOpts.angleList = mod(angleList - 45, -180) + 45;
    end

%%
    function unwrapAndFlag(obj)
        %UNWRAPANDFLAG Unwrap the current phase map and flag anomalies.
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

        tempDeltaPhi = obj.phaseData.curr;
        rows = margin + 1:ny - margin;
        cols = margin + 1:nx - margin;

        tempDeltaPhi(rows, cols) = ...
            unwrap(tempDeltaPhi(rows, cols), [], 1);
        tempDeltaPhi(rows, cols) = ...
            unwrap(tempDeltaPhi(rows, cols), [], 2);

        if strcmp(obj.prcOpts.unwrapMethod, '2D')
            gMag = imgradient(tempDeltaPhi(rows, cols));
            gMagSmooth = medfilt2(gMag, [8 8]);
            residual = abs(gMag - gMagSmooth);
            resMed = median(residual, 'all');
            outlier_mask = residual > 1000*resMed;
            sumOutliers = sum(outlier_mask, 'all');

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
            gMagSmooth = medfilt2(gMag, [8 8]);
            residual = abs(gMag - gMagSmooth);
            resMed = median(residual, 'all');
            outlier_mask = residual > 1000*resMed;
            sumOutliers = sum(outlier_mask, 'all');
         
            if sumOutliers > 3
                outlier_mask = imdilate(outlier_mask, strel("disk", 8));
                outlier_mask = imclose(outlier_mask, strel("disk", 10));
                regions = regionprops(outlier_mask);


                rectArr = cat(1, regions([regions.Area] > 3).BoundingBox);

                if ~isempty(rectArr)
                    dispRect = round(obj.prcOpts.cropRectDisplay); 
                    rectArr(:, 1:2) = (rectArr(:, 1:2) + margin - 0.5)*rf_o + 0.5;
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
                                'frame', obj.loopState.surfLoopInd, ...
                                'rects', rectArr);
                    end
                end
            end
        end

        obj.phaseData.curr = imresize(obj.phaseData.curr, ...
                                      obj.prcOpts.resizeFactorDisplay, ...
                                      'bilinear', Antialiasing=false);
    end
    %%
    function correctPhase(obj)
        %CORRECTPHASE Remove 2*pi ambiguities from the current phase map.
        %   'temporal': adds the modal 2*pi jump between this frame and the
        %   previous one, measured at three interior points. 'spatial': uses
        %   the tracked fringe peak on each correction line to predict the
        %   true phase there and adds the modal 2*pi offset. Frames where the
        %   estimate is unreliable are listed in pCorrOpts.warningFrames and
        %   left uncorrected.

        if strcmp(obj.pCorrOpts.method, 'none')
            return
        end
 
        lineInds = obj.pCorrOpts.lineInds;
 
        if ~obj.loopState.firstTimestep && strcmpi(obj.pCorrOpts.method, 'temporal')
            [Ny, Nx] = size(obj.phaseData.curr);
            i_ind = floor(linspace(1, Ny, 5));
            i_ind = i_ind(2:4);
            j_ind = floor(linspace(1, Nx, 5));
            j_ind = j_ind(2:4);
            linInd = sub2ind([Ny, Nx], i_ind, j_ind);
 
            temporalUnwrap = unwrap([obj.phaseData.old(linInd); obj.phaseData.curr(linInd)]);
            phaseDiff = temporalUnwrap(2,:) - obj.phaseData.curr(linInd);
            tol = 1e-08;
            phaseDiff = round(phaseDiff/tol)*tol;
            [temporalJump, F] = mode(phaseDiff);
            if F == 1
                obj.pCorrOpts.warningFrames(end + 1, :) = [obj.loopState.surfLoopInd, temporalJump];
                temporalJump = 0;
            end
        end
 
        if strcmp(obj.pCorrOpts.method, 'temporal')
            if obj.loopState.firstTimestep
                return
            end
            % [jump, frequency] = mode(phaseDiffTemporal);
            if temporalJump ~= 0
                obj.phaseData.curr = obj.phaseData.curr + temporalJump;
            end
        else                        % spatial or manual phase correction
            if ~isempty(obj.pCorrOpts.manualVals)
                obj.phaseData.curr = obj.phaseData.curr + obj.pCorrOpts.manualVals(obj.loopState.surfLoopInd);
                return
            end
 
            peaksCurr = obj.pCorrOpts.peaksCurr;
            peaksRef = obj.pCorrOpts.peaksRef;
            wvLength = obj.pCorrOpts.wvLength;
 
            % keep only the lines on which both the reference and the
            % current peak were detected (failures are marked with NaN)
            validLine = ~isnan(peaksCurr) & ~isnan(peaksRef) & ~isnan(wvLength);
 
            if ~any(validLine)
                obj.pCorrOpts.warningFrames(end + 1, :) = ...
                    [obj.loopState.surfLoopInd, nan(1, length(lineInds))];
                return
            end
 
            peaksCurr = peaksCurr(validLine);
            peaksRef = peaksRef(validLine);
            wvLength = wvLength(validLine);
            lineIndsUsed = lineInds(validLine);
 
            if strcmp(obj.prcOpts.patNormAxis, 'X')
                targetPoint = obj.RAW2STACK([lineIndsUsed, peaksCurr']);
            elseif  strcmp(obj.prcOpts.patNormAxis, 'Y')
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
                    obj.loopState.surfLoopInd)
                obj.pCorrOpts.warningFrames(end + 1, :) = ...
                    [obj.loopState.surfLoopInd, nan(1, length(lineInds))];
                return
            end
 
            targetPoint = targetPoint(inBounds, :);
            peaksCurr = peaksCurr(inBounds);
            peaksRef = peaksRef(inBounds);
            wvLength = wvLength(inBounds);
 
            realDeltaPhi = ( peaksRef - peaksCurr ) ...
                                ./ wvLength ...
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
                obj.pCorrOpts.warningFrames(end + 1, :) = [obj.loopState.surfLoopInd, corrRow];
                correction = 0;
            else
                correction = M;
            end
 
            obj.phaseData.curr = obj.phaseData.curr + correction;
        end
    end

%%
    function calculateSurface(obj)
        if ~obj.outputConfig.saveSurf
            return
        end

        if strcmp(obj.surfParams.profModel, 'takeda')
            period = obj.surfParams.period;
            pixelPitch = obj.worldCoords.scaling.mmPerPixel;
            d = obj.surfParams.d;
            L = obj.surfParams.L;
            obj.surfData.curr = obj.phaseData.curr*L ./ ...
                            (2*pi/period/pixelPitch*d + obj.phaseData.curr);
        elseif strcmp(obj.surfParams.profModel, 'poly')
            ord = size(obj.surfParams.profPolynomial.coeffMat, 3) - 1;
            n = reshape(0:ord, 1, 1, []);
            obj.surfData.curr = sum(obj.surfParams.profPolynomial.coeffMat .* ...
                                    obj.phaseData.curr.^n, 3);
        end
        
        if obj.prcOpts.lateralShiftCorrection
            obj.resampleAtTrueCoords()
        end
    end
    %%
    function resampleAtTrueCoords(obj)
    %RESAMPLEATTRUECOORDS Resample surfData.curr onto the nominal world grid.
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

        surfCrop = FtpSolver.cropOutMargin(obj.surfData.curr, mrg);
        u        = FtpSolver.cropOutMargin(obj.worldCoords.mesh.xPixel, mrg);
        v        = FtpSolver.cropOutMargin(obj.worldCoords.mesh.yPixel, mrg);
    
        % Back-project every pixel to a world-frame camera ray
        K_inv = inv(obj.cameraCalib.intrinsicsMatlab.K);
        R     = obj.cameraCalib.extrinsics.R;
        t     = obj.cameraCalib.extrinsics.Translation';
    
        pix   = [u(:), v(:), ones(numel(u), 1)]';   % 3 x N homogeneous pixels
        RKuv  = R' * (K_inv * pix);                 % 3 x N ray directions
        rinvT = R' * t;                             % camera-position term
    
        % surfData is positive when the surface is raised with respect to the
        % reference plane. If the world z-axis points away from the camera,
        % positive world z corresponds to a surface depression; -sign(rinvT(3))
        % resolves the convention (positive when z points towards the camera)
        zSign = -sign(rinvT(3));

        kField = (zSign * surfCrop(:)' + rinvT(3)) ./ RKuv(3, :);
        res    = RKuv .* kField - rinvT;    % 3 x N true world positions

        % handle 1-D data
        if min(size(obj.imgState.ref)) < 2            
            if strcmpi(obj.prcOpts.patNormAxis, 'X')
                obj.surfData.curr = interp1(res(1, :), surfCrop(1, :), ...
                                            obj.worldCoords.mesh.x, 'linear', 'extrap');
            else
                obj.surfData.curr = interp1(res(2, :), surfCrop(:, 1), ...
                                            obj.worldCoords.mesh.y, 'linear', 'extrap');
            end
            return
        end
    
        % 2-D resampling
        % Compute the true world position of every sample explicitly and
        % rebuild a Delaunay triangulation each frame     
        X_corr = res(1, :);
        Y_corr = res(2, :);

        F = scatteredInterpolant(X_corr(:), Y_corr(:), surfCrop(:), ...
                                 'linear', 'none');
        obj.surfData.curr = F(obj.worldCoords.mesh.x, obj.worldCoords.mesh.y);
    end

%%
    function storeData(obj)
        if obj.outputConfig.savePhase
            if ~obj.loopState.discardCurr
                obj.phaseData.stack(:,:, obj.loopState.surfDataInd) = obj.phaseData.curr;
            else
                obj.phaseData.stack(:,:,obj.loopState.surfDataInd) = nan;
            end
        end

        if obj.outputConfig.saveSurf
            if ~obj.loopState.discardCurr
                obj.surfData.stack(:,:,obj.loopState.surfDataInd) = obj.surfData.curr;
            else
                obj.surfData.stack(:,:, obj.loopState.surfDataInd) = nan;
            end
        end          

        if obj.outputConfig.saveImages
            obj.imgState.stack(:,:,obj.loopState.surfDataInd) = obj.imgState.curr;
        end
    end

%%
    function calibrateTakeda(obj, heightVec, opts)
    %CALIBRATETAKEDA Fit L and d in Takeda's phase-to-height formula.
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
    %   Stores obj.surfParams.L and obj.surfParams.d and sets
    %   obj.surfParams.profModel = 'takeda' - the exact form evaluated by
    %   calculateSurface, whose carrier term 2*pi/period/mmPerPixel is
    %   resizeFactor-invariant, so the fit remains valid for solves at any
    %   resizeFactor.
    
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
        f0 = 1 / (obj.surfParams.period * obj.worldCoords.scaling.mmPerPixel);
    
        % Assemble fit data
        [rows, cols, nPlanes] = size(phaseCrop);
        phaseVec = reshape(double(phaseCrop), rows*cols, nPlanes);
    
        % Two global parameters gain nothing from millions of samples;
        % stride the pixels to a manageable count. All planes are kept.
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
    
        % Fall back to geometry-scale guesses when the linear estimate is
        % degenerate or violates the bounds (midpoint of finite bounds,
        % otherwise a multiple of the height range).
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
    
        % Bounded nonlinear least-squares fit
        % Coefficient order [L, d] follows the argument order of the handle.
        % f0 is captured from the workspace. Weights are squared so the
        % objective sum((w*r)^2) matches calibratePoly's weighting semantics
        % (fit() itself minimizes sum(w*r^2))
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
        obj.surfParams.L = L;
        obj.surfParams.d = d;
        obj.surfParams.profModel = 'takeda';
    end
%%
    function calibratePoly(obj, polyOrder, heightVec, opts)
    %CALIBRATEPOLY Fit a per-pixel polynomial phase-to-height model.
    %
    %   For every pixel in the calibration region, fits
    %       h = a0 + a1*dPhi + a2*dPhi^2 + ... + a_polyOrder*dPhi^polyOrder
    %   (a0 fitted only when constOffset = true, zero otherwise) by weighted
    %   least squares against the calibration stack in phaseData.stack.
    %
    %   Stores in obj.surfParams.profPolynomial:
    %       coeffMatOrg  - pixelwise coefficients, RAW frame, zero-padded
    %       fittedPlanes - poly22 fits of each coefficient over the region
    %       pixelwise    - 0 (initParams selects smooth planes by default)
     
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
     
        % Embed pixelwise coefficients into the RAW frame (zero-padded)
        coeffMatOrg = zeros(geom.rawSize(1), geom.rawSize(2), size(coeffMat, 3));
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
        obj.surfParams.profPolynomial.coeffMatOrg  = coeffMatOrg;
        obj.surfParams.profPolynomial.fittedPlanes = fittedPlanes;
        obj.surfParams.profPolynomial.pixelwise = false;
        obj.surfParams.profModel = 'poly';
    end
%%
    function setProfModelTakeda(obj, L, d)
        %SETPROFMODELTAKEDA Select Takeda's phase-to-height model.
        %   In:  L - camera height above the reference plane (mm)
        %        d - camera-projector separation (mm)

        obj.surfParams.profModel = 'takeda';
        obj.surfParams.L = L;
        obj.surfParams.d = d;
    end

%%
    function setProfModelPoly(obj, addr)
        %SETPROFMODELPOLY Select a polynomial phase-to-height model.
        % calibration from a .mat file.
        %   The path is remembered so initImages can reload it for a saved case.
        %   In:  addr - .mat file written by writePolyCalibration

        obj.surfParams.profModel = 'poly';
        obj.surfParams.profPolynomial = load(addr);
        obj.surfParams.profPolynomialAddr = addr;
    end

%%
    function writePolyCalibration(obj, addr)
        %WRITEPOLYCALIBRATION Save the polynomial calibration to a .mat file.
        %   The resized working copy (coeffMat) is dropped; only the RAW-frame
        %   coefficients and the fitted planes are written.
        %   In:  addr - output .mat file path

        temp = obj.surfParams.profPolynomial;
        if isfield(temp, 'coeffMat')
            temp = rmfield(temp, 'coeffMat');
        end
        save(addr, '-struct', 'temp')
    end

%%
    function setPolyPixelwise(obj, state)
        %SETPOLYPIXELWISE Choose pixelwise or smoothed polynomial coefficients.
        %   In:  state - true for the raw pixelwise fit, false for the poly22
        %                smoothed coefficient planes

        if isfield(obj.surfParams, 'profPolynomial')
            obj.surfParams.profPolynomial.pixelwise = logical(state);
        else
            warning("Polynomial model not found. Run readPolyCalibration() first.")
        end

    end

%%
    function setROI(obj)
        % Define the computational domain interactively.
        % Shows the raw reference image and waits for a rectangle to be
        % drawn. The result becomes cropRectOrg, the display rectangle is
        % set one unwrap margin inside it, and the geometry is refreshed.

        if isempty(obj.imgState.refScaled)
            error("Missing reference image.")
        end
        mg = obj.prcOpts.unwrapMarginOrg;

        figure()
        tempHandle = imshow(obj.imgState.refRaw, []);

        while true
            rect = drawrectangle(tempHandle.Parent);
            rect = round(rect);
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

        if isvalid(tempHandle)
            close(tempHandle.Parent.Parent)
        end

        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry();
    end

%%
    function showROI(obj)
        % Draw the domain rectangles on the reference image.
        % Overlays the computational domain (red), the unwrap margin (blue)
        % and the display domain (green) for visual checking.

        figure;
        image = obj.imgState.refRaw;
        image = rescale(image);

        if ~isempty(obj.prcOpts.cropRectOrg)
            image = insertShape(image, 'rectangle', obj.prcOpts.cropRectOrg, 'LineWidth', 5, 'ShapeColor', 'red');
            image = insertText(image, obj.prcOpts.cropRectOrg(1:2), 'Computational domain', 'FontColor', 'red', 'FontSize', 20);
            % draw unwrap margin
            mg = obj.prcOpts.unwrapMarginOrg;
            unwrapRect = obj.prcOpts.cropRectOrg + [mg, mg, -2*mg, -2*mg];
            image = insertShape(image, 'rectangle', unwrapRect, 'LineWidth', 5, 'ShapeColor', 'blue');
            image = insertText(image, unwrapRect(1:2) + [unwrapRect(3) 0], 'Unwrap margin', 'FontColor', 'blue', 'FontSize', 20, 'AnchorPoint', 'RightTop');

        end

        if ~isempty(obj.prcOpts.cropRectDisplayOrg)
            absoluteCropDisplay = obj.prcOpts.cropRectDisplayOrg;
            absoluteCropDisplay(1:2) = absoluteCropDisplay(1:2) + obj.prcOpts.cropRectOrg(1:2);
            image = insertShape(image, 'rectangle', absoluteCropDisplay, 'LineWidth', 2, 'ShapeColor', 'green');
            image = insertText(image, [absoluteCropDisplay(1), ...
                                       absoluteCropDisplay(2) ...
                                       + absoluteCropDisplay(4)], ...
                                       'Display domain', 'FontColor', ...
                                       'green', 'FontSize', 20, ...
                                       'AnchorPoint', 'LeftBottom');
        end
        imshow(image, [])
    end
    
%%
    function tseries = probeElev(obj, location, shape, operation)
        %PROBEELEV Elevation time series at a world location.
        %   Finds the grid point nearest to the requested position; with a
        %   size given, aggregates over a square or circular patch instead.
        %   In:  location  - [x y] or [x y size] in world units (mm)
        %        shape     - 'square' or 'circle' (only used with a size)
        %        operation - 'min', 'mean' or 'max' over the patch
        %   Out: tseries   - elevation versus frame

        if ~isequal(size(location), [1 2]) && ~isequal(size(location), [1 3])
            error("Probe target must be either a 1x2 or 1x3 vector.")
        end

        if isempty(obj.surfData.stack)
            error("No surface elevation data.")
        end

        pitch = hypot(obj.worldCoords.mesh.x(1,2) - obj.worldCoords.mesh.x(1,1), ...
            obj.worldCoords.mesh.y(1,2) - obj.worldCoords.mesh.y(1,1));

        [~, ind] = min((obj.worldCoords.mesh.x - location(1)).^2 + ...
            (obj.worldCoords.mesh.y - location(2)).^2, [], 'all');
        [y_ind, x_ind] = ind2sub(size(obj.worldCoords.mesh.x), ind);

        n = size(location, 2);
        if n == 2
            tseries = squeeze(obj.surfData.stack(y_ind, x_ind, :));
            return
        end

        r = location(3)/2;
        r_pixel = max(1, round(r / pitch));

        switch shape
            case 'square'
                % square centered at x = location(1), y = location(2) with side length of location(3)
                probeData = obj.surfData.stack(y_ind - r_pixel:y_ind + r_pixel, x_ind - r_pixel:x_ind + r_pixel, :);
            case 'circle'
                angles = linspace(0, 2*pi, 10000);
                x = cos(angles) * r_pixel + x_ind;
                y = sin(angles) * r_pixel + y_ind;
                mask = poly2mask(x, y, size(obj.worldCoords.mesh.x, 1), size(obj.worldCoords.mesh.x, 2));
                mask = double(mask);
                mask(mask == 0) = nan;
                probeData = obj.surfData.stack .* mask;
            otherwise
                error("Invalid probe shape. Valid shapes are 'square' and 'circle'.")
        end

        switch operation
                case 'min'
                    processedData = min(probeData, [], [1 2]);
                case 'mean'
                    processedData = mean(probeData, [1 2], "omitnan");
                case 'max'
                    processedData = max(probeData, [], [1 2]);
                otherwise
                    error("Invalid operation. Valid operations are 'min', 'mean' and 'max'")
        end

        tseries = squeeze(processedData);
    end

%%
    function tseries = probePhase(obj, location, shape, operation)
        %PROBEPHASE Phase time series at a world location.
        %   Same addressing as probeElev, applied to phaseData.stack, so the
        %   case must have been solved with savePhase on.
        %   In:  location  - [x y] or [x y size] in world units (mm)
        %        shape     - 'square' or 'circle' (only used with a size)
        %        operation - 'min', 'mean' or 'max' over the patch
        %   Out: tseries   - phase versus frame

        if ~isequal(size(location), [1 2]) && ~isequal(size(location), [1 3])
            error("Probe target must be either a 1x2 or 1x3 vector.")
        end

        if ~obj.outputConfig.savePhase
            error("No phase data. Switch savePhase property on before running solve().")
        end

        pitch = hypot(obj.worldCoords.mesh.x(1,2) - obj.worldCoords.mesh.x(1,1), ...
            obj.worldCoords.mesh.y(1,2) - obj.worldCoords.mesh.y(1,1));

        [~, ind] = min((obj.worldCoords.mesh.x - location(1)).^2 + ...
            (obj.worldCoords.mesh.y - location(2)).^2, [], 'all');
        [y_ind, x_ind] = ind2sub(size(obj.worldCoords.mesh.x), ind);
        
        n = length(location);
        if n == 2
            tseries = squeeze(obj.phaseData.stack(y_ind, x_ind, :));
            return
        end

        % assumes x and y scaling are equal
        r = location(3)/2;
        r_pixel = max(1, round(r / pitch));

        switch shape
            case 'square'
                % square centered at x = location(1), y = location(2) with side length of location(3)
                probeData = obj.phaseData.stack(y_ind - r_pixel:y_ind + r_pixel, x_ind - r_pixel:x_ind + r_pixel, :);
            case 'circle'
                angles = linspace(0, 2*pi, 10000);
                x = cos(angles) * r_pixel + x_ind;
                y = sin(angles) * r_pixel + y_ind;
                mask = poly2mask(x, y, size(obj.worldCoords.mesh.x, 1), size(obj.worldCoords.mesh.x, 2));
                mask = double(mask);
                mask(mask == 0) = nan;
                probeData = obj.phaseData.stack .* mask;
            otherwise
                error("Invalid probe shape. Valid shapes are 'square' and 'circle'.")
        end

        switch operation
                case 'min'
                    processedData = min(probeData, [], [1 2]);
                case 'mean'
                    processedData = mean(probeData, [1 2], "omitnan");
                case 'max'
                    processedData = max(probeData, [], [1 2]);
                otherwise
                    error("Invalid operation. Valid operations are 'min', 'mean' and 'max'")
        end

        tseries = squeeze(processedData);
    end

%%
    function animate(obj, startTime, endTime, framerate, opts)
        %ANIMATE Animate the reconstructed surface or phase.
        %   Plots the surface elevation or phase (or an array passed in opts.data) frame
        %   by frame, with optional plane/mean/constant subtraction, NaN-frame
        %   interpolation, smoothing and a time or frame annotation.
        %   In:  startTime, endTime - frame range ([] = to the end)
        %        framerate          - playback rate (Hz)
        %        opts               - name-value display options (data, ZLim,
        %                             view, colormap, camFPS, subtractPlane,
        %                             subtractMean, smoothing, ...)

        arguments
            obj
            startTime = [];
            endTime = [];
            framerate (1,1) double = 30;
            opts.target string = "surf";
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
            opts.anotTStart (1,1) double = 0;
            opts.nanInterp (1,1) logical = false;
            opts.cleanup (1,1) logical = false;
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

            if isempty(startTime)
                startTime = 1;
            end

            if isempty(endTime)
                endTime = size(dataArray, 3);
            end
            dataArray = dataArray(:,:,startTime:endTime);
            if length(opts.data) > 1
                x_crop = opts.data{2};
                y_crop = opts.data{3};
            else
                [x_crop, y_crop] = meshgrid(1:size(dataArray, 2), 1:size(dataArray, 1));
            end
        else
            if isempty(startTime)
                startTime = obj.inputData.solveRange(1);
            end

            if isempty(endTime)
                endTime = obj.inputData.solveRange(1) ...
                            + size(obj.surfData.stack, 3) - 1;
            end
            if strcmpi(opts.target, "surf")
                [dataArray, x_crop, y_crop] = returnDisplaySubarray(obj, startTime, endTime);
            elseif strcmpi(opts.target, "phase")
                [dataArray, x_crop, y_crop] = returnPhaseDisplaySubarray(obj, startTime, endTime);
            else
                error("Valid inputs to opts.target are 'surf' and 'phase'")
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
            zMin = prctile(dataArray(:,:,1:4:end), 0.01, 'all');
            zMax = prctile(dataArray(:,:,1:4:end), 99.99, 'all');
        else
            zMin = opts.ZLim(1);
            zMax = opts.ZLim(2);
        end

        midHeight = (zMin + zMax)/2;
        amp = (zMax - zMin)/2;    

        figH = figure('Units', 'normalized', 'WindowStyle', 'normal', ...
                        'Position', opts.figPosition);
        axHndl = axes();

        sl = surfl(axHndl, x_crop, y_crop, dataArray(:,:,1), 'light', 'EdgeColor', 'none');
        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % sl(1).AmbientStrength = 0.4;
        sl(2).Style = 'infinite';
        sl(2).Position = opts.lightPosition;

        xlim([min(x_crop, [], 'all'), max(x_crop, [], 'all')])
        ylim([min(y_crop, [], 'all'), max(y_crop, [], 'all')])
        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        c_hndl = colorbar;
        c_hndl.Label.String = 'Elevation (mm)';
        clim(axHndl, [0.9*zMin 0.9*zMax])

        xlabel(axHndl, 'x (mm)')
        ylabel(axHndl, 'y (mm)')
        zlabel(axHndl, 'z (mm)')

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(figH, 'Position', opts.figPosition)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        c_hndl.Position = [0.917,0.11,0.02,0.815];

        if opts.dispTime
            h_anot = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.tFormatString) ~= 0
                tFormatString = opts.tFormatString;
            else
                if opts.camFPS == 1
                    tFormatString = "Frame = %05g";
                else
                    tFormatString = "t = %.3f s";
                end
            end
            timeString = sprintf(tFormatString, opts.anotTStart);
            set(h_anot, 'String', timeString)
        end

        set(gca, 'DataAspectRatio', opts.DataAspectRatio);

        pause(pauseTime)
    
        for i = 2:size(dataArray, 3)
            if ~ishghandle(figH)
                break
            end
            sl(1).ZData = dataArray(:,:,i);

            if opts.dispTime
                if opts.camFPS == 1                    
                    timeString = sprintf(tFormatString, opts.anotTStart + (i - 1)*opts.nominalSkip);
                else
                    timeString = sprintf(tFormatString, opts.anotTStart + (i - 1)*opts.nominalSkip/opts.camFPS);
                end
                set(h_anot, 'String', timeString)
            end

            pause(pauseTime)
        end
    end

%%
    function animateImages(obj, startTime, endTime, framerate)
        arguments
            obj
            startTime (1,1) double = obj.inputData.solveRange(1)
            endTime (1,1) double = obj.inputData.solveRange(2)
            framerate (1,1) double = 5
        end
        figH = figure(120);
        
        pauseTime = 1/framerate;

        obj.loopState.surfLoopInd = startTime;
        obj.loadNextImage();
        obj.preprocessCurrentImage();

        [values, edges] = histcounts(obj.imgState.curr, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));
        minThresh = edges(find(values > 0.001, 1, 'first'));
        
        imHndl = imshow(obj.imgState.curr, [minThresh maxThresh]);
        title(num2str(startTime))                
        pause(pauseTime)

        for i = startTime + 1:endTime
            if ~ishghandle(figH)
                break
            end

            obj.loopState.surfLoopInd = i;
            obj.loadNextImage();
            obj.preprocessCurrentImage();

            set(imHndl, 'CData', obj.imgState.curr)
            title(num2str(i))                
            pause(pauseTime)
        end
    end

%%
    function animateRawImages(obj, startTime, endTime, framerate)
        arguments
            obj
            startTime (1,1) double = obj.inputData.solveRange(1)
            endTime (1,1) double = obj.inputData.solveRange(2)
            framerate (1,1) double = 5
        end
        figH = figure(120);
        
        pauseTime = 1/framerate;

        obj.loopState.surfLoopInd = startTime;
        obj.loadNextImage();

        [values, edges] = histcounts(obj.imgState.curr, 'Normalization','cdf');
        maxThresh = 1.1*edges(find(values > 0.999, 1, 'first'));

        for i = startTime:endTime
            if ~ishghandle(figH)
                break
            end

            obj.loopState.surfLoopInd = i;
            obj.loadNextImage();
    
            if isfield(obj.cameraCalib, 'GSx')
                obj.imgState.curr = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, double(obj.imgState.curr), obj.cameraCalib.GSx, obj.cameraCalib.GSy, 'linear');
                obj.imgState.curr(isnan(obj.imgState.curr)) = 0;
            elseif isfield(obj.cameraCalib, 'intrinsicsMatlab')
                obj.imgState.curr = undistortImage(obj.imgState.curr, obj.cameraCalib.intrinsicsMatlab);
                obj.imgState.curr = imwarp(obj.imgState.curr, obj.cameraCalib.pTransform);
            end
   
            if i == startTime
                imHndl = imshow(obj.imgState.curr, [0 maxThresh]);
            else
                set(imHndl, 'CData', obj.imgState.curr)
            end
            title(num2str(i))                
            pause(pauseTime)
        end
    end

%%
    function imageCell = returnRawImages(obj, startTime, endTime)
        if ~exist('startTime', 'var')
            startTime = obj.inputData.solveRange(1);
            endTime = obj.inputData.solveRange(2);
        end

        imageCell = cell(1, endTime - startTime + 1);

        for i = startTime:endTime
            obj.loopState.surfLoopInd = i;
            obj.loadNextImage();
            imageCell{i - startTime + 1} = obj.imgState.curr;
        end            
    end
%%
    function imageCell = returnImages(obj, startTime, endTime)
        if ~exist('startTime', 'var')
            startTime = obj.inputData.solveRange(1);
            endTime = obj.inputData.solveRange(2);
        end

        imageCell = cell(1, endTime - startTime + 1);

        for i = startTime:endTime
            obj.loopState.surfLoopInd = i;
            obj.loadNextImage();
            obj.preprocessCurrentImage();
            imageCell{i - startTime + 1} = obj.imgState.curr;
        end            
    end
%%
    function setCropRect(obj, cropRect)
        obj.prcOpts.cropRectOrg = cropRect;
        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry(); 
     end
%%
    function exportVideo(obj, startTime, endTime, framerate, opts)
        arguments
            obj
            startTime = [];
            endTime = [];
            framerate (1,1) double = 30;
            opts.target string = "surf";
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
            opts.anotTStart (1,1) double = 0;
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
            if isempty(startTime)
                startTime = 1;
            end
            if isempty(endTime)
                endTime = size(dataArray, 3);
            end
            dataArray = dataArray(:,:,startTime:endTime);
            if length(opts.data) > 1
                x_crop = opts.data{2};
                y_crop = opts.data{3};
            else
                [x_crop, y_crop] = meshgrid(1:size(dataArray, 2), 1:size(dataArray, 1));
            end
        else
            if isempty(startTime)
                startTime = obj.inputData.solveRange(1);
            end
            if isempty(endTime)
                endTime = obj.inputData.solveRange(1) ...
                            + size(obj.surfData.stack, 3) - 1;
            end

            if strcmpi(opts.target, "surf")
                [dataArray, x_crop, y_crop] = returnDisplaySubarray(obj, startTime, endTime);
            elseif strcmpi(opts.target, "phase")
                [dataArray, x_crop, y_crop] = returnPhaseDisplaySubarray(obj, startTime, endTime);
            else
                error("Valid inputs to opts.target are 'surf' and 'phase'")
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

        figH = figure('Units', 'normalized', 'WindowStyle', 'normal', ...
                      'Position', opts.figPosition);
        axHndl = axes;
        
        if ~opts.fast
            mkdir(fullfile(opts.saveAddr, "tempImages_" + obj.caseID));
        end

        if isequal(opts.ZLim, [0 0])
            zMin = prctile(dataArray(:,:,1:4:end), 0.01, 'all');
            zMax = prctile(dataArray(:,:,1:4:end), 99.99, 'all');
        else
            zMin = opts.ZLim(1);
            zMax = opts.ZLim(2);
        end

        midHeight = (zMin + zMax)/2;
        amp = (zMax - zMin)/2;

        sl = surfl(axHndl, x_crop, y_crop, dataArray(:,:,1), 'light', ...
                    'EdgeColor', 'none');

        material dull
        if ~isempty(opts.colormap)
            colormap(opts.colormap)
        end

        % sl(1).AmbientStrength = 0.4;
        sl(2).Style = 'infinite';
        sl(2).Position = opts.lightPosition;

        xlim([min(x_crop, [], 'all'), max(x_crop, [], 'all')])
        ylim([min(y_crop, [], 'all'), max(y_crop, [], 'all')])

        zlim(axHndl, [midHeight - 1.5*amp, midHeight + 2*amp])

        c_hndl = colorbar;
        c_hndl.Label.String = 'Elevation (mm)';
        clim(axHndl, [0.9*zMin 0.9*zMax])

        xlabel(axHndl, 'x (mm)')
        ylabel(axHndl, 'y (mm)')
        zlabel(axHndl, 'z (mm)')

        set(axHndl, 'fontsize', opts.FontSize)
        set(axHndl, 'view' , opts.view)
        set(axHndl, 'Position', [0.09,0.10,0.775,0.815])
        c_hndl.Position = [0.917,0.11,0.02,0.815];

        if opts.dispTime
            h_anot = annotation('textbox', [0.09, 0.8, 0.15, 0.0550], 'String', '', 'FontSize', opts.FontSize + 2, 'EdgeColor', 'none');
            if strlength(opts.tFormatString) ~= 0
                tFormatString = opts.tFormatString;
            else
                if opts.camFPS == 1
                    tFormatString = "Frame = %05g";
                else
                    tFormatString = "t = %.3f s";
                end
            end
            timeString = sprintf(tFormatString, opts.anotTStart);
            set(h_anot, 'String', timeString)
        end

        set(gca, 'DataAspectRatio', opts.DataAspectRatio);

        if opts.fast
            frame = getframe(figH);
            writeVideo(v, frame)
        else
            exportgraphics(figH, ...
                fullfile(opts.saveAddr, "tempImages_" + obj.caseID, ...
                           sprintf("tempImage_%04g.png", 1)), ...
                           'Resolution', opts.imageResolution)
        end

            
        if size(dataArray, 3) < 2
            close(v);
            return
        end

        for i = 2:size(dataArray, 3)
            if ~ishghandle(figH)
                break
            end
            sl(1).ZData = dataArray(:,:,i);
            
            if opts.dispTime
                if opts.camFPS == 1                    
                    timeString = sprintf(tFormatString, opts.anotTStart + (i - 1)*opts.nominalSkip);
                else
                    timeString = sprintf(tFormatString, opts.anotTStart + (i - 1)*opts.nominalSkip/opts.camFPS);
                end
                set(h_anot, 'String', timeString)
            end

            if opts.fast
                frame = getframe(figH);
                writeVideo(v, frame)
            else
                exportgraphics(figH, ...
                    fullfile(opts.saveAddr, "tempImages_" + obj.caseID, ...
                               sprintf("tempImage_%04g.png", i)), ...
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
%%
    function makeVideo(obj, framerate)            
        addr = uigetdir();
        if addr == 0
            return
        end

        filenames = dir(fullfile(addr, "*.png"));

        videoFilename = fullfile(addr, (obj.caseID + ".mp4"));
        videoSuffix = 2;
        while isfile(videoFilename)
            videoFilename = sprintf('%s_%g.mp4', fullfile(addr, obj.caseID), videoSuffix);
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
%%
    function loadNextImage(obj)
        if strcmp(obj.inputData.dataFiletype, 'im7') % read from Davis im7 images
            addr = fullfile(obj.inputData.dataAddr(obj.loopState.surfLoopInd).folder, ...
                                        obj.inputData.dataAddr(obj.loopState.surfLoopInd).name);
            obj.imgState.curr = obj.getDavisFrame(addr, obj.inputData.imgFrameNum);
        elseif  strcmp(obj.inputData.dataFiletype, 'set')  % read from Davis set file
            obj.imgState.curr = obj.getDavisFrame(obj.inputData.dataAddr, obj.inputData.imgFrameNum, obj.loopState.surfLoopInd);
        else
            obj.imgState.curr = imread(   fullfile(obj.inputData.dataAddr(obj.loopState.surfLoopInd).folder, ...
                                    obj.inputData.dataAddr(obj.loopState.surfLoopInd).name) ...
                                );
        end

        obj.imgState.curr = double(obj.imgState.curr);
     end 
%%
    function initPhaseCorrArr(obj)
        if ~isempty(obj.imgState.ref)
            obj.pCorrOpts.arr = false(size(obj.imgState.refRaw));
        end
        
        cropRectOrg = obj.prcOpts.cropRectOrg;
        if strcmpi(obj.prcOpts.patNormAxis, 'Y')
            targetInds = round(linspace(cropRectOrg(1), cropRectOrg(1) + cropRectOrg(3), 5));
            obj.pCorrOpts.lineInds = unique(targetInds(2:end - 1))';
            obj.pCorrOpts.arr(:, obj.pCorrOpts.lineInds) = 1;
        else                                                    % normal vector to pattern is parallel to X axis
            targetInds = round(linspace(cropRectOrg(2), cropRectOrg(2) + cropRectOrg(4), 5));
            obj.pCorrOpts.lineInds = unique(targetInds(2:end - 1))';
            obj.pCorrOpts.arr(obj.pCorrOpts.lineInds, :) = 1;
        end
    end

%%
    function setPhaseCorrTemporal(obj)
        obj.pCorrOpts.method = 'temporal';
    end
%%
    function setPhaseCorrSpatial(obj)
        obj.pCorrOpts.method = 'spatial';
    end
%%
    function setDisplayToMargin(obj)
        mrg = obj.prcOpts.unwrapMarginOrg;
        obj.prcOpts.cropRectDisplayOrg = [mrg + 1, mrg + 1, ... 
                                    obj.prcOpts.cropRectOrg(3) - 2*mrg, ...
                                    obj.prcOpts.cropRectOrg(4) - 2*mrg ...
                                    ];
        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry();
    end
%%
    function setDisplayFromMargin(obj, offset)
        mrg = obj.prcOpts.unwrapMarginOrg;
        obj.prcOpts.cropRectDisplayOrg = [mrg + 1 + offset, mrg + 1 + offset, ... 
                                obj.prcOpts.cropRectOrg(3) - 2*mrg - 2*offset, ...
                                obj.prcOpts.cropRectOrg(4) - 2*mrg - 2*offset];
        obj.pCorrOpts.peakInd = [];
        obj.updateGeometry();
    end
%%
    function setResizeFactor(obj, factor) 
        obj.prcOpts.resizeFactor = factor;
        obj.updateGeometry();
    end

%%
    function [rowInds, colInds] = outputInds(obj)
        % Row and column indices of the OUTPUT frame
        r = round(obj.prcOpts.cropRectDisplay);
        [Ny, Nx] = size(obj.worldCoords.mesh.x);
        rowInds = r(2) : min(Ny, r(2) + r(4));
        colInds = r(1) : min(Nx, r(1) + r(3));
    end
%%
    function [subarray, x_crop, y_crop] = returnDisplaySubarray(obj, startTime, endTime)
        solveRange = obj.inputData.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = startTime - solveRange(1) + 1:endTime - solveRange(1) + 1; 

        if ~isempty(obj.surfData.stack)
            subarray = obj.surfData.stack(rowInds, colInds, tRange);
        else
            subarray = 0;
        end
        x_crop = obj.worldCoords.mesh.x(rowInds, colInds);
        y_crop = obj.worldCoords.mesh.y(rowInds, colInds);
    end
%%
    function [subarray, x_crop, y_crop] = returnFilteredDisplaySubarray(obj, startTime, endTime)
        solveRange = obj.inputData.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = startTime - solveRange(1) + 1:endTime - solveRange(1) + 1; 

        if ~isempty(obj.surfData.stack)
            subarray = obj.surfData.stack(rowInds, colInds, tRange);
        else
            subarray = 0;
        end
        x_crop = obj.worldCoords.mesh.x(rowInds, colInds);
        y_crop = obj.worldCoords.mesh.y(rowInds, colInds);
       
        for i = 1:size(subarray, 3)
            % fit plane to data and subtract
            temp = subarray(:,:,i);
            sf = fit([x_crop(:), y_crop(:)], temp(:), 'poly11');
            fittedPlane = sf.p00 + x_crop*sf.p10 + y_crop*sf.p01;
            subarray(:,:,i) = subarray(:,:,i) - fittedPlane;
        end

        if size(subarray, 3) > 99
            subarray = subarray - mean(subarray, 'all', 'omitmissing');
        end
    end

%%
    function [subarray, x_crop, y_crop] = returnPhaseDisplaySubarray(obj, startTime, endTime)
        solveRange = obj.inputData.solveRange;
        [rowInds, colInds] = obj.outputInds();
        tRange = startTime - solveRange(1) + 1:endTime - solveRange(1) + 1; 

        if ~isempty(obj.phaseData.stack)
            subarray = obj.phaseData.stack(rowInds, colInds, tRange);
        else
            subarray = 0;
        end 
        x_crop = obj.worldCoords.mesh.x(rowInds, colInds);
        y_crop = obj.worldCoords.mesh.y(rowInds, colInds);
    end

%%
    function importData(obj, data)
        obj.inputData.importedDataCell{end + 1} = data;
    end

%%
    function reportProgress(obj)
        if ~obj.outputConfig.reportProgressEnabled
            return
        end
        textstring = sprintf('Processing image %5d...', obj.loopState.surfLoopInd);
        if obj.loopState.surfLoopInd > obj.inputData.solveRange(1)
            textstringOld = sprintf('Processing image %5d...', obj.loopState.surfLoopInd - 1);
            fprintf(repmat('\b', 1, numel(textstringOld)))
            fprintf(textstring)
        else
            % first timestep
            fprintf(textstring)
        end
    end

%%
    function userFunction(obj)
        if obj.prcOpts.userFncEnabled
            % store or do something ...
            % obj.postData.postArray(obj.loopState.surfDataInd,1) = ... ;
        end
    end

%%
    function clbModeOn(obj)
        obj.prcOpts.clbMode = true;
        obj.prcOpts.resizeFactor = 1;
        obj.prcOpts.resizeFactorDisplay = 1;
        obj.outputConfig.savePhase = true;
        obj.outputConfig.saveSurf = false;
    end

%%
    function clbModeOff(obj)
        obj.prcOpts.clbMode = false;
        obj.outputConfig.savePhase = false;
        obj.outputConfig.saveSurf = true;
    end

%%
    function readInterpPixelMask(obj, mask)
        obj.prcOpts.interpPixelMask = mask;
        [obj.prcOpts.interpPixelRows, ...
            obj.prcOpts.interpPixelCols] = find(mask);
    end
%%
    function drawPeaks(obj, startTime, endTime, skip)
        %DRAWPEAKS Diagnostic plot of phase-correction peak detection.
        %   Shows the reference line with its detected peaks and the tracked
        %   peak highlighted, then the same line in the requested images.
        if ~exist('skip', 'var')
            skip = 1;
        end

        lineNo = min(2, length(obj.pCorrOpts.lineInds));
        pcInd = obj.pCorrOpts.lineInds(lineNo);

        % reference peaks from the shared detection/selection machinery
        [optInd, ~, peaksLineCell] = obj.findOptimalPeakInd(obj.pCorrOpts.peakInd);

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

        peaksRef = peaksLineCell{lineNo};
        if isempty(peaksRef) || peakInd > length(peaksRef)
            error("Peak %g was not detected on phase-correction line %g of the reference image.", ...
                peakInd, pcInd)
        end

        if strcmp(obj.prcOpts.patNormAxis, 'X')
            refLine = double(obj.imgState.refRaw(pcInd, :));
        else
            refLine = double(obj.imgState.refRaw(:, pcInd));
        end
        peakValsRef = refLine(peaksRef);

        % detect the same line in the requested images
        targetImages = startTime:skip:endTime;
        N = length(targetImages);

        peaksCurrCell = cell(1, N);
        corrLineCell = cell(1, N);

        for n = 1:N
            obj.loopState.surfLoopInd = targetImages(n);
            obj.loadNextImage();
            if isfield(obj.cameraCalib, 'GSx')
                obj.imgState.curr = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, double(obj.imgState.curr), obj.cameraCalib.GSx, obj.cameraCalib.GSy, 'linear');
                obj.imgState.curr(isnan(obj.imgState.curr)) = 0;
            elseif isfield(obj.cameraCalib, 'intrinsicsMatlab')
                obj.imgState.curr = undistortImage(obj.imgState.curr, obj.cameraCalib.intrinsicsMatlab);
                obj.imgState.curr = imwarp(obj.imgState.curr, obj.cameraCalib.pTransform);
            end

            if strcmp(obj.prcOpts.patNormAxis, 'X')
                phaseCorrLine = double(obj.imgState.curr(pcInd, :));
            else
                phaseCorrLine = double(obj.imgState.curr(:, pcInd));
            end
            [peakVals, peaksCurr] = findpeaks(phaseCorrLine, ...
                'MinPeakDistance', 0.6*obj.surfParams.periodOrg, ...
                'MinPeakProminence', obj.pCorrOpts.peakProm, ...
                'MinPeakHeight', obj.pCorrOpts.peakMin(lineNo));

            if strcmp(obj.pCorrOpts.startEdge, 'right') || strcmp(obj.pCorrOpts.startEdge, 'bottom')
                peaksCurr = flip(peaksCurr);
                peakVals = flip(peakVals);
            end

            peaksCurrCell{n} = [peaksCurr(:), peakVals(:)];
            corrLineCell{n} = phaseCorrLine;
        end

        % plot
        figure;
        tiledlayout(3, ceil(N/3) + 1)
        nexttile([3 1])
        plot(refLine)
        hold on
        scatter(peaksRef, 1.01*peakValsRef, 'v', 'filled')
        scatter(peaksRef(peakInd), 1.01*peakValsRef(peakInd), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
        hold off
        if strcmp(obj.pCorrOpts.startEdge, 'left') || strcmp(obj.pCorrOpts.startEdge, 'top')
            xLimit = [1 round(1.5*peaksRef(peakInd))];
        else
            xLimit = [round(0.9*peaksRef(peakInd)) length(refLine)];
        end
        yLimit = [0.8*min(refLine(xLimit(1):xLimit(2))) 1.2*max(refLine(xLimit(1):xLimit(2)))];

        xlim(xLimit)
        ylim(yLimit)
        title("Reference: x_p = " + num2str(peaksRef(peakInd)) )
        for i = 1:N
            nexttile
            plot(corrLineCell{i})
            hold on
            scatter(peaksCurrCell{i}(:,1), 1.05*peaksCurrCell{i}(:,2), 'v', 'filled')
            if size(peaksCurrCell{i}, 1) >= peakInd
                scatter(peaksCurrCell{i}(peakInd, 1), 1.05*peaksCurrCell{i}(peakInd,2), 'v', 'filled', 'MarkerFaceColor', [0.47 0.9 0.19])
                title(num2str(targetImages(i)) + ": x_p = " + num2str(peaksCurrCell{i}(peakInd, 1)) )
            else
                title(num2str(targetImages(i)) + ": peak " + num2str(peakInd) + " not detected")
            end
            hold off
            xlim(xLimit)
            ylim(yLimit)
        end
    end

%%
    function readCameraCalibration(obj, calibAddr, camNumber)
        [~, ~, ext] = fileparts(calibAddr);

        if strcmpi(ext, '.mat')
            obj.cameraCalib = load(calibAddr);
            obj.cameraCalib.type = "Pinhole";

            imSize = obj.cameraCalib.intrinsics.ImageSize;
            [~, pTransform, scaling, xMesh, yMesh] = rectifyImagePinhole(zeros(imSize), obj.cameraCalib, 0, obj.prcOpts.imgRotAngle);
            obj.cameraCalib.pTransform = pTransform;
            obj.cameraCalib.scaling = scaling;

            obj.worldCoords.mesh.xOrg = xMesh;
            obj.worldCoords.mesh.yOrg = yMesh;

            obj.storePinholeScaling(scaling);

        elseif strcmpi(ext, '.xml')
            calibCell = readDavisCalibration(calibAddr);
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
                
                mmPerPixel_x = polyCalib.scaling.X.Slope;
                mmPerPixel_y = polyCalib.scaling.Y.Slope;
                offSetX = polyCalib.scaling.X.Offset;
                offSetY = polyCalib.scaling.Y.Offset;

                obj.worldCoords.scaling.mmPerPixelOrg = abs(mmPerPixel_x);
                obj.worldCoords.scaling.X.Slope = mmPerPixel_x;
                obj.worldCoords.scaling.X.SlopeOrg = mmPerPixel_x;
                obj.worldCoords.scaling.X.Unit = 'mm';
                obj.worldCoords.scaling.Y.Slope = mmPerPixel_y;
                obj.worldCoords.scaling.Y.SlopeOrg = mmPerPixel_y;
                obj.worldCoords.scaling.Y.Unit = 'mm';
                obj.worldCoords.scaling.X.Offset = offSetX;
                obj.worldCoords.scaling.Y.Offset = offSetY;
            end
            obj.cameraCalib.type = calibType;
        end

        obj.cameraCalib.addr = calibAddr;
        obj.cameraCalib.camNum = camNumber;
    end

%%
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
%%
    function writeSurfData(obj, addr)
        if ~exist('addr', 'var')
            addr = uigetdir();
        end                
        
        if obj.outputConfig.saveSurf
            fID = fopen(fullfile(addr, obj.caseID + "_surfData.bin"), 'w');
        end
        if obj.outputConfig.savePhase
            fID_phase = fopen(fullfile(addr, obj.caseID + "_phaseData.bin"), 'w');
        end

        if ~obj.prcOpts.useSegmentation
            if obj.outputConfig.savePhase
                [displayArr, ~, ~] = obj.returnPhaseDisplaySubarray(obj.inputData.solveRange(1), obj.inputData.solveRange(2));
                fwrite(fID_phase, numel(size(displayArr)), 'uint32');
                fwrite(fID_phase, size(displayArr), 'uint32');
                fwrite(fID_phase, displayArr, 'single');
                fclose(fID_phase);
            end

            if obj.outputConfig.saveSurf
                displayArr = obj.returnDisplaySubarray(obj.inputData.solveRange(1), obj.inputData.solveRange(2));
    
                fwrite(fID, numel(size(displayArr)), 'uint32');
                fwrite(fID, size(displayArr), 'uint32');
                fwrite(fID, displayArr, 'single');
                fclose(fID);
            end
        else
            tempFolderAddr = fullfile(obj.outputConfig.tempWriteAddr, "TEMP_" + obj.caseID);
            if obj.outputConfig.saveSurf
                % read temporary files and write contents to binary one by one
                blockFilenames = dir(fullfile(tempFolderAddr, "*surfData.bin"));
                % read first block and array size
                filename = fullfile(blockFilenames(1).folder, blockFilenames(1).name);
                [blockArr, arrSize] = obj.readData(filename);
    
                if numel(arrSize) > 2
                    arrSize(end) = diff(obj.inputData.solveRange) + 1;
                end
    
                fwrite(fID, numel(arrSize), 'uint32');
                fwrite(fID, arrSize, 'uint32');
                fwrite(fID, blockArr, 'single');
    
                for i = 2:length(blockFilenames)
                    filename = fullfile(blockFilenames(i).folder, blockFilenames(i).name);
                    blockArr = obj.readData(filename);
                    fwrite(fID, blockArr, 'single');
                end
                fclose(fID);
            end

            if obj.outputConfig.savePhase
                % read temporary phase files and write contents to binary one by one
                blockFilenames = dir(fullfile(tempFolderAddr, "*phaseData.bin"));
                % read first block and array size
                filename = fullfile(blockFilenames(1).folder, blockFilenames(1).name);
                [blockArr, arrSize] = obj.readData(filename);

                if numel(arrSize) > 2
                    arrSize(end) = diff(obj.inputData.solveRange) + 1;
                end

                fwrite(fID_phase, numel(arrSize), 'uint32');
                fwrite(fID_phase, arrSize, 'uint32');
                fwrite(fID_phase, blockArr, 'single');

                for i = 2:length(blockFilenames)
                    filename = fullfile(blockFilenames(i).folder, blockFilenames(i).name);
                    blockArr = obj.readData(filename);
                    fwrite(fID_phase, blockArr, 'single');
                end
                fclose(fID_phase);
            end
        end

        [rowInds, colInds] = obj.outputInds();
        xMesh = obj.worldCoords.mesh.x(rowInds,colInds);
        yMesh = obj.worldCoords.mesh.y(rowInds,colInds);

        temp.xMesh = xMesh;
        temp.yMesh = yMesh;
        save(fullfile(addr, "surfMesh.mat"), '-struct', 'temp')

        if obj.prcOpts.useSegmentation
            rmdir(tempFolderAddr, 's')
        end
    end

%%
    function writeInParts(obj, blocksize, writeAddr)
        obj.prcOpts.blockSize = blocksize;
        obj.outputConfig.tempWriteAddr = writeAddr;
        obj.prcOpts.useSegmentation = true;
    end
%%
    function writeTempBlock(obj)
        if ~obj.prcOpts.useSegmentation
            return
        end
        tempFolderAddr = fullfile(obj.outputConfig.tempWriteAddr, "TEMP_" + obj.caseID);
        if ~isfolder(tempFolderAddr)
            mkdir(tempFolderAddr);
        end

        solveRange = obj.inputData.solveRange;
        [Ny, Nx] = size(obj.worldCoords.mesh.x);
        Nt = diff(solveRange) + 1;
        totalBlocks = ceil(Nt / obj.prcOpts.blockSize);

        if obj.loopState.blockNumber == totalBlocks - 1
            Nt_sub = Nt - obj.loopState.blockNumber*obj.prcOpts.blockSize;
        else
            Nt_sub = obj.prcOpts.blockSize;
        end

        if obj.outputConfig.saveSurf
            filename = sprintf("tempFile_%03g_surfData.bin", obj.loopState.blockNumber);
            fID = fopen(fullfile(tempFolderAddr, filename), 'w');
    
            tempArr = obj.returnDisplaySubarray(solveRange(1), ...
                            solveRange(1) + size(obj.surfData.stack, 3) - 1);
    
            fwrite(fID, numel(size(tempArr)), 'uint32');
            fwrite(fID, size(tempArr), 'uint32');
            fwrite(fID, tempArr, 'single');
            fclose(fID);
        end

        if obj.outputConfig.savePhase
            filenamePhase = sprintf("tempFile_%03g_phaseData.bin", obj.loopState.blockNumber);
            fID = fopen(fullfile(tempFolderAddr, filenamePhase), 'w');
            tempArr = obj.returnPhaseDisplaySubarray(solveRange(1), ...
                        solveRange(1) + size(obj.phaseData.stack, 3) - 1);
            
            fwrite(fID, numel(size(tempArr)), 'uint32');
            fwrite(fID, size(tempArr), 'uint32');
            fwrite(fID, tempArr, 'single');
            fclose(fID);
        end

        if obj.loopState.blockNumber == totalBlocks
            return
        end

        if obj.outputConfig.savePhase
            obj.phaseData.stack = single(nan(Ny, Nx, Nt_sub));
        end

        if obj.outputConfig.saveSurf
            obj.surfData.stack = single(nan(Ny, Nx, Nt_sub));
        end
        
    end
%%
    function writeCase(obj, addr)
      if ~exist('addr', 'var')
        addr = uigetdir();
      end

      if ~isfolder(fullfile(addr, obj.caseID))
        mkdir(fullfile(addr, obj.caseID))
      end

      temp.(obj.caseID) = obj;
      save(fullfile(addr, obj.caseID, obj.caseID + ".mat"), '-struct', 'temp')
      obj.writeSurfData(fullfile(addr, obj.caseID))
    end

%%
    function imgOut = interpBadPixels(obj, img)
        if strcmpi(obj.prcOpts.interpImgMethod, 'fastLinear')
            imgOut = obj.fastFill(img, obj.prcOpts.interpPixelRows, ...
                                    obj.prcOpts.interpPixelCols);
        elseif strcmpi(obj.prcOpts.interpImgMethod, 'smooth')
            imgOut = regionfill(img, obj.prcOpts.interpPixelMask);
        end
    end

%%
    function updateGeometry(obj)
        %UPDATEGEOMETRY Recompute all resizeFactor-derived quantities
        rf = obj.prcOpts.resizeFactor;
        rf_o = obj.prcOpts.resizeFactorDisplay;

        obj.surfParams.period = obj.surfParams.periodOrg * rf;

        if ~isempty(obj.worldCoords.scaling)
            obj.worldCoords.scaling.X.Slope = obj.worldCoords.scaling.X.SlopeOrg / rf;
            obj.worldCoords.scaling.Y.Slope = obj.worldCoords.scaling.Y.SlopeOrg / rf;
            obj.worldCoords.scaling.mmPerPixel = obj.worldCoords.scaling.mmPerPixelOrg / rf;
        end

        obj.prcOpts.unwrapMargin = round(obj.prcOpts.unwrapMarginOrg * rf * rf_o);

        if isempty(obj.prcOpts.cropRectOrg)
            [Ny, Nx] = size(obj.imgState.refRaw);
            obj.prcOpts.cropRectOrg = [1, 1, Nx - 1, Ny - 1];
        end

        if ~isempty(obj.imgState.refRaw)
            obj.imgState.refScaled = imresize(obj.imgState.refRaw, rf);
        end

        cropRectOrg = round(obj.prcOpts.cropRectOrg);
        if isfield(obj.imgState, 'refRaw') && ~isempty(obj.imgState.refRaw)
            [NyRaw, NxRaw] = size(obj.imgState.refRaw);
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

        [NyScaled, NxScaled] = size(obj.imgState.refScaled);
        cropRect = round(cropRectOrg * rf);
        cropRect(1:2) = max(1, cropRect(1:2));
        cropRect(3) = min(NxScaled - cropRect(1), cropRect(3));
        cropRect(4) = min(NyScaled - cropRect(2), cropRect(4));
        obj.prcOpts.cropRect = cropRect;

        % the dimensions of the resized and cropped reference image is
        % required to initialize arrays
        % obj.imgState.ref is reinitialized in preprocessRefImage()
        obj.imgState.ref = imcrop(obj.imgState.refScaled, obj.prcOpts.cropRect);
        stackSize = ceil(rf_o*size(obj.imgState.ref));
        
        if isempty(obj.prcOpts.cropRectDisplayOrg)
            obj.prcOpts.cropRectDisplayOrg = [1 1 cropRectOrg(3) cropRectOrg(4)];
        end
        cropRectDisplay = round(obj.prcOpts.cropRectDisplayOrg * rf * rf_o);

        cropRectDisplay(1:2) = max(1, cropRectDisplay(1:2));
        cropRectDisplay(3) = min(stackSize(2) - cropRectDisplay(1), cropRectDisplay(3));
        cropRectDisplay(4) = min(stackSize(1) - cropRectDisplay(2), cropRectDisplay(4));
        obj.prcOpts.cropRectDisplay = cropRectDisplay;

        obj.initPhaseCorrArr();

        % initialize coordinate mesh
        [ny, nx, ~] = size(obj.imgState.refScaled);
        [Ny, Nx] = size(obj.imgState.refRaw);

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
%%
    function cLims = plotCalibrationDiagnostics(obj, phaseCrop, elevCrop, ...
                fitPhase, fitHeight, heightVec, calibRect, figTitle, cLims)
        %PLOTCALIBRATIONDIAGNOSTICS One diagnostics figure for a coefficient set.
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
            cHndl = colorbar;
            cHndl.Label.String = 'Error (mm)';
            if applyLims
                clim(cLims{k});
            else
                cLims{k} = cHndl.Limits;
            end
            xlabel('X (mm)')
            ylabel('Y (mm)')
            title(sprintf('H = %.2f mm', heightVec(showErrInd(k))))
        end
    end
%%
    function output = RAW2STACK(obj, input)
        %RAW2STACK Convert (row, col) values from RAW to STACK frame
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
end

methods (Access = private)
    function scaleAndTransformRef(obj)
        if ~isempty(obj.prcOpts.interpPixelMask)
            obj.imgState.refRaw = obj.interpBadPixels(obj.imgState.refRaw);
        end

        if isfield(obj.cameraCalib, 'GSx')
            obj.imgState.refRaw = interp2(obj.cameraCalib.Gx, obj.cameraCalib.Gy, double(obj.imgState.refRaw), obj.cameraCalib.GSx, obj.cameraCalib.GSy, 'linear');
            obj.imgState.refRaw(isnan(obj.imgState.refRaw)) = 0;
        elseif isfield(obj.cameraCalib, 'intrinsicsMatlab')
            obj.imgState.refRaw = undistortImage(obj.imgState.refRaw, obj.cameraCalib.intrinsicsMatlab);
            obj.imgState.refRaw = imwarp(obj.imgState.refRaw, obj.cameraCalib.pTransform);
        end

        obj.updateGeometry();
    end
%%
    function [phaseCrop, heightVec, weights, geom] = prepareCalibrationData( ...
            obj, heightVec, excludeHeights, weights)
        %PREPARECALIBRATIONDATA Shared validation, trimming, and geometry.
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
    
        [Ny, Nx] = size(obj.imgState.refRaw);
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
%%
    function savedObj = saveobj(obj)
        savedObj.caseID = obj.caseID;
        savedObj.inputData = obj.inputData;
        savedObj.imgState = obj.defaultImgState();         
        savedObj.phaseData = obj.defaultPhaseData();    
        savedObj.surfData = obj.defaultSurfData();
        savedObj.pCorrOpts = obj.pCorrOpts;           
        savedObj.demodOpts = obj.demodOpts;
        savedObj.prcOpts = obj.prcOpts;
        savedObj.loopState = obj.loopState;
        savedObj.cameraCalib = obj.cameraCalib;         
        savedObj.worldCoords = obj.worldCoords;
        savedObj.worldCoords.mesh = [];
        savedObj.surfParams = obj.surfParams;
        savedObj.outputConfig = obj.outputConfig;
        savedObj.postData = obj.postData;
        if isfield(savedObj.surfParams, 'profPolynomial')
            savedObj.surfParams = rmfield(savedObj.surfParams, 'profPolynomial');
        end
    end
end

methods (Access = private, Static)
    function output = defaultInput()
        output.dataAddr = "";
        output.refAddr = "";
        output.refFiletype = "";
        output.dataFiletype = "";
        output.fullRange = [];
        output.solveRange = [];
        output.imgFrameNum = 1;
        output.refImgFrameNum = 1;
        output.importedDataCell = {};
    end

    function output = defaultImgState()
        output.refRaw = [];
        output.refScaled = [];
        output.ref = [];
        output.curr = [];
    end

    function output = defaultPhaseData()
        output.curr = [];
        output.old = [];
        output.stack = [];
        output.refComplexCoeffs = [];
    end

    function output = defaultPCorrOpts()
        output.method = "spatial";
        output.lineInds = [];
        output.arr = [];
        output.startEdge = "left";
        output.edgeSafetyFactor = 2.5;
        output.peakInd = [];
        output.manualVals = [];
        output.peakProm = 60;
        output.peakMin = [];
        output.warningFrames = [];
    end

    function output = defaultDemodOpts()
        output.method = "FT";
        output.scaleList = [];
        output.angleList = 0;
        output.sigma = 0.6;
        output.gamma = 1;
        output.filtWidthFrac = 0.6;
    end

    function output = defaultSurfData()
        output.stack = [];
        output.curr = [];
    end

    function output = defaultPrcOpts()
        output.unwrapEnabled = true;
        output.unwrapMethod = "1D";
        output.unwrapMarginOrg = 20;
        output.unwrapMargin = [];
        output.patNormVec = [];
        output.patNormAxis = "";
        output.imgRotAngle = 0; % in degrees and positive clockwise
        output.lateralShiftCorrection = false;
        output.useSegmentation = false;
        output.discardThreshold = [];
        output.blockSize = 1000;
        output.userFncEnabled = false;
        output.clbMode = false;
        output.dcSubtraction = true;
        output.equalizeExposure = false;
        output.smoothRefImg = false;
        output.interpImgMethod = 'fastLinear';
        output.interpPixelMask = [];
        output.interpPixelRows = [];
        output.interpPixelCols = [];
        output.cropRectOrg = [];
        output.cropRect = [];
        output.cropRectDisplayOrg = [];
        output.cropRectDisplay = [];
        output.resizeFactor = 1;
        output.resizeFactorDisplay = 0.5;
    end

    function output = defaultLoopState()
        output.surfDataInd = 0;
        output.surfLoopInd = [];
        output.firstTimestep = false;
        output.blockNumber = 1;
        output.discardCurr = false;
        output.prcTime = [];
    end

    function output = defaultCameraCalib()
        output.type = "none";
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

    function output = defaultOutputConfig()
        output.savePhase = false;
        output.saveImages = false;
        output.saveSurf = true;
        output.reportProgressEnabled = true;
        output.tempWriteAddr = "";
    end
end
    
methods(Static)
    function newObj = loadobj(fileObj)
        if isstruct(fileObj)
            newObj = FtpSolver("temp");
            newObj.caseID = fileObj.caseID;
            newObj.inputData = fileObj.inputData;
            newObj.imgState = fileObj.imgState;         
            newObj.phaseData = fileObj.phaseData;    
            newObj.surfData = fileObj.surfData;
            newObj.pCorrOpts = fileObj.pCorrOpts;           
            newObj.demodOpts = fileObj.demodOpts;
            newObj.prcOpts = fileObj.prcOpts;
            newObj.loopState = fileObj.loopState;
            newObj.cameraCalib = fileObj.cameraCalib;         
            newObj.worldCoords = fileObj.worldCoords;             
            newObj.surfParams = fileObj.surfParams;
            newObj.outputConfig = fileObj.outputConfig;
            newObj.postData = fileObj.postData;
        else
            newObj = fileObj;
        end
    end
    %%
    function [data, dataSizeAll] = readData(addr, seekTime)
            if ~exist('addr', 'var')
                [file, path] = uigetfile('*.bin');
                addr = fullfile(path, file);
            end
            fID = fopen(addr, 'rb');

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
                % skip first three numbers which are matrix dimensions
                fseek(fID, skipLength, 'cof');
                data = fread(fID, readElements, '*single');
            else
                data = fread(fID, '*single');
            end

            data = reshape(data, dataSize);            
            fclose(fID);
    end
%%
    function data = subtractPlane(data)
        downscale = 0.2;
        [Ny, Nx, Nt] = size(data);
        [xGrid, yGrid] = meshgrid(1:Nx, 1:Ny);
        xGrid_coarse = imresize(xGrid, downscale, 'bilinear', Antialiasing=false);
        yGrid_coarse = imresize(yGrid, downscale, 'bilinear', Antialiasing=false);
    
        for i = 1:Nt
            % fit plane to data and subtract
            temp = double(imresize(data(:,:,i), downscale));
            sf = fit([xGrid_coarse(:), yGrid_coarse(:)], temp(:), 'poly11');
            fittedPlane = sf.p00 + xGrid*sf.p10 + yGrid*sf.p01;
            data(:,:,i) = data(:,:,i) - fittedPlane;
        end
    end
%%
    function img = subtractDC(img, period)
            %SUBTRACTDC Remove the low-frequency (DC) component from an image.
            %   Applies a double box filter of width ~4 fringe periods (rounded
            %   up to an odd size) and subtracts the result. Shared by
            %   preprocessRefImage and preprocessCurrentImage.
            filterSize = 4*period;
            filterSize = filterSize + 1 - mod(filterSize, 2);
            dc = imboxfilt(img, filterSize, 'padding', 'symmetric');
            dc = imboxfilt(dc, filterSize, 'padding', 'symmetric');
            img = double(img) - double(dc);
    end
%%
    function [period, patNormAxis, tiltAngle, normVec] = analyzeFringe(img)
    %ANALYZEFRINGE Estimate fringe period and orientation via FFT
    %
    %   [period, patNormAxis, tiltAngle, kCarrier] = analyzeFringe(img)
    %
    %   Applies a Hann window to the image and locates the dominant
    %   spectral peak outside the low-frequency core.
    %
    %   period      - fringe period in pixels
    %   patNormAxis - 'X' if the pattern normal is closer to the x axis,
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
            patNormAxis = 'X';       % fx > 0 by construction of the half-plane
            tiltAngle = atand(fy/fx);
        else
            patNormAxis = 'Y';
            if fy < 0                % report the conjugate with fy > 0
                fx = -fx;
                fy = -fy;
            end
            tiltAngle = atand(fx/fy);
        end
        normVec = [fx, fy]/norm([fx, fy]);
    end
%%
    function delta = quadPeakOffset(fm, f0, fp)
        %QUADPEAKOFFSET Sub-bin offset of a spectral peak from its neighbors
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
%%
    function writeBinary(data, filename, addr)
        filename = string(filename);
        if ~exist('addr', 'var')
            addr = uigetdir();
        end                
        fID = fopen(fullfile(addr, filename + "_surfData.bin"), 'w');

        fwrite(fID, numel(size(data)), 'uint32');
        fwrite(fID, size(data), 'uint32');
        fwrite(fID, data, 'single');
        fclose(fID);
    end
%%
    function imgOut = fastFill(img, interpRows, interpCols)
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

%%
    function W4D = morletCWT(img, scales, anglesDeg, sigma, epsilon, fringeAxis)
        %MORLETCWT Continuous Morlet wavelet transform of a 1-D or 2-D input.
        %   W4D = MORLETCWT(IMG, SCALES, ANGLESDEG, SIGMA, EPSILON) returns the
        %   complex CWT coefficients of IMG as an array of size
        %   [size(IMG,1), size(IMG,2), numel(SCALES), nAngles], L1-normalized.
        %
        %   2-D input: Wavelet Toolbox CWTFT2 with the anisotropic Morlet
        %   wavelet {Omega0 = 6, SIGMA, EPSILON}
        %
        %   1-D input (row or column vector): self-contained Fourier-domain
        %   analytic Morlet (k0 = 6, tunable SIGMA); ANGLESDEG and EPSILON are
        %   ignored and nAngles = 1.
        if ~exist('fringeAxis', 'var')
            fringeAxis = '';
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
                % Analytic Morlet in the frequency domain, L2-normalized
                psiHat = exp(-sigma^2 * (a*k - k0).^2 / 2);
                W4D(1, :, iS) = ifft(fSig .* psiHat);
            end
            W4D = W4D(1, padSize+1:end-padSize, :);

            if isCol
                W4D = permute(W4D, [2 1 3 4]);
            end
        else
            % 2-D image ----------------
            if strcmpi(fringeAxis, 'X')
                paddedImg = padarray(img, [0 padSize], 0);
                paddedImg = padarray(paddedImg, [padSize 0], 'symmetric');
            elseif strcmpi(fringeAxis, 'Y')
                paddedImg = padarray(img, [padSize 0], 0);
                paddedImg = padarray(paddedImg, [0 padSize], 'symmetric');
            else
                paddedImg = padarray(img, [padSize padSize], 0);
            end
            cwtOut = cwtft2(paddedImg, ...
                wavelet = {"morlet", {6, sigma, epsilon}}, ...
                scales  = scales, ...
                angles  = deg2rad(anglesDeg), ...
                norm    = "L1");

            % cfs is [rows x cols x 1 x nScales x nAngles]; drop the
            % singleton plane dimension and crop the padding.
            W4D = permute(cwtOut.cfs, [1 2 4 5 3]);
            W4D = W4D(padSize+1:end-padSize, padSize+1:end-padSize, :, :);
        end
    end
%%
    function elev = evalPolyElev(phase, coeffs)
        %EVALPOLYELEV h = sum_i coeffs(:,:,i) .* phase.^(i-1)
        %   coeffs planes are 2-D (or scalar) and broadcast over the time
        %   dimension of phase.
        elev = zeros(size(phase));
        for i = 1:size(coeffs, 3)
            elev = elev + phase.^(i - 1) .* coeffs(:,:,i);
        end
    end
%%
    function elev = evalTakedaElev(phase, L, d, f0)
        %EVALTAKEDAELEV h = L*dPhi/ (2*pi*f0*d + dPhi)
        elev = L*phase ./ (2*pi*f0*d + phase);
    end
%%
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
%%
    function out = rescaleCoord(x, scale)
        %RESCALECOORD Map a coordinate between frames related by imresize.
        % scale about the image edge
        out = (x - 0.5)*scale + 0.5;
    end
end
end
