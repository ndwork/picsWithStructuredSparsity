
function run_structuredSparsity( varargin )
  close all; rng(1);
  addpath( './dworkLib' );

  showScale = 3;
  %datacases = [ 9 11 7 1 0 2:6 ];
  datacases = [ 4 10 9 11 7 1 0 2:6 ];
  %sampleFractions = 0.08 : 0.02 : 0.4;
  sampleFractions = 0.10 : 0.02 : 0.34;
  lambdas = [ 0.001:0.001:0.01 0.02:0.01:0.1 0.2:0.1:1.0 2:1:10 ];
  logFilename = 'log.csv';
  nReweightIter = 1;
  mainOut = './out_vdpd/';
  [~, nCores] = getNumCores();
  nWorkers = nCores - 5;

  p = inputParser;
  p.addParameter( 'datacases', datacases, @isnonnegative );
  p.addParameter( 'lambdas', lambdas, @isnumeric );
  p.addParameter( 'logFilename', logFilename, @(x) true );
  p.addParameter( 'nWorkers', nWorkers, @ispositive );
  p.addParameter( 'sampleFractions', sampleFractions, @ispositive );
  p.parse( varargin{:} );
  datacases = p.Results.datacases;
  lambdas = p.Results.lambdas ;
  logFilename = p.Results.logFilename;
  nWorkers = p.Results.nWorkers;
  sampleFractions = p.Results.sampleFractions;

  nDatacases = numel( datacases );
  nLambdas = numel( lambdas );


  if ~parpoolExists() && nWorkers > 1, parpool( nWorkers ); end

  for datacaseIndx = 1 : nDatacases
    close all;
    datacase = datacases( datacaseIndx );

    [ data, noiseCoords, ~, trueRecon ] = loadDatacase( datacase );
    data = squeeze( data ) / max( abs( data(:) ) );  % assumes one slice.

    datacaseOut = [ mainOut, filesep, 'datacase_', indx2str( datacase, max( max(datacases), 10 ) ) ];
    logFile = [ datacaseOut, filesep, logFilename ];

    if ~exist( datacaseOut, 'dir' ), mkdir( datacaseOut ); end
    if numel( trueRecon ) > 0
      showAndSaveThisImg( datacaseOut, abs(trueRecon), 'origImg', 'showScale', showScale );
    end

    sData = size( data );
    sImg = sData(1:2);
    nCoils = sData(end);
    noiseCov = [];
    if numel( trueRecon ) == 0
      coilRecons = mri_reconIFFT( data );
      roemerRecon = mri_reconRoemer( coilRecons );
      showAndSaveThisImg( datacaseOut, abs(roemerRecon), 'roemerReconTrue', 'showScale', showScale );
      trueRecon = abs( roemerRecon );

      noiseVecs = coilRecons( noiseCoords(1):noiseCoords(3), noiseCoords(2):noiseCoords(4), : );
      noiseVecs = reshape( noiseVecs, [], nCoils );
      noiseCov = cov( noiseVecs );
      noiseCov = noiseCov ./ max( abs( noiseCov(:) ) );
    end

    logID = fopen( logFile, 'w' );
    fprintf( logID, 'datacase, sampleFraction, nSamples, Algorithm, mse, mae, relErr, ' );
    fprintf( logID, 'ssim, mdmScore, mdmScore2, niqeScore, piqeScore, ms_ssimValue, correlation, ' );
    fprintf( logID, 'angleErr, mutualInfo, autofocusValue\n' );
    fclose( logID );

    wavSplit = makeWavSplit( sImg, 'minSplitSize', 24 );
    [ wavACR, sWavACR ]   = makeWavAutoCalRegion( sImg, wavSplit );
    [ curvACR, sCurvACR ] = makeCurvAutoCalRegion( sImg );

    nNSampleFracs = numel( sampleFractions );
    maxNSamples = max( round( sampleFractions * prod( sImg ) ) );

    sampleMasks     = cell( 3, nNSampleFracs );
    parfor sampleMaskIndx = 1 : nNSampleFracs * 3
      [ maskTypeIndx, nSamplesIndx ] = ind2sub( [ 3 nNSampleFracs ], sampleMaskIndx );

      sampleFraction = sampleFractions( nSamplesIndx );   %#ok<PFBNS>
      nSamples = round( sampleFraction * prod( sImg ) );

      if maskTypeIndx == 1
        sampleMask = mri_makeSampleMask( sImg, nSamples, 'maskType', 'VDPD' );
        sampleMasks{sampleMaskIndx} = sampleMask;
      elseif maskTypeIndx == 2
        wavMaskACR = mri_makeSampleMask( sImg, nSamples, 'maskType', 'VDPD', 'startMask', wavACR>0 );
        sampleMasks{sampleMaskIndx} = wavMaskACR;
      elseif maskTypeIndx == 3
        curvMaskACR = mri_makeSampleMask( sImg, nSamples, 'maskType', 'VDPD', 'startMask', curvACR>0 );
        sampleMasks{sampleMaskIndx} = curvMaskACR;
      end
    end
    sampleMasksWav = sampleMasks(2,:)';
    sampleMasksCurv = sampleMasks(3,:)';
    sampleMasks = sampleMasks(1,:)';

    parfor nSamplesIndx = 1 : nNSampleFracs
      sampleFraction = sampleFractions( nSamplesIndx );
      nSamples = round( sampleFraction * prod( sImg ) );

      sampleMask = sampleMasks{nSamplesIndx};
      wavMaskACR = sampleMasksWav{nSamplesIndx};
      curvMaskACR = sampleMasksCurv{nSamplesIndx};

      outDir = [ datacaseOut, filesep, 'nSamples_', indx2str( nSamples, maxNSamples ) ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end

      fftSamples = bsxfun( @times, data, sampleMask );
      fftSamples_wavACR = bsxfun( @times, data, wavMaskACR );
      fftSamples_curvACR = bsxfun( @times, data, curvMaskACR );

      if sum( abs( wavACR(:) ) > 0 ) > sum( abs( curvACR(:) ) > 0 )   %#ok<PFBNS>
        wavCurvMaskACR = wavMaskACR;
        fftSamples_wavCurvACR = fftSamples_wavACR;
        sWavCurvACR = sWavACR;
      else
        wavCurvMaskACR = curvMaskACR;
        fftSamples_wavCurvACR = fftSamples_curvACR;
        sWavCurvACR = sCurvACR;
      end

      maskDir = [ outDir, filesep, 'sampleMasks' ];
      if ~exist( maskDir, 'dir' ), mkdir( maskDir ); end
      imwrite( sampleMask, [ maskDir, filesep, 'sampleMask.png' ] );
      imwrite( wavMaskACR, [ maskDir, filesep, 'wavMask.png' ] );
      imwrite( curvMaskACR, [ maskDir, filesep, 'curvMask.png' ] );
      imwrite( wavCurvMaskACR, [ outDir, filesep, 'wavCurvMask.png' ] );


      %--- --- --- --- --- --- --- --- ---
      %-- Roemer

      saveImages = false;
      reconRoemer = loadSavedResult( outDir, 'roemer' );
      if numel( reconRoemer ) == 0
        coilRecons = mri_reconIFFT( fftSamples, 'multiSlice', true );
        reconRoemer = mri_reconRoemer( coilRecons );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconRoemer, ...
        'roemer', outDir, trueRecon, 'saveImages', saveImages );

      saveImages = false;
      reconRoemerWav = loadSavedResult( outDir, 'roemerWav' );
      if numel( reconRoemerWav ) == 0
        coilRecons = mri_reconIFFT( fftSamples_wavACR, 'multiSlice', true );
        reconRoemerWav = mri_reconRoemer( coilRecons );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconRoemerWav, ...
        'roemerWav', outDir, trueRecon, 'saveImages', saveImages );

      saveImages = false;
      reconRoemerCurv = loadSavedResult( outDir, 'roemerCurv' );
      if numel( reconRoemerCurv ) == 0
        coilRecons = mri_reconIFFT( fftSamples_curvACR, 'multiSlice', true );
        reconRoemerCurv = mri_reconRoemer( coilRecons );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconRoemerCurv, ...
        'roemerCurv', outDir, trueRecon, 'saveImages', saveImages );

      saveImages = false;
      reconRoemerWavCurv = loadSavedResult( outDir, 'roemerWavCurv' );
      if numel( reconRoemerWavCurv ) == 0
        coilRecons = mri_reconIFFT( fftSamples_wavCurvACR, 'multiSlice', true );
        reconRoemerWavCurv = mri_reconRoemer( coilRecons );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconRoemerWavCurv, ...
        'roemerWavCurv', outDir, trueRecon, 'saveImages', saveImages );

    end
    clear reconRoemerWavCurv reconRoemerWav reconRoemer

    parfor outIndx = 1 : nNSampleFracs * nLambdas
      [ nSamplesIndx, indxLambda ] = ind2sub( [ nNSampleFracs nLambdas ], outIndx );
      lambda = lambdas( indxLambda );   %#ok<PFBNS>

      sampleFraction = sampleFractions( nSamplesIndx );   %#ok<PFBNS>
      nSamples = round( sampleFraction * prod( sImg ) );

      sampleMask = sampleMasks{nSamplesIndx};   %#ok<PFBNS>
      wavMaskACR = sampleMasksWav{nSamplesIndx};   %#ok<PFBNS>
      curvMaskACR = sampleMasksCurv{nSamplesIndx};   %#ok<PFBNS>

      outDir = [ datacaseOut, filesep, filesep, 'nSamples_', indx2str( nSamples, maxNSamples ) ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end

      fftSamples = bsxfun( @times, data, sampleMask );
      fftSamples_wavACR = bsxfun( @times, data, wavMaskACR );
      fftSamples_curvACR = bsxfun( @times, data, curvMaskACR );

      if sum( abs( wavACR(:) ) > 0 ) > sum( abs( curvACR(:) ) > 0 )   %#ok<PFBNS>
        fftSamples_wavCurvACR = fftSamples_wavACR;
        sWavCurvACR = sWavACR;
      else
        fftSamples_wavCurvACR = fftSamples_curvACR;
        sWavCurvACR = sCurvACR;
      end


      %--- --- --- --- --- --- --- --- ---
      %--- Structured Sparse SENSE

      saveImages = false;
      SSp_Dir = [ outDir, '/strucSpars/wavACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSSP_wav,reconSSP_wav_sMaps] = loadSavedResult( SSp_Dir, 'sspWav' );
      if numel( reconSSP_wav ) == 0
        disp( [ 'Working on ', SSp_Dir, '/sspWav' ] );
        reconSSP_wav_sMaps = mri_makeSensitivityMaps( fftSamples_wavACR, 'alg', 'sortaSimple' );
        if ~exist( SSp_Dir, 'dir' ), mkdir( SSp_Dir ); end
        reconSSP_wav = mri_reconStructuredSparseSENSE( fftSamples_wavACR, reconSSP_wav_sMaps, ...
          lambda, 'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'wavelet', 'wavSplit', wavSplit );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSSP_wav, ...
        'sspWav', SSp_Dir, trueRecon, 'senseMaps', reconSSP_wav_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

      saveImages = false;
      SSp_Dir = [ outDir, '/strucSpars/curvACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSSP_curv,reconSSP_curv_sMaps] = loadSavedResult( SSp_Dir, 'sspCurv' );
      if numel( reconSSP_curv ) == 0
        disp( [ 'Working on ', SSp_Dir, '/sspCurv' ] );
        reconSSP_curv_sMaps = mri_makeSensitivityMaps( fftSamples_curvACR, 'alg', 'sortaSimple' );
        if ~exist( SSp_Dir, 'dir' ), mkdir( SSp_Dir ); end
        reconSSP_curv = mri_reconStructuredSparseSENSE( fftSamples_curvACR, reconSSP_curv_sMaps, ...
          lambda, 'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'curvelet' );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSSP_curv, ...
        'sspCurv', SSp_Dir, trueRecon, 'senseMaps', reconSSP_curv_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

      saveImages = false;
      SSp_Dir = [ outDir, '/strucSpars/wavCurvACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSSP_wavCurv,reconSSP_wavCurv_sMaps] = loadSavedResult( SSp_Dir, 'sspWavCurv' );
      if numel( reconSSP_wavCurv ) == 0
        disp( [ 'Working on ', SSp_Dir, '/sspWavCurv' ] );
        reconSSP_wavCurv_sMaps = mri_makeSensitivityMaps( fftSamples_wavCurvACR, 'alg', 'sortaSimple' );
        if ~exist( SSp_Dir, 'dir' ), mkdir( SSp_Dir ); end
        reconSSP_wavCurv = mri_reconStructuredSparseSENSE( fftSamples_wavCurvACR, reconSSP_wavCurv_sMaps, ...
          lambda, 'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'wavCurv', 'wavSplit', wavSplit );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSSP_wavCurv, ...
        'sspWavCurv', SSp_Dir, trueRecon, 'senseMaps', reconSSP_wavCurv_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

    end
    clear reconSSP_wav_sMaps reconSSP_wav
    clear reconSSP_curv_sMaps reconSSP_curv
    clear reconSSP_wavCurv_sMaps reconSSP_wavCurv

    parfor outIndx = 1 : nNSampleFracs * nLambdas
      [ nSamplesIndx, indxLambda ] = ind2sub( [ nNSampleFracs nLambdas ], outIndx );
      lambda = lambdas( indxLambda );   %#ok<PFBNS>

      sampleFraction = sampleFractions( nSamplesIndx );   %#ok<PFBNS>
      nSamples = round( sampleFraction * prod( sImg ) );

      sampleMask = sampleMasks{nSamplesIndx};   %#ok<PFBNS>
      wavMaskACR = sampleMasksWav{nSamplesIndx};   %#ok<PFBNS>
      curvMaskACR = sampleMasksCurv{nSamplesIndx};   %#ok<PFBNS>

      outDir = [ datacaseOut, filesep, filesep, 'nSamples_', indx2str( nSamples, maxNSamples ) ];
      if ~exist( outDir, 'dir' ), mkdir( outDir ); end

      fftSamples = bsxfun( @times, data, sampleMask );
      fftSamples_wavACR = bsxfun( @times, data, wavMaskACR );
      fftSamples_curvACR = bsxfun( @times, data, curvMaskACR );

      if sum( abs( wavACR(:) ) > 0 ) > sum( abs( curvACR(:) ) > 0 )   %#ok<PFBNS>
        fftSamples_wavCurvACR = fftSamples_wavACR;
        sWavCurvACR = sWavACR;
      else
        fftSamples_wavCurvACR = fftSamples_curvACR;
        sWavCurvACR = sCurvACR;
      end


      %--- --- --- --- --- --- --- --- ---
      %--- Sparse SENSE

      saveImages = false;
      SparseSense_Dir = [ outDir, '/sparseSENSE/noACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSparseSense,reconSparseSense_sMaps] = loadSavedResult( SparseSense_Dir, 'sparseSense' );
      if numel( reconSparseSense ) == 0
        reconSparseSense_sMaps = mri_makeSensitivityMaps( fftSamples, 'alg', 'sortaSimple' );
        if ~exist( SparseSense_Dir, 'dir' ), mkdir( SparseSense_Dir ); end
        reconSparseSense = mri_reconSparseSENSE( fftSamples, reconSparseSense_sMaps, lambda, ...
          'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'wavelet' );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSparseSense, ...
        'sparseSense', SparseSense_Dir, trueRecon, 'senseMaps', reconSparseSense_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

      saveImages = false;
      SparseSense_Dir = [ outDir, '/sparseSENSE/wavACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSparseSense_wav,reconSparseSense_wav_sMaps] = loadSavedResult( SparseSense_Dir, 'sparseSenseWav' );
      if numel( reconSparseSense_wav ) == 0
        reconSparseSense_wav_sMaps = mri_makeSensitivityMaps( fftSamples_wavACR, 'alg', 'sortaSimple' );
        if ~exist( SparseSense_Dir, 'dir' ), mkdir( SparseSense_Dir ); end
        reconSparseSense_wav = mri_reconSparseSENSE( fftSamples_wavACR, reconSparseSense_wav_sMaps, lambda, ...
          'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'wavelet', 'wavSplit', wavSplit );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSparseSense_wav, ...
        'sparseSenseWav', SparseSense_Dir, trueRecon, 'senseMaps', reconSparseSense_wav_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

      saveImages = false;
      SparseSense_Dir = [ outDir, '/sparseSENSE/curvACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSparseSense_curv,reconSparseSense_curv_sMaps] = loadSavedResult( SparseSense_Dir, 'sparseSenseCurv' );
      if numel( reconSparseSense_curv ) == 0
        reconSparseSense_curv_sMaps = mri_makeSensitivityMaps( fftSamples_curvACR, 'alg', 'sortaSimple' );
        if ~exist( SparseSense_Dir, 'dir' ), mkdir( SparseSense_Dir ); end
        reconSparseSense_curv = mri_reconSparseSENSE( fftSamples_curvACR, reconSparseSense_curv_sMaps, lambda, ...
          'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'curvelet' );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSparseSense_curv, ...
        'sparseSenseCurv', SparseSense_Dir, trueRecon, 'senseMaps', reconSparseSense_curv_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

      saveImages = false;
      SparseSense_Dir = [ outDir, '/sparseSENSE/wavCurvACR/lam_', num2str( lambda, '%4.3f' ) ];
      [reconSparseSense_wavCurv,reconSparseSense_wavCurv_sMaps] = loadSavedResult( SparseSense_Dir, 'sparseSenseWavCurv' );
      if numel( reconSparseSense_wavCurv ) == 0
        reconSparseSense_wavCurv_sMaps = mri_makeSensitivityMaps( fftSamples_wavCurvACR, 'alg', 'sortaSimple' );
        if ~exist( SparseSense_Dir, 'dir' ), mkdir( SparseSense_Dir ); end
        reconSparseSense_wavCurv = mri_reconSparseSENSE( fftSamples_wavCurvACR, reconSparseSense_wavCurv_sMaps, lambda, ...
          'noiseCov', noiseCov, 'nReweightIter', nReweightIter, 'optAlg', 'fista_wLS', ...
          'transformType', 'wavCurv', 'wavSplit', wavSplit );
        saveImages = true;
      end
      printMetrics( logFile, datacase, sampleFraction, nSamples, reconSparseSense_wavCurv, ...
        'sparseSenseWavCurv', SparseSense_Dir, trueRecon, 'senseMaps', reconSparseSense_wavCurv_sMaps, 'saveImages', saveImages, ...
        'optLogNames', {'lambda'}, 'optLogValues', lambda );

    end
    clear reconSparseSense_wav_sMaps reconSparseSense_wav
    clear reconSparseSense_curv_sMaps reconSparseSense_curv
    clear reconSparseSense_wavCurv_sMaps reconSparseSense_wavCurv

    disp('Going on to next set of samples.')
  end

end

