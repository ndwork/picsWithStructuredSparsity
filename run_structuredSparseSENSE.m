
function run_structuredSparseSENSE
  close all;  rng(1);

  sampleFraction = 0.16;
  lambda = 0.04;
  noiseCoords = [ 1 1 10 10 ];

  load( './data/brain_8ch.mat', 'DATA' );
  data = DATA;  clear DATA;
  data = squeeze( data ) / max( abs( data(:) ) );  % assumes one slice.
  data = padData( data, [ 224 224 8] );  % zero pad data for wavelet transform

  sImg = size( data, [1 2] );
  nCoils = size( data, 3 );

  coilRecons = mri_reconIFFT( data );
  trueRecon = mri_reconRoemer( coilRecons );
  absTrueRecon = abs( trueRecon );

  noiseVecs = coilRecons( noiseCoords(1):noiseCoords(3), noiseCoords(2):noiseCoords(4), : );
  noiseVecs = reshape( noiseVecs, [], nCoils );
  noiseCov = cov( noiseVecs );
  noiseCov = noiseCov ./ max( abs( noiseCov(:) ) );

  nSamples = round( sampleFraction * prod( sImg ) );

  lapMask = mri_makeSampleMask( sImg, nSamples, round( 0.3*sImg ) );
  fftSamples_lap = bsxfun( @times, data, lapMask );

  wavSplit = makeWavSplit( sImg, 'minSplitSize', 24 );
  wavACR   = makeWavAutoCalRegion( sImg, wavSplit );
  wavMaskACR = mri_makeSampleMask( sImg, nSamples, 'maskType', 'VDPD', 'startMask', wavACR>0 );
  fftSamples_wavACR = bsxfun( @times, data, wavMaskACR );

  lap_sMaps = mri_makeSensitivityMaps( fftSamples_lap, 'alg', 'sortaSimple' );
  reconSparseSense = mri_reconSparseSENSE( fftSamples_lap, lap_sMaps, lambda, ...
    'noiseCov', noiseCov, 'optAlg', 'fista_wLS', 'transformType', 'wavelet', 'wavSplit', wavSplit );
  pccSparseSense = dotP( abs( reconSparseSense ), absTrueRecon ) / ...
    norm( abs( reconSparseSense(:) ) ) / norm( absTrueRecon(:) );

  wavACR_sMaps = mri_makeSensitivityMaps( fftSamples_wavACR, 'alg', 'sortaSimple' );

  reconSparseSense_wav = mri_reconSparseSENSE( fftSamples_wavACR, wavACR_sMaps, lambda, ...
    'noiseCov', noiseCov, 'optAlg', 'fista_wLS', 'transformType', 'wavelet', 'wavSplit', wavSplit );
  pccSparseSense_wav = dotP( abs( reconSparseSense_wav ), absTrueRecon ) / ...
    norm( abs( reconSparseSense_wav(:) ) ) / norm( absTrueRecon(:) );

  reconSSP_wav = mri_reconStructuredSparseSENSE( fftSamples_wavACR, wavACR_sMaps, lambda, ...
    'noiseCov', noiseCov, 'optAlg', 'fista_wLS', 'transformType', 'wavelet', 'wavSplit', wavSplit );
  pccSSP_wav = dotP( abs( reconSSP_wav ), absTrueRecon ) / ...
    norm( abs( reconSSP_wav(:) ) ) / norm( absTrueRecon(:) );

  figure;  imshowscale( absTrueRecon, 3 );  titlenice( 'True Recon' );  drawnow;
  figure;
  showImageCube( abs( cat( 3, reconSparseSense, reconSparseSense_wav, reconSSP_wav ) ), 3 );
  title( 'Structured Sparsity / Sparse SENSE' );

  disp( 'PCC of sparseSense, sparseSense_wav structuredSparseSense' );
  disp( [ pccSparseSense, pccSparseSense_wav, pccSSP_wav ] );
end

