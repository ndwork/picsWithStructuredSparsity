
function run_structuredSparseSENSE

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
  noiseVecs = coilRecons( noiseCoords(1):noiseCoords(3), noiseCoords(2):noiseCoords(4), : );
  noiseVecs = reshape( noiseVecs, [], nCoils );
  noiseCov = cov( noiseVecs );
  noiseCov = noiseCov ./ max( abs( noiseCov(:) ) );

  nSamples = round( sampleFraction * prod( sImg ) );
  wavSplit = makeWavSplit( sImg, 'minSplitSize', 24 );
  wavACR   = makeWavAutoCalRegion( sImg, wavSplit );
  wavMaskACR = mri_makeSampleMask( sImg, nSamples, 'maskType', 'VDPD', 'startMask', wavACR>0 );
  fftSamples_wavACR = bsxfun( @times, data, wavMaskACR );

  reconSSP_wav_sMaps = mri_makeSensitivityMaps( fftSamples_wavACR, 'alg', 'sortaSimple' );
  reconSSP_wav = mri_reconStructuredSparseSENSE( fftSamples_wavACR, reconSSP_wav_sMaps, ...
    lambda, 'noiseCov', noiseCov, 'optAlg', 'fista_wLS', ...
    'transformType', 'wavelet', 'wavSplit', wavSplit );

  reconSparseSense_wav_sMaps = mri_makeSensitivityMaps( fftSamples_wavACR, 'alg', 'sortaSimple' );
  reconSparseSense_wav = mri_reconSparseSENSE( fftSamples_wavACR, reconSparseSense_wav_sMaps, lambda, ...
    'noiseCov', noiseCov, 'optAlg', 'fista_wLS', ...
    'transformType', 'wavelet', 'wavSplit', wavSplit );

  figure;
  showImageCube( cat( 3, abs(reconSSP_wav), abs(reconSparseSense_wav) ) );
  title( 'Structured Sparsity / Sparse SENSE' );
end

