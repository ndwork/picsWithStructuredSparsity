
function printMetrics( logFile, sampleFraction, nSamples, recon, imgTitle, outDir, varargin )

  p = inputParser;
  p.addOptional( 'trueRecon', [], @isnumeric );
  p.addParameter( 'saveImages', true );
  p.addParameter( 'senseMaps', [], @isnumeric );
  p.addParameter( 'showScale', 3, @ispositive );
  p.addParameter( 'optLogNames', [] );
  p.addParameter( 'optLogValues', [] );  % must be numeric
  p.parse( varargin{:} );
  saveImages = p.Results.saveImages;
  trueRecon = p.Results.trueRecon;
  senseMaps = p.Results.senseMaps;
  showScale = p.Results.showScale;
  optLogNames = p.Results.optLogNames;
  optLogValues = p.Results.optLogValues;

  if max( abs( recon(:) ) ) == 0
    normAbsRecon = recon;
  else
    normAbsRecon = abs( recon ) / max( abs( recon(:) ) );
  end

  mdmScore = calcMetricMDM( normAbsRecon );
  mdmScore2 = calcMetricMDM( 1-normAbsRecon );
  if max( normAbsRecon(:) ) == 0
    niqeScore = 0;
    piqeScore = 0;
    autofocusValue = 0;
  else
    niqeScore = niqe( normAbsRecon );
    piqeScore = piqe( normAbsRecon );
    autofocusValue = calcAutofocusMetric( normAbsRecon );
  end

  if numel( trueRecon ) > 0
    absRecon = abs( recon );
    absTrueRecon = abs( trueRecon );
    errImg = absRecon - absTrueRecon;
    if saveImages == true
      errFig = figure;  imshowscale( abs(errImg), showScale );
      title( [ imgTitle, ' absErr' ] );
      saveas( errFig, [ outDir, '/err_', imgTitle, '.png' ] ); close( errFig );
    end
    mse = norm( abs( errImg(:) ), 2 ).^2 / numel( recon );  % mean squared error
    mae = sum( abs( errImg(:) ) ) / numel( recon );  % mean absolute error
    relErr = norm( errImg(:) ) / norm( trueRecon(:) );
    ssimValue = ssim( min( absRecon / max( abs( trueRecon(:) ) ), 1 ), ...
                      absTrueRecon / max( abs( trueRecon(:) ) ) );
    ms_ssimValue = multiScaleSSIM( min( absRecon / max( abs( trueRecon(:) ) ), 1 ), ...
                                   absTrueRecon / max( abs( trueRecon(:) ) ) );
    correlation = dotP( absRecon, absTrueRecon ) / ...
      norm( absRecon(:) ) / norm( absTrueRecon(:) );
    angleErr = acos( correlation );
    mutualInfo = mi( absRecon * 255, absTrueRecon * 255 );
  else
    mse = -1;
    relErr = -1;
    mae = -1;
    ssimValue = -1;
    ms_ssimValue = -1;
    correlation = -999;
    angleErr = -1;
    mutualInfo = -1;
  end

  logID = fopen( logFile, 'a' );

  if numel( optLogValues ) > 0
    logValuesString = sprintf( '%d, ' , optLogValues );
    disp([ 'Logging ', imgTitle, ...
      ' / ', num2str(sampleFraction), ' / ', strjoin( optLogNames, ', ' ), ...
      ', ', logValuesString(1:end-1) ]);
    fprintf( logID, ['%f, %d, ', imgTitle, ', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f ', ...
      repmat( ', %d', [1 numel(optLogValues)] ), '\n' ], ...
      sampleFraction, nSamples, mse, mae, relErr, ssimValue, mdmScore, mdmScore2, niqeScore, ...
      piqeScore, ms_ssimValue, correlation, angleErr, mutualInfo, autofocusValue, optLogValues(:)' );
  else
    disp([ 'Logging ', imgTitle, ...
      ' / ', num2str(sampleFraction) ]);
    fprintf( logID, ['%f, %d, ', imgTitle, ', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n'], ...
      sampleFraction, nSamples, mse, mae, relErr, ssimValue, mdmScore, mdmScore2, niqeScore, ...
      piqeScore, ms_ssimValue, correlation, angleErr, mutualInfo, autofocusValue );
  end
  fclose( logID );

  if saveImages == true
    figH = figure; imshowscale( abs(recon), showScale );  title( imgTitle );
    saveas( figH, [ outDir, '/recon_', imgTitle, '.png' ] );  close( figH );

    figH = figure; imshowscale( abs(recon), showScale, 'range', 'nice', 'thresh', 0.02 );  title( imgTitle );
    saveas( figH, [ outDir, '/recon_', imgTitle, '_nice.png' ] );  close( figH );

    if numel( senseMaps ) > 0
      mapsFig = figure;  showImageCube( abs( senseMaps ), showScale );
      saveas( mapsFig, [outDir, '/senseMaps_', imgTitle, '.png'] );  close( mapsFig );
  
      senseRecons = bsxfun( @times, senseMaps, recon );
      senseReconsFig = figure;  showImageCube( abs(senseRecons), showScale );
      saveas( senseReconsFig, [outDir, '/senseRecons_', imgTitle, '.png'] );
      close( senseReconsFig );
  
      save( [ outDir, '/mat_recon_', imgTitle, '.mat' ], 'recon', 'senseMaps' );
    else
  
      save( [ outDir, '/mat_recon_', imgTitle, '.mat' ], 'recon' );
    end
  end

  close all;
end
