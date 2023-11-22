
function figH = showAndSaveThisImg( outDir, img, name, varargin )

  p = inputParser;
  p.addOptional( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.addParameter( 'range', [], @isnumeric );
  p.addParameter( 'showScale', 1, @ispositive );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  verbose = p.Results.verbose;
  range = p.Results.range;
  showScale = p.Results.showScale;
  wavSplit = p.Results.wavSplit;

  if verbose == true
    figH = figure; imshowscale( img, showScale );
  end

  if numel( wavSplit ) > 0
    outImg = wavScale( img, wavSplit );
  else
    outImg = scaleImg( img, [ 0 1 ], range );
  end
  imwrite( imresize( outImg, showScale ), [ outDir, filesep(), name, '.png' ] );

  thresh = 0.05;
  lowScalingLevel = findValueBelowFraction( outImg(:), 1-thresh );
  highScalingLevel = findValueBelowFraction( outImg(:), thresh );
  scaling = [ lowScalingLevel highScalingLevel ];
  outImgNice = scaleImg( outImg, [ 0 1 ], scaling );
  imwrite( outImgNice, [ outDir, filesep(), name, '_nice.png' ] );  
end
