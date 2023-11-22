
function sensitivities = simulateSliceSensitivities( img, varargin )

  p = inputParser;
  p.addParameter( 'H', 0.3, @ispositive );
  p.addParameter( 'type', 'axial', @(x) true );
  p.parse( varargin{:} );
  H = p.Results.H;
  type = p.Results.type;

  sImg = size( img );

  coils = makeBirdcageCoils( );
  nCoils = numel( coils );

  sensitivities = zeros( [ sImg nCoils ] );

  coords = size2fftCoordinates( sImg );

  if strcmp( type, 'axial' )
    W = sImg(2) / sImg(1) * H;  % width of slice of interest
    [ xs, ys ] = meshgrid( 0.5*H * coords{2}, 0.5*W * coords{1} );
    locs = [ xs(:)  ys(:)  zeros(numel(ys),1) ];
  elseif strcmp( type, 'sagittal' )
    D = sImg(2) / sImg(1) * H;  % width of slice of interest
    [ ys, zs ] = meshgrid( 0.5*H * coords{2}, 0.5*D * coords{1} );
    locs = [ zeros(numel(ys),1)  ys(:)  zs(:)  ];
  end


  for coilIndx = 1 : nCoils
    thisCoil = coils{ coilIndx }';

    tmp = mri_computeSensitivityBiotSavart( thisCoil, locs );
    sensitivities( :, :, coilIndx ) = reshape( tmp, sImg );
  end

  sensitivities = reshape( sensitivities, [ prod(sImg) nCoils ] );
  [u,s,v] = svd( sensitivities, 'econ' );
  for coilIndx = nCoils-3 : nCoils
    s( coilIndx, coilIndx ) = 0;
  end
  sensitivities = reshape( u * s * v', [ sImg nCoils ] );
  sensitivities = max( sensitivities, 0 );

  sensitivities = bsxfun( @rdivide, sensitivities, sqrt( sum( abs( sensitivities ).^2, 3 ) ) );
end

