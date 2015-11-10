function bowtieThickness = loadBowtieFilter( p, geom )
% Load bowtie filter parameters and return resulting bowtie filter
% thickness for each pixel
%
% Meng Wu at Stanford University
% 2012 - 2013


if length( geom.detSize ) == 1
    
    nu = geom.detSize(1);
    du = geom.detSpacing(1);
    ou = geom.detOffset(1);
    uu = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du ;
    
elseif length( geom.detSize ) == 2

    nu = geom.detSize(1);
    nv = geom.detSize(2);
    du = geom.detSpacing(1);
    dv = geom.detSpacing(2);
    ou = geom.detOffset(1);
    ov = geom.detOffset(2);
    
    u = ( ( -(nu-1)/2:(nu-1)/2) + ou ) * du ;
    v = ( ( -(nv-1)/2:(nv-1)/2) + ov ) * dv ;
    uu = meshgrid( u, v);

    
end

if strcmpi(p.Bowtie.shapeType, 'cosine')
    
    angles = p.Bowtie.alpha * atan( uu ./ geom.SDD );
    
    angles( angles < -pi  ) = -pi ;
    angles( angles >  pi  ) =  pi ;
    
    a = p.Bowtie.minimumThickness;
    b = p.Bowtie.maximumThickness;
    
    bowtieThickness = a + ( b - a ) / 2 * ( 1 - cos(  angles ) ) ;
    
    
elseif strcmpi(p.Bowtie.shapeType, 'constant')
    
    bowtieThickness = p.Bowtie.maximumThickness * ones( size(uu), 'single');
    
elseif strcmpi(p.Bowtie.shapeType, 'quadratic')
    
    angles = p.Bowtie.alpha * atan( uu ./ geom.SDD );
    
    a = p.Bowtie.minimumThickness;
    b = p.Bowtie.maximumThickness;
    
    bowtieThickness = b * angles.^2 +  a ;
    
else
    error('unknown bowtie shape. \n');
end



%bowtieThickness = bowtieThickness';


end