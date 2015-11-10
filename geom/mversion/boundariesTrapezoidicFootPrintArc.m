function [taus] = boundariesTrapezoidicFootPrintArc(dsd_ds, tproj0, dsxy0, wu, dx_cos_2, dx_sin_2, dy_sin_2, n_dy_cos_2)

% define ratio(x,y) (tproj0 x dx_cos_2 y dy_sin_2 ) / (dsxy0 x dx_sin_2 y n_dy_cos_2)
% define taufun(x,y) dsd_ds * ratio(x,y)

taus = zeros(1,4);

taus(1) = dsd_ds * atan(( tproj0 + dx_cos_2 + dy_sin_2 ) / ( dsxy0 + dx_sin_2 + n_dy_cos_2 )) - wu;
taus(2) = dsd_ds * atan(( tproj0 + dx_cos_2 - dy_sin_2 ) / ( dsxy0 + dx_sin_2 - n_dy_cos_2 )) - wu;
taus(3) = dsd_ds * atan(( tproj0 - dx_cos_2 + dy_sin_2 ) / ( dsxy0 - dx_sin_2 + n_dy_cos_2 )) - wu;
taus(4) = dsd_ds * atan(( tproj0 - dx_cos_2 - dy_sin_2 ) / ( dsxy0 - dx_sin_2 - n_dy_cos_2 )) - wu;

taus = sort(taus);

end