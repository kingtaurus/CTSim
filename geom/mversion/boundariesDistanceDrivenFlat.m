function [umin, umax] = boundariesDistanceDrivenFlat(dsd_ds, ...
    tproj0, dsxy0, wu, dx_cos_2, dx_sin_2, dy_sin_2, n_dy_cos_2, alongx)

if alongx
    umin = dsd_ds * ( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ) - wu;
    umax = dsd_ds * ( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 ) - wu;
else
    umin = dsd_ds * ( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ) - wu;
    umax = dsd_ds * ( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 ) - wu;
end

%swap if umin > umax
if umin > umax
    temp = umax;
    umax = umin;
    umin = temp;
end

end