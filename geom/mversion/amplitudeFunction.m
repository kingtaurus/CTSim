function [ lu, alongx ] = amplitudeFunction(dsxy0, tproj0, tproj0_dsxy0, sin_a, cos_a)

cos_gamma = dsxy0 / tproj0_dsxy0;
sin_gamma = tproj0 / tproj0_dsxy0;
cos_phi = abs ( cos_a * cos_gamma - sin_a * sin_gamma );
sin_phi = abs ( sin_a * cos_gamma + cos_a * sin_gamma );


if cos_phi > sin_phi
    alongx = true;
    lu = 1 / cos_phi;
else
    alongx = false;
    lu = 1 / sin_phi;
end

end