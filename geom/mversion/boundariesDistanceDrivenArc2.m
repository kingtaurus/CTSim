function [umin, umax] = boundariesDistanceDrivenArc2(sdd, add, du, ...
    tproj0, dsxy0, wu, dx_cos_2, dx_sin_2, dy_sin_2, n_dy_cos_2, alongx)
% For special curved detector by He Yang, Nov 2014
if alongx
    gam_min = atan(( tproj0 + dx_cos_2 ) / ( dsxy0 + dx_sin_2 ));
    xi_min = atan( sdd*sin(gam_min)*cos(gam_min)/(ADD-SDD*sin(gam_min)*sin(gam_min)) );
    umin = xi_min*add/du - wu;
    
    gam_max = atan(( tproj0 - dx_cos_2 ) / ( dsxy0 - dx_sin_2 )); 
    xi_max = atan( sdd*sin(gam_max)*cos(gam_max)/(ADD-SDD*sin(gam_max)*sin(gam_max)) );
    umax = xi_max*add/du - wu;
else
    gam_min = atan(( tproj0 + dy_sin_2 ) / ( dsxy0 + n_dy_cos_2 ));
    xi_min = atan( sdd*sin(gam_min)*cos(gam_min)/(ADD-SDD*sin(gam_min)*sin(gam_min)) );
    umin = xi_min*add/du - wu;
    
    gam_max = atan(( tproj0 - dy_sin_2 ) / ( dsxy0 - n_dy_cos_2 )); 
    xi_max = atan( sdd*sin(gam_max)*cos(gam_max)/(ADD-SDD*sin(gam_max)*sin(gam_max)) );
    umax = xi_max*add/du - wu;
end

%swap if umin > umax
if umin > umax
    temp = umax;
    umax = umin;
    umin = temp;
end

end