function [ weights, isproj, iumin, iumax ] = weightTrapezoidicFootPrint(taus, nu)


iumin = max (floor(taus(1) + 0.5),0);
iumax = min (ceil(taus(4) + 0.5), nu-1);

if iumin > iumax
    isproj = false;
    weights = 0;
    return;
else
    isproj = true;
end

weights = zeros(1, iumax-iumin+1);

i = 1;
%compute the trapeziod footprint
for iu = iumin : iumax
    
    temp = 0;
    bl = iu - 0.5;
    br = iu + 0.5;
    % left piece
    xl = max(bl, taus(1));
    xr = min(br, taus(2));
    if (xr > xl)
        temp = temp + ((xr - taus(1))^2 - (xl - taus(1))^2) / (taus(2) - taus(1)) / 2;
    end
    % middle piece
    xl = max(bl, taus(2)); 
    xr = min(br, taus(3));
    if (xr > xl)
        temp = temp + (xr - xl);
    end
    % right piece
    xl = max(bl, taus(3)); 
    xr = min(br, taus(4));
    if (xr > xl)
        temp = temp + ((xl - taus(4))^2 - (xr - taus(4))^2) / (taus(4) - taus(3)) / 2;
    end
    
    weights(i) = temp ;
    i = i+1;
end

end