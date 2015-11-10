function [ weights, isproj, iumin, iumax ] = weightDistanceDriven(umin, umax, nu)


umin = max( umin + 0.5 , 0);
umax = min( umax + 0.5 , nu);

iumin = floor(umin);
iumax = floor(umax);
iumax = min(iumax, nu-1);

if umin >= umax
    isproj = false;
    weights = 0;
    return;
else
    isproj = true;
end

iuLength = iumax-iumin+1;

weights = ones(1, iuLength);

if (iuLength == 1)
    weights(1) = umax - umin;
else
    weights(1) = iumin + 1 - umin;
    weights(iuLength) =  umax - iumax;
end

end