function proj = tilteProjction(  proj, v1, v2 )
% function proj = tilteProjction(  proj, v1, v2 )
%   Tilte the projection from vertical postion v1 to v2 in helical scan
%   
% Meng Wu, 2014.3

for i = 1:size(proj,1)
    proj(:,i) = interp1( v1(:,i), proj(:,i), v2(:,i) );
end
proj( isnan(proj) ) = 0;

end