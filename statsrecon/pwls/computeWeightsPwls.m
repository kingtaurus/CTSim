function w = computeWeightsPwls( sino, r , v)
% compute the weights for pwls algorithm
%
% Meng Wu at Stanford University
% 2013


w = single( ( sino - r ).^2 ./ (sino + v )) ;
w( sino <= r) = 0;

end