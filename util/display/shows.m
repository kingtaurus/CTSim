function  shows( imgAtt )
% Evaluted image reconstruction results
%
% Meng Wu
% 2012 - 2013

if min( imgAtt(:) ) < -500
    window = [-800 1500];
else
    window = [0 0.4];
end


if length( size(imgAtt) ) == 3
    imgAtt = squeeze( imgAtt(:,:, ceil(size(imgAtt,3)/2) ) );
end

imshow( imgAtt', window);
