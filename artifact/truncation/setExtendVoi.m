function img = setExtendVoi( img, k, value )

if length( size(img) ) == 2
    return;
end


for i = 1:k
   img(:,:,i) = value ;
   img(:,:,end-i+1) = value ;
end



end