function img = extendVoi( img, k )

if length( size(img) ) == 2
    return;
end


for i = 1:k
   img(:,:,i) = ( img(:,:,i) + img(:,:,k+1) * ( k + 1 - i ) ) /  ( k + 2 - i );
   img(:,:,end-i+1) = ( img(:,:,end-i+1) + img(:,:,end-k) * ( k + 1 - i ) ) / ( k + 2 - i );
end



end