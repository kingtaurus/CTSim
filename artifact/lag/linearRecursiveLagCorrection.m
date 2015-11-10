function x = linearRecursiveLagCorrection( y, b, a )

N = length(a);

S = zeros( size(y,1), size(y,2), N);
x = zeros( size(y), 'single');


for iv = 1:size(y, 3)
    
    for n = 1 : N
       x(:,:,iv) = y(:,:,iv) - b(n) * exp( - a(n) ) * S(:,:,n); 
    end
        
    x(:,:,iv) = x(:,:,iv) / ( 1 + sum(b) );
    
    for n = 1 : N
       S(:,:,n) = x(:,:,iv) +  exp( - a(n) ) * S(:,:,n); 
    end

end


end
