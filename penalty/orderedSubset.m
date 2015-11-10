function y =  orderedSubset( x, M, isub )

if M == 1
    y = x;
elseif length( size(x) ) == 2
    y = x(:,1+mod(isub,M):M:end);
elseif length( size(x) ) == 3
    y = x(:,:,1+mod(isub,M):M:end);
end