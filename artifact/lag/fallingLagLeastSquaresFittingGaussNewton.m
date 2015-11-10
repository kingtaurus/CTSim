function [b, a, f] = fallingLagLeastSquaresFittingGaussNewton( y, w, N )

y = y(:);

K = length(y);
%initialization 
b = logspace( -1, -N, N );
a = logspace( -1, -N, N );
x = [ b(:) ; a(:)];



for itr = 1:10
    
    f = fallingLagResponse( a, b, w, K, N );
    r = y - f;
    J = fallingLagJacobian( a, b, w, K, N );
    
    if cond(  J' * J ) > 1e8
        break;
    end
    x = x  + ( J' * J ) \ ( J' *  r ) ;
    b = x(1:N);
    a = x(N+1:end);
    w = (1-sum(b.*exp(-a))) * y(1);
    %norm(r)
    
end

figure;
plot( y );
hold on; plot( f, '-.' ); 

end



function f = fallingLagResponse( a, b, w, K, N )

k = [1:K]';


for i = 1 : N

f = b(i) * exp( - a(i) * k ) ;
f(1) = f(1) + w;
f(2:end) = f(2:end) + w * b(i) * exp( - a(i) *  k(1:end-1) ) ;

end

end

function J = fallingLagJacobian( a, b, w, K, N )

k = [1:K]';
J = zeros( K, N * 2 );

for i = 1 : N
   J(:,i) =  exp( - a(i) * k )  ;
   J(2:end, i) = J(2:end, i) + w * exp( - a(i) * k(1:end-1) ) ;
   
   J(:,N+i) = - b(i) * k .* exp( - a(i) * k );
   J(2:end,N+i) = J(2:end,N+i) - b(i) * w * k(1:end-1) .*  exp( - a(i) * k(1:end-1) ) ;
   
end



end

