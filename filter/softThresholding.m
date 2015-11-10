function x = softThresholding( x, ctr, step )


x( x >= ctr - step & x <= ctr + step ) = ctr;
x( x > ctr + step ) = x( x > ctr + step ) - step;
x( x < ctr - step ) = x( x < ctr - step ) + step;


end