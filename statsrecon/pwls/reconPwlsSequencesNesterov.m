function [xs, phis, regs, Betas ] = reconPwlsSequencesNesterov( y, w, geom, betas, pfun, delta, numos, nitn, noFrames )

Betas = logspace( log10(betas(1)), log10(betas(2)), noFrames );

xs = zeros( [ geom.reconSize(1) geom.reconSize(2) noFrames], 'single' );
phis = zeros( noFrames, 1);
regs = zeros( noFrames, 1);

R  = loadPenaltyOperator( pfun, delta );

for i = 1 : noFrames 
    
    if i == 1
       [x, ps] = reconPwlsSeNesterovSqs( y, w, geom, Betas(i), pfun, nitn, delta, numos );
    else
       [x, ps] =  reconPwlsSeNesterovSqs( y, w, geom, Betas(i), pfun, round(nitn/2), delta, numos, x );
    end
    
    xs(:,:,i) = x(:,:,ceil(end/2));
    phis(i) = ps(end);
    regs(i) = R(x);


end

