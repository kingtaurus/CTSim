
close all;

N = 10000;
M = 100;

Lambda = logspace(2,8,M);

electronicNoiseStd1 = 5;

electronicNoiseStd2 = 20;

mu = zeros(1,M);
va = zeros(1,M);
mue1 = zeros(1,M);
vae1 = zeros(1,M);
mue2 = zeros(1,M);
vae2 = zeros(1,M);

for i = 1:M
    
    r =  Lambda(i);
    
    if r > 1e6
        mu(i) = log(r);
        va(i) = 1/ r;
        mue1(i) = log(r);
        vae1(i) = 1/ r;
        mue2(i) = log(r);
        vae2(i) = 1/ r;
        continue;
    end
    
    
    x = poissrnd(r, [N 1]);
    x( x < 1) = 1;
    mu(i) = mean( log(x));
    va(i) = var( log(x));
    
    x = poissrnd(r, [N 1]) + electronicNoiseStd1 * randn([N, 1]);
    x( x <= 1) =  1;
    mue1(i) = mean( log(x));
    vae1(i) = var( log(x));
    
    x = poissrnd(r, [N 1]) + electronicNoiseStd2 * randn([N, 1]);
    x( x <= 1) =   1;
    mue2(i) = mean( log(x));
    vae2(i) = var( log(x));
    
    
end

%%
figure(1);
semilogx( Lambda, mu); hold on;
semilogx( Lambda, mue1, '--');
semilogx( Lambda, mue2, '-.');
xlabel( 'Expected photon counts ', 'fontSize', 14);
ylabel( 'Means of log signals g', 'fontSize', 14);
legend( '\sigma_e = 0', '\sigma_e = 10', '\sigma_e = 20');
set( gca, 'fontSize', 12);


figure(2);
semilogx( Lambda, mu - log(Lambda)); hold on;
semilogx( Lambda, mue1 - log(Lambda), '--');
semilogx( Lambda, mue2 - log(Lambda), '-.');
xlabel( 'Expected photon counts ', 'fontSize', 14);
ylabel( 'Bias of log signals g', 'fontSize', 14);
legend( '\sigma_e = 0', '\sigma_e = 10', '\sigma_e = 20');
set( gca, 'fontSize', 12);

figure(3);
loglog( Lambda, va  ); hold on;
loglog( Lambda,  vae1 , '--');
loglog( Lambda,  vae2 , '-.');
xlabel ('Expected photon counts ', 'fontSize', 14);
ylabel ('Variances of log signals', 'fontSize', 14);
legend( '\sigma_e = 0', '\sigma_e = 10', '\sigma_e = 20');
set( gca, 'fontSize', 12);