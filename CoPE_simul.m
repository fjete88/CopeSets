%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%        Covering rate of CoPE sets
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% Set working path
cd '/home/drtea/Research/MatlabPackages/CopeSets'

% Field parameters
NOISE_TYPE = 'isotropic';
% NOISE_TYPE = 'scale-space';
L0 = 1;     % EC of domain
D = 2;      % Dimension of domain (including scale)
L = 50;    % size of domain
N = 200;     % replicates (sample size)
nsim = 500;  % number of simulations
nu = 5;     % parameter for isotropic
gamma = 4:.2:40;  % bandwiths for scale-space
b = 5;      % size of filter domain in std
lvls = [0.85 0.90 0.95];

% Signal
SIGNAL_TYPE = 'SNR'; % 'signal';%
SIGNAL_SHAPE = 'linear';
aa = 4;      % Slope of signal


c = 2;      % Target SNR level


%% Generate fields
switch(NOISE_TYPE),
    case 'isotropic', params = nu;
    case 'scale-space', params = gamma;
end

% Noise
tic
[Y, delta, LKC] = generateProcess(NOISE_TYPE, N, nsim, D, L, params, b, SIGNAL_TYPE, SIGNAL_SHAPE);
toc
%%
%%%%% Save error processes for later use to accelarate simulations (files might get large (N=200, nsim=1000 => 4gB!))
%save('N50Nsim1000LN50_isotropic_nu5.mat', '-v7.3')
%% Compute Signal and Observations if precomputed error process is used
%load('N20Nsim1000L50_isotropic_nu5.mat')
load 'N200Nsim1000L50_isotropic_nu5.mat';

t = (0.5:L)/L;
[xx, yy] = meshgrid(t, t);
switch(SIGNAL_SHAPE),
    case 'linear', mu = aa*xx;
    case 'quadratic', mu = aa*(1 - xx.^2 - yy.^2);
end

% construct signal
y = repmat(mu, [1 1 N nsim]) + eps;

% Observed
switch(SIGNAL_TYPE),
    case 'signal', delta = mu;
    case 'SNR', delta = mu/ 1;%sigma;
end
clear xx yy t
Y=y;
%% % estimate empirical variance and mean of the fields
switch(SIGNAL_TYPE),
    case 'signal'
        deltahat = permute(mean(Y,3), [1 2 4 3]);
        deltahat2 = permute(mean(Y.^2,3), [1 2 4 3]);
        tic
        hatsigma = sqrt( N/(N-1)*(deltahat2 - deltahat.^2));
        toc
    case 'SNR'
        tic
        mY       = squeeze(mean(Y,3));
        sdY      = squeeze(std(Y, 0, 3));
        deltahat = mY ./sdY;
        hatsigma = sqrt(1+deltahat.^2/2 );
        toc
end

%% Simulation of covering rate
% Initialize vector for threshold a
a_truebdry  = zeros(3,nsim);
a_estimbdry = zeros(3,nsim);
a_estimbdry2 = zeros(3,nsim);

% Simulation of estimated quantiles
switch(SIGNAL_TYPE)
    case 'signal'
        tic
        for k = 1:nsim
            % Obtain quantile estimate    
            [threshtmp, a_truebdry(:,k), hatdeltatmp, hatsigmatmp] = CopeSets( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'true', 1, 1, delta );
            [threshtmp, a_estimbdry(:,k), hatdeltatmp, hatsigmatmp] = CopeSets( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'linear', 1, 1 );
            [threshtmp, a_estimbdry2(:,k), hatdeltatmp, hatsigmatmp] = CopeSets( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'erodilation', 1, 1 );
        end
        toc
    case 'SNR'
        tic
        for k = 1:nsim
            % Obtain quantile estimate    
            [threshtmp, a_truebdry(:,k), hatdeltatmp, hatsigmatmp] = CopeSets_SNR( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'true', delta );
            [threshtmp, a_estimbdry(:,k), hatdeltatmp, hatsigmatmp] = CopeSets_SNR( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'linear' );
            [threshtmp, a_estimbdry2(:,k), hatdeltatmp, hatsigmatmp] = CopeSets_SNR( Y(:,:,:,k),c,...
                1-lvls, 5e3, 'erodilation' );
        end
        toc
end
%%
% compute upper and lower threshold for CoPE sets
thresh_truebdry = ones([size(delta) 3 nsim 2]);
thresh_truebdry(:,:,:,:,1) = c - repmat(shiftdim(a_truebdry, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
thresh_truebdry(:,:,:,:,2) = c + repmat(shiftdim(a_truebdry, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
% compute upper and lower threshold for CoPE sets
thresh_estimbdry = ones([size(delta) 3 nsim 2]);
thresh_estimbdry(:,:,:,:,1) = c - repmat(shiftdim(a_estimbdry, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
thresh_estimbdry(:,:,:,:,2) = c + repmat(shiftdim(a_estimbdry, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
% compute upper and lower threshold for CoPE sets
thresh_estimbdry2 = ones([size(delta) 3 nsim 2]);
thresh_estimbdry2(:,:,:,:,1) = c - repmat(shiftdim(a_estimbdry2, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
thresh_estimbdry2(:,:,:,:,2) = c + repmat(shiftdim(a_estimbdry2, -2), [L L 1 1])...
    .*repmat(permute(shiftdim(hatsigma, -1), [2 3 1 4]), [1 1 3 1]) / sqrt(N);
% Compute the covering rate
tic
[covRate_truebdry]  = CovRateLvlSets( delta, deltahat, thresh_truebdry, c );
[covRate_estimbdry] = CovRateLvlSets( delta, deltahat, thresh_estimbdry, c );
[covRate_estimbdry2] = CovRateLvlSets( delta, deltahat, thresh_estimbdry2, c );
toc
%%
[[covRate_truebdry];[covRate_estimbdry];[covRate_estimbdry2]]
%% %%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%
[thresh2, quantiles2, hatdelta2, hatsigma2] = CopeSets( Y(:,:,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );
thresh2 = squeeze(thresh2(:,:,1,:) );

figure(1)
subplot(1,2,1), hold on
imagesc(hatdelta2(:,:)), colorbar
contour(delta(1:L, 1:L), [1 1]*c, 'r', 'Linewidth', 2)
contour(hatdelta2(:,:), [1 1]*c, 'k', 'Linewidth', 2)
contour(hatdelta2(:,:)-thresh2(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
contour(hatdelta2(:,:)-thresh2(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
title('Excursion sets L=100, alpha=0.05'), hold off

[thresh3, quantiles3, hatdelta3, hatsigma3] = CopeSets( Y(1:2:L,1:2:L,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );
thresh3 = squeeze(thresh3(:,:,1,:) );

subplot(1,2,2), hold on
imagesc(hatdelta3(:,:)), colorbar
contour(delta(1:2:L, 1:2:L), [1 1]*c, 'r', 'Linewidth', 2)
contour(hatdelta3(:,:), [1 1]*c, 'k', 'Linewidth', 2)
contour(hatdelta3(:,:)-thresh3(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
contour(hatdelta3(:,:)-thresh3(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
contour(hatdelta2(1:2:L,1:2:L)-thresh2(1:2:L,1:2:L,1), [1 1]*0, 'g--', 'Linewidth', 2)
contour(hatdelta2(1:2:L,1:2:L)-thresh2(1:2:L,1:2:L,2), [1 1]*0, 'g--', 'Linewidth', 2)
text(8, 10, 'L=50');
text(5.1, 10, '----', 'Color', 'k', 'Linewidth', 2);
text(8, 8, 'L=100');
text(5.1, 8.6, '___', 'Color', 'g', 'Linewidth', 2);
%legend( 'L=50', 'k--', 'L=100', 'g--' )
title('Excursion sets L=50, alpha=0.05'), hold off,
%% %%%% Plot confidence bands for different alpha
[thresh2, quantiles2, hatdelta2, hatsigma2] = CopeSets( Y(:,:,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );

figure(2)
hold on
imagesc(hatdelta2(:,:)), colorbar
contour(delta(1:L, 1:L), [1 1]*c, 'r', 'Linewidth', 2)
contour(hatdelta2(:,:), [1 1]*c, 'k', 'Linewidth', 2)
thresh = squeeze(thresh2(:,:,1,:) );
contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
thresh = squeeze(thresh2(:,:,2,:) );
contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'g--', 'Linewidth', 2)
contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'g--', 'Linewidth', 2)
thresh = squeeze(thresh2(:,:,3,:) );
contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'b--', 'Linewidth', 2)
contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'b--', 'Linewidth', 2)
text(8, 12, 'alpha=0.15');
text(5.1, 12, '----', 'Color', 'k', 'Linewidth', 2);
text(8, 10, 'alpha=0.1');
text(5.1, 10, '----', 'Color', 'g', 'Linewidth', 2);
text(8, 8, 'alpha=0.05');
text(5.1, 8, '----', 'Color', 'b', 'Linewidth', 2);
title('Excursion sets L=100, different alpha'), hold off