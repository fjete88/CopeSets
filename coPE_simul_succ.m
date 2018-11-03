%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%        Covering rate of CoPE sets
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% Set working path
%cd '/Users/maullz/Desktop/CopeSets'
% Set working path
cd '/home/drtea/Research/MatlabPackages/CopeSets'

% Precompute fields y/n
precomp = 1;

% Sample size and simulation size
n    = 1e2;
nsim = 1e3;

% Considered confidence levels
lvls  = [0.85 0.90 0.95];
nlvls = length(lvls);

% Number of GPUs used for parallel computing
pool_num = 1;

% Error Field parameters
stddev   = [4, 4];
dim      = [50, 50];
noise    = 'normal'; % 't'; 'uniform'; %
nu       = '';
kernel   = 'gauss'; % 'quartic'; %
bin      = 0;
D        = length(dim);

% Signal parameters
SIGNAL_TYPE  = 'signal'; % 'SNR'; %
SIGNAL_SHAPE = 'linear'; %'circle';
SIGNAL_SD    = ones(dim);
param        = [1 3];%[5, 3, 3];

% Target SNR or signal strength by CoPe sets
c = 2;

% index for variable access on fields with different dimensions
index  = repmat( {':'}, 1, D );
index1 = repmat( {':'}, 1, D+1 );
index2 = repmat( {':'}, 1, D+2 );

%% Generate fields, if not precomputed
tic
if precomp == 1
    [Y, delta] = generateProcess( n, nsim, stddev, dim, noise, nu,...
                    kernel, bin, SIGNAL_SHAPE, param, SIGNAL_TYPE, SIGNAL_SD, pool_num );
end
toc
%
%%%%% Save error processes for later use to accelarate simulations (files might get large (N=200, nsim=1000 => 4gB!))
%save('N50Nsim1000LN50_isotropic_nu5.mat', '-v7.3')

%% Simulation of covering rate
% Initialize vector for threshold a
a_truebdry   = zeros(nlvls,nsim);
a_estimbdry  = zeros(nlvls,nsim);
a_estimbdry2 = zeros(nlvls,nsim);

% Initialize fields to save the thresholds and estimators
hatdelta = zeros([dim nsim]);
hatsigma = zeros([dim nsim]);
thresh_truebdry   = zeros([dim nlvls nsim 2]);
thresh_estimbdry  = zeros([dim nlvls nsim 2]);
thresh_estimbdry2 = zeros([dim nlvls nsim 2]);

quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );

% Simulation of estimated quantiles
switch(SIGNAL_TYPE)
    case 'signal'
        tic
        for k = 1:nsim
            % Obtain quantile estimate    
            [thresh_truebdry(index1{:},k,:), a_truebdry(:,k), ~, ~]...
                    = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'true', 1, 1, delta );
            [thresh_estimbdry(index1{:},k,:), a_estimbdry(:,k), ~, ~]...
                    = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'linear' );
            [thresh_estimbdry2(index1{:},k,:), a_estimbdry2(:,k), hatdelta(index{:},k), hatsigma(index{:},k)]...
                    = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'erodilation' );
        end
        toc
    case 'SNR'
        tic
        for k = 1:nsim
            % Obtain quantile estimate    
            [thresh_truebdry(index1{:},k,:), a_truebdry(:,k), ~, ~]...
                    = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'true', 't', delta );
            [thresh_estimbdry(index1{:},k,:), a_estimbdry(:,k), ~, ~]...
                    = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'linear', 't' );
            [thresh_estimbdry2(index1{:},k,:), a_estimbdry2(:,k), hatdelta(index{:},k), hatsigma(index{:},k)]...
                    = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'erodilation', 't' );

        end
        toc
end
%% % Compute the covering rate
tic
[covRate_truebdry]   = CovRateLvlSets( delta, hatdelta, thresh_truebdry, c, 0 );
[covRate_estimbdry]  = CovRateLvlSets( delta, hatdelta, thresh_estimbdry, c, 0 );
[covRate_estimbdry2] = CovRateLvlSets( delta, hatdelta, thresh_estimbdry2, c, 0 );
[covRate_truebdry_new]   = CovRateLvlSets( delta, hatdelta, thresh_truebdry, c, 1 );
[covRate_estimbdry_new]  = CovRateLvlSets( delta, hatdelta, thresh_estimbdry, c, 1 );
[covRate_estimbdry2_new] = CovRateLvlSets( delta, hatdelta, thresh_estimbdry2, c, 1 );
toc

%
[[covRate_truebdry];[covRate_truebdry_new];[covRate_estimbdry];[covRate_estimbdry_new];...
    [covRate_estimbdry2];[covRate_estimbdry2_new]]

% %% %%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%
% [thresh2, quantiles2, hatdelta2, hatsigma2] = CopeSets( Y(:,:,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );
% thresh2 = squeeze(thresh2(:,:,1,:) );
% 
% figure(1)
% subplot(1,2,1), hold on
% imagesc(hatdelta2(:,:)), colorbar
% contour(delta(1:L, 1:L), [1 1]*c, 'r', 'Linewidth', 2)
% contour(hatdelta2(:,:), [1 1]*c, 'k', 'Linewidth', 2)
% contour(hatdelta2(:,:)-thresh2(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
% contour(hatdelta2(:,:)-thresh2(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
% title('Excursion sets L=100, alpha=0.05'), hold off
% 
% [thresh3, quantiles3, hatdelta3, hatsigma3] = CopeSets( Y(1:2:L,1:2:L,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );
% thresh3 = squeeze(thresh3(:,:,1,:) );
% 
% subplot(1,2,2), hold on
% imagesc(hatdelta3(:,:)), colorbar
% contour(delta(1:2:L, 1:2:L), [1 1]*c, 'r', 'Linewidth', 2)
% contour(hatdelta3(:,:), [1 1]*c, 'k', 'Linewidth', 2)
% contour(hatdelta3(:,:)-thresh3(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
% contour(hatdelta3(:,:)-thresh3(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
% contour(hatdelta2(1:2:L,1:2:L)-thresh2(1:2:L,1:2:L,1), [1 1]*0, 'g--', 'Linewidth', 2)
% contour(hatdelta2(1:2:L,1:2:L)-thresh2(1:2:L,1:2:L,2), [1 1]*0, 'g--', 'Linewidth', 2)
% text(8, 10, 'L=50');
% text(5.1, 10, '----', 'Color', 'k', 'Linewidth', 2);
% text(8, 8, 'L=100');
% text(5.1, 8.6, '___', 'Color', 'g', 'Linewidth', 2);
% %legend( 'L=50', 'k--', 'L=100', 'g--' )
% title('Excursion sets L=50, alpha=0.05'), hold off,
% %% %%%% Plot confidence bands for different alpha
% [thresh2, quantiles2, hatdelta2, hatsigma2] = CopeSets( Y(:,:,:,1), c, 1-lvls, 5e3, 'linear', 1, 1 );
% 
% figure(2)
% hold on
% imagesc(hatdelta2(:,:)), colorbar
% contour(delta(1:L, 1:L), [1 1]*c, 'r', 'Linewidth', 2)
% contour(hatdelta2(:,:), [1 1]*c, 'k', 'Linewidth', 2)
% thresh = squeeze(thresh2(:,:,1,:) );
% contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'k--', 'Linewidth', 2)
% contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'k--', 'Linewidth', 2)
% thresh = squeeze(thresh2(:,:,2,:) );
% contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'g--', 'Linewidth', 2)
% contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'g--', 'Linewidth', 2)
% thresh = squeeze(thresh2(:,:,3,:) );
% contour(hatdelta2(:,:)-thresh(:,:,1), [1 1]*0, 'b--', 'Linewidth', 2)
% contour(hatdelta2(:,:)-thresh(:,:,2), [1 1]*0, 'b--', 'Linewidth', 2)
% text(8, 12, 'alpha=0.15');
% text(5.1, 12, '----', 'Color', 'k', 'Linewidth', 2);
% text(8, 10, 'alpha=0.1');
% text(5.1, 10, '----', 'Color', 'g', 'Linewidth', 2);
% text(8, 8, 'alpha=0.05');
% text(5.1, 8, '----', 'Color', 'b', 'Linewidth', 2);
% title('Excursion sets L=100, different alpha'), hold off