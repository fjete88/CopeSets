%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Simulates the covering rate in 2D and 3D for SNR of signal plus noise
%   models
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear workspace
clear all
close all

%unit = 'personal';
unit = 'wias-server';


% add path to Hermite projecton estimator and Cope Set repository
if strcmp(unit, 'personal')
    addpath('../../CopeSets')
    addpath('../../HermiteProjector')
    addpath('/home/drtea/Documents/MATLAB/spm12')
elseif strcmp(unit, 'wias-server')
    addpath('~/projects/CopeSets')
    addpath('~/projects/HermiteProjector')
    addpath('~/projects/spm12')
end

clear unit

%%%%%% Define parameters of the simulation
% Simulation parameters
Msim  = 1%1e3;
Mboot = 5%e3; 
Nvec  = [30, 60, 120, 240, 400];
lvls  = [0.85 0.9, 0.95];

% parameters: error process
dim = [75 75];
[~, sdField ] = generateProcess( 1, 1, [4 4], dim, 'normal', '',...
    'gauss', 0, 'linear', [1, 2], 'signal', ones(dim), 0 );

paramsNoise = struct( 'FWHMvec', [3 5 8 12],...
                      'dim', dim,...
                      'noise', 'normal',...
                      'nu', '',...
                      'kernel', 'gauss',...
                      'bin', 0,...                        
                      'sddev', sdField);

% parameters: signal
rng(666);
[params, ~] = generateProcess( 1, 1, [14 18], dim, 'normal', '',...
                               'gauss', 0, 'none', '', 'signal', ones(dim), 0 );
paramsSignal = struct( 'signal', params / max(params(:)),...
                       'Amplitudevec', [0.5 1 2 4]);

% % Plot mean function
% figure(1)
% imagesc(params./sdField); colorbar

%%%% Define analyzed methods for estimation of the quantile and put them
%%%% into a cell
quantEstim = cell([1 5]);
quantEstim1 = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', Mboot,...
                                                     'weights', 'rademacher',...   
                                                     'method', 'regular',...
                                                     'mask', ones(dim)...
                                                     ));
quantEstim{1} = quantEstim1;
                                                 
quantEstim1 = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', Mboot,...
                                                     'weights', 'rademacher',...   
                                                     'method', 't',...
                                                     'mask', ones(dim)...
                                                     ));
quantEstim{2} = quantEstim1;

quantEstim1 = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', Mboot,...
                                                     'weights', 'gaussian',...   
                                                     'method', 'regular',...
                                                     'mask', ones(dim)...
                                                     ));
quantEstim{3} = quantEstim1;

quantEstim1 = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', Mboot,...
                                                     'weights', 'gaussian',...   
                                                     'method', 't',...
                                                     'mask', ones(dim)...
                                                     ));
quantEstim{4} = quantEstim1;

quantEstim1 = struct('name', "GKF",...
                    'params', struct('Mboot', Mboot,...
                                     'weights', "gaussian",...   
                                     'u', -6:0.01:6,...
                                     'mask', ones(dim),...
                                     'pool_num', 1)...
                                        );                                                 
quantEstim{5} = quantEstim1;
clear params sdField dim quantEstim1

%% %%%% simulate the covering rate and quantiles
[covRates, quantiles] = Sim_SNRSCB( Msim, Nvec, lvls, quantEstim, ...
                                             paramsSignal, paramsNoise );
save('Sim_SNRSCB_isotropicGauss.mat')