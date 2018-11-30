%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This script provides simulations of the approximation of the max
%%% distribution and the approximation using different quantile estimation
%%% methods
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% clear workspace
clear
close all

%%%%%% set path for the files
cd  /home/drtea/Research/MatlabPackages/CopeSets

%% %%%%%%%%%%%%%%% Parameters of the simulation
%%% General simulation parameters
% number of realisations used to approximate the max distribution
maxDistnsim = 25000;
% degrees of freedom for the t-process
nvec = [0,29,59,119, 239, 299];
%%% parameters of the error process
% dimensions
dim =  [100 50]; % [50, 50]; % [500 50];
% smoothing parameters of the noise
FWHM = [5 10];

%%% technical parameters
% number of parallel computing units
pool_num = 1;
% number of small batches used
batchnumber = 20;

% define cell array to store the results of the simulations
simResults = cell(1, length(FWHM));

%% %%%%%%%%%%%%%%% Simulate and save the different values for the maximum
% of the different used error process or its derived t-processes
for f = FWHM
    % construct the noise parameter structure
    paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                 'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );
    % construct boundary mask
    bdry       = zeros(dim);
    bdry(:,25) = 1;
    bdry       = logical(bdry);

    % Simulate the process and its t-version on the boundary
    tic
    simResults{1,f==FWHM} = maxDist_sim( maxDistnsim, nvec, paramNoise, bdry, batchnumber, pool_num );
    toc
    % name for the file containing the simulations
    simname = strcat( 'ResultMaxDistSim', num2str(dim(1)), num2str(dim(2)),...
                      'Nsim',  num2str(maxDistnsim), ...
                      'maxDF', num2str(max(nvec)), 'isotropicFWHM510_1','.mat');
    % Just in case it breaks saving the results
    save(strcat('simulations/',simname), 'simResults', '-v7.3')
end