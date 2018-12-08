%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This script computes and saves some interesting simulations of the
%%% covering rates of Cope Sets
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% clear workspace
clear
close all

%%%%%% set path for the files
cd  /home/drtea/Research/MatlabPackages/CopeSets
[base, ~, ~] = fileparts(mfilename('fullpath'));
cd base;

%%%%%%%%%%%%%%%%%%%%%%%%%% Linear Ramp simulations %%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%% 50x50 domain simulation
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240 300];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
nsim = 3e3;
dim  = [50 50];
FWHM = 10; %[3 5]; %[5 10];
a    = [[1.5 2.5]; [0 4]; [-1 5]; [-2 6]]; % [[1 3];[0 4]]
% saving working memory by breaking down the simulation in batches of size
batchnumber = 20;

% Considered quantile estimation parameters
weightsVec = "rademacher"; %["rademacher", "gaussian"]
methodVec  = "t"; %["regular", "t"]


% define cell array to store the results of the simulations
simResults = cell(length(FWHM), size(a,1), length(nvec));

%% %%%%%% simulations of covering rates
% loop over bootstrap method
for weights = weightsVec
    for method = methodVec
        % quantile estimation parameters
        quantEstim = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', 5e3,...
                                                     'weights', weights,...   
                                                     'method', method)...
                                                        );
        % name for the file containing the simulations
        simname = strcat('ResultSim', num2str(dim(1)), num2str(dim(2)), 'Nsim',  num2str(nsim), ...
                    'maxN', num2str(max(nvec)), 'isotropicFWM',sprintf('%i_', FWHM), 'slopes',sprintf('%i_', diff(a')),'_boot_', method, quantEstim.params.weights, '.mat');
        % loop over the smoothing parameter
        % Initialize counter for filling the results cell structure
        countf = 0;
        for f = FWHM
            countf = countf+1;
            paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                                 'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );
            % loop over hte slope of the ramp signal
            for counta = 1:size(a,1)
                % parameters for the linear ramp
                paramSignal = struct('shape', "linear", 'shapeparam', a(counta,:), 'type', "signal");
                % loop over the sample size
                countn = 0;
                for n = nvec
                    countn = countn+1;
                    tic
                    simResults{countf, counta, countn} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                                paramNoise, quantEstim, batchnumber, pool_num );
                    toc
                    % Just in case it breaks saving the results
                    save(strcat('simulations/',simname), 'simResults', '-v7.3')
                end % end loop over sample size
            end % end loop over slope of the ramp
        end % end loop over the smoothing parameter
    end
end
%%%% save final results
save(strcat('simulations/',simname), 'simResults', '-v7.3')

%% %%%%%%%%%%% 100x50 domain simulation Gauss
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240 300];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
nsim = 3e3;
dim  = [100 50];
FWHM = [5 10];
a    = [[1.5 2.5]; [1.75 2.25]]; % [[1 3];[0 4]]
% saving working memory by breaking down the simulation in batches of size
batchnumber = 20;

% Considered quantile estimation parameters
weightsVec = "rademacher"; %["rademacher", "gaussian"]
methodVec  = "t"; %["regular", "t"]

% define cell array to store the results of the simulations
simResults = cell(length(FWHM), size(a,1), length(nvec));

%%%%%% simulations of covering rates
% loop over bootstrap method
for weights = weightsVec
    for method = methodVec
        % quantile estimation parameters
        quantEstim = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', 5e3,...
                                                     'weights', weights,...   
                                                     'method', method)...
                                                        );
        % name for the file containing the simulations
        simname = strcat('ResultSim', num2str(dim(1)), num2str(dim(2)), 'Nsim',  num2str(nsim), ...
                    'maxN', num2str(max(nvec)), 'isotropicFWM',sprintf('%i', FWHM),'_boot_', method, quantEstim.params.weights, '.mat');
        % loop over the smoothing parameter
        % Initialize counter for filling the results cell structure
        countf = 0;
        for f = FWHM
            countf = countf+1;
            paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                                 'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );
            % loop over hte slope of the ramp signal
            for counta = 1:size(a,1)
                % parameters for the linear ramp
                paramSignal = struct( 'shape', "linear", 'shapeparam', a(counta,:),...
                                      'type', "signal");
                % loop over the sample size
                countn = 0;
                for n = nvec
                    countn = countn+1;
                    tic
                    simResults{countf, counta, countn} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                                paramNoise, quantEstim, batchnumber, pool_num );
                    toc
                    % Just in case it breaks saving the results
                    save(strcat('simulations/',simname), 'simResults', '-v7.3')
                end % end loop over sample size
            end % end loop over slope of the ramp
        end % end loop over the smoothing parameter
    end
end
%%%% save final results
save(strcat('simulations/',simname), 'simResults', '-v7.3')

%% %%%%%%%%%%% 500x50 domain simulation Rademacher
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240];% 300];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
nsim = 3e3;
dim  = [500 50];
FWHM = [5 10];
a    = [[1.5 2.5]; [1.75 2.25]]; % [[1 3];[0 4]]
% saving working memory by breaking down the simulation in batches of size
batchnumber = 20;

% Considered quantile estimation parameters
weightsVec = "rademacher"; %["rademacher", "gaussian"]
methodVec  = "t"; %["regular", "t"]

% define cell array to store the results of the simulations
simResults = cell(length(FWHM), size(a,1), length(nvec));

%%%%%% simulations of covering rates
% loop over bootstrap method 
for weights = weightsVec
    for method = methodVec
        % quantile estimation parameters
        quantEstim = struct('name', 'MultiplierBootstrap',...
                                    'params', struct('Mboot', 5e3,...
                                                     'weights', 'rademacher',...   
                                                     'method', method)...
                                                        );
        % name for the file containing the simulations
        simname = strcat('ResultSim', num2str(dim(1)), num2str(dim(2)), 'Nsim',  num2str(nsim), ...
                    'maxN', num2str(max(nvec)), 'isotropicFWM',sprintf('%i', FWHM),'_boot_', method, quantEstim.params.weights, '.mat');
        % loop over the smoothing parameter
        % Initialize counter for filling the results cell structure
        countf = 0;
        for f = FWHM
            countf = countf+1;
            paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                                 'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );
            % loop over hte slope of the ramp signal
            for counta = 1:size(a,1)
                % parameters for the linear ramp
                paramSignal = struct( 'shape', "linear", 'shapeparam', a(counta,:),...
                                      'type', "signal");
                % loop over the sample size
                countn = 0;
                for n = nvec
                    countn = countn+1;
                    tic
                    simResults{countf, counta, countn} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                                paramNoise, quantEstim, batchnumber, pool_num );
                    toc
                    % Just in case it breaks saving the results
                    save(strcat('simulations/',simname), 'simResults', '-v7.3')
                end % end loop over sample size
            end % end loop over slope of the ramp
        end % end loop over the smoothing parameter
    end
end
%%%% save final results
save(strcat('simulations/',simname), 'simResults', '-v7.3')