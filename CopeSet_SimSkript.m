%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This script computes and saves some interesting simulations
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% clear workspace
clear
close all

%%%%%% set path for the files
cd  /home/drtea/Research/MatlabPackages/CopeSets

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% precompute data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% paramters for the 50x50 simulation
dim  = [50 50];
FWHM = [5 10];

%%%%%% precompute the data as pieces
for f = FWHM
    for i = 1:3
        paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                         'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );

        eps = SmoothField2D( 240, 1e3, paramNoise.FWHM, paramNoise.dim, paramNoise.noise, paramNoise.nu,...
                                paramNoise.kernel, paramNoise.bin, 1 );
        save(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM', num2str(f),'_', num2str(i), '.mat'], ...
                'eps', 'paramNoise', '-v7.3')
        clear eps
    end
end
% clear the workspace from unneccesary variables
clear i f FWHM dim paramNoise

%%%%%% paramters for the 100x50 simulation
dim  = [100 50];
FWHM = [5 10];

%%%%%% precompute the data as pieces
for f = FWHM
    for i = 1:3
        paramNoise = struct( 'FWHM', [f, f], 'dim', dim, 'noise', "normal", 'nu', '',...
                         'kernel', "gauss", 'bin', 0, 'sd', ones(dim) );

        eps = SmoothField2D( 240, 1e3, paramNoise.FWHM, paramNoise.dim, paramNoise.noise, paramNoise.nu,...
                                paramNoise.kernel, paramNoise.bin, 1 );
        save(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM', num2str(f),'_', num2str(i), '.mat'], ...
                'eps', 'paramNoise', '-v7.3')
        clear eps
    end
end
% clear the workspace from unneccesary variables
clear i f FWHM dim paramNoise

%% %%%%%%%%%%%%%%%%%%%%%%% Linear Ramp simulations %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 50x50 domain simulation quantile estim t multiplier bootstrap
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
% quantile estimation parameters
quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );
nsim = 1e3;
dim  = [50 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim5050tResults = cell(length(FWHM), size(a,2), length(nvec), 3);

%%%%%% simulations of covering rates
% loop over the smoothing parameter
% Initialize counter for filling the results cell structure
countf = 0;
for f = FWHM
    countf = countf+1;
    % loop over precomputed errorfields
    for i = 1:3
        % load precomputed data
        load(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM',...
                num2str(f),'_', num2str(i), '.mat'])
        % loop over hte slope of the ramp signal
        counta = 0;
        for a = [1 4]
            counta = counta+1;
            % specify the paramter for the slope of the ramp
            if a==1
                aa = [1 3];
            else
                aa = a;
            end
            % parameters for the linear ramp
            paramSignal = struct('shape', "linear", 'shapeparam', aa, 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim5050tResults{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim5050tResults)
        save('simulations/sim5050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050tResults', '-v7.3')
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter

save('simulations/sim5050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050tResults', '-v7.3')
            
% compute the concatinated results and save
sim5050tResult = ConcatResults(sim5050tResults);
save('simulations/ResultSim5050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050tResult', '-v7.3');
            
%%%%%%%%%%% 50x50 domain simulation quantile estim regular multiplier
%%%%%%%%%%% bootstrap
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
% quantile estimation parameters
quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 'regular')...
                                                );
nsim = 1e3;
dim  = [50 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim5050regResults = cell(length(FWHM), size(a,2), length(nvec), 3);

%%%%%% simulations of covering rates
% loop over the smoothing parameter
% Initialize counter for filling the results cell structure
countf = 0;
for f = FWHM
    countf = countf+1;
    % loop over precomputed errorfields
    for i = 1:3
        % load precomputed data
        load(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM',...
                num2str(f),'_', num2str(i), '.mat'])
        % loop over hte slope of the ramp signal
        counta = 0;
        for a = [1 4]
            counta = counta+1;
            % specify the paramter for the slope of the ramp
            if a==1
                aa = [1 3];
            else
                aa = a;
            end
            % parameters for the linear ramp
            paramSignal = struct('shape', "linear", 'shapeparam', aa, 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim5050regResults{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim5050regResults)
        save('simulations/sim5050regNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050regResults', '-v7.3')
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter
save('simulations/sim5050regNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050regResults', '-v7.3')
% concatinate simulation results and save these results
sim5050regResult = ConcatResults(sim5050regResults);
save('simulations/ResultSim5050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim5050regResult', '-v7.3');         
%% %%%%%%%%%%% 100x50 domain simulation
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
% quantile estimation parameters
quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );
nsim = 1e3;
dim  = [100 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim10050tResults = cell(length(FWHM), size(a,2), length(nvec), 3);

%%%%%% simulations of covering rates
% loop over the smoothing parameter
% Initialize counter for filling the results cell structure
countf = 0;
for f = FWHM
    countf = countf+1;
    % loop over precomputed errorfields
    for i = 1:3
        % load precomputed data
        load(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM',...
                num2str(f),'_', num2str(i), '.mat'])
        % loop over hte slope of the ramp signal
        counta = 0;
        for a = [1 4]
            counta = counta+1;
            % specify the paramter for the slope of the ramp
            if a==1
                aa = [1 3];
            else
                aa = a;
            end
            % parameters for the linear ramp
            paramSignal = struct('shape', "linear", 'shapeparam', aa, 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim10050tResults{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim10050tResults)
        save('simulations/sim10050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050tResults', '-v7.3')
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter

save('simulations/sim10050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050tResults', '-v7.3')
            
% compute the concatinated results and save
sim10050tResult = ConcatResults(sim10050tResults);
save('simulations/ResultSim10050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050tResult', '-v7.3');
            
%%%%%%%%%%% 100x50 domain simulation quantile estim regular multiplier
%%%%%%%%%%% bootstrap
%%%%%% simulation paramters
clear all
close all
% Sample size and simulation size
nvec = [30 60 120 240];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% number of parallel computing units
pool_num = 1;
% quantile estimation parameters
quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 'regular')...
                                                );
nsim = 1e3;
dim  = [100 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim10050regResults = cell(length(FWHM), size(a,2), length(nvec), 3);

%%%%%% simulations of covering rates
% loop over the smoothing parameter
% Initialize counter for filling the results cell structure
countf = 0;
for f = FWHM
    countf = countf+1;
    % loop over precomputed errorfields
    for i = 1:3
        % load precomputed data
        load(['errorfields/dim',num2str(dim(1)), num2str(dim(2)),'Nsim1000LN240_isotropic_FWHM',...
                num2str(f),'_', num2str(i), '.mat'])
        % loop over hte slope of the ramp signal
        counta = 0;
        for a = [1 4]
            counta = counta+1;
            % specify the paramter for the slope of the ramp
            if a==1
                aa = [1 3];
            else
                aa = a;
            end
            % parameters for the linear ramp
            paramSignal = struct('shape', "linear", 'shapeparam', aa, 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim10050regResults{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim10050regResults)
        save('simulations/sim10050regNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050regResults', '-v7.3')
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter
save('simulations/sim10050regNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050regResults', '-v7.3')
% concatinate simulation results and save these results
sim10050regResult = ConcatResults(sim5050regResults);
save('simulations/ResultSim10050tNsim3000maxN240_isotropic_FWHM.mat', ...
                'sim10050regResult', '-v7.3');       