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

%%%%%% general simulation paramters
% Sample size and simulation size
nvec = [30 60 120 240];
nsim = 3e3;
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% threshold for the simulation
c    = 2;
% quantile estim structure
quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );
% number of parallel computing units
pool_num = 1;
 
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%% Linear Ramp simulations %%%%%%%%%%%%%%%%%%%%%%%%%
nsim = 1e3;

%%%%%%%%%%% 50x50 domain simulation
dim  = [50 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim5050_results = cell(length(FWHM), size(a,2), length(nvec), 3);

% Initialize counter for filling the results cell structure
countn = 0;
countf = 0;
counta = 0;

%%%%%% simulations of covering rates
% loop over the smoothing parameter
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
            paramSignal = struct('shape', "linear", 'shapeparam', [1 3], 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim5050_results{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim5050_results)
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter

save(['simulations/dim5050Nsim3000maxN240_isotropic_FWHM.mat'], ...
                'sim5050_results', '-v7.3')

%% %%%%%%%%%%% 100x50 domain simulation
dim  = [100 50];
FWHM = [5 10];
a    = [[1 3]; [0 5]];

% define index to access spatial domain of the fields
index  = repmat( {':'}, 1, length(dim) );

% define cell array to store the results of the simulations
sim10050_results = cell(length(FWHM), size(a,2), length(nvec), 3);

% Initialize counter for filling the results cell structure
countn = 0;
countf = 0;
counta = 0;

%%%%%% simulations of covering rates
% loop over the smoothing parameter
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
            paramSignal = struct('shape', "linear", 'shapeparam', [1 3], 'type', "signal");
            % loop over the sample size
            countn = 0;
            for n = nvec
                countn = countn+1;
                sim10050_results{countf, counta, countn, i} = CopeSets_sim( nsim, n, lvls, paramSignal, c,...
                                            paramNoise, quantEstim, eps(index{:},1:n, :), pool_num );    
            end % end loop over sample size
        end % end loop over slope of the ramp
        size(sim10050_results)
    end % end over loop of precomputed errorfields
end % end loop over the smoothing parameter

save(['simulations/dim10050Nsim3000maxN240_isotropic_FWHM.mat'], ...
                'sim10050_results', '-v7.3')