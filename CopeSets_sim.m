function results = CopeSets_sim( nsim, n, lvls, paramSignal, c, paramNoise, quantEstim, precomp, pool_num )

% Simulates the covering rate of Cope sets for various signals and error processes.
% Input:
%  nsim:      number of simulations used to estimate the covering rate
%  n:         sample size in each simulation. ("number of subjects")
%  lvls:      vector containing the confidence levels. Must be between 0 and 1.
%  paramSignal:  structure containing the parameters for the Signal
%  c:          targeted levelset
%  paramNoise: structure containing the parameters for the error process
%  quantEstim: structure containing the parameters for the quantile
%              estimation
%  precomp:   If '1' the random fields are precomputed for all nsim
%             simulations. Note that it might run into storage problems.
%             If '0' random fields are separetly computed in a for loop.
%  pool_num:  number of GPUs used for parallizing, must be greater than 1
%             to enable. (Default is 1)
%
% Output:
%  - results is a structure containg the results and parameter of the
%  simulation
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/25/2018
%__________________________________________________________________________

%%%%%% Fill in unset optional values.
switch nargin
    case 6
        quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );
        precomp = 1;
        pool_num = 1;
    case 7
        quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "gaussian",...   
                                             'method', 't')...
                                                );
        precomp  = 1;
        pool_num = 1;
    case 8
        pool_num  = 1;
end

%%%%%% Compute some useful constants for the simulation
% Get the number of considered confidence levels
nlvls  = length(lvls);
% Compute the dimension of the domain of the data
dim    = paramNoise.dim;
D      = length(dim);
% compute indices for variable access on fields with different dimensions
index  = repmat( {':'}, 1, D );
index1 = repmat( {':'}, 1, D+1 );

%%%%%% Initialize various variables
% Initialize vector for threshold a
a_truebdry       = zeros(nlvls,nsim);
a_linbdry        = zeros(nlvls,nsim);
a_erodbdry       = zeros(nlvls,nsim);

% Initialize fields to save the thresholds and estimators
hatdelta          = zeros([dim nsim]);
hatsigma          = zeros([dim nsim]);
thresh_truebdry   = zeros([dim nlvls nsim 2]);
thresh_linbdry    = zeros([dim nlvls nsim 2]);
thresh_erodbdry   = zeros([dim nlvls nsim 2]);

if precomp(1) == 1 || length(size(precomp)) > 2
    %%%%%% Simulate the fields for the simulation
    if precomp == 1
        tic
            [Y, delta] = generateProcess( n, nsim, paramNoise.FWHM, paramNoise.dim,...
                                          paramNoise.noise, paramNoise.nu, paramNoise.kernel,...
                                          paramNoise.bin, paramSignal.shape, paramSignal.shapeparam,...
                                          paramSignal.type, paramNoise.sd, pool_num );
        toc
    %%%%%% Save the precomputed fields into Y
    else
        tic
            [~, delta] = generateProcess( 1, 1, paramNoise.FWHM, paramNoise.dim,...
                                          paramNoise.noise, paramNoise.nu, paramNoise.kernel,...
                                          paramNoise.bin, paramSignal.shape, paramSignal.shapeparam,...
                                          paramSignal.type, paramNoise.sd, pool_num );
             Y = precomp;
             clear precomp;
        toc
    end
    %%%%%% Compute the quantiles and threshold fields from the precomputed
    %%%%%% fields
    % s
    switch( paramSignal.type )
        case 'signal'
            tic
            for k = 1:nsim
                % Obtain quantile estimate    
                [thresh_truebdry(index1{:},k,:), a_truebdry(:,k), ~, ~]...
                        = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'true', 1, 1, delta );
                [thresh_linbdry(index1{:},k,:), a_linbdry(:,k), ~, ~]...
                        = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'linear' );
                [thresh_erodbdry(index1{:},k,:), a_erodbdry(:,k), hatdelta(index{:},k), hatsigma(index{:},k)]...
                        = CopeSets( Y(index1{:},k), c, lvls, quantEstim, 'erodilation' );
            end
            toc
        case 'SNR'
            tic
            for k = 1:nsim
                % Obtain quantile estimate    
                [thresh_truebdry(index1{:},k,:), a_truebdry(:,k), ~, ~]...
                        = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'true', 't', delta );
                [thresh_linbdry(index1{:},k,:), a_linbdry(:,k), ~, ~]...
                        = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'linear', 't' );
                [thresh_erodbdry(index1{:},k,:), a_erodbdry(:,k), hatdelta(index{:},k), hatsigma(index{:},k)]...
                        = CopeSets_SNR( Y(index1{:},k), c, lvls, 5e3, 'erodilation', 't' );

            end
            toc
    end
    %%%%%% Compute the covering rate
    tic
    covRate_truebdry     = CovRateLvlSets( delta, hatdelta, thresh_truebdry, c, 0 );
    covRate_linbdry      = CovRateLvlSets( delta, hatdelta, thresh_linbdry,  c, 0 );
    covRate_erodbdry     = CovRateLvlSets( delta, hatdelta, thresh_erodbdry, c, 0 );
    covRate_truebdry_new = CovRateLvlSets( delta, hatdelta, thresh_truebdry, c, 1 );
    covRate_linbdry_new  = CovRateLvlSets( delta, hatdelta, thresh_linbdry,  c, 1 );
    covRate_erodbdry_new = CovRateLvlSets( delta, hatdelta, thresh_erodbdry, c, 1 );
    toc
end

%%%%%%% report the results of the simulation
% save the covering rate results
results.covRate.truebdry = [covRate_truebdry; covRate_truebdry_new];
results.covRate.linbdry = [covRate_linbdry; covRate_linbdry_new];
results.covRate.erodbdry = [covRate_erodbdry; covRate_erodbdry_new];
% standard error of the simulation
results.stdErr.rough = sqrt(lvls.*(1-lvls)/nsim);
results.stdErr.truebdry = sqrt(results.covRate.truebdry .* (1-results.covRate.truebdry)/nsim);
results.stdErr.linbdry = sqrt(results.covRate.linbdry .* (1-results.covRate.linbdry)/nsim);
results.stdErr.erodbdry = sqrt(results.covRate.erodbdry .* (1-results.covRate.erodbdry)/nsim);
% estimates of the quantile
results.quant.truebdry     = a_truebdry;
results.quant.linbdry      = a_linbdry;
results.quant.erodbdrybdry = a_erodbdry;
% params of the simulation
results.lvls = lvls;
results.c = c;
results.n = n;
results.nsim = nsim;
results.quantEstim  = quantEstim;
results.paramSignal = paramSignal;
results.paramNoise  = paramNoise;

