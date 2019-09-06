function [hatCohen, SCB, quantiles, asymptSD] = SCB_cohen( F, lvls, quantEstim )
% Computes simultaneous confidence bands for one sample cohen's d.
% Input:
%  F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%        where the last dimension enumerates the samples
%  lvls: vector containing the required confidence levels. Must be
%        between 0 and 1.
%  quantEstim: structure containing the name and the parameters for the
%              quantile estimation method. Choices:
%               {
%                quantEstim.name = 'multiplierbootstrap'
%                quantEstim.params:
%                   Mboot:  amount of bootstrap replicates (default=5e3)
%                   method: option for the bootstrap estimator (default='t')
%               }
%               {
%                quantEstim.name = 'GKF'
%                quantEstim.params:
%                   Mboot:  amount of bootstrap replicates for LKC estimator (default=5e3)
%                   method: option for the bootstrap estimator (default='t')
%               }
% Output:
%  - hatCohen is the the estimate of Cohen's d from the data
%  - SCB is a structure containing the simultaneous confidence bands for
%    different quantiles
%  - quantiles are the quantiles estimated for the different confidence levels
%  - hatStd is the estimate of the asymptotic variance of hatCohen
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 08/22/2019
%__________________________________________________________________________
%%%%%% Check user specified input
if any(lvls >=1) || any(lvls<=0)
    error("The vector lvls needs to have entries between 0 and 1!")
end
sF = size(F);
%%%%%% Fill in unset optional values.
switch nargin
   case 1
       lvls = [0.9, 0.95, 0.99];
       quantEstim = struct('name', "GKF",...
                            'params', struct('Mboot', 3e3,...
                                             'weights', "gaussian",...   
                                             'u', -6:0.01:6,...
                                             'mask', ones(sF(1:end-1)),...
                                             'pool_num', 1)...
                                                );
    case 2
       quantEstim = struct('name', "GKF",...
                            'params', struct('Mboot', 3e3,...
                                             'weights', "gaussian",...   
                                             'u', -6:0.01:6,...
                                             'mask', ones(sF(1:end-1)),...
                                             'pool_num', 1)...
                                                );
end

%%%% compute the dimension of the domain of the field
D = length(sF)-1;

%%%% convert F to residuals
[residuals, hatCohen, asymptSD] = SNR_residuals( F );

%%%% Estimate the quantiles of the Gaussian process on the boundary
if strcmp( quantEstim.name, 'MultiplierBootstrap' )
    % estimate using the multiplier bootstrap
    quantiles = MultiplierBoots( residuals, 2*lvls, ...
                    quantEstim.params.Mboot,...
                    quantEstim.params.mask,...
                    quantEstim.params.weights,...
                    quantEstim.params.method );
elseif strcmp( quantEstim.name, 'GKF' )
    % Estimate LKCs of asymptotic process
    [~, EEC] = LKCestim_HermProjExact( residuals, D,...
                                       quantEstim.params.mask,...
                                       quantEstim.params.Mboot,...
                                       quantEstim.params.weights,...
                                       quantEstim.params.u,...
                                       quantEstim.params.pool_num );
    quantiles = repmat(-Inf, [1 length(lvls)]);
    for l = 1:length(lvls)
        % find quantiles from estimated expected Euler characteristic curve
        % of asymptotic process
        [~, q] = crossing(EEC.hatn, quantEstim.params.u, (1-lvls(l))/2, 'linear');
        quantiles(l) = q(end);
    end
else
    error("Please specify a valid method for quantile estimation")
end

%%%% Compute upper and lower threshold for CoPE sets
sF = size(F);
N  = sF(end);

%%%% compute bias correction
if N<250
    biasfac = gamma( (N-1)/2) / gamma((N-2)/2)*sqrt(2/(N-1));
else
    biasfac = 1;
end

SCB = cell([1 3]);

for l = 1:length(lvls)
    SCB{l} = struct( 'level',  lvls(l),...
                     'SCBlow', biasfac*hatCohen - quantiles(l)*asymptSD/sqrt(N),...
                     'SCBup',  biasfac*hatCohen + quantiles(l)*asymptSD/sqrt(N) );
end