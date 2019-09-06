function [thresh, quantiles, hatdelta, hatsigma, len_bdry] = CopeSets_simultaneous( F, c, lvls, quantEstim,...
                                                        bdry_type, center, normalize, delta )
% Computes CoPe all ingredients for CoPe sets.
% Input:
%  F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%        where the last dimension enumerates the samples
%  c:    threshold for excursions
%  lvls:       vector containing the required confidence levels. Must be
%              between 0 and 1.
%  quantEstim: structure containing the name and the parameters for the
%              quantile estimation method. Choices:
%               {
%                quantEstim.name = 'multiplierbootstrap'
%                quantEstim.params:
%                   Mboot:     amount of bootstrap replicates (default=5e3)
%                   method:    option for the bootstrap estimator (default='t')
%               }
%  bdry_type: currently 'interval' and 'points' are supported.
%  center:    option to center the field using the sample mean (default=1)
%  normalize: option to normalize the field by sample variance (default=1)
%  delta:     required, if bdry_type is equal to 'true'. This is the true
%             population mean function given on a D-dimensional array
% Output:
%  - thresh is the threshold lower and upper for the sample mean in order to
%    be in the estimated lower and upper excursion sets 
%  - quantile is the bootstrapped quantile of the maximum distribution of the 
%    input processes
%  - hatdelta is the sample mean of the fields
%  - hatsigma is the sample variance of the fields
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 11/30/2018
%__________________________________________________________________________
%%%%%% Check user specified input
if any(lvls >=1) || any(lvls<=0)
    error("The vector lvls needs to have entries between 0 and 1!")
end

%%%%%% Fill in unset optional values.
switch nargin
    case 3
        quantEstim = struct('name', "MultiplierBootstrap",...
                            'params', struct('Mboot', 5e3,...
                                             'weights', "rademacher",...   
                                             'method', 't')...
                                                );
        bdry_type = 'linear';
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 4
        bdry_type = 'linear';
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 5
        center    = 1;
        normalize = 1;
        delta     = 666;
    case 6
        normalize = 1;
        delta     = 666;
    case 7
        delta     = 666;
end

%%%% compute the dimension of the domain of the field
D = length(size(F))-1;

%%%% Compute mean curve and variance
hatdelta = mean( F, D+1 ); 
hatsigma = std( F, 0, D+1 );

%%%% ensure there are values above the threshold!
if sum(hatdelta(:) >= c) == 0
    error("The threshold c is to high. There are now values exceeding it!")
end

%%%% Compute the process on the boundary and its mask
switch(bdry_type)
    case 'points'
        F_bdry = [];
        for cc = c
            F_bdry = [F_bdry; linBdryEstim( F, cc, hatdelta )];
        end
        mask   = ones([size(F_bdry,1) 1] );
        DD = 2;
    case 'interval'
        F_bdry = F;
        mask = hatdelta >= min(c) & hatdelta <= max(c);
        DD = D+1;
end

%%%% Compute length of estimated boundary
len_bdry = sum(mask);

%%%% convert F_bdry to residuals if neccessary
% Center/normalize, if required
if center
    F_bdry = F_bdry - mean(F_bdry,DD);
end
if normalize
   F_bdry = F_bdry ./ std( F_bdry, 0, DD );
end

%%%% Estimate the quantiles of the Gaussian process on the boundary
if strcmp( quantEstim.name, 'MultiplierBootstrap' )
    quantiles = MultiplierBoots( F_bdry, lvls, ...
                    quantEstim.params.Mboot, mask, quantEstim.params.weights,...
                    quantEstim.params.method );
elseif strcmp( quantEstim.name, 'GKF' ) && strcmp( bdry_type, 'interval' )
    LKCestim_HermProjExact( F_bdry, D, mask, quantEstim.params.Mboot,...
                            quantEstim.params.weights, quantEstim.params.u,...
                            quantEstim.params.pool_num )
    quantiles = repmat(-Inf, [1 length(lvls)]);
    for l = 1:length(lvls)
        [~, q] = crossing(EEC.hatn, u, lvls(l), 'linear');
        quantiles(l) = q(end);
    end
else
    error("Please specify a valid method for quantile estimation")
end

%%%% Compute upper and lower threshold for CoPE sets
sF = size(F);
N  = sF(end);
index = repmat( {':'}, 1, length(sF) );
sF = sF(1:end-1);

thresh = ones([ sF(1:end) length(lvls) 2 length(c)]);
countc = 0;
for cc = c
    countc = countc + 1;
    thresh(index{:},1,countc) = cc - repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                         .*repmat(hatsigma,[1 1 length(lvls)]) / sqrt(N);
    thresh(index{:},2,countc) = cc + repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                         .*repmat(hatsigma,[1 1 length(lvls)]) / sqrt(N);
end