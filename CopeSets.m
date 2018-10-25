function [thresh, quantiles, hatdelta, hatsigma] = CopeSets( F, c, lvls, quantEstim,...
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
%  bdry_type: currently 'linear' or 'true' are supported
%  center:    option to center the field using the sample mean (default=1)
%  normalize: option to normalize the field by sample variance (default=1)
%  delta:     required, if bdry_type is equal to 'true'. This is the true
%            population mean function given on a D-dimensional array
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
% Last changes: 10/25/2018
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
                                             'weights', "gaussian",...   
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

%%%% Compute the process on the boundary and its mask
switch(bdry_type)
    case 'linear'
        F_bdry = linBdryEstim( F, c );
        mask   = ones([size(F_bdry,1) 1] );
    case 'erodilation'
        F_bdry = erodedilateBdryEstim( F, c, delta );
        mask   = ones([size(F_bdry,1) 1] );
    case 'true'
        F_bdry = linBdryEstim( F, c, delta );
        mask   = ones([size(F_bdry,1) 1] );
end

%%%% convert F_bdry to residuals if neccessary
% Center/normalize, if required
if center
    F_bdry = F_bdry - mean(F_bdry,2);
end
if normalize
   F_bdry = F_bdry ./ std( F_bdry, 0, 2 );
end

%%%% Estimate the quantiles of the Gaussian process on the boundary
if strcmp( quantEstim.name, 'MultiplierBootstrap' )
    quantiles = MultiplierBoots( F_bdry, lvls, ...
                    quantEstim.params.Mboot, mask, quantEstim.params.weights, quantEstim.params.method );
else
    error("Please specify a valid method for quantile estimation")
end

%%%% Compute upper and lower threshold for CoPE sets
sF = size(F);
N  = sF(end);
index = repmat( {':'}, 1, length(sF) );
sF = sF(1:end-1);

thresh = ones([ sF(1:end) length(lvls) 2]);
thresh(index{:},1) = c - repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                     .*repmat(hatsigma,[1 1 length(lvls)]) / sqrt(N);
thresh(index{:},2) = c + repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                     .*repmat(hatsigma,[1 1 length(lvls)]) / sqrt(N);