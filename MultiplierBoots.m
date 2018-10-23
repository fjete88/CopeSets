function [quantiles] = MultiplierBoots( R, alpha, Mboot, mask, method )

% Bootstraps the alpha-quantile of the maximum of a Gaussian field
% Input:
%  R:    residual field over a domain in R^D, it is an (D+1)-dimensional array,
%        where the last dimension enumerates the samples
%  alpha:     vector for which the quantile needs to be estimated
%  Mboot:     amount of bootstrap replicates (default=5e3)
%  mask:      array containing the values which are inside the domain (default=ones(size(R)))
%  method:    Either 't' or 'regular'. (Default 't')
% Output:
%  quantile is the bootstrapped quantile of the maximum distribution of the 
%  input processes
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/05/2018
%__________________________________________________________________________

% Check number of inputs.
if nargin > 5
    error('MultiplierBoots requires at most 3 optional inputs');
end

% Compute the dimension of the field
dimR  = size(R);

% Fill in unset optional values.
switch nargin
    case 2
        Mboot     = 5e3;
        mask      = ones(dimR(1:end-1));
        method    = 't';
    case 3
        mask      = ones(dimR(1:end-1));
        method    = 't';
    case 4
        method    = 't';
end

%%%%%% Get parameters for simulation
N    = dimR(end) ;

% Combine values of R in the mask to an matrix for faster computation
R       = reshape(R, prod(dimR(1:end-1)), N);
mask    = repmat( reshape(logical(mask), prod(dimR(1:end-1)), 1), [1 N] );

R = R(mask) ;
R = reshape(R, [length(R)/N N] );
      
% % Center/normalize, if required
% if center
%     R = R - repmat( mean(R,2), [1 N] );
% end
% if normalize
%    R = R ./ repmat( transpose(sqrt(var( transpose(R)))), [1 N] );
% end
    

%%%%% compute bootstrap replicates
% compute multipliers
multiplier = normrnd( 0, 1, [N,Mboot] );

% compute bootstrapped means
bootMeans      = R * multiplier / N;

if method == 't'
    bootSecMoments = R.^2 * multiplier.^2 / N;
    % we put an abs here to make sure that no NaNs are produced due to machine precision error.
    bootSigma = sqrt( abs(bootSecMoments - bootMeans.^2) / (N-1) * N );
else
    bootSigma = 1;
end

%compute bootstrapped values of maximum
bootMax = max(abs( sqrt(N)*bootMeans./bootSigma ));

%compute quantiles from the bootstrapp distribution of maximum
quantiles = quantile( bootMax, alpha );