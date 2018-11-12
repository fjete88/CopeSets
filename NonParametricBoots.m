function [quantiles, bootMax] = NonParametricBoots( R, alpha, Mboot, mask, method )

% Bootstraps alpha-quantiles of the maximum of a random field using a
% classical non parametric bootstrap.
% Input:
%  R:         random field over a domain in R^D, it is an (D+1)-dimensional array,
%             where the last dimension enumerates the samples
%  alpha:     vector of quantiles
%  Mboot:     amount of bootstrap replicates (default=5e3)
%  mask:      option of specifying an boolean array of size of the first D
%             components of size(R) containing the locations used for the
%             multplier bootstrap.
%             Default = ones(size(R)[1:D]))
%  method:    options are 't' for normalizing by the bootstrapped variance
%             or 'regular' for not normalizing.
%             Default = 't'.
% Output:
%  quantile is the bootstrapped quantile of the maximum distribution of the 
%  input processes
%  bootMax the bootstrap distribution of the maximum of the process
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/25/2018
%__________________________________________________________________________

% Check number of inputs.
if nargin > 5
    error('MultiplierBoots requires at most 3 optional inputs');
end

%%%%%% Get parameters for simulation
dimR = size(R) ;
N    = dimR(end) ;

%%%%%% Fill in unset optional values.
switch nargin
    case 2
        Mboot     = 5e3;
        if( length(dimR) == 2 )
           mask   = ones([dimR(1) 1]);
        else
           mask   = ones(dimR(1:end-1));
        end
        method    = 't';
    case 3
        if( length(dimR) == 2 )
           mask   = ones([dimR(1) 1]);
        else
           mask   = ones(dimR(1:end-1));
        end
        method    = 't';
    case 4
        method    = 't';
end

%%%%%% Combine values of R in the mask to an matrix for faster computation
R       = reshape(R, prod(dimR(1:end-1)), N);
mask    = repmat( reshape(logical(mask), prod(dimR(1:end-1)), 1), [1 N] );

R = R(mask) ;
R = reshape(R, [length(R)/N N] );        

%%%%%% compute bootstrap replicates
% compute multinomial multipliers from standard bootstrap
multiplier = mnrnd(N, 1/N*ones([1,N]), Mboot)';

% compute bootstrapped means
bootMeans = R * multiplier / N;

% Compute the sample mean of R
meanR = mean(R,2);

% Compute the variance depending on the method
if method == 't'
    % note that here is no square of the multiplier as in the multiplier
    % bootstrap! It is because the multinomial does not has variance 1 and
    % this version matches the varaiance of the sample with replacement.
    bootSecMoments = R.^2 * multiplier / N;
    % we put an abs here to make sure that no NaNs are produced due to machine precision error.
    bootSigma = sqrt( abs(bootSecMoments - bootMeans.^2) / (N-1) * N );
else
    bootSigma = 1;
end

%%%%%% compute bootstrapped values of maximum
bootMax = max(abs( sqrt(N)*(meanR - bootMeans)./bootSigma ));

%%%%%% compute quantiles from the bootstrapp distribution of maximum
quantiles = quantile( bootMax, alpha );