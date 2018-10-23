function [thresh, quantiles, SNR, asymptSD] = CopeSets_SNR( F, c, lvl, Mboot, bdry_type, delta )

% Computes all ingredients for CoPe sets of the signal-plus-noise-ratio.
% Input:
% F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%       where the last dimension enumerates the samples
% c:    threshold for excursions
% lvl:       vector for which the quantile needs to be estimated
% Mboot:     amount of bootstrap replicates (default=1e3)
% bdry_type: currently 'linear', 'erodilation' or 'true' are supported
% delta:     required, if bdry_type is equal to 'true'. This is the true
%            population SNR function given on a D-dimensional array
%Output:
% thresh is the threshold lower and upper for the SNR in order to
% be in the estimated lower and upper excursion sets 
% quantile is the bootstrapped quantile of the maximum distribution of the 
% input processes
% hatdelta is the sample mean of the fields
% hatsigma is the sample variance of the fields
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/05/2018
%__________________________________________________________________________

% Fill in unset optional values.
switch nargin
    case 7
        delta     = 666;
end

%%%% Compute the SNR residuals and the SNR
[F, SNR, asymptSD] = SNR_residuals(F);

%%%% Compute the process on the boundary and its mask
switch(bdry_type)
    case 'linear'
        F_bdry = linBdryEstim( F, c, 1, SNR );
        mask   = ones([size(F_bdry,1) 1] );
    case 'erodilation'
        F_bdry = erodedilateBdryEstim( F, c, SNR );
        mask   = ones([size(F_bdry,1) 1] );
    case 'true'
        F_bdry = linBdryEstim( F, c, 1, delta );
        mask   = ones([size(F_bdry,1) 1] );
end

%%%% Estimate the quantiles of the Gaussian process on the boundary
quantiles = MultiplierBoots( F_bdry, lvl, Mboot, mask, 0, 0 );

%%%% Compute upper and lower threshold for CoPE sets
sF    = size(F);
N     = sF(end);
index = repmat( {':'}, 1, length(sF) );
sF    = sF(1:end-1);

thresh = ones([ sF(1:end) length(lvl) 2]);
thresh(index{:},1) = c - repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                     .*repmat(asymptSD,[1 1 length(lvl)]) / sqrt(N);
thresh(index{:},2) = c + repmat(shiftdim(quantiles, -length(sF)+1), [sF 1])...
                     .*repmat(asymptSD,[1 1 length(lvl)]) / sqrt(N);