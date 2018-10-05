function [Y, delta, LKC] = generateProcess( NOISE_TYPE, N, nsim, D, L, params, b, SIGNAL_TYPE, SIGNAL_SHAPE )

% Interpolate values of random field to boundary of an excursion set using
% 4 -connectivity
% Input:
% field:     random field over a domain in R^2, it is an (1+1)-dimensional array,
%            where the last dimension enumerates realisations
% c:         threshold value for excursion
%Output:
% bdryValues is a vector containing the values of the random field on the
% boundary estimated using linear interpolation of a 4-connectivity grid

%%%%% Generate noise fields

% Noise
[eps, LKC] = generateField(NOISE_TYPE, N, nsim, D, L, params, b);

%%%%%% Compute Signal plus Observations
aa = 4;      % Slope of signal
t = (0.5:L)/L;
[xx, yy] = meshgrid(t, t);
switch(SIGNAL_SHAPE),
    case 'linear', mu = aa*xx;
    case 'quadratic', mu = aa*(1 - xx.^2 - yy.^2);
end

% construct signal
Y = repmat(mu, [1 1 N nsim]) + eps;

% Observed
switch(SIGNAL_TYPE),
    case 'signal', delta = mu;
    case 'SNR', delta = mu;%/sigma;
end
