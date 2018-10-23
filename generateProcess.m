function [Y, delta] = generateProcess( n, nsim, stddev, dim, noise, nu,...
    kernel, bin, SIGNAL_SHAPE, SIGNAL_TYPE, SIGNAL_SD, pool_num )

% Generates an error process
% Input:
%   n        -  Sample size
%   nsim     -  Number of simulations
%   stddev   -  2-D or 3-D vector containing the std in the different dimensions for
%               smoothing
%   dim      -  dimensions of the field
%   noise    -  options are 'normal', 't', 'uniform' (default: 'normal')
%   nu       -  parameters for noise,
%                 't'       = degrees of freedom
%                 'uniform' = half length of interval
%   kernel   -  options 'gauss' and 'quartic' and 't-density'
%   bin      -  2x3 matrix allowing to bin parts of the cube to produce
%               non-stationary noise. Note that we need bin(1,i)*bin(2,i) < dim(i)
%   SIGNAL_SHAPE - Shape of the mean function of noise. Currently, 'linear'
%                  and 'quadratic' are supported.
%   SIGNAL_TYPE  - gives the signal type which is aimed for. Currently, 'signal'
%                  and 'SNR' are possible.
%   SIGNAL_SD    - Array of size dim specifying the standard deviation of
%                  the error field
%   pool_num -  number of GPUs used for parallizing, must be greater than 1
%               to enable.
% Output:
%   f        -  Array of size dim x n x nsim
%   delta    -  True signal or True SNR of the simulated process
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/22/2018
%__________________________________________________________________________

%%%%% dimension of the domain
D = length(dim);
if D ~= length(stddev)
    error('stddev and dim must have the same size!')
end
if size(SIGNAL_SD) ~= dim
    error('SIGNAL_SD needs to have dimension "dim"!')
end

%%%%% Generate noise fields and change its variance according to SIGNAL_SD
if D == 2
    eps = SmoothField2D( n*nsim, 1, stddev, dim, noise, nu,...
                                    kernel, bin, pool_num ) .* SIGNAL_SD;
elseif D == 3
    eps = SmoothField3D( n*nsim, 1, stddev, dim, noise, nu,...
                                    kernel, bin, pool_num ) .* SIGNAL_SD;
else
    error('Currently, only fields up to 3D domain are supported!')
end


%%%%%% Compute Signal plus Observations
aa = 4; % Slope of signal
if D==2
    [xx, yy] = meshgrid( (0.5:dim(1))/dim(1), (0.5:dim(2))/dim(2) );
    switch(SIGNAL_SHAPE)
        case 'linear',    mu = aa*xx;
        case 'quadratic', mu = aa*( 1 - (xx-dim(1)/2).^2 + (yy-dim(2)/2).^2);
    end    
else
    [xx, yy, zz] = meshgrid( (0.5:dim(1))/dim(1), (0.5:dim(2))/dim(2),...
                             (0.5:dim(3))/dim(3) );
    switch(SIGNAL_SHAPE)
        case 'linear',    mu = aa*xx;
        case 'quadratic', mu = aa*( 1 - (xx-dim(1)/2).^2 + (yy-dim(2)/2).^2 ...
                                      + (zz-dim(3)/2).^2 );
    end
end

% construct signal and reshape
Y = reshape(mu + eps, [dim n nsim]) ;

% Observed
switch(SIGNAL_TYPE)
    case 'signal', delta = mu;
    case 'SNR',    delta = mu ./ SIGNAL_SD;
end