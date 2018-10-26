function [Y, delta] = generateProcess( n, nsim, FWHM, dim, noise, nu,...
    kernel, bin, SIGNAL_SHAPE, param, SIGNAL_TYPE, SIGNAL_SD, pool_num )

% Generates an error process
% Input:
%   n        -  Sample size
%   nsim     -  Number of simulations
%   FWHM   -  2-D or 3-D vector containing the std in the different dimensions for
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
%                  'quadratic' and 'circle' are supported.
%   param    -  If SIGNAL_SHAPE is 'linear', param is either a number Y in
%               which case the signal is a ramp from 0 to X, or a vector
%               [X,Y] which generates a signal ramp from X to Y. If
%               SIGNAL_SHAPE is 'circle', param is a vector (rad, mag, smo)
%               where rad is radius of signal, mag is magniture and smo is
%               smoothing FWHM.
%   SIGNAL_TYPE  - gives the signal type which is aimed for. Currently, 'signal'
%                  and 'SNR' are possible.
%   SIGNAL_SD    - Array of size dim specifying the standard deviation of
%                  the error field
%   pool_num -  number of GPUs used for parallizing, must be greater than 1
%               to enable.
%
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
if D ~= length(FWHM)
    error('FWHM and dim must have the same size!')
end
if size(SIGNAL_SD) ~= dim
    error('SIGNAL_SD needs to have dimension "dim"!')
end

%%%%% Generate noise fields and change its variance according to SIGNAL_SD
if D == 2
    eps = SmoothField2D( n*nsim, 1, FWHM, dim, noise, nu,...
                                    kernel, bin, pool_num ) .* SIGNAL_SD;
elseif D == 3
    eps = SmoothField3D( n*nsim, 1, FWHM, dim, noise, nu,...
                                    kernel, bin, pool_num ) .* SIGNAL_SD;
else
    error('Currently, only fields up to 3D domain are supported!')
end


%%%%%% Compute Signal plus Observations
if D==2
    [xx, yy] = ndgrid( (0.5:dim(1))/dim(1), (0.5:dim(2))/dim(2) );
    switch(SIGNAL_SHAPE)
        case 'linear'
            aa = param;
            if size(aa,2) == 1
                mu = aa*xx;
            elseif size(aa,2) == 2
                mu = (aa(2) - aa(1))*xx + aa(1);
            end 
        case 'quadratic', mu = aa*( 1 - (xx-dim(1)/2).^2 + (yy-dim(2)/2).^2);
        case 'circle'
            rad   = param(1);
            mag   = param(2);
            smo   = param(3);
            
            cent0 = dim/2 + 1/2;
            Comm  = sprintf('Circle r=%02d sm=%1d',rad,smo);
            
            [x,y] = ndgrid([1:dim(1)], [1:dim(2)]);
            Rmap  = sqrt((x-cent0(1)).^2 + (y-cent0(2)).^2);
            mu    = mySmooth(Rmap<=rad,smo)*mag;
    end    
else
    [xx, yy, zz] =ndgrid( (0.5:dim(1))/dim(1), (0.5:dim(2))/dim(2),...
                             (0.5:dim(3))/dim(3) );
    switch(SIGNAL_SHAPE)
        case 'linear'    
            aa = param;             
            if size(aa,2) == 1
                mu = aa*xx;
            elseif size(aa,2) == 2
                mu = (aa(2) - aa(1))*xx + aa(1);
            end 
            mu = aa*xx;
        case 'quadratic', mu = aa*( 1 - (xx-dim(1)/2).^2 + (yy-dim(2)/2).^2 ...
                                      + (zz-dim(3)/2).^2 );
        case 'circle'
            rad = param(1);
            mag = param(2);
            smo = param(3);
            
            cent0 = dim/2 + 1/2;
            Comm  = sprintf('Sphere r=%02d sm=%1d',rad,smo);
            
            [x, y, z] = ndgrid([1:dim(1)],[1:dim(2)],[1:dim(3)]);
            Rmap = sqrt((x-cent0(1)).^2 + (y-cent0(2)).^2 + (z-cent0(3)).^2);
            mu    = mySmooth(Rmap<=rad,smo)*mag;        
    end
end

% construct signal and reshape
Y = reshape(mu + eps, [dim n nsim]) ;

% Observed
switch(SIGNAL_TYPE)
    case 'signal', delta = mu;
    case 'SNR',    delta = mu ./ SIGNAL_SD;
end