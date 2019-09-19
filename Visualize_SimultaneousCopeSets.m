%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Short script visualising the simultaneous CoPe sets.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

nsim = 100;
FWHM = [3 3];
aa   = [-3, 7];
dim  = [100 100];
Nvec = [30, 60, 120, 240, 400];
lvl  = 0.95;

strengthArray = [[0.5 0.9]; [1 1.5]; [3 6] ];

% signal = 'linear';
% params = aa;

signal = 'dcircle';
center    = 1;
normalize = 1;
weights = 'rademacher';
method  = 't';

% quantile estimation parameters
quantEstim = struct('name', 'MultiplierBootstrap',...
                            'params', struct('Mboot', 5e3,...
                                             'weights', weights,...   
                                             'method', method)...
                                                );
       
bdry_type = ["interval", "linear", "points"];

for nstrength = 1:size(strengthArray,1)
    strength = strengthArray(nstrength,:);
    params   = [9 strength(1) 9 9 strength(2) 9];
    [Y, delta] = generateProcess( Nvec(end), 1, FWHM, dim, 'normal', '',...
    'gauss', 0, signal, params, 'signal', ones(dim), 0 );

    for N = Nvec
        c = strength*0.75;

        location = strcat('pics/SimultaneousCopeSets_signal', signal, '_lowBd',...
                           num2str(10*strength(1)), '_upBd', num2str(10*strength(2)),...
                           '_N', num2str(N) );

        visualize_CopeSets( Y(:,:,1:N), c, lvl, quantEstim, bdry_type,...
                                          center, normalize, location, delta )
    end
end