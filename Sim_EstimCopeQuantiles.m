function [] = Sim_EstimCopeQuantiles( outpostfix, msim, Nsubj, Lvec, FWHM, SNR, lvls, Mboot)

% This function simulates the quantile estimator of the SNR processes and
% is written to parallize using CPUs
% Input:
%  outpostfix: string attached to the end of the simulation name
%  msim:       number of monte carlo simulations
%  Nsubj:      vector containing the considered number of subjects of the
%              estimator
%  Lvec:       vector containing the considered sizes of the 2D array
%  FWHM:       vector containg the FWHM of the gaussian smoothing kernel
%  SNR:        vector containing the signal to noise ratios
%  lvls:       vector containing the values of the quantiles
%  Mboot:      number of bootstrap samples used in the estimators
%  
% Output:
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 03/21/2019
%__________________________________________________________________________
%%%%% Fill default parameters
if ~exist('msim', 'var'), msim = 63; end
if ~exist('Nsubj', 'var'), Nsubj = [30 60 120]; end
if ~exist('Lvec', 'var'), Lvec = [10 60 124]; end
if ~exist('FWHM', 'var'), FWHM = 3; end
if ~exist('SNR', 'var'), SNR = [0.2 0.7 1 2 5]; end
if ~exist('mboot', 'var'), Mboot = 3e3; end
if ~exist('lvls', 'var'), lvls = [0.85, 0.9, 0.95]; end

%% %%%%%%%%%%%%%%% Bootstrap quantile estimator simulations %%%%%%%%%%%%%%%

% Initialize buckets for the estimated quantiles
TrueVarMGauss  = zeros([msim length(lvls) length(SNR) length(Lvec) length(Nsubj)]);
TrueVarMRadem  = TrueVarMGauss;
CohenVarMGauss = TrueVarMGauss;
CohenVarMRadem = TrueVarMGauss;
CohenVarMtGauss = TrueVarMGauss;
CohenVarMtRadem = TrueVarMGauss;
SNRVarMGauss = TrueVarMGauss;
SNRVarMRadem = TrueVarMGauss;
SNRVarMtGauss = TrueVarMGauss;
SNRVarMtRadem = TrueVarMGauss;

tic
for f = FWHM
    for m=1:msim
        Y = SmoothField2D( max(Nsubj), 1, [f f], [max(Lvec) max(Lvec)] );
        for c = SNR
            countc = find(c==SNR);
            % Get the data with mean c, meaning the SNR is also equal c
            Yc = Y + c;
            for n = Nsubj
                countn = find(n==Nsubj);
                
                % Get the subsampled data to subject size
                Ycn = Yc(:,:,1:n);
                % Get the SNR residuals and the estimated Cohen's-d
                % variance
                [SNRresYcn, ~, CohenVar] = SNR_residuals( Ycn );

                % remove the standardization of SNR residuals
                trueVar  = sqrt( 1 + c^2/2 );
                SNRVar   = std(SNRresYcn, 1, 3);


                for L = Lvec
                    countL = find(L==Lvec);
                    mask = ones([L L]);
                    SNRresYcnL = SNRresYcn(1:L,1:L,:);
                    SNRVarL    = SNRVar(1:L,1:L,:);
                    CohenVarL  = CohenVar(1:L,1:L,:);
                    
                    % Cope quantile estimates using different methods
                    TrueVarMGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ trueVar, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    TrueVarMRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ trueVar, 2*lvls, Mboot, ones([L L]), 'rademacher', 'regular' );

                    CohenVarMGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ CohenVarL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    CohenVarMRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ CohenVarL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );

                    CohenVarMtGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ CohenVarL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    CohenVarMtRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ CohenVarL, 2*lvls, Mboot, mask, 'rademacher', 't' );

                    SNRVarMGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    SNRVarMRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );

                    SNRVarMtGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    SNRVarMtRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'rademacher', 't' );
                end
            end
        end
    end
end
toc

sim_name = strcat('simulations/estimQuantile_SNRCopeSet_processes', outpostfix, '.mat');
save(sim_name, 'TrueVarMGauss', 'TrueVarMRadem', 'CohenVarMGauss', 'CohenVarMRadem', 'CohenVarMtGauss',...
    'CohenVarMtRadem', 'SNRVarMGauss', 'SNRVarMRadem', 'SNRVarMtGauss', 'SNRVarMtRadem',...
    'msim', 'Nsubj', 'Lvec', 'FWHM', 'SNR', 'lvls', 'Mboot')

end
