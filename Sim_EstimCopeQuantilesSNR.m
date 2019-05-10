function [] = Sim_EstimCopeQuantilesSNR( outpostfix, msim, Nsubj, Lvec, FWHM, SNR, lvls, Mboot)

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
if ~exist('msim', 'var'), msim = 125; end
if ~exist('Nsubj', 'var'), Nsubj = [30 60 120]; end
if ~exist('Lvec', 'var'), Lvec = [10 60 124]; end
if ~exist('FWHM', 'var'), FWHM = 3; end
if ~exist('SNR', 'var'), SNR = [0.2 0.7 2]; end
if ~exist('Mboot', 'var'), Mboot = 5e3; end
if ~exist('lvls', 'var'), lvls = [0.85, 0.9, 0.95]; end

%% %%%%%%%%%%%%%%% Bootstrap quantile estimator simulations %%%%%%%%%%%%%%%

% Initialize buckets for the estimated quantiles
TrueStdMGauss   = zeros([msim length(lvls) length(SNR) length(Lvec) length(Nsubj)]);
TrueStdMRadem   = TrueStdMGauss;
TrueStdMtGauss  = TrueStdMGauss;
TrueStdMtRadem  = TrueStdMGauss;

AsymStdMGauss   = TrueStdMGauss;
AsymStdMRadem   = TrueStdMGauss;
AsymStdMtGauss  = TrueStdMGauss;
AsymStdMtRadem  = TrueStdMGauss;

CohenStdMGauss  = TrueStdMGauss;
CohenStdMRadem  = TrueStdMGauss;
CohenStdMtGauss = TrueStdMGauss;
CohenStdMtRadem = TrueStdMGauss;

SNRStdMGauss    = TrueStdMGauss;
SNRStdMRadem    = TrueStdMGauss;
SNRStdMtGauss   = TrueStdMGauss;
SNRStdMtRadem   = TrueStdMGauss;

StabStdMGauss   = TrueStdMGauss;
StabStdMRadem   = TrueStdMGauss;
StabStdMtGauss  = TrueStdMGauss;
StabStdMtRadem  = TrueStdMGauss;

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
                
                % constant which might cause numerical failure, so
                % approximation for large n is used
                if n < 200
                    fac = sqrt((n-1)/2)*(gamma((n-2)/2)/gamma((n-1)/2));
                else
                    fac = 1 / ( 1 - (3/(4*(n-1)-1)) );
                end

                % true variance at d=c derived from the variance of
                % non-central t
                TrueStd  = sqrt( (n-1)^2/n/(n-3) + ( (n-1)^2/(n-3)-fac^2*(n-1) )*c^2 );
                
                % get constants for variance stabilisation
                a     = (n-1)/n/sqrt(n-3);
                b     = sqrt( (n-1)/fac^2/(n-3)-1 ) * (n-1)/n;
                alpha = 1/b/sqrt(n);
                beta  = b/a;
                
                % Get the SNR residuals
                [SNRresYcn, etaYcn, CohenStd] = SNR_residuals( Yc(:,:,1:n) );
                
                % further possibilities for variance normalisation
                AsymStd  = sqrt( 1 + c^2/2 );
                SNRStd   = std(SNRresYcn, 0, 3);
                
                % stabilising factor for residuals
                StabSNRres = (alpha*beta / sqrt( beta*etaYcn^2+1 )) .* SNRresYcn;

                for L = Lvec
                    countL      = find(L==Lvec);
                    mask        = ones([L L]);
                    
                    SNRresYcnL  = SNRresYcn(1:L,1:L,:);
                    StabSNRresL = StabSNRres(1:L,1:L,:);
                    SNRStdL     = SNRStd(1:L,1:L);
                    CohenStdL   = CohenStd(1:L,1:L);
                    
                    % Cope quantile estimates using different methods
                    TrueStdMGauss( m, :, countc, countL, countn)   = MultiplierBoots( SNRresYcnL ./ TrueStd, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    TrueStdMRadem( m, :, countc, countL, countn)   = MultiplierBoots( SNRresYcnL ./ TrueStd, 2*lvls, Mboot, ones([L L]), 'rademacher', 'regular' );
                    TrueStdMtGauss( m, :, countc, countL, countn)  = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    TrueStdMtRadem( m, :, countc, countL, countn)  = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, ones([L L]), 'rademacher', 't' );
                    
                    AsymStdMGauss( m, :, countc, countL, countn)   = MultiplierBoots( SNRresYcnL ./ AsymStd, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    AsymStdMRadem( m, :, countc, countL, countn)   = MultiplierBoots( SNRresYcnL ./ AsymStd, 2*lvls, Mboot, ones([L L]), 'rademacher', 'regular' );
                    AsymStdMtGauss( m, :, countc, countL, countn)  = TrueStdMtGauss( m, :, countc, countL, countn);
                    AsymStdMtRadem( m, :, countc, countL, countn)  = TrueStdMtRadem( m, :, countc, countL, countn);
                    
                    CohenStdMGauss( m, :, countc, countL, countn)  = MultiplierBoots( SNRresYcnL ./ CohenStdL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    CohenStdMRadem( m, :, countc, countL, countn)  = MultiplierBoots( SNRresYcnL ./ CohenStdL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );
                    CohenStdMtGauss( m, :, countc, countL, countn) = TrueStdMtGauss( m, :, countc, countL, countn);
                    CohenStdMtRadem( m, :, countc, countL, countn) = TrueStdMtRadem( m, :, countc, countL, countn);

                    SNRStdMGauss( m, :, countc, countL, countn)    = MultiplierBoots( SNRresYcnL ./ SNRStdL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    SNRStdMRadem( m, :, countc, countL, countn)    = MultiplierBoots( SNRresYcnL ./ SNRStdL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );
                    SNRStdMtGauss( m, :, countc, countL, countn)   = TrueStdMtGauss( m, :, countc, countL, countn);
                    SNRStdMtRadem( m, :, countc, countL, countn)   = TrueStdMtRadem( m, :, countc, countL, countn);

                    StabStdMGauss( m, :, countc, countL, countn)   = MultiplierBoots( StabSNRresL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    StabStdMRadem( m, :, countc, countL, countn)   = MultiplierBoots( StabSNRresL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );
                    StabStdMtGauss( m, :, countc, countL, countn)  = MultiplierBoots( StabSNRresL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    StabStdMtRadem( m, :, countc, countL, countn)  = MultiplierBoots( StabSNRresL, 2*lvls, Mboot, mask, 'rademacher', 't' );
                end
            end
        end
    end
end
toc

sim_name = strcat('simulations/estimQuantile_SNRCopeSet_processes', outpostfix, '.mat');
save(sim_name, 'TrueStdMGauss', 'TrueStdMRadem', 'TrueStdMtGauss', 'TrueStdMtRadem',...
               'AsymStdMGauss', 'AsymStdMRadem', 'AsymStdMtGauss', 'AsymStdMtRadem',...
               'CohenStdMGauss', 'CohenStdMRadem', 'CohenStdMtGauss','CohenStdMtRadem', ...
               'SNRStdMGauss', 'SNRStdMRadem', 'SNRStdMtGauss', 'SNRStdMtRadem',...
               'StabStdMGauss', 'StabStdMRadem', 'StabStdMtGauss', 'StabStdMtRadem',...
               'msim', 'Nsubj', 'Lvec', 'FWHM', 'SNR', 'lvls', 'Mboot')

end
