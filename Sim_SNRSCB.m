function [covRates, quantiles] = Sim_SNRSCB( Msim, Nvec, lvls, quantEstim, ...
                                             paramsSignal, paramsNoise )

% Simulates simultaneous confidence bands for Cohen's-d/signal-to-noise
% ratio for different methods.
% Input:
% F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%       where the last dimension enumerates the samples
% Output:
% residuals is an array of the same size as F containing the SNR residuals
% hatdelta is the sample mean of the fields
% hatsigma is the sample variance of the fields
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
%__________________________________________________________________________

% Get static parameters from input
Namplitudes = length(paramsSignal.Amplitudevec);
NFWHM       = length(paramsNoise.FWHMvec);
Nq = length(quantEstim);
index = repmat( {':'}, 1, length(paramsNoise.dim) );

% Initialize the covering rate output array
covRates  = zeros([length(Nvec), Nq, Namplitudes, NFWHM, length(lvls)]);
quantiles = zeros([length(Nvec), Nq, Namplitudes, NFWHM, length(lvls)]);

tic
% loop for MC simulation of covering rate                                                 
for m = 1:Msim
    for f = 1:NFWHM % loop FWHM 
        % generate random field
        [Y, ~] = generateProcess( Nvec(end), 1, paramsNoise.FWHM(f), paramsNoise.dim, paramsNoise.noise, paramsNoise.nu,...
                                  paramsNoise.kernel, paramsNoise.bin, 'predefined', zeros(paramsNoise.dim), 'signal',...
                                  paramsNoise.sddev, 0 );
        for n = 1:length(Nvec) % loop over sample size
            for a = 1:Namplitudes % loop over number of amplitudes
                signal_tmp = paramsSignal.Amplitudevec(a)*paramsSignal.signal;
                for q = 1:Nq % loop over number of quantile estimation methods
                    % compute simultaneous confidence bands and quantiles
                    [~, SCB, quantiles_tmp] = SCB_cohen( Y(index{:},1:Nvec(n))...
                                                         + signal_tmp,...
                                                           lvls, quantEstim{q});
                    % save quantiles
                    quantiles(n,q,a,f,:) = quantiles_tmp;

                    for l = 1:length(lvls) % loop over number of quantile levels
                        tmp1 = signal_tmp./paramsNoise.sddev > SCB{l}.SCBlow;
                        tmp2 = signal_tmp./paramsNoise.sddev < SCB{l}.SCBup;
                        % update covering, if SNR is contained in SCBs
                        % everywhere
                        if all(tmp1(:)) && all(tmp2(:))
                            covRates(n,q,a,f,l) = covRates(n,q,a,f,l)+1;
                        end
                    end % loop over number of quantile levels
                end % loop over number of quantile estimation methods
            end % loop over number of amplitudes
        end % loop over sample size
    end % loop FWHM 
end % loop for MC simulation of covering rate
toc

covRates = covRates / Msim;