function [residuals, SNR, asymptSD] = SNR_residuals( F )

% Computes the residuals of the Gaussian limit of the SNR estimator and the
% SNR estimate.
% Input:
% F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%       where the last dimension enumerates the samples
% Output:
% residuals is an array of the same size as F containing the SNR residuals
% hatdelta is the sample mean of the fields
% hatsigma is the sample variance of the fields
%
% Author: Dr. Fabian J.E. Telschow
% Last Changes: Oct. 5 2018

%%%% Compute mean curve and variance
hatdelta = mean( F, 3 ); 
hatsigma = std( F, 0, 3 );

%%%% Compute the estimated SNR
SNR      = hatdelta ./ hatsigma;
%%%% Compute the asmptotic variance of \sqrt(N) (\hat SNR - SNR)
asymptSD = sqrt(1+SNR.^2/2 );

%%%% Compute the residuals (R2016>)
residuals = ( (F-hatdelta)./ hatsigma - SNR/2.*(((F-hatdelta)./hatsigma).^2-1) )...
                ./ asymptSD;
end

