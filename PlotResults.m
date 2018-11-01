function[] = PlotResults( result, location, type, titlename, stdfac, colorVec, dotVec)
% Plots results from a simulation and outputs figures for the different lvls.
% Input:
%  F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%        where the last dimension enumerates the samples
%  c:    threshold for excursions
%  lvls:       vector containing the required confidence levels. Must be
%              between 0 and 1.
%  quantEstim: structure containing the name and the parameters for the
%              quantile estimation method. Choices:
%               {
%                quantEstim.name = 'multiplierbootstrap'
%                quantEstim.params:
%                   Mboot:     amount of bootstrap replicates (default=5e3)
%                   method:    option for the bootstrap estimator (default='t')
%               }
%  bdry_type: currently 'linear' or 'true' are supported
%  center:    option to center the field using the sample mean (default=1)
%  normalize: option to normalize the field by sample variance (default=1)
%  delta:     required, if bdry_type is equal to 'true'. This is the true
%            population mean function given on a D-dimensional array
% Output:
%  - thresh is the threshold lower and upper for the sample mean in order to
%    be in the estimated lower and upper excursion sets 
%  - quantile is the bootstrapped quantile of the maximum distribution of the 
%    input processes
%  - hatdelta is the sample mean of the fields
%  - hatsigma is the sample variance of the fields
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/25/2018
%__________________________________________________________________________

switch nargin
    case 4
        stdfac   = 1.96;
        colorVec = 'rbmc';
        dotVec   = 'xxoo';
    case 5
        colorVec = 'rbmc';
        dotVec   = 'xxoo';

    case 6
        dotVec   = 'xxoo';
end

% Set font to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(gca, 'fontsize', 20);

dim    = size(result);
result = reshape(result, [prod(dim(1:end-1)), dim(end)]);
dim    = size(result);
% initialize FWHMnames
FWHMnames = [];
slopenames = [];
% get the vector with the subjects
nVec   = 1:dim(2);
for l = 1:dim(2)
    nVec(l) = result{1,l}.n;
end
% get vector of levels
lvls   = result{1,1}.lvls;

% restructure to get matrix with covering rates and with the scaled stdErr
covRates = zeros([ dim ,length(lvls)] );
scStdErr = stdfac*result{1,1}.stdErr.rough;

for i = 1:dim(1)
    % fill FWHMnames and slopenames
    FWHMnames  = [FWHMnames result{i,1}.paramNoise.FWHM(1)];
    slope      = result{i,1}.paramSignal.shapeparam
    if length(slope)==2
        slope = diff(slope)
    end
    slopenames = [slopenames slope]
    for j = 1:dim(2)
        for l = 1:length(lvls)
            switch type
                case 'truebdry'
                     covRates(i,j,l) = result{i,j}.covRate.truebdry(2,l);
                case 'linbdry'
                     covRates(i,j,l) = result{i,j}.covRate.linbdry(2,l);
                case 'erodbdry'
                     covRates(i,j,l) = result{i,j}.covRate.erodbdry(2,l);    
            end

        end
    end
end
% initialize counter on the covering levels
countl = 0;
for l = lvls
    % update couner on the covering levels
    countl = countl + 1;
    % open a new figure
    figure(countl), clf, hold on
    % plot the simulation results with the errorbars
    for i = 1:dim(1)
        plot( nVec, covRates(i,:,countl), [colorVec(i) '-' dotVec(i)], 'linewidth', 1.5)
    end
    % plot the nominal lvl
    plot([nVec(1)-5, nVec(end)+5], [lvls(countl) lvls(countl)], 'k', 'linewidth', 1.5)
    plot([nVec(1)-5, nVec(end)+5], [lvls(countl)-scStdErr(countl) lvls(countl)-scStdErr(countl)], 'k--')
    plot([nVec(1)-5, nVec(end)+5], [lvls(countl)+scStdErr(countl) lvls(countl)+scStdErr(countl)], 'k--')
    % set the range for the axis
    xlim([nVec(1)-5, nVec(end)+5])
    ylim([lvls(countl)-0.05 lvls(countl) + 0.05])
    % Add title
    title(titlename)
    % Add a legend
    legend( ['FWHM ', num2str(FWHMnames(1)), ' slope ', num2str(slopenames(1))],...
            ['FWHM ', num2str(FWHMnames(2)), ' slope ', num2str(slopenames(2))],...
            ['FWHM ', num2str(FWHMnames(3)), ' slope ', num2str(slopenames(3))],...
            ['FWHM ', num2str(FWHMnames(4)), ' slope ', num2str(slopenames(4))],...
            'nominal level', '1.96 x stdErr', 'Location', 'southeast' );
    set(gca, 'fontsize', 14);
    axis square;
    hold off
    saveas( gcf, strcat(location, '_level', num2str(100*lvls(countl)), '.png' ) )
end
end