function[] = PlotResults( result, location, type, titlename, linbdryCov, ylimType, stdfac, colorVec, dotVec)
% Plots results from a simulation and outputs figures for the different lvls.
% Input:
%  result:
%  location:
%  type:
%  titlename:
%  linbdryCov:
%  ylimType: 'static', 'variable'
%  stdfac:
%  colorVec:
%  dotVec:
% Output:
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 11/03/2018
%__________________________________________________________________________
switch nargin
    case 4
        linbdryCov = 'y';
        ylimType   = [0.05, 0.05];
        stdfac     = 1.96;
        colorVec   = 'rbmc';
        dotVec     = 'xxoo';
    case 5
        ylimType   = 'static';
        stdfac     = 1.96;
        colorVec   = 'rbmc';
        dotVec     = 'xxoo';
    case 6
        stdfac     = 1.96;
        colorVec   = 'rbmc';
        dotVec     = 'xxoo';
    case 7
        colorVec   = 'rbmc';
        dotVec     = 'xxoo';
    case 8
        dotVec   = 'xxoo';
end

% Set font to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(gca, 'fontsize', 20);
% Get constants from the results cell structure
dim    = size(result);
result = reshape(result, [prod(dim(1:end-1)), dim(end)]);
dim    = size(result);

% Get the index for the correct covering method. Old (without checking on boundary 'n', new checking along boundary 'y')
if linbdryCov == 'y'
    CovIndex = 2;
else
    CovIndex = 1;
end

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
scStdErr = stdfac*result{1}.stdErr.rough;

for i = 1:dim(1)
    % fill FWHMnames and slopenames
    FWHMnames  = [FWHMnames result{i,1}.paramNoise.FWHM(1)];
    slope      = result{i,1}.paramSignal.shapeparam;
    if length(slope)==2
        slope = diff(slope);
    end
    slopenames = [slopenames slope];
    for j = 1:dim(2)
        for l = 1:length(lvls)
            switch type
                case 'truebdry'
                     covRates(i,j,l) = result{i,j}.covRate.truebdry(CovIndex,l);
                case 'linbdry'
                     covRates(i,j,l) = result{i,j}.covRate.linbdry(CovIndex,l);
                case 'erodbdry'
                     covRates(i,j,l) = result{i,j}.covRate.erodbdry(CovIndex,l);    
            end
        end
    end
end

dynamiclegend = cell([1 dim(1)+2]);
for j = 1:dim(1)
    dynamiclegend{j} = ['FWHM ', num2str(FWHMnames(j)), ' slope ', num2str(slopenames(j))];
end
dynamiclegend{dim(1)+1} = 'nominal level';
dynamiclegend{dim(1)+2} = '1.96 x stdErr';

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
    if strcmp(ylimType, 'variable')
        tmp = covRates(:,:,countl);
        ylim( [max(0,0.95*min(tmp(:))) min(1,1.05*max(tmp(:)))] );
    else
        ylim([lvls(countl)-ylimType(1) lvls(countl) + ylimType(2)]);
    end
    % specify tiks
    xticks( nVec )
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel('emp. Covering Rate')
    % Add title
    title(titlename)
    % Add a legend
    legend( dynamiclegend{:}, 'Location', 'southeast' )

    set(gca, 'fontsize', 14);
    axis square;
    hold off
    saveas( gcf, strcat(location, '_level', num2str(100*lvls(countl)), '.png' ) )
end
end