function[] = PlotQuantiles( result, location, titlename, colorVec, dotVec)
% Plots average estimate of the bootstrap quantiles.
% Input:
%  result:
%  location:
%  titlename:
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
    case 3
        colorVec = 'rbmc';
        dotVec   = 'xxoo';
    case 4
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
a_true = zeros( [ dim ,length(lvls)] );
a_lin  = a_true; 

for i = 1:dim(1)
    % fill FWHMnames and slopenames
    FWHMnames  = [FWHMnames result{i,1}.paramNoise.FWHM(1)];
    slope      = result{i,1}.paramSignal.shapeparam;
    if length(slope)==2
        slope = diff(slope);
    end
    slopenames = [slopenames slope];
    for j = 1:dim(2)
         a_true(i,j,:) = mean( result{i,j}.quant.truebdry,2 );
         a_lin(i,j,:)  = mean( result{i,j}.quant.linbdry,2  );
    end
end

dynamiclegend = cell([1 dim(1)+2]);
for j = 1:dim(1)
    dynamiclegend{j} = ['FWHM ', num2str(FWHMnames(j)), ' slope ', num2str(slopenames(j))];
end
dynamiclegend{dim(1)+1} = 'Lin. Bdry';
dynamiclegend{dim(1)+2} = 'True.Bdry';

% initialize counter on the covering levels
countl = 0;
% open new figure
figure('pos',[10 10 2*900 2*600]), clf, hold on
for l = lvls
    % update couner on the covering levels
    countl = countl + 1;
    h = zeros([1 2*dim(1)+2]);
    % open a new subfigure
    subplot(1, length(lvls), countl),  hold on
    % plot the simulation results of the mean quantile estimate for lin
    counti = 0;
    for i = 1:dim(1)
        h(counti+1)  = plot( nVec, a_lin(i,:,countl), [colorVec(i) '-' dotVec(i)], 'linewidth', 1.5);
        h(counti+2)  = plot( nVec, a_true(i,:,countl), [colorVec(i) '--' dotVec(i)], 'linewidth', 1.5);      
        counti = counti+2;
    end
    % fake plots for legend
    h(counti+1) = plot([nVec(1)-5, nVec(end)+5], [-100 -100], 'k', 'linewidth', 1.5);
    h(counti+2) = plot([nVec(1)-5, nVec(end)+5], [-100 -100], 'k--', 'linewidth', 1.5);
    % plot([nVec(1)-5, nVec(end)+5], [lvls(countl)+scStdErr(countl) lvls(countl)+scStdErr(countl)], 'k--')
    
    % set the range for the axis
    xlim([nVec(1)-5, nVec(end)+5])
    tmp = [a_true(:,:,countl) a_lin(:,:,countl)];
    ylim( [0.9*min(tmp(:)) 1.1*max(tmp(:))] );
    % specify tiks
    xticks( nVec )
    % put label onto the axis
    xlabel('Sample Size [N]');
    ylabel(strcat(num2str(lvls(countl)),' - quantile value'))
    % Add title
    title(titlename)
    % Add a legend
    if(countl==1)
        legend( h([2*(1:dim(1))-1 2*dim(1)+1 2*dim(1)+2]), dynamiclegend{:}, 'Location', 'northeast' );
    end
    set(gca, 'fontsize', 14);
    axis square;
    hold off
end
hold off
saveas( gcf, strcat(location, '.png' ) )
end