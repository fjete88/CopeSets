close all
clear all
cd /home/drtea/Research/MatlabPackages/CopeSets

lvls = [0.85 0.9 0.95];
res  = 0;
% load quantile simulation results
load(strcat('simulations/ResultBootstrappedQuantilesRes',num2str(res),'.mat'))
% mean estimated quantiles
mnonpareg = squeeze(mean(nonpareg,3));
mrademreg = squeeze(mean(rademreg,3));
mgaussreg = squeeze(mean(gaussreg,3));
mnonpat  = squeeze(mean(nonpat,3));
mrademt  = squeeze(mean(rademt,3));
mgausst  = squeeze(mean(gausst,3));
stdnonpareg = squeeze(std(nonpareg,0,3));
stdrademreg = squeeze(std(rademreg,0,3));
stdgaussreg = squeeze(std(gaussreg,0,3));
stdnonpat  = squeeze(std(nonpat,0,3));
stdrademt  = squeeze(std(rademt,0,3));
stdgausst  = squeeze(std(gausst,0,3));

%% Get the true quantiles from the simulations
% Set font to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(gca, 'fontsize', 20);

for d=1:3
    load(strcat('simulations/ResultMaxDistSim50000maxDF299isotropicFWHM510Bdry',...
                 num2str(dim(d,1)), '.mat'))
    for f=1:2             
        simMax    = simResults{1,f}.maxDist;
        trueQuant = [quantile( simMax(2,:), lvls )', quantile( simMax(3,:), lvls )',...
                     quantile( simMax(4,:), lvls )', quantile( simMax(5,:), lvls )',...
                     quantile( simMax(6,:), lvls )', quantile( simMax(1,:), lvls )'];
        figure('pos',[10 10 2*900 2*600]); clf;
        for l=1:length(lvls)
            subplot(2,3,l); hold on;
            plot( nvec-1, trueQuant(l,1:5), '-k', 'LineWidth', 1.5)
            plot( [0 400], [trueQuant(l,6) trueQuant(l,6)], '--k', 'LineWidth', 1.5 );
            plot( nvec-1, mnonpat(l,:,f,d), '-g', 'LineWidth', 1.2 )
            plot( nvec-1, mrademt(l,:,f,d), '-r', 'LineWidth', 1.2 )
            plot( nvec-1, mgausst(l,:,f,d), '-b', 'LineWidth', 1.2 )
            plot( nvec-1, mnonpareg(l,:,f,d), '--g', 'LineWidth', 1.2 )
            plot( nvec-1, mrademreg(l,:,f,d), '--r', 'LineWidth', 1.2 )
            plot( nvec-1, mgaussreg(l,:,f,d), '--b', 'LineWidth', 1.2 )

            % specify tiks
            xticks( nvec-1 )
            % put label onto the axis
            xlabel('Degrees of freedom [df]');
            ylabel('Quantile Value')
            set(gca, 'fontsize', 14);
            if(l==2)
              % Add title
              if(res==1)
                title(strcat('Residuals BdryLength ', num2str(dim(d,1)), 'FWHM',...
                              num2str(FWHM(f)) ))
              else
                title(strcat('Direct BdryLength ', num2str(dim(d,1)), 'FWHM',...
                              num2str(FWHM(f)) ))
              end
                legend( 't-quantile MC', 'Gauss-quantile MC','NonPa-t','Radem-t', 'Gauss-t','NonPa-reg','Radem-reg',...
                       'Gauss-reg', 'Location', 'northeast' )
            end
        end
        for l=1:length(lvls)
            subplot(2,3,l+3); hold on;
            plot( nvec-1, stdnonpat(l,:,f,d), '-g', 'LineWidth', 1.2 )
            plot( nvec-1, stdrademt(l,:,f,d), '-r', 'LineWidth', 1.2 )
            plot( nvec-1, stdgausst(l,:,f,d), '-b', 'LineWidth', 1.2 )
            plot( nvec-1, stdnonpareg(l,:,f,d), '--g', 'LineWidth', 1.2 )
            plot( nvec-1, stdrademreg(l,:,f,d), '--r', 'LineWidth', 1.2 )
            plot( nvec-1, stdgaussreg(l,:,f,d), '--b', 'LineWidth', 1.2 )


            % specify tiks
            xticks( nvec-1 )
            % put label onto the axis
            xlabel('Degrees of freedom [df]');
            ylabel('Std of Quantile Value')
            set(gca, 'fontsize', 14);
         end

        hold off
        saveas( gcf, strcat('pics/simResultsQuantileFWHM', num2str(FWHM(f)),'Bdry', num2str(dim(d,1)),'Res', num2str(res), '.png' ) )
    end
end
close all