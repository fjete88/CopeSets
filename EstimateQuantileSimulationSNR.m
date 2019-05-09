%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     This simulation estimates the maximum quantile of Z and T processes
%%%     with different bootstrap estimators
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% clear workspace
clear
close all

printBool    = 0;
trueSimBool  = 1;
estimSimBool = 0;

%%%%%% set path for the files
%cd  /home/drtea/Research/MatlabPackages/CopeSets

%% %%%%%  Obtain true quantile for different definitions of CoPe sets %%%%%
%%%  Simulation parameters
% number of simulations to approximate the maximum of the absolute value
msim = 5e4;
mreport = [1e4 2.5e4 4e4];
% Considered FWHM
FWHM = 3;
% Size of the field
Lvec = [10 60 124];
% number subjects
nsubj = [30 60 120];
% SNR
SNR = [0.2 0.7 2];
% quantile levels
lvls = [0.85 0.9 0.95];
%% % for loop generating the distribution of the maximum along the boundary of Cope sets
if(trueSimBool)
% Initialize array to save the results of the simulation
maxDistr = zeros([msim length(SNR) length(Lvec) length(nsubj) 6]);

tic
for f = FWHM
    for m=1:msim
        Y = SmoothField2D( max(nsubj), 1, [f f], [max(Lvec) max(Lvec)] );
        for c = SNR
            countc = find(c==SNR);
            % Get the data with mean c, meaning the SNR is also equal c
            Yc = Y + c;
            for n = nsubj
                countn   = find(n==nsubj);
                % constant which might cause numerical failure
                if n < 200
                    fac = sqrt((n-1)/2)*(gamma((n-2)/2)/gamma((n-1)/2));
                else
                    fac = 1 / ( 1 - (3/(4*(n-1)-1)) );
                end
                % Factor in the bias, i.e. E[^d] = biasfac*d
                biasfac  = sqrt((n-1)/n) * fac;

                % true variance at d=c derived from the variance of
                % non-central t
                trueStd  = sqrt( (n-1)^2/n/(n-3) + ( (n-1)^2/(n-3)-fac^2*(n-1) )*c^2 );
                
                % get constants for variance stabilisation
                a     = (n-1)/n/sqrt(n-3);
                b     = sqrt( (n-1)/fac^2/(n-3)-1 ) * (n-1)/n;
                alpha = 1/b/sqrt(n);
                beta  = b/a;
                
                
                % Get the SNR residuals
                [SNRresYcn, etaYcn, CohenStd] = SNR_residuals( Yc(:,:,1:n) );
                
                for L = Lvec
                    countL = find(L==Lvec);
                    % asymptotic variance of cohen's d
                    asymStd  = sqrt( 1 + c^2/2 );
                    % empirical variance of the SNR residuals
                    SNRStd   = std(SNRresYcn(1:L,1:L,:), 0, 3);
                    
                    % Different possibilities of true Cope set processes
                    % Note that we need to include the biasfac to make the
                    % process mean zero!
                    % normalizing by true known variance on the ^d=c if the
                    % data is Gaussian
                    maxDistr( m, countc, countL, countn, 1 ) = max(max(abs( sqrt(n)*(etaYcn(1:L,1:L,:) - c*biasfac) ./ trueStd  )));
                    % asymptotic variance on the contour line d=c
                    maxDistr( m, countc, countL, countn, 2 ) = max(max(abs( sqrt(n)*(etaYcn(1:L,1:L,:) - c*biasfac) ./ asymStd  )));
                    % estimated asymptotic variance on contour line d=c
                    maxDistr( m, countc, countL, countn, 3 ) = max(max(abs( sqrt(n)*(etaYcn(1:L,1:L,:) - c*biasfac) ./ CohenStd(1:L,1:L,:) )));
                    % Estimated empirical variance of the process from SNR
                    % residuals
                    maxDistr( m, countc, countL, countn, 4 ) = max(max(abs( sqrt(n)*(etaYcn(1:L,1:L,:) - c*biasfac) ./ SNRStd   )));
                    % variance stabilization with first order mean approximation
                    maxDistr( m, countc, countL, countn, 5 ) = max(max(abs( sqrt(n)*( alpha*asinh(beta*etaYcn(1:L,1:L,:)) - alpha*asinh(beta*c*biasfac) )   )));
                    % variance stabilization with second order mean approximation
                    maxDistr( m, countc, countL, countn, 6 ) = max(max(abs( sqrt(n)*( alpha*asinh(beta*etaYcn(1:L,1:L,:)) ...
                                                                                      - alpha*asinh(beta*c*biasfac) + b^2*c*biasfac*trueStd^(-1)/2 )  )));
                end
            end
        end
        if m == mreport(1) || m == mreport(2) || m == mreport(3)
            m/msim
        end
    end
end
toc

clear f m L Y Yc Ycn countc countL countn trueAsymVar EstimAsymVar ResVariance etaYcn stdYcn meanYcn SNRresYcn n
%save('simulations/maxDistr_SNRCopeSet_processes')
end

%% %%%%%%%%%%%%%%% Bootstrap quantile estimator simulations %%%%%%%%%%%%%%%
% Note that there is also an ssh script for parallel computing in the /scripts
% folder 
if(estimSimBool)
    
Mboot= 2.5e3;
msim = 1e3;
mreport = [250 500 750];

% Initialize buckets for the estimated quantiles
TrueVarMGauss   = zeros([msim length(lvls) length(SNR) length(Lvec) length(nsubj)]);
TrueVarMRadem   = TrueVarMGauss;
CohenVarMGauss  = TrueVarMGauss;
CohenVarMRadem  = TrueVarMGauss;
CohenVarMtGauss = TrueVarMGauss;
CohenVarMtRadem = TrueVarMGauss;
SNRVarMGauss    = TrueVarMGauss;
SNRVarMRadem    = TrueVarMGauss;
SNRVarMtGauss   = TrueVarMGauss;
SNRVarMtRadem   = TrueVarMGauss;

tic
for f = FWHM
    for m=1:msim
        Y = SmoothField2D( max(nsubj), 1, [f f], [max(Lvec) max(Lvec)] );
        for c = SNR
            countc = find(c==SNR);
            % Get the data with mean c, meaning the SNR is also equal c
            Yc = Y + c;
            for n = nsubj
                countn = find(n==nsubj);
                
                % Get the subsampled data to subject size
                Ycn = Yc(:,:,1:n);
                % Get the SNR residuals and the estimated Cohen's-d
                % variance
                [SNRresYcn, ~, CohenVar] = SNR_residuals( Ycn );

                % remove the standardization of SNR residuals
                trueVar  = sqrt( 1 + c^2/2 );
                SNRVar   = std(SNRresYcn, 0, 3);


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

                    CohenVarMtGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    CohenVarMtRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, mask, 'rademacher', 't' );

                    SNRVarMGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'gaussian', 'regular' );
                    SNRVarMRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL ./ SNRVarL, 2*lvls, Mboot, mask, 'rademacher', 'regular' );

                    SNRVarMtGauss( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, mask, 'gaussian', 't' );
                    SNRVarMtRadem( m, :, countc, countL, countn) = MultiplierBoots( SNRresYcnL, 2*lvls, Mboot, mask, 'rademacher', 't' );
                end
            end
        end
        if m == mreport(1) || m == mreport(2) || m == mreport(3)
            m/msim
        end
    end
end
toc

clear f m L Y Yc Ycn countc countL countn trueVar CohenVar SNRVar etaYcn stdYcn meanYcn SNRresYcn n c
save('simulations/estimQuantile_SNRCopeSet_processes')

end
%%%%%% Get the summarized simulation results and plot them
%% true Quantiles from maximum simulation
if(printBool)

load('simulations/maxDistr_SNRCopeSet_processes')

trueQuantTrue  = zeros([length(lvls) length(SNR) length(Lvec) length(nsubj)]);
trueQuantCohen = trueQuantTrue;
trueQuantSNR   = trueQuantTrue;

for f = FWHM
    for c = SNR
        countc = find(c==SNR);
        for n = nsubj
            countn = find(n==nsubj);
            for L = Lvec
                countL = find(L==Lvec);
                % Get true quantiles
                trueQuantTrue( :, countc, countL, countn)  = quantile( maxDistr( :, countc, countL, countn, 1 ), lvls );
                trueQuantCohen( :, countc, countL, countn) = quantile( maxDistr( :, countc, countL, countn, 2 ), lvls );
                trueQuantSNR( :, countc, countL, countn)   = quantile( maxDistr( :, countc, countL, countn, 3 ), lvls );
            end
        end
    end
end


%% Estimated Quantiles
load('simulations/estimQuantile_SNRCopeSet_processes')

QuantCohenMGauss  = zeros([length(lvls) length(SNR) length(Lvec) length(Nsubj)]);
QuantCohenMRadem  = QuantCohenMGauss;
QuantCohenMtGauss = QuantCohenMGauss;
QuantCohenMtRadem = QuantCohenMGauss;

QuantSNRMGauss  = QuantCohenMGauss;
QuantSNRMRadem  = QuantCohenMGauss;
QuantSNRMtGauss = QuantCohenMGauss;
QuantSNRMtRadem = QuantCohenMGauss;

QuantTrueMGauss  = QuantCohenMGauss;
QuantTrueMRadem  = QuantCohenMGauss;

for f = FWHM
    for c = SNR
        countc = find(c==SNR);
        % Get the data with mean c, meaning the SNR is also equal c
        for n = Nsubj
            countn = find(n==Nsubj);
            for L = Lvec
                countL = find(L==Lvec);
                % Get true quantiles
                QuantTrueMGauss( :, countc, countL, countn)   = mean(TrueVarMGauss( :, :, countc, countL, countn));
                QuantTrueMRadem( :, countc, countL, countn)   = mean(TrueVarMRadem( :, :, countc, countL, countn));

                QuantSNRMGauss( :, countc, countL, countn)    = mean(SNRVarMGauss( :, :, countc, countL, countn));
                QuantSNRMRadem( :, countc, countL, countn)    = mean(SNRVarMRadem( :, :, countc, countL, countn));
                QuantSNRMtGauss( :, countc, countL, countn)   = mean(SNRVarMtGauss( :, :, countc, countL, countn));
                QuantSNRMtRadem( :, countc, countL, countn)   = mean(SNRVarMtRadem( :, :, countc, countL, countn));       
                
                QuantCohenMGauss( :, countc, countL, countn)  = mean(CohenVarMGauss( :, :, countc, countL, countn));
                QuantCohenMRadem( :, countc, countL, countn)  = mean(CohenVarMRadem( :, :, countc, countL, countn));
                QuantCohenMtGauss( :, countc, countL, countn) = mean(CohenVarMtGauss( :, :, countc, countL, countn));
                QuantCohenMtRadem( :, countc, countL, countn) = mean(CohenVarMtRadem( :, :, countc, countL, countn));  
            end
        end
    end
end

%% Plot the results
Vibrant    = [[0 119 187];... % blue
              [51 187 238];...% cyan
              [0 153 136];... % teal
              [238 119 51];...% orange
              [204 51 17];... % red
              [238 51 119];...% magenta
              [187 187 187]...% grey
              ]/255;

colMat = Vibrant([1 3 4 5],:);
path_pics = '/home/drtea/Research/MatlabPackages/CopeSets/pics';

% Global figure settings
sfont = 20;
addf  = 5;
scale = 10/12;
WidthFig   = 800;
HeightFig  = WidthFig * scale;
xvec       = [10 20 30]; %nvec;
xtickcell  = {'10', '60', '124'};
%yvec1      = [12 13 14 15 16];
%yvec2      = [40 45 50 55];
%ytickcell1 = {'12' '13' '14' '15' '16'};
%ytickcell2 = {'40' '45' '50' '55'};

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

for a = 1:length(lvls)
for n = Nsubj
    countn = find(n==Nsubj);
    for c = SNR
        countc = find(c==SNR);
        
        figure('pos',[10 10 WidthFig HeightFig]), clf, hold on
        set(gca, 'fontsize', sfont);
        plot([10 20 30], squeeze(trueQuantSNR(a,countc, :, countn)), '-k', 'LineWidth', 2 )
        plot([10 20 30], squeeze(QuantSNRMGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(1,:) )
        plot([10 20 30], squeeze(QuantSNRMtGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(2,:) )
        plot([10 20 30], squeeze(QuantSNRMRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(3,:) )
        plot([10 20 30], squeeze(QuantSNRMtRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(4,:) )
        

        % Modify gloabal font size for this plot
        set(gca,'FontSize', sfont)

        % Change axis style
        xlim([xvec(1)-5 xvec(end)+5])
        xticks(xvec)
        xticklabels(xtickcell)

        h = xlabel('Size of image [LxL]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        h = ylabel('quantile value', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        
        % add legend
        legend( 'True Quantile', 'gMult', 'gMult-t', 'rMult', 'rMult-t',...
                                'Location', 'southeast' );
        set(legend, 'fontsize', sfont);
        legend boxoff
        
        h=title(strcat('SNR=',num2str(c), ' Nsubj=',num2str(n), ' Quant= ',num2str(lvls(a)) ));
        set(h, 'Interpreter', 'latex');
        set(h, 'fontsize', sfont+addf);

        saveas( gcf, [path_pics,'/ResultsQuantileSimulation_VarSNR_SNR',num2str(c), '_Nsubj',num2str(n), '_Quant',num2str(100*lvls(a)),'.png'] )

        figure('pos',[10 10 WidthFig HeightFig]), clf, hold on
        set(gca, 'fontsize', sfont);
        plot([10 20 30], squeeze(trueQuantCohen(a,countc, :, countn)), '-k', 'LineWidth', 2 )
        plot([10 20 30], squeeze(QuantCohenMGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(1,:) )
        plot([10 20 30], squeeze(QuantCohenMtGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(2,:) )
        plot([10 20 30], squeeze(QuantCohenMRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(3,:) )
        plot([10 20 30], squeeze(QuantCohenMtRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(4,:) )
        

        % Modify gloabal font size for this plot
        set(gca,'FontSize', sfont)

        % Change axis style
        xlim([xvec(1)-5 xvec(end)+5])
        xticks(xvec)
        xticklabels(xtickcell)

        h = xlabel('Size of image [LxL]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        h = ylabel('quantile value', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        
        % add legend
        legend( 'True Quantile', 'gMult', 'gMult-t', 'rMult', 'rMult-t',...
                                'Location', 'southeast' );
        set(legend, 'fontsize', sfont);
        legend boxoff
        
        h=title(strcat('SNR=',num2str(c), ' Nsubj=',num2str(n), ' Quant= ',num2str(lvls(a)) ));
        set(h, 'Interpreter', 'latex');
        set(h, 'fontsize', sfont+addf);

        saveas( gcf, [path_pics,'/ResultsQuantileSimulation_VarCohen_SNR',num2str(c), '_Nsubj',num2str(n), '_Quant',num2str(100*lvls(a)),'.png'] )

        figure('pos',[10 10 WidthFig HeightFig]), clf, hold on
        set(gca, 'fontsize', sfont);
        plot([10 20 30], squeeze(trueQuantTrue(a,countc, :, countn)), '-k', 'LineWidth', 2 )
        plot([10 20 30], squeeze(QuantTrueMGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(1,:) )
%        plot([10 20 30], squeeze(QuantCohenMtGauss(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(2,:) )
        plot([10 20 30], squeeze(QuantTrueMRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(3,:) )
%        plot([10 20 30], squeeze(QuantCohenMtRadem(a,countc, :, countn)), '-', 'LineWidth', 1.5,'Color',colMat(4,:) )
        

        % Modify gloabal font size for this plot
        set(gca,'FontSize', sfont)

        % Change axis style
        xlim([xvec(1)-5 xvec(end)+5])
        xticks(xvec)
        xticklabels(xtickcell)

        h = xlabel('Size of image [LxL]', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        h = ylabel('quantile value', 'fontsize', sfont+addf); set(h, 'Interpreter', 'latex');
        
        % add legend
        legend( 'True Quantile', 'gMult', 'gMult-t', 'rMult', 'rMult-t',...
                                'Location', 'southeast' );
        set(legend, 'fontsize', sfont);
        legend boxoff
        
        h=title(strcat('SNR=',num2str(c), ' Nsubj=',num2str(n), ' Quant= ',num2str(lvls(a)) ));
        set(h, 'Interpreter', 'latex');
        set(h, 'fontsize', sfont+addf);

        saveas( gcf, [path_pics,'/ResultsQuantileSimulation_VarTrueAsympt_SNR',num2str(c), '_Nsubj',num2str(n), '_Quant',num2str(100*lvls(a)),'.png'] )
        close all
    end
end
end
end