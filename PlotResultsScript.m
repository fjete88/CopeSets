%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%         This script plots the results from CopeSets_SimSkript.m
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

cd /home/drtea/Research/MatlabPackages/CopeSets
mkdir pics

simpathvector = ["simulations/ResultSim50050Nsim3000maxN240isotropic_boot_regulargaussian.mat",...
                 "simulations/ResultSim50050Nsim3000maxN240isotropic_boot_tgaussian.mat",...
                 "simulations/ResultSim50050Nsim3000maxN240isotropic_boot_regularrademacher.mat",...
                 "simulations/ResultSim10050Nsim3000maxN300isotropic_boot_regulargaussian.mat",...
                 "simulations/ResultSim10050Nsim3000maxN300isotropic_boot_tgaussian.mat",...
                 "simulations/ResultSim5050Nsim3000maxN300isotropic_boot_regularrademacher.mat",...
                 "simulations/ResultSim5050Nsim3000maxN300isotropic_boot_trademacher.mat",...
                 "simulations/ResultSim5050Nsim3000maxN300isotropic_boot_regulargaussian.mat",...
                 "simulations/ResultSim5050Nsim3000maxN300isotropic_boot_tgaussian.mat"];


%% Visualize 500x50 simulations Rademacher regular
for simpath = simpathvector
    % Plot regular results
    load(simpath)

    dim     = simResults{1}.paramNoise.dim;
    method  = simResults{1}.quantEstim.params.method;
    weights = simResults{1}.quantEstim.params.weights;

    % remove spaces in string
    tmp = num2str(dim);
    tmp = tmp(find(~isspace(tmp)));
    % define path identifier for plots
    picpathID = strcat( 'pics/', tmp, method, weights );

    if strcmp( method, 'regular' )
        methodflag = ' Reg. ';
    else
        methodflag = ' t ';
    end    
    if strcmp( weights, 'rademacher' )
        weightsflag = 'Rad.';
    else
        weightsflag = 'Gauss';
    end

    titleSnippet = strcat(tmp,methodflag, weightsflag);
    clear weightsflag methodflag tmp;

    %%% plot results for regular bootstrap
    % cov rates
    PlotResults( simResults, strcat(picpathID, 'linResult'), 'linbdry', ...
                 strcat(titleSnippet,'/linbdry'), 'y', 'variable')
    PlotResults( simResults, strcat(picpathID, 'trueResult')', 'truebdry', ...
                strcat(titleSnippet,'/truedbdry'), 'y', 'variable')
    %PlotResults( sim5050regResult, 'pics/5050regGausserodResult', 'erodbdry', 'CovRate: Reg. boots/erodbdry Estim', 'static')
    % quantile estimates
    PlotQuantiles( simResults, strcat(picpathID, 'meanQuantiles'), titleSnippet)

    close all
end
