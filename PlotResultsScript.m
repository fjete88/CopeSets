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
             
simpathvector = "simulations/ResultSim5050Nsim3000maxN300isotropicFWHM35_boot_trademacher.mat"

%% Visualize simulation results
% Loop over the simulation results
for simpath = simpathvector(1)
    % load the simulation results
    load(simpath)

    % get the parameters the simulation used
    dim     = simResults{1}.paramNoise.dim;
    method  = simResults{1}.quantEstim.params.method;
    weights = simResults{1}.quantEstim.params.weights;
    FWHM    = [simResults{1}.paramNoise.FWHM(1) simResults{2}.paramNoise.FWHM(1)];

    % remove spaces in string for the dimension
    tmp = num2str(dim);
    tmp = tmp(find(~isspace(tmp)));

    % define path identifier for plots
    picpathID = strcat( 'pics/', tmp, method, weights,'FWHM', sprintf('%i', FWHM));

    % get shortcuts for the method and weights for the title and labels
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
    
    % Create a title string for the plot
    titleSnippet = strcat(tmp,methodflag, weightsflag);
    clear weightsflag methodflag tmp;

    %%% plot results for the simulation
    % plot cov rates linear bdry estimation 
    PlotResults( simResults, strcat(picpathID, 'linResult'), 'linbdry', ...
                 strcat(titleSnippet,'/linbdry'), 'y', 'variable')
    % plot cov rates true bdry
    PlotResults( simResults, strcat(picpathID, 'trueResult')', 'truebdry', ...
                strcat(titleSnippet,'/truedbdry'), 'y', 'variable')
    % plot cov rates erod/dilation bdry estimation
    PlotResults( simResults, strcat(picpathID, 'erodResult')', 'erodbdry', ...
                strcat(titleSnippet,'/erodbdry'), 'y', 'variable')
    % plot mean quantile estimates
    PlotQuantiles( simResults, strcat(picpathID, 'meanQuantiles'), titleSnippet)

    close all
end
