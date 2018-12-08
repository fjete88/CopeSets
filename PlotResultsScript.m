%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%         This script plots the results from CopeSets_SimSkript.m
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

cd /home/drtea/Research/MatlabPackages/CopeSets
mkdir pics

simpathvector = ["ResultSim50050Nsim3000maxN240isotropic_boot_regulargaussian.mat",...
                 "ResultSim50050Nsim3000maxN240isotropic_boot_tgaussian.mat",...
                 "ResultSim50050Nsim3000maxN240isotropic_boot_regularrademacher.mat",...
                 "ResultSim10050Nsim3000maxN300isotropic_boot_regulargaussian.mat",...
                 "ResultSim10050Nsim3000maxN300isotropic_boot_tgaussian.mat",...
                 "ResultSim5050Nsim3000maxN300isotropic_boot_regularrademacher.mat",...
                 "ResultSim5050Nsim3000maxN300isotropic_boot_trademacher.mat",...
                 "ResultSim5050Nsim3000maxN300isotropic_boot_regulargaussian.mat",...
                 "ResultSim5050Nsim3000maxN300isotropic_boot_tgaussian.mat"];
             
simpathvector = "ResultSim5050Nsim3000maxN300isotropicFWM35_boot_trademacher.mat"
simpathvector = "ResultSim5050Nsim3000maxN300isotropicFWM3_boot_trademacher.mat"
simpathvector = "ResultSim5050Nsim3000maxN300isotropicFWM10_slopes1_4_6_8__boot_trademacher.mat"

%% Visualize simulation results
% Loop over the simulation results
for simpath = simpathvector
    % load the simulation results
    load(strcat('simulations/',simpath))

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
