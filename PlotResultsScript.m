clear all
close all

cd /home/drtea/Research/MatlabPackages/CopeSets
%% Visualize 50x50 simulations
load('simulations/ResultSim5050regNsim3000maxN240_isotropic_FWHM.mat')
load('simulations/ResultSim5050tNsim3000maxN240_isotropic_FWHM.mat')

PlotResults( sim5050regResult, 'pics/5050reglinResult', 'linbdry', 'CovRate: Reg. boots/linbdry Estim')
PlotResults( sim5050regResult, 'pics/5050regerodResult', 'erodbdry', 'CovRate: Reg. boots/erodbdry Estim')
PlotResults( sim5050regResult, 'pics/5050regtrueResult', 'truebdry', 'CovRate: Reg. boots/truedbdry Estim')

PlotResults( sim5050tResult, 'pics/5050tlinResult', 'linbdry', 'CovRate: t-boots/linbdry Estim')
PlotResults( sim5050tResult, 'pics/5050terodResult', 'erodbdry', 'CovRate: t-boots/erodbdry Estim')
PlotResults( sim5050tResult, 'pics/5050ttrueResult', 'truebdry', 'CovRate: t-boots/truedbdry Estim')

%% Visualize 100x50 simulations
load('simulations/ResultSim10050regNsim3000maxN240_isotropic_FWHM.mat')
load('simulations/ResultSim10050tNsim3000maxN240_isotropic_FWHM.mat')

PlotResults( sim10050regResults, 'pics/10050reglinResult', 'linbdry', 'CovRate: Reg. boots/linbdry Estim')
PlotResults( sim10050regResults, 'pics/10050regerodResult', 'erodbdry', 'CovRate: Reg. boots/erodbdry Estim')
PlotResults( sim10050regResults, 'pics/10050regtrueResult', 'truebdry', 'CovRate: Reg. boots/truedbdry Estim')

PlotResults( sim10050tResults, 'pics/10050tlinResult', 'linbdry', 'CovRate: t-boots/linbdry Estim')
PlotResults( sim10050tResults, 'pics/10050terodResult', 'erodbdry', 'CovRate: t-boots/erodbdry Estim')
PlotResults( sim10050tResults, 'pics/10050ttrueResult', 'truebdry', 'CovRate: t-boots/truedbdry Estim')

%% Visualize 100x50 simulations
%load('simulations/ResultSim10050regNsim3000maxN240_isotropic_FWHM.mat')
load('simulations/ResultSim10050tNsim3000maxN240_isotropic_FWHM.mat')

%PlotResults( sim10050regResults, 'pics/10050reglinResult', 'linbdry', 'CovRate: Reg. boots/linbdry Estim')
%PlotResults( sim10050regResults, 'pics/10050regerodResult', 'erodbdry', 'CovRate: Reg. boots/erodbdry Estim')
%PlotResults( sim10050regResults, 'pics/10050regtrueResult', 'truebdry', 'CovRate: Reg. boots/truedbdry Estim')

PlotResults( sim10050tResults, 'pics/10050tlinResult', 'linbdry', 'CovRate: t-boots/linbdry Estim')
PlotResults( sim10050tResults, 'pics/10050terodResult', 'erodbdry', 'CovRate: t-boots/erodbdry Estim')
PlotResults( sim10050tResults, 'pics/10050ttrueResult', 'truebdry', 'CovRate: t-boots/truedbdry Estim')
