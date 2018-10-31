clear all
close all

load('simulations/ResultSim5050regNsim3000maxN240_isotropic_FWHM.mat')
load('simulations/ResultSim5050tNsim3000maxN240_isotropic_FWHM.mat')

PlotResults( sim5050regResult, 'pics/5050reglinResult', 'linbdry', 'CovRate: Reg. boots/linbdry Estim')
PlotResults( sim5050regResult, 'pics/5050regerodResult', 'erodbdry', 'CovRate: Reg. boots/erodbdry Estim')
PlotResults( sim5050regResult, 'pics/5050regtrueResult', 'truebdry', 'CovRate: Reg. boots/truedbdry Estim')

PlotResults( sim5050tResult, 'pics/5050tlinResult', 'linbdry', 'CovRate: t-boots/linbdry Estim')
PlotResults( sim5050tResult, 'pics/5050tregResult', 'erodbdry', 'CovRate: t-boots/erodbdry Estim')
PlotResults( sim5050tResult, 'pics/5050terodResult', 'truebdry', 'CovRate: t-boots/truedbdry Estim')

