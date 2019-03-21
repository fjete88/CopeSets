%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Script concatinating the simulations from the server
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..

parallel_count = 16;

for k = 1:16
    load( strcat('simulations/estimQuantile_SNRCopeSet_processes', num2str(k), '.mat'));
    if k==1
        TrueVarMGauss_tmp   = TrueVarMGauss;
        TrueVarMRadem_tmp   = TrueVarMRadem;
        CohenVarMGauss_tmp  = CohenVarMGauss;
        CohenVarMRadem_tmp  = CohenVarMRadem;
        CohenVarMtGauss_tmp = CohenVarMtGauss;
        CohenVarMtRadem_tmp = CohenVarMtRadem;
        SNRVarMGauss_tmp    = SNRVarMGauss;
        SNRVarMRadem_tmp    = SNRVarMRadem;
        SNRVarMtGauss_tmp   = SNRVarMtGauss;
        SNRVarMtRadem_tmp   = SNRVarMtRadem;        
    else
        TrueVarMGauss_tmp   = [ TrueVarMGauss_tmp; TrueVarMGauss ];
        TrueVarMRadem_tmp   = [ TrueVarMRadem_tmp; TrueVarMRadem ];
        CohenVarMGauss_tmp  = [ CohenVarMGauss_tmp; CohenVarMGauss ];
        CohenVarMRadem_tmp  = [ CohenVarMRadem_tmp; CohenVarMRadem ];
        CohenVarMtGauss_tmp = [ CohenVarMtGauss_tmp; CohenVarMtGauss ];
        CohenVarMtRadem_tmp = [ CohenVarMtRadem_tmp; CohenVarMtRadem ];
        SNRVarMGauss_tmp    = [ SNRVarMGauss_tmp; SNRVarMGauss ];
        SNRVarMRadem_tmp    = [ SNRVarMRadem_tmp; SNRVarMRadem ];
        SNRVarMtGauss_tmp   = [ SNRVarMtGauss_tmp; SNRVarMtGauss ];
        SNRVarMtRadem_tmp   = [ SNRVarMtRadem_tmp; SNRVarMtRadem ];        
    end
end

TrueVarMGauss   = TrueVarMGauss_tmp;
TrueVarMRadem   = TrueVarMRadem_tmp;
CohenVarMGauss  = CohenVarMGauss_tmp;
CohenVarMRadem  = CohenVarMRadem_tmp;
CohenVarMtGauss = CohenVarMtGauss_tmp;
CohenVarMtRadem = CohenVarMtRadem_tmp;
SNRVarMGauss    = SNRVarMGauss_tmp;
SNRVarMRadem    = SNRVarMRadem_tmp;
SNRVarMtGauss   = SNRVarMtGauss_tmp;
SNRVarMtRadem   = SNRVarMtRadem_tmp;       

% Total number of simulations
msim = size(SNRVarMtRadem,1);

save('simulations/estimQuantile_SNRCopeSet_processes.mat', 'TrueVarMGauss', 'TrueVarMRadem', 'CohenVarMGauss', 'CohenVarMRadem', 'CohenVarMtGauss',...
    'CohenVarMtRadem', 'SNRVarMGauss', 'SNRVarMRadem', 'SNRVarMtGauss', 'SNRVarMtRadem',...
    'msim', 'Nsubj', 'Lvec', 'FWHM', 'SNR', 'lvls', 'Mboot')