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
        TrueStdMGauss_tmp  = TrueStdMGauss;
        TrueStdMRadem_tmp  = TrueStdMRadem;
        TrueStdMtGauss_tmp = TrueStdMtGauss;
        TrueStdMtRadem_tmp = TrueStdMtRadem;

        AsymStdMGauss_tmp  = AsymStdMGauss;
        AsymStdMRadem_tmp  = AsymStdMRadem;
        AsymStdMtGauss_tmp = AsymStdMtGauss;
        AsymStdMtRadem_tmp = AsymStdMtRadem;        
        
        CohenStdMGauss_tmp  = CohenStdMGauss;
        CohenStdMRadem_tmp  = CohenStdMRadem;
        CohenStdMtGauss_tmp = CohenStdMtGauss;
        CohenStdMtRadem_tmp = CohenStdMtRadem;
        
        SNRStdMGauss_tmp  = SNRStdMGauss;
        SNRStdMRadem_tmp  = SNRStdMRadem;
        SNRStdMtGauss_tmp = SNRStdMtGauss;
        SNRStdMtRadem_tmp = SNRStdMtRadem;
        

        StabStdMGauss_tmp  = StabStdMGauss;
        StabStdMRadem_tmp  = StabStdMRadem;
        StabStdMtGauss_tmp = StabStdMtGauss;
        StabStdMtRadem_tmp = StabStdMtRadem;   
    else
        TrueStdMGauss_tmp   = [ TrueStdMGauss_tmp; TrueStdMGauss ];
        TrueStdMRadem_tmp   = [ TrueStdMRadem_tmp; TrueStdMRadem ];
        TrueStdMtGauss_tmp  = [ TrueStdMtGauss_tmp; TrueStdMtGauss ];
        TrueStdMtRadem_tmp  = [ TrueStdMtRadem_tmp; TrueStdMtRadem ];
        
        AsymStdMGauss_tmp   = [ AsymStdMGauss_tmp; AsymStdMGauss ];
        AsymStdMRadem_tmp   = [ AsymStdMRadem_tmp; AsymStdMRadem ];
        AsymStdMtGauss_tmp  = [ AsymStdMtGauss_tmp; AsymStdMtGauss ];
        AsymStdMtRadem_tmp  = [ AsymStdMtRadem_tmp; AsymStdMtRadem ];
        
        CohenStdMGauss_tmp  = [ CohenStdMGauss_tmp; CohenStdMGauss ];
        CohenStdMRadem_tmp  = [ CohenStdMRadem_tmp; CohenStdMRadem ];
        CohenStdMtGauss_tmp = [ CohenStdMtGauss_tmp; CohenStdMtGauss ];
        CohenStdMtRadem_tmp = [ CohenStdMtRadem_tmp; CohenStdMtRadem ];
        
        SNRStdMGauss_tmp    = [ SNRStdMGauss_tmp; SNRStdMGauss ];
        SNRStdMRadem_tmp    = [ SNRStdMRadem_tmp; SNRStdMRadem ];
        SNRStdMtGauss_tmp   = [ SNRStdMtGauss_tmp; SNRStdMtGauss ];
        SNRStdMtRadem_tmp   = [ SNRStdMtRadem_tmp; SNRStdMtRadem ];
        
        StabStdMGauss_tmp    = [ StabStdMGauss_tmp; StabStdMGauss ];
        StabStdMRadem_tmp    = [ StabStdMRadem_tmp; StabStdMRadem ];
        StabStdMtGauss_tmp   = [ StabStdMtGauss_tmp; StabStdMtGauss ];
        StabStdMtRadem_tmp   = [ StabStdMtRadem_tmp; StabStdMtRadem ];
    end
end

TrueStdMGauss  = TrueStdMGauss_tmp;
TrueStdMRadem  = TrueStdMRadem_tmp;
TrueStdMtGauss = TrueStdMtGauss_tmp;
TrueStdMtRadem = TrueStdMtRadem_tmp;

AsymStdMGauss  = AsymStdMGauss_tmp;
AsymStdMRadem  = AsymStdMRadem_tmp;
AsymStdMtGauss = AsymStdMtGauss_tmp;
AsymStdMtRadem = AsymStdMtRadem_tmp;

CohenStdMGauss  = CohenStdMGauss_tmp;
CohenStdMRadem  = CohenStdMRadem_tmp;
CohenStdMtGauss = CohenStdMtGauss_tmp;
CohenStdMtRadem = CohenStdMtRadem_tmp;

SNRStdMGauss    = SNRStdMGauss_tmp;
SNRStdMRadem    = SNRStdMRadem_tmp;
SNRStdMtGauss   = SNRStdMtGauss_tmp;
SNRStdMtRadem   = SNRStdMtRadem_tmp;

StabStdMGauss    = StabStdMGauss_tmp;
StabStdMRadem    = StabStdMRadem_tmp;
StabStdMtGauss   = StabStdMtGauss_tmp;
StabStdMtRadem   = StabStdMtRadem_tmp;   

% Total number of simulations
msim = size(SNRStdMtRadem,1);

save('simulations/estimQuantile_SNRCopeSet_processes.mat', 'TrueStdMGauss', 'TrueStdMRadem',...
     'TrueStdMtGauss', 'TrueStdMtRadem', 'CohenStdMGauss', 'CohenStdMRadem', 'CohenStdMtGauss',...
    'CohenStdMtRadem', 'SNRStdMGauss', 'SNRStdMRadem', 'SNRStdMtGauss', 'SNRStdMtRadem',...
    'AsymStdMGauss', 'AsymStdMtGauss', 'AsymStdMRadem', 'AsymStdMtRadem',...
    'StabStdMGauss', 'StabStdMtGauss', 'StabStdMRadem', 'StabStdMtRadem',...
    'msim', 'Nsubj', 'Lvec', 'FWHM', 'SNR', 'lvls', 'Mboot')