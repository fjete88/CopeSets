%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     This simulation estimates the maximum quantile of Z and T processes
%%%     with different bootstrap estimators
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% clear workspace
clear
close all

%%%%%% set path for the files
cd  /home/drtea/Research/MatlabPackages/CopeSets

%% %%%%%%%%%%%%%%%%%%%% Bootstrap estimator simulations %%%%%%%%%%%%%%%%%%%
%%%%%% simulation paramters
% Dimension of the field
dim  = [[50 50];[100 50];[500,50]];
% smoothing FWHM for the errorfield
FWHM = [5 10];
% Sample sizes
nvec = [30 60 120 240 300];
% Considered confidence levels
lvls = [0.85 0.90 0.95];
% number of simulations
nsim = 3e3;
% number of bootstrap replicates used
Mboot = 5e3;
% compute residuals or not!
Res = 1;

%%%%%% Initialize the output from the simulations
gaussreg   = zeros([length(lvls), length(nvec), nsim, length(FWHM), size(dim,1)]);
gausst     = gaussreg;
rademreg   = gaussreg;
rademt     = gaussreg;
nonpareg   = gaussreg;
nonpat     = gaussreg;

%%%%%% simulations of maximum distribution estimates using different bootstrap
%%%%%% estimators
for d = 1:3
    mask = zeros(dim(d,:));
    mask(:,25) = ones([dim(d,1) 1]);
    mask = boolean(mask);
    for f = 1:length(FWHM)
        paramNoise = struct( 'FWHM', [FWHM(f), FWHM(f)], 'dim', dim(d,:), 'noise', "normal", 'nu', '',...
                                     'kernel', "gauss", 'bin', 0, 'sd', ones(dim(d,:)) );
        tic
        for j = 1:nsim
           % Generate errorfields
           R = SmoothField2D( max(nvec), 1, paramNoise.FWHM, paramNoise.dim,...
                              paramNoise.noise, paramNoise.nu,...
                              paramNoise.kernel, paramNoise.bin );
           if Res == 1
                   R = R - mean(R,3);
                   R = R ./ std(R, 0, 3);
           end
           for n = nvec
               nind = find(n==nvec);
               [gaussreg(:,nind,j,f,d), ~] = MultiplierBoots( R(:,:,1:n), lvls, ...
                                                Mboot, mask, 'gaussian', 'regular' );
               [gausst(:,nind,j,f,d), ~]   = MultiplierBoots( R(:,:,1:n), lvls, ...
                                                Mboot, mask, 'gaussian', 't' );
               [rademreg(:,nind,j,f,d), ~] = MultiplierBoots( R(:,:,1:n), lvls,...
                                                Mboot, mask, 'rademacher', 'regular' ); 
               [rademt(:,nind,j,f,d), ~]   = MultiplierBoots( R(:,:,1:n), lvls, ...
                                                Mboot, mask, 'rademacher', 't' );
               [nonpareg(:,nind,j,f,d), ~] = NonParametricBoots( R(:,:,1:n), lvls, ...
                                                Mboot, mask, 'regular' );
               [nonpat(:,nind,j,f,d), ~]   = NonParametricBoots( R(:,:,1:n), lvls, ...
                                                Mboot, mask, 't' );
           end
        end
        toc
    end
end
%%%% save final results
save(strcat('simulations/ResultBootstrappedQuantilesRes', num2str(Res)), 'gaussreg', 'gausst',...
     'rademreg', 'rademt', 'nonpareg', 'nonpat', 'nvec', 'Mboot', 'FWHM', 'dim',...
     'nsim', 'lvls')