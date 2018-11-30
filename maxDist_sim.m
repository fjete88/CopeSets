function results = maxDist_sim( nsim, df, paramNoise, bdry_mask, batchnumber, pool_num )
% Simulates the maximum distribution of an error process or an t-field constructed from an error process.
% Input:
%  nsim:        number of simulations used to estimate the covering rate
%  df:          integer vector containing the degrees of freedom of the considered
%               t-process. If df=0 the generated process is Gaussian.
%  paramNoise:  structure containing the parameters for the error process
%  bdry_mask:   boolean array of size paramNoise.dim indicating which
%               voxels belong to the boundary
%  batchnumber: amount of simultaneously constructed fields
%  pool_num:    number of parallel computing units (default=1)
%
% Output:
%  - results is a structure containg the results and parameter of the
%  simulation
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 11/08/2018
%__________________________________________________________________________

%%%%%% Add missing input by using default values
if nargin ==4
    pool_num = 1;
end

%%%%%% Check whether batch number divides nsim
if( mod(nsim,batchnumber)~=0 )
    error("Choose the batch number such that nsim/batchnumber is an integer!")
end
if ~all(size(bdry_mask)==paramNoise.dim)
    error("Error the boundary mask must have the same dimension as the noise field!")
end

% open connection to GPUs, if not already established
if( pool_num > 1 )
    % save the state of the GPU's are open already
    state_gcp = isempty(gcp('nocreate'));
    if( state_gcp )
        parpool( pool_num );
        state_gcp = 42;
    end
end

%%%%%% Compute some useful constants for the simulation
% Compute the dimension of the domain of the data
dim    = paramNoise.dim;
D      = length(dim);
% compute indices for variable access on fields with different dimensions
index  = repmat( {':'}, 1, D );
% initialize the variable storing the maximum distribution
maxDist = zeros([length(df) nsim]);

%%% Loop over simulations/batch number
for nn = 1:(nsim/batchnumber)
    %%%%%% Set the batch range
    batchrange = (batchnumber*(nn-1)+1):nn*batchnumber;
    %%%%%% Simulate a small batch of fields along the boundary
    % generate random field
    if D==2
        f = SmoothField2D( max(df)+1, batchnumber, paramNoise.FWHM, paramNoise.dim,...
                           paramNoise.noise, paramNoise.nu, paramNoise.kernel,...
                           paramNoise.bin, pool_num );
    elseif D==3
        f = SmoothField3D( max(df)+1, batchnumber, paramNoise.FWHM, paramNoise.dim,...
                           paramNoise.noise, paramNoise.nu, paramNoise.kernel,...
                           paramNoise.bin, pool_num );
    else
        error("Currently only 2D and 3D fields are supported.")
    end
    % loop over the different df
    for countdf = 1:length(df)
        if df(countdf) > 0
            % construct t-processes with df degrees of freedom
            tmpf      = f(index{:}, df(countdf)+1,:) ./ ...
                        sqrt( sum(f(index{:}, 1:df(countdf),:).^2, D+1) / df(countdf) );
            % restrict t-field to the boundary
            tmpf = reshape( tmpf(repmat(bdry_mask, [ones([1 D]) batchnumber])),...
                            [sum(bdry_mask(:)), batchnumber] );
        else
            % 'construct' the Gaussian field
            tmpf      = f(index{:}, df(countdf)+1,:);
            % restrict field to the boundary
            tmpf = reshape( tmpf(repmat(bdry_mask, [ones([1 D]) batchnumber])),...
                            [sum(bdry_mask(:)), batchnumber] );
        end
        % store maximum over the boundary
        maxDist(countdf, batchrange) = max(abs(tmpf),[],1);
    end
end % Loop over simulations/batch number

% close connection to GPUs
if( pool_num > 1 ) 
    if ( state_gcp == 42 )
        delete(gcp)
    end
end

%%%%%%% report the results of the simulation
% save the maxStat
results.maxDist = maxDist;
% params of the simulation
results.df = df;
results.nsim = nsim;
results.paramNoise  = paramNoise;
results.bdry  = bdry_mask;