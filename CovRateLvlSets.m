function [coveringRate] = CovRateLvlSets( truef, hatf, thresh, c, Bdry_test )

% Compute covering rate of level sets
% Input:
% truef:  underlying true function on an N-dimensional array, it is saved as
%         an N-dimensional array, where the values correspond to the heights
% hatf:   array of size [size(truef), M] containing M estimates of truef
% thresh: array of size [size(truef), M, 2] corresponding to the thresholds for the
%         hatf, [ , , 1] are lower bounds, [ , , 2] are upper bounds
% c:      targeted levelset
% Bdry_test: if 0, applies old boundary test comparing just binarized sets.
%        if 1, also used lintear interpolation boundary test. 
%Output:
% coveringRate is computed from all available hatf and computes the
% frequency that the true c-levelset of truef is completly contained in the
% thresholded hatfs.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/05/2018
%__________________________________________________________________________

sf    = size( truef );
N     = ndims( truef );
M     = size( hatf, N+1 );
index = repmat( {':'}, 1, N+1 );

truefm = reshape( repmat(reshape( truef, [prod(sf) 1] ), 1, M ) , [sf M] ) ;

% Finding the edges and weights for interpolating true field = c
if ( Bdry_test )
    bdry_params = getBdryparams(truef, c);
    bdry_length = length(bdry_params.lshift.w1);
    low_thresh_bdry_values  = zeros(bdry_length, M);
    high_thresh_bdry_values = zeros(bdry_length, M);    
    hatf_bdry_values        = zeros(bdry_length, M);   
end

nquantiles = size(thresh,N+1);
coveringRate = zeros(1,nquantiles);
index2 = repmat( {':'}, 1, N );
index3 = repmat( {':'}, 1, 2 );

for i = 1:nquantiles
    thresh_tmp = squeeze(thresh(index2{:}, i, index3{:}));

    if ( Bdry_test )
        % finding values of estimated field at locations where true
        % field = c 
        for j = 1:M          
            low_thresh_bdry_values(:,j)  = getBdryvalues(squeeze(thresh_tmp(index2{:},j,1)), bdry_params);
            high_thresh_bdry_values(:,j) = getBdryvalues(squeeze(thresh_tmp(index2{:},j,2)), bdry_params);
            hatf_bdry_values(:,j)        = getBdryvalues(hatf(index2{:},j), bdry_params);
        end
        violations = sum( ...
            any(reshape(...
                ( truefm < c & (hatf >= squeeze(thresh_tmp(index{:},2)))) |...
                ( truefm >= c & (hatf < squeeze(thresh_tmp(index{:},1))))...
            ,[prod(sf) M]), 1) | ...
            any(...
                hatf_bdry_values >= high_thresh_bdry_values | ...
                hatf_bdry_values < low_thresh_bdry_values...
            , 1 )) / M;
    else 
       violations = sum(any(reshape(...
        ( truefm < c & (hatf >= squeeze(thresh_tmp(index{:},2))) ) |...
        ( truefm >= c & (hatf < squeeze(thresh_tmp(index{:},1))))...
       ,[prod(sf) M]), 1 )) / M ; 
    end
    coveringRate(i) = 1-violations; 
end
