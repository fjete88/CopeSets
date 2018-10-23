function [coveringRate] = CovRateLvlSets( truef, hatf, thresh, c )

% Compute covering rate of level sets
% Input:
% truef:  underlying true function on an N-dimensional array, it is saved as
%         an N-dimensional array, where the values correspond to the heights
% hatf:   array of size [size(truef), M] containing M estimates of truef
% thresh: array of size [size(truef), M, 2] corresponding to the thresholds for the
%         hatf, [ , , 1] are lower bounds, [ , , 2] are upper bounds
% c:      targeted levelset
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

if ndims( thresh ) == N+2
    coveringRate = 1-sum(any(reshape(...
                ( truefm < c & (hatf >= squeeze(thresh(index{:},2))) ) |...
                ( truefm >= c & (hatf < squeeze(thresh(index{:},1))))...
               ,[prod(sf) M]), 1 )) / M ;
else
    nquantiles = size(thresh,N+1);
    coveringRate = zeros(1,nquantiles);
    index2 = repmat( {':'}, 1, N );
    index3 = repmat( {':'}, 1, 2 );
    
    for i = 1:nquantiles
        thresh_tmp = squeeze(thresh(index2{:}, i, index3{:}));
        coveringRate(i) = 1-sum(any(reshape(...
                ( truefm < c & (hatf >= squeeze(thresh_tmp(index{:},2))) ) |...
                ( truefm >= c & (hatf < squeeze(thresh_tmp(index{:},1))))...
               ,[prod(sf) M]), 1 )) / M ;
    end
end
