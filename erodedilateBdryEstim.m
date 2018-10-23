function [bdryValues] = erodedilateBdryEstim( field, c, delta )

% Erode/dilate values of random field to boundary of an excursion set using
% 4 -connectivity
% Input:
% field:     random field over a domain in R^2, it is an (2+1)-dimensional array,
%            where the last dimension enumerates realisations
% c:         threshold value for excursion
% erode:     if 1 include eroded points into boundary estimation (either erode or dilate has to be 1)
% dilate:    if 1 include dilated points into boundary estimation
%Output:
% bdryValues is a vector containing the values of the random field on the
% estimated boundary from the mean. The values are interpolations using
% a 4-connectivity grid and simple averages of the points inside and
% outside the excursion set.
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/05/2018
%__________________________________________________________________________

sf = size(field) ; 

% Fill in unset optional values.
switch nargin
    case 2
        delta     = 666;
end

if sum(delta(:))==666
    % Compute the excursion set estimated from the sample mean
    A_c = ( mean(field, length(sf))  >= c);   
else
    % Compute the excursion set estimated from the sample mean
    A_c = ( delta  >= c);   
end

sA = size(A_c);

switch length(sA)
    case 2
        mask       = dilateSet(A_c) | dilateSet(~A_c);
        bdryValues = reshape(field(repmat(mask, [1,1,sf(end)]) ), [sum(mask(:)) sf(end)]);
    case 3
        
end