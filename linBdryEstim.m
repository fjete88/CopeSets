function [bdryValues] = linBdryEstim( field, c, delta )

% Interpolate values of random field to boundary of an excursion set using
% 4 -connectivity
% Input:
%  field:     random field over a domain in R^2, it is an (2+1)-dimensional array,
%            where the last dimension enumerates realisations
%  c:         threshold value for excursion
%  delta:     field used to estimate the excursion set. Default value is
%            sample mean of the Input 'field'
% Output:
%  bdryValues is a 2D array containing the values of the random fields on the
%  estimated boundary from the mean. Each column is a realisation. The values
%  are interpolations using a 4-connectivity grid and simple averages of the points inside and
%  outside the excursion set.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 11/30/2018
%__________________________________________________________________________

sf = size(field) ; 
index = repmat( {':'}, 1, length(sf)-1);

% Fill in unset optional values.
switch nargin
    case 1
        delta     = mean(field, length(sf));
end

bdry_param  = getBdryparams(delta, c);
bdryValues = zeros(bdry_param.length, sf(end));

for n = 1:sf(end)
    bdryValues(:,n) = getBdryvalues(field(index{:},n), bdry_param);
end
end