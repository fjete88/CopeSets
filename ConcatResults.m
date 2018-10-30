function result = ConcatResults( results )
% This function concatinates the results from CopeSet_SimSkript into the
% results with the final number of simulations. 
% Input:
%  results: cell array containing the results of different simulations. The last
%           dimension is assumed to enumerate identical simulations
% Output:
%  - result cell array containg the simulation results with the total
%    number of simulations.
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/29/2018
%__________________________________________________________________________

% get the dimension of the cell array
dimR  = size(results);
index = repmat( {':'}, 1, length(dimR)-1 );

% allocate variable for the output
result = reshape(results(index{:},1), [prod(dimR(1:end-1)) 1]);

% reshape to loop over the different simulations
tmp = reshape(results, [prod(dimR(1:end-1)), dimR(end)]);

% loop over the different simulations
for k = 1:prod(dimR(1:end-1))
    % loop over the results of a single experiment
    for l = 1:dimR(end)
        % add the covering rates from the same experiments and save the
        % total number of simulations
        if l==1
            % initialize the variable counting the number of simulations
            nsim = tmp{k,l}.nsim;
            % weight covering rate by number of simulations
            result{k}.covRate.truebdry = nsim*tmp{k}.covRate.truebdry;
            result{k}.covRate.linbdry  = nsim*tmp{k}.covRate.linbdry;
            result{k}.covRate.erodbdry = nsim*tmp{k}.covRate.erodbdry;
        elseif l==dimR(end)
            % update the number of simulations
            nsim = nsim + tmp{k,l}.nsim;
            % concatinate the estimates of the quantiles
            result{k}.quant.truebdry     = [result{k}.quant.truebdry tmp{k,l}.quant.truebdry];
            result{k}.quant.linbdry      = [result{k}.quant.linbdry tmp{k,l}.quant.linbdry];
            result{k}.quant.erodbdrybdry = [result{k}.quant.erodbdry tmp{k,l}.quant.erodbdry];
            % save the final number of simualtions to the output
            result{k}.nsim = nsim;
            % transform back to covering percentages
            result{k}.covRate.truebdry = result{k}.covRate.truebdry / nsim;
            result{k}.covRate.linbdry  = result{k}.covRate.linbdry  / nsim;
            result{k}.covRate.erodbdry = result{k}.covRate.erodbdry / nsim;
            % compute the new standard errors of the simulation
            result{k}.stdErr.rough = sqrt(result{k}.lvls .* (1-result{k}.lvls)/nsim);
            result{k}.stdErr.truebdry = sqrt(result{k}.covRate.truebdry .* (1-result{k}.covRate.truebdry) / nsim);
            result{k}.stdErr.linbdry  = sqrt(result{k}.covRate.linbdry  .* (1-result{k}.covRate.linbdry)  / nsim);
            result{k}.stdErr.erodbdry = sqrt(result{k}.covRate.erodbdry .* (1-result{k}.covRate.erodbdry) / nsim);
        else
            % update the number of simulations
            nsim = nsim + tmp{k,l}.nsim;
            % add the weighted by number of simulations to the already weighted covering rate
            result{k}.covRate.truebdry = result{k}.covRate.truebdry + tmp{k,l}.nsim * tmp{k}.covRate.truebdry;
            result{k}.covRate.linbdry  = result{k}.covRate.linbdry  + tmp{k,l}.nsim * tmp{k}.covRate.linbdry;
            result{k}.covRate.erodbdry = result{k}.covRate.erodbdry + tmp{k,l}.nsim * tmp{k}.covRate.erodbdry;
            % concatinate the estimates of the quantiles
            result{k}.quant.truebdry     = [result{k}.quant.truebdry tmp{k,l}.quant.truebdry];
            result{k}.quant.linbdry      = [result{k}.quant.linbdry tmp{k,l}.quant.linbdry];
            result{k}.quant.erodbdrybdry = [result{k}.quant.erodbdry tmp{k,l}.quant.erodbdry];
        end
    end % loop over the results of a single experiment
end % loop over the different simulations
% reshape results to match the input size
result = reshape(result, dimR(1:end-1));
end % end of function