function [bdryValues] = linBdryEstim( field, c, weights, delta )

% Interpolate values of random field to boundary of an excursion set using
% 4 -connectivity
% Input:
% field:     random field over a domain in R^2, it is an (2+1)-dimensional array,
%            where the last dimension enumerates realisations
% c:         threshold value for excursion
% weights:   weights for weightening the values in the estimation of the
%            boundary. Options are:
%               0: outside and inside the excursion set has the same weight
%                  (default)
%               1: outside and inside weighted by expected value of fields
% delta:     field used to estimate the excursion set. Default value is
%            sample mean of the Input 'field'
%Output:
% bdryValues is a 2D array containing the values of the random fields on the
% estimated boundary from the mean. Each column is a realisation. The values
% are interpolations using a 4-connectivity grid and simple averages of the points inside and
% outside the excursion set.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 10/22/2018
%__________________________________________________________________________

sf = size(field) ; 

% Fill in unset optional values.
switch nargin
    case 2
        delta     = 666;
        weights   = 0;
    case 3
        delta     = 666;
end

if sum(delta(:))==666
    % Compute the excursion set estimated from the sample mean
    delta = mean(field, length(sf));
    A_c = ( delta  >= c);   
else
    % Compute the excursion set estimated from the sample mean
    A_c = ( delta  >= c);   
end

sA  = size(A_c);

switch length(sA)
    case 2
        %%%%%%%%%%%%%%%% Case 2D random field with 4-connectivity %%%%%%%%%%%%%%%%%
        %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the second component might be called vertical)
        horz = A_c(:,2:end) | A_c(:,1:end-1);

        %%% Compute the left shifted horizontal edges
        lshift            = A_c; % initialize
        lshift(:,1:end-1) = horz;
        lshift            = lshift & ~A_c;
        %%% Values of random field on left shifted horizontal edges of the
        %%% estimated value on boundary, note that the output is a concatinated
        %%% vector
        if( weights )
            % weight the points by its fraction in the mean sample
            hl_edges = ( repmat(delta(lshift(:,[sA(2) 1:sA(2)-1]))-c, [sf(end),1]) .* field(repmat(lshift, [1, 1, sf(end)]))...
                          + repmat(c-delta(lshift), [sf(end),1]) .* field(repmat(lshift(:,[sA(2) 1:sA(2)-1]), [1, 1, sf(end)])) )...
                          ./ repmat( -delta(lshift) + delta(lshift(:,[sA(2) 1:sA(2)-1])), [sf(end),1] ) ;  
        else
            % equally weight the points
            hl_edges = ( field(repmat(lshift, [1, 1, sf(end)]) ) + ...
                         field(repmat(lshift(:,[sA(2) 1:sA(2)-1]), [1, 1, sf(end)])) ) / 2;
        end
        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        hl_edges = reshape(hl_edges, [sum(lshift(:)) sf(end)]);

        %%% Compute the right shifted horizontal edges
        rshift          = A_c; % initialize
        rshift(:,2:end) = horz;
        rshift          = rshift & ~A_c;
        %%% Values of random field on right shifted horizontal edges
        if( weights )
            % weight the points by its fraction in the mean sample
            hr_edges = ( repmat(delta(rshift(:,[2:sA(2) 1]))-c, [sf(end),1]) .* field(repmat(rshift, [1, 1, sf(end)]))...
                          + repmat(c-delta(rshift), [sf(end),1]) .* field(repmat(rshift(:,[2:sA(2) 1]), [1, 1, sf(end)])) )...
                          ./ repmat( -delta(rshift) + delta(rshift(:,[sA(2) 1:sA(2)-1])), [sf(end),1] ) ;  
        else
            % equally weight the points
            hr_edges = ( field(repmat(rshift, [1,1,sf(end)]) ) + ...
                         field(repmat(rshift(:,[2:sA(2) 1]), [1,1,sf(end)])) ) / 2;
        end
        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        hr_edges = reshape(hr_edges, [sum(rshift(:)) sf(end)]);

        clear horz

        %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the first component might be called horizontal)
        vert = A_c(1:end-1,:) | A_c(2:end,:);
        %%% Compute the right shifted horizontal edges
        ushift = A_c;
        ushift(1:end-1,:) = vert;
        ushift = ushift & ~A_c;
        %%% Values of random field on up shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            vu_edges = ( repmat(delta(ushift([sA(1) 1:sA(1)-1],:))-c, [sf(end),1]) .* field(repmat(ushift, [1, 1, sf(end)]))...
                          + repmat(c-delta(ushift), [sf(end),1]) .* field(repmat(ushift([sA(1) 1:sA(1)-1],:), [1, 1, sf(end)])) )...
                          ./ repmat( -delta(ushift) + delta(ushift([sA(1) 1:sA(1)-1],:)), [sf(end),1] ) ;  
        else
            % equally weight the points
            vu_edges = ( field(repmat(ushift, [1,1,sf(end)]) ) + ...
                         field(repmat(ushift([sA(1) 1:sA(1)-1],:), [1,1,sf(end)])) ) / 2;
        end

        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        vu_edges = reshape(vu_edges, [sum(ushift(:)) sf(end)]);

        %%% Compute the down shifted vertical edges
        dshift = A_c;
        dshift(2:end,:)   = vert;
        dshift = dshift & ~A_c;
        %%% Values of random field on down shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            vd_edges = ( repmat(delta(dshift([2:sA(1) 1],:))-c, [sf(end),1]) .* field(repmat(dshift, [1, 1, sf(end)]))...
                          + repmat(c-delta(dshift), [sf(end),1]) .* field(repmat(dshift([2:sA(1) 1],:), [1, 1, sf(end)])) )...
                          ./ repmat( -delta(dshift) + delta(dshift([2:sA(1) 1],:)), [sf(end),1] ) ;  
        else
            % equally weight the points
            vd_edges = ( field(repmat(dshift, [1,1,sf(end)]) ) + ...
                         field(repmat(dshift([2:sA(1) 1],:), [1,1,sf(end)])) ) / 2;
        end
        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        vd_edges = reshape(vd_edges, [sum(dshift(:)) sf(end)]);
        % concatinated values of field on the linear edges 
        bdryValues = [ hr_edges; hl_edges; vu_edges; vd_edges ] ;
      case 3
        %%%%%%%%%%%%%%%% Case 3D random field with 6-connectivity %%%%%%%%%%%%%%%%%
        %%%%%% Horizontal edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the second component might be called vertical)
        horz = A_c(:,2:end,:) | A_c(:,1:end-1,:);

        %%% Compute the left shifted horizontal edges
        lshift              = A_c; % initialize
        lshift(:,1:end-1,:) = horz;
        lshift              = lshift & ~A_c;
        %%% Values of random field on left shifted horizontal edges of the
        %%% estimated value on boundary, note that the output is a concatinated
        %%% vector
        if( weights )
            % weight the points by its fraction in the mean sample
            hl_edges = ( repmat(delta(lshift(:,[sA(2) 1:sA(2)-1],:))-c, [sf(end),1]) .* field(repmat(lshift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(lshift), [sf(end),1]) .* field(repmat(lshift(:,[sA(2) 1:sA(2)-1],:), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(lshift) + delta(lshift(:,[sA(2) 1:sA(2)-1],:)), [sf(end),1] ) ;  
        else
            % equally weight the points
            hl_edges = ( field(repmat(lshift, [1, 1, 1, sf(end)]) ) + ...
                         field(repmat(lshift(:,[sA(2) 1:sA(2)-1],:), [1, 1, 1, sf(end)])) ) / 2;
        end

        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        hl_edges = reshape(hl_edges, [sum(lshift(:)) sf(end)]);

        %%% Compute the right shifted horizontal edges
        rshift            = A_c; % initialize
        rshift(:,2:end,:) = horz;
        rshift            = rshift & ~A_c;
        %%% Values of random field on right shifted horizontal edges
        if( weights )
            % weight the points by its fraction in the mean sample
            hr_edges = ( repmat(delta(rshift(:,[2:sA(2) 1],:))-c, [sf(end),1]) .* field(repmat(rshift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(rshift), [sf(end),1]) .* field(repmat(rshift(:,[2:sA(2) 1],:), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(rshift) + delta(rshift(:,[2:sA(2) 1],:)), [sf(end),1] ) ;  
        else
            % equally weight the points
            hr_edges = ( field(repmat(rshift, [1, 1, 1, sf(end)]) ) + ...
                         field(repmat(rshift(:,[2:sA(2) 1],:), [1, 1, 1, sf(end)])) ) / 2;
        end

        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        hr_edges = reshape(hr_edges, [sum(rshift(:)) sf(end)]);

        clear horz

        %%%%%% Vertical edges (note that this uses the matlab image nomenclature,
        %%%%%% usually the first component might be called horizontal)
        vert = A_c(1:end-1,:,:) | A_c(2:end,:,:);
        %%% Compute the right shifted horizontal edges
        ushift = A_c;
        ushift(1:end-1,:,:) = vert;
        ushift = ushift & ~A_c;
        %%% Values of random field on up shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            vu_edges = ( repmat(delta(ushift([sA(1) 1:sA(1)-1],:,:))-c, [sf(end),1]) .* field(repmat(ushift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(ushift), [sf(end),1]) .* field(repmat(ushift([sA(1) 1:sA(1)-1],:,:), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(ushift) + delta(ushift([sA(1) 1:sA(1)-1],:,:)), [sf(end),1]) ;  
        else
            % equally weight the points
            vu_edges = ( field(repmat(ushift, [1, 1, 1, sf(end)]) ) + ...
                         field(repmat(ushift([sA(1) 1:sA(1)-1],:,:), [1, 1, 1, sf(end)])) ) / 2;
        end

        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        vu_edges = reshape(vu_edges, [sum(ushift(:)) sf(end)]);

        %%% Compute the down shifted vertical edges
        dshift = A_c;
        dshift(2:end,:,:)   = vert;
        dshift = dshift & ~A_c;
                %%% Values of random field on down shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            vd_edges = ( repmat(delta(dshift([2:sA(1) 1],:,:))-c, [sf(end),1]) .* field(repmat(dshift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(dshift), [sf(end),1]) .* field(repmat(dshift([2:sA(1) 1],:,:), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(dshift) + delta(dshift([2:sA(1) 1],:,:)), [sf(end),1] ) ;  
        else
            % equally weight the points
            vd_edges = (    field(repmat(dshift, [1, 1, 1, sf(end)])) + ...
                            field(repmat(dshift([2:sA(1) 1],:,:), [1, 1, 1, sf(end)])) ) / 2;
        end
        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        vd_edges = reshape(vd_edges, [sum(dshift(:)) sf(end)]);
        
        %%%%%% depth edges
        depth = A_c(:,:,1:end-1) | A_c(:,:,2:end);
        %%% Compute the back shifted depth edges
        bshift = A_c;
        bshift(:,:,1:end-1) = depth;
        bshift = bshift & ~A_c;
        %%% Values of random field on up shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            db_edges = ( repmat(delta(bshift(:,:,[sA(3) 1:sA(3)-1]))-c, [sf(end),1]) .* field(repmat(bshift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(bshift), [sf(end),1]) .* field(repmat(bshift(:,:,[sA(3) 1:sA(3)-1]), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(bshift) + delta(bshift(:,:,[sA(3) 1:sA(3)-1])), [sf(end),1] ) ;  
        else
            % equally weight the points
            db_edges = (    field(repmat(bshift, [1, 1, 1, sf(end)]) ) + ...
                            field(repmat(bshift(:,:,[sA(3) 1:sA(3)-1]), [1, 1, 1, sf(end)])) ) / 2;
        end

        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        db_edges = reshape(db_edges, [sum(bshift(:)) sf(end)]);

        %%% Compute the fron shifted depth edges
        fshift = A_c;
        fshift(:,:,2:end)   = depth;
        fshift = fshift & ~A_c;
        %%% Values of random field on down shifted vertical edges
        if( weights )
            % weight the points by its fraction in the mean sample
            df_edges = ( repmat(delta(fshift(:,:,[2:sA(3) 1]))-c, [sf(end),1]) .* field(repmat(fshift, [1, 1, 1, sf(end)]))...
                          + repmat(c-delta(fshift), [sf(end),1]) .* field(repmat(fshift(:,:,[2:sA(3) 1]), [1, 1, 1, sf(end)])) )...
                          ./ repmat( -delta(fshift) + delta(fshift(:,:,[2:sA(3) 1])), [sf(end),1] ) ;  
        else
            % equally weight the points
            df_edges = ( field(repmat(fshift, [1, 1, 1, sf(end)]) ) + ...
                     field(repmat(fshift(:,:,[2:sA(3) 1]), [1, 1, 1, sf(end)])) ) / 2;
        end
        
        %%% Form the boundary value estimates back to a matrix with second
        %%% dimension equal to the number of realisations
        df_edges = reshape(df_edges, [sum(fshift(:)) sf(end)]);
        
        % concatinated values of field on the linear edges 
        bdryValues = [ hr_edges; hl_edges; vu_edges; vd_edges; db_edges; df_edges ] ;
end