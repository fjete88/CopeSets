function [dilation] = dilateSet( F )

% dilation of a 2D logical field along a dimension using 4 -connectivity
% Input:
% F:         logical field over a domain in R^2
%Output:
% dilation are the voxel which are in the complement of F and its dilation


vert = F(1:end-1,:) | F(2:end,:);
%%% Compute the right shifted horizontal edges
ushift = F;
ushift(1:end-1,:) = vert;
ushift = ushift & ~F;

%%% Compute the down shifted vertical edges
dshift = F;
dshift(2:end,:)   = vert;
dshift   = dshift & ~F;


horz = F(:,2:end) | F(:,1:end-1);
%%% Compute the left shifted horizontal edges
lshift            = F; % initialize
lshift(:,1:end-1) = horz;
lshift            = lshift & ~F;

%%% Compute the right shifted horizontal edges
rshift          = F; % initialize
rshift(:,2:end) = horz;
rshift          = rshift & ~F;

%%% Produce the output field
dilation        = rshift | lshift | dshift | ushift; 