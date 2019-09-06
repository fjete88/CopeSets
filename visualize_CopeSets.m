function [] = visualize_CopeSets( F, c, lvl, quantEstim, bdry_type,...
                                  center, normalize, location, delta )
% Computes CoPe all ingredients for CoPe sets.
% Input:
%  F:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%        where the last dimension enumerates the samples
%  c:    threshold for excursions
%  lvls:       vector containing the required confidence levels. Must be
%              between 0 and 1.
%  quantEstim: structure containing the name and the parameters for the
%              quantile estimation method. Choices:
%               {
%                quantEstim.name = 'multiplierbootstrap'
%                quantEstim.params:
%                   Mboot:     amount of bootstrap replicates (default=5e3)
%                   method:    option for the bootstrap estimator (default='t')
%               }
%  bdry_type: currently 'interval' and 'points' are supported.
%  center:    option to center the field using the sample mean (default=1)
%  normalize: option to normalize the field by sample variance (default=1)
%  delta:     required, if bdry_type is equal to 'true'. This is the true
%             population mean function given on a D-dimensional array
% Output:
%  - thresh is the threshold lower and upper for the sample mean in order to
%    be in the estimated lower and upper excursion sets 
%  - quantile is the bootstrapped quantile of the maximum distribution of the 
%    input processes
%  - hatdelta is the sample mean of the fields
%  - hatsigma is the sample variance of the fields
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 11/30/2018
%__________________________________________________________________________

sF  = size(F);
dim = sF(1:end-1);
scale = 4/9;
WidthFig   = 1300;
HeightFig  = WidthFig * scale;

[intervalCopeSets, quantInt, hatdelta] = ...
        CopeSets_simultaneous(  F, c, lvl, quantEstim,...
                       bdry_type(1), center, normalize, delta );
[CopeSets_c1, quantIntc1] = ...
        CopeSets( F, c(1), lvl, quantEstim,...
                       bdry_type(2), center, normalize, delta );
[CopeSets_c2, quantIntc2] = ...
        CopeSets( F, c(2), lvl, quantEstim,...
                       bdry_type(2), center, normalize, delta );
                   
close all
ncol = 256;
colMap = jet(ncol);
I_c = [min(hatdelta(:)) max(hatdelta(:))] .* [0.975, 1.025];  %[min(Y(:)), max(Y(:))] %      

a = (1-ncol)/(I_c(1)-I_c(2));
b = ncol - I_c(2)*a;

col_c1 = round(c(1) * [0.8, 1.1],2);
col_c2 = round(c(2) * [0.8, 1.1],2);

copeImage = ones(dim)*I_c(1);
copeImage((hatdelta >= intervalCopeSets(:,:,1,1,1))) = col_c1(1);
copeImage((hatdelta >= intervalCopeSets(:,:,1,2,1))) = col_c1(2);
copeImage((hatdelta >= intervalCopeSets(:,:,1,1,2))) = col_c2(1);
copeImage((hatdelta >= intervalCopeSets(:,:,1,2,2))) = col_c2(2);

copeImage_c1 = ones(dim)*I_c(1);
copeImage_c1((hatdelta >= CopeSets_c1(:,:,1))) = col_c1(1);
copeImage_c1((hatdelta >= CopeSets_c1(:,:,2))) = col_c1(2);

copeImage_c2 = ones(dim)*I_c(1);
copeImage_c2((hatdelta >= CopeSets_c2(:,:,1))) = col_c2(1);
copeImage_c2((hatdelta >= CopeSets_c2(:,:,2))) = col_c2(2);

% true support
true_c = zeros(dim);
true_c(delta(:,:)>=c(1)) = c(1);
true_c(delta(:,:)>=c(2)) = c(2);

%%
figure, clf
set(gcf, 'Position', [ 300 300 WidthFig HeightFig]);
set(gcf,'PaperPosition', [ 300 300 WidthFig HeightFig])
colormap(colMap)
if delta ~= -666
    subplot(2, 3, 1)
    imagesc(true_c)
    title("True Level Sets")
    caxis(I_c)
    colorbar
    axis square;

    subplot(2, 3, 2)
    imagesc(hatdelta)
    title("Sample Mean")
    colorbar
    axis square;

    subplot(2, 3, 3)
    imagesc(F(:,:,1))
    title("Sample")
    colorbar
    axis square;
else
    subplot(2, 3, 1)
    imagesc(hatdelta)
    title("Sample Mean")
    colorbar
    axis square;
    
    subplot(2, 3, 2)
    imagesc(F(:,:,1))
    title("Sample")
    colorbar
    axis square;

    subplot(2, 3, 3)
    imagesc(F(:,:,2))
    title("Sample")
    colorbar
    axis square;
end

subplot(2, 3, 4);
imagesc(copeImage)
title([strcat("Thresholds ", num2str(c(1)), " and ", num2str(c(2))),...
      strcat("Quantile ", num2str(round(quantInt,2)))])
caxis(I_c)
colorbar
axis square;

subplot(2, 3, 5);
imagesc(copeImage_c1)
title([strcat("Thresholds ", num2str(c(1))),...
      strcat("Quantile ", num2str(round(quantIntc1,2)))])
caxis(I_c)
colorbar
axis square;

subplot(2, 3, 6);
imagesc(copeImage_c2)
title([strcat("Thresholds ", num2str(c(2))),...
      strcat("Quantile ", num2str(round(quantIntc2,2)))])
caxis(I_c)
colorbar
axis square;

hold off

set(gcf,'papersize',[12 12*scale])
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(strcat(location, '.png' ), '-dpng')