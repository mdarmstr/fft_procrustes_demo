clear all; clc;

%Data from online repository - see manuscript.
c29 = load('C29_1.mat');
c29r = load('C29R_1.mat');

numPerm = 10000;

set(groot, 'DefaultAxesFontSize', 12);
set(groot, 'DefaultTextFontSize', 12);

ticc29 = c29.tempstruct.tic(1:843000);
ticc29r = c29r.tempstruct.tic;

maxVal = max(ticc29r(:));
newMax = 0.66 * maxVal;
newMin = min(ticc29r(:));

%correlation before analysis
cosCorr_1 = dot(ticc29(:), ticc29r(:)) / (norm(ticc29(:)) * norm(ticc29r(:)));
cosCorrP = zeros(numPerm,1);

tic_matc29 = reshape(ticc29,[200*2.5,1686]);
tic_matc29r = reshape(ticc29r,[200*2.5,1686]);

[imgH,imgW] = size(tic_matc29r);

%% 3. Frequency-Domain Alignment Procedure
% Compute the 2D FFT of the original image and center it.
F_original = fft2(tic_matc29);
F_orig_centered = fftshift(F_original);

% Define down-sampling parameters.
downsampleFactor = 1;  % Do downsampling
imgH_small = imgH / downsampleFactor;
imgW_small = imgW / downsampleFactor;
idxH = (imgH/2 - imgH_small/2 + 1):(imgH/2 + imgH_small/2);
idxW = (imgW/2 - imgW_small/2 + 1):(imgW/2 + imgW_small/2);

% Extract and vectorize the central block from the original FFT.
F_orig_small = F_orig_centered(idxH, idxW);
Fvec_orig_small = F_orig_small(:);

% Define a helper function for frequency-domain alignment.
alignVariant = @(imgVar) frequencyAlign(imgVar, idxH, idxW, imgH, imgW, Fvec_orig_small);

imgAlign = alignVariant(tic_matc29r);

cosCorr_2 = dot(ticc29(:), imgAlign(:)) / (norm(ticc29(:)) * norm(imgAlign(:)));

ticc29rA = imgAlign(:);

% Create figure
figure;

% Top-left subplot: original matrix1
subplot(2,2,1);
imagesc(tic_matc29);
title(sprintf('C29: fNorm=%.4e', norm(tic_matc29, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
clim([newMin, newMax]);


% Top-right subplot: shifted matrix2
subplot(2,2,2);
imagesc(tic_matc29r);
title(sprintf('C29R: fNorm=%.4e', norm(tic_matc29r, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
clim([newMin, newMax]);


% Bottom-left subplot: aligned image
subplot(2,2,3);
imagesc(imgAlign);
title(sprintf('C29RAlign: fNorm=%.4e', norm(imgAlign, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
clim([newMin, newMax]);


% Export figure
exportgraphics(gcf, 'chrom_analysis2.pdf', 'Resolution', 500);



