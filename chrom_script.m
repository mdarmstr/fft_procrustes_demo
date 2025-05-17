function chrom_script(matrix1,matrix2,shift)

% Set default font size for all text and axes
set(groot, 'DefaultAxesFontSize', 12);
set(groot, 'DefaultTextFontSize', 12);

matrix2 = shift_operator(matrix2,shift);

maxVal = max(matrix2(:));
newMax = 0.66 * maxVal;
newMin = min(matrix2(:));

% apply as z‑limits

[imgH,imgW] = size(matrix1);

% Compute FFT of shifted matrix2
F_original = fft2(matrix2);
F_orig_centered = fftshift(F_original);

% Define down-sampling parameters.
downsampleFactor = 1;  % For a 256x256 image, this gives a 256x256 block.
imgH_small = floor(imgH / downsampleFactor);
imgW_small = floor(imgW / downsampleFactor);
idxH = (imgH/2 - imgH_small/2 + 1):(imgH/2 + imgH_small/2);
idxW = (imgW/2 - imgW_small/2 + 1):(imgW/2 + imgW_small/2);

% Down-sampling
F_orig_small = F_orig_centered(idxH, idxW);
Fvec_orig_small = F_orig_small(:);

alignVariant = @(imgVar) frequencyAlign(imgVar, idxH, idxW, imgH, imgW, Fvec_orig_small);

imgAlign = alignVariant(matrix1);


% Create figure
figure;

% Top-left subplot: original matrix1
subplot(2,2,1);
imagesc(matrix1);
title(sprintf('Chrom1: fNorm=%.4e', norm(matrix1, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
yline(shift(1),'r-','LineWidth',1)
xline(shift(2),'r-','LineWidth',1)
clim([newMin, newMax]);


% Top-right subplot: shifted matrix2
subplot(2,2,2);
imagesc(matrix2);
title(sprintf('Chrom1Shift: fNorm=%.4e', norm(matrix2, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
clim([newMin, newMax]);


% Bottom-left subplot: aligned image
subplot(2,2,3);
imagesc(imgAlign);
title(sprintf('Chrom1Align: fNorm=%.4e', norm(imgAlign, 'fro')), 'FontSize', 12);
xlabel('Relative ^1t_R', 'FontSize', 12);
ylabel('Relative ^2t_R', 'FontSize', 12);
set(gca, 'YDir', 'normal');
clim([newMin, newMax]);


% Export figure
exportgraphics(gcf, 'chrom_analysis.pdf', 'Resolution', 500);

end
