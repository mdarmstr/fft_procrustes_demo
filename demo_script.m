%% Set the random seed for reproducibility
rng(1);

%% 1. Generate the Original and Basic Shifted Image
% Parameters
imgH   = 256;            % image height (rows)
imgW   = 256;            % image width (columns)
nBlobs = 20;             % number of Gaussian blobs
alpha  = 0.5;            % blending factor for the shift

% Create coordinate grid
[xGrid, yGrid] = meshgrid(1:imgW, 1:imgH);

% Preallocate images and parameter arrays
imgOriginal = zeros(imgH, imgW);
imgShifted  = zeros(imgH, imgW);
centers     = zeros(nBlobs, 2);
amplitudes  = zeros(nBlobs, 1);
sigmas      = zeros(nBlobs, 1);

for k = 1:nBlobs
    % Choose random blob parameters (avoid edges)
    centerX   = randi([10, imgW-10]);
    centerY   = randi([10, imgH-10]);
    centers(k,:) = [centerX, centerY];
    
    amplitude = 0.5 + rand();      % between 0.5 and 1.5
    sigma     = 3 + 7*rand();       % between 3 and 10 pixels
    amplitudes(k) = amplitude;
    sigmas(k) = sigma;
    
    % Create a symmetric Gaussian blob for the original image.
    blob = amplitude * exp( -((xGrid - centerX).^2 + (yGrid - centerY).^2) / (2*sigma^2) );
    imgOriginal = imgOriginal + blob;
    
    % Compute shifted coordinates using a logarithmic mapping.
    fullShiftX = 1 + (imgW - 1) * log(centerX) / log(imgW);
    fullShiftY = 1 + (imgH - 1) * log(centerY) / log(imgH);
    newCenterX = (1 - alpha)*centerX + alpha*fullShiftX;
    newCenterY = (1 - alpha)*centerY + alpha*fullShiftY;
    
    % Create a symmetric (circular) Gaussian blob for the shifted image.
    blobShifted = amplitude * exp( -((xGrid - newCenterX).^2 + (yGrid - newCenterY).^2) / (2*sigma^2) );
    imgShifted = imgShifted + blobShifted;
end

%% 2. Create Variants of the Shifted Image
% (A) Shifted (No Noise): already in imgShifted.

% (B) Shifted (With Noise)
noiseLevel = 0.1 * max(imgShifted(:));  % 10% of maximum intensity
imgShifted_noisy = imgShifted + noiseLevel * randn(size(imgShifted));

% (C) Shifted (+ Extra Peaks): add three extra Gaussian peaks.
imgShifted_extra = imgShifted;  % start with the basic shifted image
extraPeaks = zeros(3,3);  % each row: [cx, cy, sigmaExtra]
for i = 1:3
    cx = randi([10, imgW-10]);
    cy = randi([10, imgH-10]);
    ampExtra = 0.5 + rand();
    sigmaExtra = 3 + 7*rand();
    extraPeak = ampExtra * exp( -((xGrid - cx).^2 + (yGrid - cy).^2) / (2*sigmaExtra^2) );
    imgShifted_extra = imgShifted_extra + extraPeak;
    extraPeaks(i,:) = [cx, cy, sigmaExtra];
end

% (D) Shifted (Wider): render each shifted blob in an elliptical shape.
imgShifted_wider = zeros(imgH, imgW);
for k = 1:nBlobs
    centerX = centers(k,1);
    centerY = centers(k,2);
    amplitude = amplitudes(k);
    sigma_val = sigmas(k);
    % Compute shifted coordinates (same as before)
    fullShiftX = 1 + (imgW - 1) * log(centerX) / log(imgW);
    fullShiftY = 1 + (imgH - 1) * log(centerY) / log(imgH);
    newCenterX = (1 - alpha)*centerX + alpha*fullShiftX;
    newCenterY = (1 - alpha)*centerY + alpha*fullShiftY;
    % Create an elliptical Gaussian blob: widen in x-direction.
    sigma_x = 1.5 * sigma_val;  % 50% wider horizontally
    sigma_y = sigma_val;        % unchanged vertically
    blobShifted_wider = amplitude * exp( - ( ((xGrid - newCenterX).^2)/(2*sigma_x^2) + ((yGrid - newCenterY).^2)/(2*sigma_y^2) ) );
    imgShifted_wider = imgShifted_wider + blobShifted_wider;
end

%% 3. Frequency-Domain Alignment Procedure
% Compute the 2D FFT of the original image and center it.
F_original = fft2(imgOriginal);
F_orig_centered = fftshift(F_original);

% Define down-sampling parameters.
downsampleFactor = 4;  % For a 256x256 image, this gives a 64x64 block.
imgH_small = imgH / downsampleFactor;
imgW_small = imgW / downsampleFactor;
idxH = (imgH/2 - imgH_small/2 + 1):(imgH/2 + imgH_small/2);
idxW = (imgW/2 - imgW_small/2 + 1):(imgW/2 + imgW_small/2);

% Extract and vectorize the central block from the original FFT.
F_orig_small = F_orig_centered(idxH, idxW);
Fvec_orig_small = F_orig_small(:);

% Define a helper function for frequency-domain alignment.
% (See the local function frequencyAlign below.)
alignVariant = @(imgVar) frequencyAlign(imgVar, idxH, idxW, imgH, imgW, Fvec_orig_small);

% Perform alignment for each variant.
imgAligned_A = alignVariant(imgShifted);         % (A) No Noise
imgAligned_B = alignVariant(imgShifted_noisy);     % (B) With Noise
imgAligned_C = alignVariant(imgShifted_extra);     % (C) + Extra Peaks
imgAligned_D = alignVariant(imgShifted_wider);     % (D) Wider

%% 4. Display the Results Using Tiled Layout with Compact Spacing
t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Row 1: Shifted (No Noise)
nexttile; imagesc(imgOriginal); axis image; colorbar; title('Original');
nexttile; imagesc(imgShifted); axis image; colorbar; title('Shifted (No Noise)');
nexttile; imagesc(imgAligned_A); axis image; colorbar; title('Aligned (No Noise)');

% Row 2: Shifted (With Noise)
nexttile; imagesc(imgOriginal); axis image; colorbar; title('Original');
nexttile; imagesc(imgShifted_noisy); axis image; colorbar; title('Shifted (With Noise)');
nexttile; imagesc(imgAligned_B); axis image; colorbar; title('Aligned (With Noise)');

% Row 3: Shifted (+ Extra Peaks)
nexttile; imagesc(imgOriginal); axis image; colorbar; title('Original');
nexttile; 
imagesc(imgShifted_extra); 
axis image; colorbar; 
title('Shifted (+ Extra Peaks)');
hold on;
for i = 1:size(extraPeaks,1)
    cx = extraPeaks(i,1);
    cy = extraPeaks(i,2);
    sigmaExtra = extraPeaks(i,3);
    r = 3 * sigmaExtra;  % choose a radius (e.g., 3*sigmaExtra)
    rectangle('Position', [cx - r, cy - r, 2*r, 2*r], 'Curvature', [1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 2, 'Clipping', 'on');
end
hold off;
% Reset the x- and y-limits after drawing the circles.
ax = gca;
set(ax, 'XLim', [1 imgW], 'YLim', [1 imgH]);
nexttile; imagesc(imgAligned_C); axis image; colorbar; title('Aligned (+ Extra Peaks)');

% Row 4: Shifted (Wider)
nexttile; imagesc(imgOriginal); axis image; colorbar; title('Original');
nexttile; imagesc(imgShifted_wider); axis image; colorbar; title('Shifted (Wider)');
nexttile; imagesc(imgAligned_D); axis image; colorbar; title('Aligned (Wider)');

%% 5. Save the Final Figure as a PNG File with Minimal Whitespace
exportgraphics(gcf, 'aligned_analysis.png', 'Resolution', 300);
vecOriginal = imgOriginal(:);
% Compute cosine correlations for the aligned variants.
cosCorr = zeros(4,1);
cosCorr(1) = dot(imgAligned_A(:), vecOriginal) / (norm(imgAligned_A(:)) * norm(vecOriginal));
cosCorr(2) = dot(imgAligned_B(:), vecOriginal) / (norm(imgAligned_B(:)) * norm(vecOriginal));
cosCorr(3) = dot(imgAligned_C(:), vecOriginal) / (norm(imgAligned_C(:)) * norm(vecOriginal));
cosCorr(4) = dot(imgAligned_D(:), vecOriginal) / (norm(imgAligned_D(:)) * norm(vecOriginal));

% The cosine correlation of the original with itself is 1.
origCorr = 1;

% Prepare LaTeX table string including the original.
latexTable = sprintf(['\\begin{tabular}{l c}\n',...
                      '\\hline\n',...
                      'Variant & Cosine Correlation \\\\\n',...
                      '\\hline\n',...
                      'Original & %.4f \\\\\n',...
                      'No Noise & %.4f \\\\\n',...
                      'With Noise & %.4f \\\\\n',...
                      '+ Extra Peaks & %.4f \\\\\n',...
                      'Wider & %.4f \\\\\n',...
                      '\\hline\n',...
                      '\\end{tabular}\n'], origCorr, cosCorr(1), cosCorr(2), cosCorr(3), cosCorr(4));

% Display the LaTeX table in the Command Window.
disp('LaTeX Table of Cosine Correlation Coefficients:');
disp(latexTable);

filename = 'cosine_correlation_table.tex';
fid = fopen(filename, 'w');  % Open the file for writing

fprintf(fid, '%s', latexTable);  % Write the LaTeX table string to the file
fclose(fid);  % Close the file

%% Local Helper Function: frequencyAlign
function imgAligned = frequencyAlign(imgShiftedVariant, idxH, idxW, imgH, imgW, Fvec_orig_small)
    % Compute the 2D FFT of the shifted variant and center it.
    F_variant = fft2(imgShiftedVariant);
    F_variant_centered = fftshift(F_variant);
    % Extract the central block.
    F_variant_small = F_variant_centered(idxH, idxW);
    Fvec_variant_small = F_variant_small(:);
    % Form the outer product with the original's vectorized block.
    S_vec = Fvec_orig_small * Fvec_variant_small';
    % Solve the SVD-based Procrustes problem.
    [U, ~, V] = svd(S_vec);
    Q = U * V';
    % Apply the correction.
    Fvec_aligned_small = Q * Fvec_variant_small;
    % Reshape back to a 2D block.
    F_aligned_small = reshape(Fvec_aligned_small, [length(idxH), length(idxW)]);
    % Up-sample by embedding the aligned block into a full-size FFT matrix.
    F_aligned_full = zeros(imgH, imgW);
    F_aligned_full(idxH, idxW) = F_aligned_small;
    % Undo fftshift and compute the inverse FFT.
    F_aligned_full = ifftshift(F_aligned_full);
    imgAligned = abs(ifft2(F_aligned_full));
end
