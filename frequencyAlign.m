function imgAligned = frequencyAlign(imgShiftedVariant, idxH, idxW, imgH, imgW, Fvec_orig_small)
    % frequencyAlign  Aligns imgShiftedVariant to the reference frequency block
    %
    % Inputs:
    %   imgShiftedVariant   – shifted image to align
    %   idxH, idxW          – vectors of row/col indices for the central block
    %   imgH, imgW          – full image dimensions
    %   Fvec_orig_small     – vectorized reference FFT block (length m)
    %
    % Output:
    %   imgAligned          – the aligned image in the spatial domain

    % 1) Compute centered FFT of the variant
    F_variant = fft2(imgShiftedVariant);
    F_variant_centered = fftshift(F_variant);

    % 2) Extract and vectorize the small central block
    F_variant_small    = F_variant_centered(idxH, idxW);
    Fvec_variant_small = F_variant_small(:);

    % 3) Rank-1 Procrustes via Gram trick (no outer product)
    a = Fvec_orig_small;   % target vector, m×1
    b = Fvec_variant_small; % moving vector, m×1

    % norms
    na = norm(a);
    nb = norm(b);

    % align b into the direction of a:
    %   u = a/na,  vTb = (b' * b) / nb
    %   Fvec_aligned_small = u * (v^T b)
    Fvec_aligned_small = (a / na) * ((b' * b) / nb);

    % 4) Reshape aligned vector back into the small 2D block
    F_aligned_small = reshape(Fvec_aligned_small, [numel(idxH), numel(idxW)]);

    % 5) Embed into the full-size FFT matrix
    F_aligned_full = zeros(imgH, imgW);
    F_aligned_full(idxH, idxW) = F_aligned_small;

    % 6) Invert the centering and do the inverse FFT
    F_aligned_full = ifftshift(F_aligned_full);
    imgAligned = abs(ifft2(F_aligned_full));
end
