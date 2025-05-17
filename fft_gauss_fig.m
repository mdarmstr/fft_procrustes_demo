%% Parameters
N      = 2048;                % number of sample points
x      = linspace(-5,5,N)';   % domain
sigma  = 1;                   % Gaussian width

%% Original data
y = exp(-x.^2/(2*sigma^2));

%% Compute FFT
Y  = fft(y);
f  = (0:N-1)'/N * (2*pi/(x(2)-x(1)));   % angular-frequency vector (unused for plotting)
mag = abs(Y);

%% Pick top 5 non-DC coefficients
%  - exclude the DC term at index 1
[~, idx_sorted] = sort(mag(2:end),'descend');
top_idx = idx_sorted(1:5) + 1;  % shift by +1 because we skipped DC

%% Build individual sinusoids and reconstruction
y_sin = zeros(N,5);
y_rec = zeros(N,1);
for k = 1:5
    iC = top_idx(k);
    
    % build a spectrum with only this coefficient and its complex conjugate
    Yk = zeros(N,1);
    Yk(iC)                    = Y(iC);
    Yk(mod(N-iC+2-1,N)+1)     = Y(mod(N-iC+2-1,N)+1);  % conjugate partner

    % get the real sinusoid
    yk = real(ifft(Yk));
    
    y_sin(:,k) = yk;
    y_rec      = y_rec + yk;
end

%% Plotting — publication quality
figure('Units','inches','Position',[1 1 6 4]);
hold on; box on;

% 1) Original Gaussian
h0 = plot(x, y, 'k-', 'LineWidth', 2);

% 2) Individual sinusoids
colors = lines(5);
for k = 1:5
    hk(k) = plot(x, y_sin(:,k), '--', 'LineWidth', 1, 'Color', colors(k,:));
end

% 3) Reconstruction
hr = plot(x, y_rec, 'r-', 'LineWidth', 2);

% Axes and labels
xlabel('x',       'FontName','Arial','FontSize',14);
ylabel('f(x)',    'FontName','Arial','FontSize',14);
set(gca, ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'LineWidth',1, ...
    'TickDir','out' ...
);

% Collect all handles into one vector
allHandles = [h0, hk, hr];

% Build a cell array of names matching the handles
legNames = cell(1, numel(allHandles));
legNames{1} = 'Original';
for k = 1:5
    legNames{1 + k} = sprintf('Coefficient %d', k);
end
legNames{end} = 'Reconstruction';

% Now call legend with the handles and the matching names cell array
leg = legend(allHandles, legNames, 'Location', 'northeast');
set(leg, 'FontName', 'Arial', 'FontSize', 10);
title('Reconstruction of a Gaussian with FFT coefficients')

% (Optional) set up for vector‑graphic export
set(gcf, 'PaperPositionMode','auto');
exportgraphics(gcf, 'gauss.pdf', 'Resolution', 500);


hold off;
