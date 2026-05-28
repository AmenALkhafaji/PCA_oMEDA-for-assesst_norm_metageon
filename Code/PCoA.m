clc; clear; close all;

% =========================================================
% PCoA visualization for phylum-level True raw data
% Raw data only
% Euclidean distance + classical PCoA
% =========================================================

dataPath = "Genus/True";

outDirPDF = fullfile(dataPath, 'Figures_PCoA_Euclidean_PDF');
outDirPNG = fullfile(dataPath, 'Figures_PCoA_Euclidean_PNG');

if ~exist(outDirPDF, 'dir')
    mkdir(outDirPDF);
end

if ~exist(outDirPNG, 'dir')
    mkdir(outDirPNG);
end

allEntries = dir(dataPath);
methodDirs = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));

R = 1;

for d = 1:length(methodDirs)

    subfolderName = methodDirs(d).name;
    methodParts   = split(subfolderName, "_");

    if numel(methodParts) >= 5
        methodName = methodParts{5};
    else
        methodName = subfolderName;
    end

    methodSaveName = regexprep(methodName, '[^\w]', '_');
    methodFullPath = fullfile(dataPath, subfolderName);

    % -----------------------------------------------------
    % Keep only RAW data folder
    % Adjust this condition if your raw folder name differs.
    % -----------------------------------------------------
    if ~contains(lower(methodName), 'raw') && ~contains(lower(subfolderName), 'raw')
        continue;
    end

    csvFiles = dir(fullfile(methodFullPath, '*.csv'));

    if isempty(csvFiles)
        fprintf('No CSV files found in: %s\n', methodFullPath);
        continue;
    end

    fprintf('Processing Euclidean PCoA for method: %s\n', methodName);

    for r = 1:min(R, length(csvFiles))

        currentFile = fullfile(methodFullPath, csvFiles(r).name);
        X = readmatrix(currentFile);

        % -------------------------------------------------
        % Group labels
        % Assumes 300 samples:
        % 50 Control + 50 Case repeated over 3 blocks.
        % Modify this if your design differs.
        % -------------------------------------------------
        y = repmat([ones(50,1)*-1; ones(50,1)], 3, 1);

        % -------------------------------------------------
        % Remove last column only if it exactly matches group coding
        % -------------------------------------------------
        if size(X,1) == numel(y) && size(X,2) > 1

            lastCol = X(:,end);

            if isequal(lastCol, y) || isequal(lastCol, -y)
                X(:,end) = [];
            end

        end

        % -------------------------------------------------
        % Checks
        % -------------------------------------------------
        if size(X,1) ~= numel(y)
            error('Mismatch: X has %d rows but y has %d labels.', ...
                size(X,1), numel(y));
        end

        if any(isnan(X(:))) || any(isinf(X(:)))
            error('X contains NaN or Inf.');
        end

        if rank(X) < 2
            error('X has rank < 2; Euclidean PCoA will be degenerate.');
        end

        % -------------------------------------------------
        % Metadata
        % -------------------------------------------------
        depth = sum(X, 2);

        % -------------------------------------------------
        % PCoA mathematics
        % Euclidean sample distance followed by classical PCoA
        % -------------------------------------------------
        D = euclideanDistance(X);

        [coords, eigvals, explained] = classicalPCoA(D);

        if size(coords,2) < 2
            error('PCoA produced fewer than two positive dimensions.');
        end

        % -------------------------------------------------
        % Final combined figure
        % Panel A = colored by sequencing depth
        % Panel B = colored by group
        % -------------------------------------------------
        comboFig = figure('Color','w', ...
            'Units','pixels', ...
            'Position',[100 100 1200 520]);

        tl = tiledlayout(comboFig, 1, 2, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        % -----------------------------
        % Panel A: sequencing depth
        % -----------------------------
        ax1 = nexttile(tl, 1);

        scatter(ax1, coords(:,1), coords(:,2), ...
            38, depth, 'filled');

        box(ax1, 'on');
        axis(ax1, 'square');

        xlabel(ax1, sprintf('PCoA 1 (%.1f%%)', explained(1)));
        ylabel(ax1, sprintf('PCoA 2 (%.1f%%)', explained(2)));

        cb = colorbar(ax1);
        ylabel(cb, 'Sequencing depth');

        title(ax1, '');
        addPanelLabel(ax1, 'A');

        % -----------------------------
        % Panel B: group
        % -----------------------------
        ax2 = nexttile(tl, 2);
        hold(ax2, 'on');

        idxControl = y == -1;
        idxCase    = y == 1;

        scatter(ax2, coords(idxControl,1), coords(idxControl,2), ...
            38, 'filled');

        scatter(ax2, coords(idxCase,1), coords(idxCase,2), ...
            38, 'filled');

        box(ax2, 'on');
        axis(ax2, 'square');

        xlabel(ax2, sprintf('PCoA 1 (%.1f%%)', explained(1)));
        ylabel(ax2, sprintf('PCoA 2 (%.1f%%)', explained(2)));

        legend(ax2, {'Control', 'Case'}, ...
            'Location', 'best');

        title(ax2, '');
        addPanelLabel(ax2, 'B');

        % -------------------------------------------------
        % Export final figure as high-quality vector PDF
        % -------------------------------------------------
        pdfFile = fullfile(outDirPDF, ...
            sprintf('%s_rep%d_PCoA_Euclidean.pdf', methodSaveName, r));

        exportgraphics(comboFig, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');

        % -------------------------------------------------
        % Optional: export high-resolution PNG
        % -------------------------------------------------
        pngFile = fullfile(outDirPNG, ...
            sprintf('%s_rep%d_PCoA_Euclidean.png', methodSaveName, r));

        exportgraphics(comboFig, pngFile, ...
            'Resolution', 600, ...
            'BackgroundColor', 'white');

        close(comboFig);

    end
end

fprintf('Done.\n');

% =========================================================
% Helper functions
% =========================================================

function D = euclideanDistance(X)

    % Euclidean distance between samples:
    % d_ij = sqrt(sum((x_i - x_j).^2))

    n = size(X,1);
    D = zeros(n,n);

    for i = 1:n

        xi = X(i,:);

        for j = i+1:n

            xj = X(j,:);

            dij = sqrt(sum((xi - xj).^2));

            D(i,j) = dij;
            D(j,i) = dij;

        end
    end

end

function [coords, eigvals, explained] = classicalPCoA(D)

    % Classical PCoA from a sample-by-sample distance matrix.

    n = size(D,1);

    if size(D,2) ~= n
        error('Distance matrix must be square.');
    end

    % Correct symmetry check.
    % Do NOT use D(:) - D(:)' because that creates a huge outer matrix.
    if max(max(abs(D - D'))) > 1e-10
        error('Distance matrix must be symmetric.');
    end

    if any(abs(diag(D)) > 1e-10)
        warning('Distance matrix diagonal is not exactly zero. Forcing diagonal to zero.');
        D(1:n+1:end) = 0;
    end

    % Squared distances
    D2 = D.^2;

    % Double-centering matrix
    J = eye(n) - ones(n,n) / n;

    % Principal coordinate matrix
    B = -0.5 * J * D2 * J;

    % Numerical symmetrization
    B = (B + B') / 2;

    % Eigen decomposition
    [V, L] = eig(B);

    eigvals = diag(L);

    % Sort eigenvalues descending
    [eigvals, idx] = sort(eigvals, 'descend');
    V = V(:,idx);

    % Keep positive eigenvalues only
    positiveIdx = eigvals > 1e-10;

    eigvalsPos = eigvals(positiveIdx);
    Vpos       = V(:,positiveIdx);

    if isempty(eigvalsPos)
        error('No positive eigenvalues found. PCoA coordinates cannot be computed.');
    end

    coords = Vpos .* sqrt(eigvalsPos)';

    % Percentage explained by positive axes only
    explained = 100 * eigvalsPos / sum(eigvalsPos);

    eigvals = eigvalsPos;

end

function addPanelLabel(ax, txt)

    text(ax, 0.02, 0.98, txt, ...
        'Units', 'normalized', ...
        'FontWeight', 'bold', ...
        'FontSize', 14, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'none', ...
        'Clipping', 'off');

end