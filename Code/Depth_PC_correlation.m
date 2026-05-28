clc; clear; close all;

% =========================================================
% Quantify association between sequencing depth and PCA axes
% Correlation between original library size and PCA score vectors
% before and after normalization
%
% IMPORTANT:
% Depth is always computed from the original RAW count table.
% Do NOT recompute depth from normalized matrices.
%
% PCA scores are obtained using scoresPca().
% =========================================================

% -----------------------------
% USER SETTINGS
% -----------------------------
dataPath = "Phylum/True";   % change to "Genus/True" when needed
outDir   = fullfile(dataPath, "Depth_PC_Correlation");

R = 100;                    % use 1 for first replicate only, 100 for all replicates
useSpearman = true;         % also compute Spearman correlation

% Create output folder
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% -----------------------------
% Find method folders
% -----------------------------
allEntries = dir(dataPath);
methodDirs = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));

% IMPORTANT FIX:
% Exclude output folder so it is not treated as a normalization method
methodDirs = methodDirs(~strcmp({methodDirs.name}, "Depth_PC_Correlation"));

if isempty(methodDirs)
    error('No method folders found in %s', dataPath);
end

% Sort method folders by name for reproducibility
methodDirs = sortFilesByName(methodDirs);

% -----------------------------
% Identify RAW folder
% -----------------------------
rawIdx = [];

for d = 1:length(methodDirs)
    folderName = methodDirs(d).name;

    if contains(lower(folderName), 'raw')
        rawIdx = d;
        break;
    end
end

if isempty(rawIdx)
    error('Could not identify RAW data folder. Folder name must contain "raw".');
end

rawFolder = fullfile(dataPath, methodDirs(rawIdx).name);
rawFiles  = dir(fullfile(rawFolder, '*.csv'));
rawFiles  = sortFilesByName(rawFiles);

if isempty(rawFiles)
    error('No CSV files found in RAW folder: %s', rawFolder);
end

fprintf('RAW folder detected: %s\n', methodDirs(rawIdx).name);

% -----------------------------
% Storage table
% -----------------------------
Result = table();

% =========================================================
% Main loop across replicates
% =========================================================
nRep = min(R, length(rawFiles));

for r = 1:nRep

    fprintf('\nProcessing replicate %d/%d\n', r, nRep);

    % -----------------------------------------------------
    % Read RAW data for this replicate
    % -----------------------------------------------------
    rawFile = fullfile(rawFolder, rawFiles(r).name);
    Xraw = readmatrix(rawFile);

    % -----------------------------------------------------
    % Create group labels
    % Assumes 300 samples:
    % 50 Control + 50 Case repeated over 3 blocks.
    % Modify if your design differs.
    % -----------------------------------------------------
    y = repmat([ones(50,1)*-1; ones(50,1)], 3, 1);

    % Remove last column from RAW only if it is a group column
    Xraw = removeGroupColumnIfPresent(Xraw, y);

    if size(Xraw,1) ~= numel(y)
        error('RAW replicate %d: Xraw has %d rows but y has %d labels.', ...
            r, size(Xraw,1), numel(y));
    end

    if any(isnan(Xraw(:))) || any(isinf(Xraw(:)))
        error('RAW replicate %d contains NaN or Inf.', r);
    end

    % -----------------------------------------------------
    % Original sequencing depth / library size
    % This depth vector must match the rows of Xraw and X.
    % -----------------------------------------------------
    depthRaw = sum(Xraw, 2);

    if any(depthRaw <= 0)
        warning('Replicate %d has samples with zero library size.', r);
    end

    % =====================================================
    % Loop across all normalization methods
    % =====================================================
    for d = 1:length(methodDirs)

        subfolderName = methodDirs(d).name;
        methodParts   = split(subfolderName, "_");

        if numel(methodParts) >= 5
            methodName = string(methodParts{5});
        else
            methodName = string(subfolderName);
        end

        methodFullPath = fullfile(dataPath, subfolderName);
        csvFiles = dir(fullfile(methodFullPath, '*.csv'));
        csvFiles = sortFilesByName(csvFiles);

        if isempty(csvFiles)
            fprintf('No CSV files found in: %s\n', methodFullPath);
            continue;
        end

        if r > length(csvFiles)
            fprintf('Method %s has fewer files than replicate %d. Skipping.\n', methodName, r);
            continue;
        end

        currentFile = fullfile(methodFullPath, csvFiles(r).name);
        X = readmatrix(currentFile);

        % Remove group column if present
        X = removeGroupColumnIfPresent(X, y);

        % Basic checks
        if size(X,1) ~= numel(depthRaw)
            error('Method %s replicate %d: X has %d rows but raw depth has %d values.', ...
                methodName, r, size(X,1), numel(depthRaw));
        end

        if any(isnan(X(:))) || any(isinf(X(:)))
            warning('Method %s replicate %d contains NaN or Inf. Skipping.', methodName, r);
            continue;
        end

        if rank(X) < 2
            warning('Method %s replicate %d has rank < 2. Skipping.', methodName, r);
            continue;
        end

        % -------------------------------------------------
        % PCA scores using scoresPca()
        % Preprocessing = 1 means mean-centering.
        % This returns the PCA score matrix directly.
        % -------------------------------------------------
        try
            [score, ~] = scoresPca(X, ...
                'PCs', 1:2, ...
                'Preprocessing', 1, ...
                'PlotType', 'Scatter');

            close all;  % close figure automatically created by scoresPca

        catch ME
            warning('scoresPca failed for method %s replicate %d: %s', ...
                methodName, r, ME.message);
            continue;
        end

        if size(score,2) < 2
            warning('Method %s replicate %d produced fewer than 2 PCs. Skipping.', methodName, r);
            continue;
        end

        % -------------------------------------------------
        % Compute variance explained consistently with scoresPca()
        % scoresPca internally mean-centers using preprocess2D.
        % -------------------------------------------------
        [Xcs, ~, ~] = preprocess2D(X, 'Preprocessing', 1);

        totalVar = sum(sum(Xcs.^2));

        if totalVar <= 0
            warning('Method %s replicate %d has zero total variance. Skipping.', methodName, r);
            continue;
        end

        PC1explained = 100 * sum(score(:,1).^2) / totalVar;
        PC2explained = 100 * sum(score(:,2).^2) / totalVar;

        % -------------------------------------------------
        % Pearson correlations:
        % original sequencing depth vs PCA score vectors
        %
        % This is the key computation:
        % correlation between the first score vector and sequencing depth
        % -------------------------------------------------
        [rPC1, pPC1] = corr(depthRaw, score(:,1), ...
            'Type', 'Pearson', ...
            'Rows', 'complete');

        [rPC2, pPC2] = corr(depthRaw, score(:,2), ...
            'Type', 'Pearson', ...
            'Rows', 'complete');

        % -------------------------------------------------
        % Spearman correlations
        % -------------------------------------------------
        if useSpearman
            [rhoPC1, psPC1] = corr(depthRaw, score(:,1), ...
                'Type', 'Spearman', ...
                'Rows', 'complete');

            [rhoPC2, psPC2] = corr(depthRaw, score(:,2), ...
                'Type', 'Spearman', ...
                'Rows', 'complete');
        else
            rhoPC1 = NaN; psPC1 = NaN;
            rhoPC2 = NaN; psPC2 = NaN;
        end

        % -------------------------------------------------
        % Store results
        % Signed correlations are stored, but manuscript should report
        % absolute correlations because PCA signs are arbitrary.
        % -------------------------------------------------
        newRow = table();

        newRow.Replicate = r;
        newRow.Method = methodName;

        newRow.PC1_Explained = PC1explained;
        newRow.PC2_Explained = PC2explained;

        newRow.Pearson_r_Depth_PC1 = rPC1;
        newRow.Pearson_p_Depth_PC1 = pPC1;

        newRow.Pearson_r_Depth_PC2 = rPC2;
        newRow.Pearson_p_Depth_PC2 = pPC2;

        newRow.Abs_Pearson_r_Depth_PC1 = abs(rPC1);
        newRow.Abs_Pearson_r_Depth_PC2 = abs(rPC2);

        newRow.Spearman_rho_Depth_PC1 = rhoPC1;
        newRow.Spearman_p_Depth_PC1 = psPC1;

        newRow.Spearman_rho_Depth_PC2 = rhoPC2;
        newRow.Spearman_p_Depth_PC2 = psPC2;

        newRow.Abs_Spearman_rho_Depth_PC1 = abs(rhoPC1);
        newRow.Abs_Spearman_rho_Depth_PC2 = abs(rhoPC2);

        Result = [Result; newRow];

    end
end

% =========================================================
% Save replicate-level results
% =========================================================
replicateFile = fullfile(outDir, "depth_pc_correlation_replicates.csv");
writetable(Result, replicateFile);

fprintf('\nSaved replicate-level results:\n%s\n', replicateFile);

% =========================================================
% Summarize across replicates
% =========================================================
methods = unique(Result.Method, 'stable');

Summary = table();

for i = 1:numel(methods)

    idx = Result.Method == methods(i);

    newRow = table();
    newRow.Method = methods(i);

    newRow.Mean_PC1_Explained = mean(Result.PC1_Explained(idx), 'omitnan');
    newRow.SD_PC1_Explained   = std(Result.PC1_Explained(idx), 'omitnan');

    newRow.Mean_PC2_Explained = mean(Result.PC2_Explained(idx), 'omitnan');
    newRow.SD_PC2_Explained   = std(Result.PC2_Explained(idx), 'omitnan');

    % Signed Pearson correlations
    newRow.Mean_r_Depth_PC1 = mean(Result.Pearson_r_Depth_PC1(idx), 'omitnan');
    newRow.SD_r_Depth_PC1   = std(Result.Pearson_r_Depth_PC1(idx), 'omitnan');

    newRow.Mean_r_Depth_PC2 = mean(Result.Pearson_r_Depth_PC2(idx), 'omitnan');
    newRow.SD_r_Depth_PC2   = std(Result.Pearson_r_Depth_PC2(idx), 'omitnan');

    % Absolute Pearson correlations
    newRow.Mean_abs_r_Depth_PC1 = mean(Result.Abs_Pearson_r_Depth_PC1(idx), 'omitnan');
    newRow.SD_abs_r_Depth_PC1   = std(Result.Abs_Pearson_r_Depth_PC1(idx), 'omitnan');

    newRow.Mean_abs_r_Depth_PC2 = mean(Result.Abs_Pearson_r_Depth_PC2(idx), 'omitnan');
    newRow.SD_abs_r_Depth_PC2   = std(Result.Abs_Pearson_r_Depth_PC2(idx), 'omitnan');

    % Signed Spearman correlations
    newRow.Mean_rho_Depth_PC1 = mean(Result.Spearman_rho_Depth_PC1(idx), 'omitnan');
    newRow.SD_rho_Depth_PC1   = std(Result.Spearman_rho_Depth_PC1(idx), 'omitnan');

    newRow.Mean_rho_Depth_PC2 = mean(Result.Spearman_rho_Depth_PC2(idx), 'omitnan');
    newRow.SD_rho_Depth_PC2   = std(Result.Spearman_rho_Depth_PC2(idx), 'omitnan');

    % Absolute Spearman correlations
    newRow.Mean_abs_rho_Depth_PC1 = mean(Result.Abs_Spearman_rho_Depth_PC1(idx), 'omitnan');
    newRow.SD_abs_rho_Depth_PC1   = std(Result.Abs_Spearman_rho_Depth_PC1(idx), 'omitnan');

    newRow.Mean_abs_rho_Depth_PC2 = mean(Result.Abs_Spearman_rho_Depth_PC2(idx), 'omitnan');
    newRow.SD_abs_rho_Depth_PC2   = std(Result.Abs_Spearman_rho_Depth_PC2(idx), 'omitnan');

    Summary = [Summary; newRow];

end

summaryFile = fullfile(outDir, "depth_pc_correlation_summary.csv");
writetable(Summary, summaryFile);

fprintf('Saved summary results:\n%s\n', summaryFile);

% =========================================================
% Display summary
% =========================================================
disp(Summary);

% =========================================================
% Print the key RAW result for quick checking
% =========================================================
rawRow = contains(lower(Summary.Method), "raw");

if any(rawRow)
    fprintf('\nKey check for RAW data:\n');

    fprintf('PC1 explained: %.2f ± %.2f%%\n', ...
        Summary.Mean_PC1_Explained(rawRow), ...
        Summary.SD_PC1_Explained(rawRow));

    fprintf('PC2 explained: %.2f ± %.2f%%\n', ...
        Summary.Mean_PC2_Explained(rawRow), ...
        Summary.SD_PC2_Explained(rawRow));

    fprintf('|Pearson r| depth-PC1: %.3f ± %.3f\n', ...
        Summary.Mean_abs_r_Depth_PC1(rawRow), ...
        Summary.SD_abs_r_Depth_PC1(rawRow));

    fprintf('|Pearson r| depth-PC2: %.3f ± %.3f\n', ...
        Summary.Mean_abs_r_Depth_PC2(rawRow), ...
        Summary.SD_abs_r_Depth_PC2(rawRow));

    fprintf('|Spearman rho| depth-PC1: %.3f ± %.3f\n', ...
        Summary.Mean_abs_rho_Depth_PC1(rawRow), ...
        Summary.SD_abs_rho_Depth_PC1(rawRow));

    fprintf('|Spearman rho| depth-PC2: %.3f ± %.3f\n', ...
        Summary.Mean_abs_rho_Depth_PC2(rawRow), ...
        Summary.SD_abs_rho_Depth_PC2(rawRow));
else
    warning('No RAW row found in Summary. Check method folder naming.');
end

fprintf('\nDone.\n');

% =========================================================
% Helper functions
% =========================================================

function X = removeGroupColumnIfPresent(X, y)

    if size(X,1) == numel(y) && size(X,2) > 1

        lastCol = X(:,end);

        if isequal(lastCol, y) || isequal(lastCol, -y)
            X(:,end) = [];
        end

    end

end

function files = sortFilesByName(files)

    if isempty(files)
        return;
    end

    names = lower({files.name});
    [~, idx] = sort(names);
    files = files(idx);

end