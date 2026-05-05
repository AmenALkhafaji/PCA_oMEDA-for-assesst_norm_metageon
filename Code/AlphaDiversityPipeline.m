
clear; close all; clc;

%% ===================== USER SETTINGS =====================
dataPath = "Genus/unblance";
R = 100;
rng(1);

alpha_sig = 0.05;

% synthetic infection_status (alternating blocks)
N = 300;
blockSize = 50;
infection_status = ones(N,1);
blockID = floor((0:N-1)'/blockSize);
infection_status(mod(blockID,2)==0) = -1;

outDir = "results_1";
if ~exist(outDir,"dir"), mkdir(outDir); end
safeName = replace(dataPath, "/", "_");

%% ===================== LIST METHODS =====================
allEntries = dir(dataPath);
methodDirs = allEntries([allEntries.isdir] & ~startsWith({allEntries.name}, '.'));

if isempty(methodDirs)
    error('No method subfolders found under %s', dataPath);
end

%% ===================== STORAGE =====================
summaryRows = {};

fprintf('Starting Shannon + Kruskal-Wallis summary analysis...\n');

%% ===================== MAIN LOOP =====================
for d = 1:numel(methodDirs)
    subfolderName = methodDirs(d).name;

    parts = split(string(subfolderName), "_");
    if numel(parts) >= 5
        methodName = parts(5);
    else
        methodName = string(subfolderName);
    end
    methodName = string(methodName);

    methodFullPath = fullfile(dataPath, subfolderName);
    csvFiles = dir(fullfile(methodFullPath, '*.csv'));

    if numel(csvFiles) < R
        warning('Method %s has only %d files. Skipping.', methodName, numel(csvFiles));
        continue;
    end

    % stable order
    [~, sortIdx] = sort({csvFiles.name});
    csvFiles = csvFiles(sortIdx);

    isCLR = contains(lower(methodName), "clr") || contains(lower(subfolderName), "clr");
    fprintf('Processing method: [%s] ', methodName);

    % Shannon per replicate
    lowH  = nan(R,1); lowU  = nan(R,1);
    highH = nan(R,1); highU = nan(R,1);
    fullH = nan(R,1); fullU = nan(R,1);

    % KW p-values per replicate
    pLow  = nan(R,1);
    pHigh = nan(R,1);
    pFull = nan(R,1);

    for r = 1:R
        currentFile = fullfile(methodFullPath, csvFiles(r).name);
        X = readmatrix(currentFile);
        y = infection_status;

        % remove label column if present
        if size(X,2) >= 1 && size(X,1) == numel(y) && all(X(:,end) == y)
            X(:, end) = [];
        end

        idxLow  = 1:100;
        idxHigh = 201:300;
        idxFull = 1:size(X,1);

        Xlow  = to_composition_for_shannon(X(idxLow,:),  isCLR);
        Xhigh = to_composition_for_shannon(X(idxHigh,:), isCLR);
        Xfull = to_composition_for_shannon(X(idxFull,:), isCLR);

        shLow  = shannon_index_rows(Xlow);
        shHigh = shannon_index_rows(Xhigh);
        shFull = shannon_index_rows(Xfull);

        yLow  = y(idxLow);
        yHigh = y(idxHigh);
        yFull = y(idxFull);

        % mean Shannon per replicate within each group
        lowH(r)  = mean(shLow(yLow==-1),   'omitnan');
        lowU(r)  = mean(shLow(yLow==1),   'omitnan');

        highH(r) = mean(shHigh(yHigh==-1), 'omitnan');
        highU(r) = mean(shHigh(yHigh==1), 'omitnan');

        fullH(r) = mean(shFull(yFull==-1), 'omitnan');
        fullU(r) = mean(shFull(yFull==1), 'omitnan');

        % KW p-value per replicate
        pLow(r)  = kruskalwallis(shLow,  yLow,  'off');
        pHigh(r) = kruskalwallis(shHigh, yHigh, 'off');
        pFull(r) = kruskalwallis(shFull, yFull, 'off');

        if mod(r,25)==0, fprintf('.'); end
    end
    fprintf(' Done.\n');

    %% ===================== SHANNON SUMMARY =====================
    lowH_str  = sprintf('%.4f (%.4f)', mean(lowH,'omitnan'),  std(lowH,'omitnan'));
    lowU_str  = sprintf('%.4f (%.4f)', mean(lowU,'omitnan'),  std(lowU,'omitnan'));

    highH_str = sprintf('%.4f (%.4f)', mean(highH,'omitnan'), std(highH,'omitnan'));
    highU_str = sprintf('%.4f (%.4f)', mean(highU,'omitnan'), std(highU,'omitnan'));

    fullH_str = sprintf('%.4f (%.4f)', mean(fullH,'omitnan'), std(fullH,'omitnan'));
    fullU_str = sprintf('%.4f (%.4f)', mean(fullU,'omitnan'), std(fullU,'omitnan'));

    %% ===================== KW COUNTS =====================
    low_sig_count    = sum(pLow  < alpha_sig, 'omitnan');
    low_nonsig_count = sum(pLow  >= alpha_sig, 'omitnan');

    high_sig_count    = sum(pHigh < alpha_sig, 'omitnan');
    high_nonsig_count = sum(pHigh >= alpha_sig, 'omitnan');

    full_sig_count    = sum(pFull < alpha_sig, 'omitnan');
    full_nonsig_count = sum(pFull >= alpha_sig, 'omitnan');

    %% ===================== FINAL COMBINED ROW =====================
    summaryRows(end+1,1:13) = { ...
        char(methodName), ...
        lowH_str,  lowU_str,  low_sig_count,  low_nonsig_count, ...
        highH_str, highU_str, high_sig_count, high_nonsig_count, ...
        fullH_str, fullU_str, full_sig_count, full_nonsig_count};
end

if isempty(summaryRows)
    error('No methods passed the file-count filter (>=R).');
end

%% ===================== BUILD FINAL TABLE =====================
SummaryTable = cell2table(summaryRows, 'VariableNames', ...
    {'Normalization_Method', ...
     'control shannon_mean(sd) at low depth',  'Low case shannon mean( sd)',  'Counts of sig p_value at low depth',  'Counts of nonsig p_value at low depth', ...
     'control shannon mean( sd) at high depth', 'case shannon mean( sd) high', 'Counts of sig p_value at high depth', 'Counts of nonsig p_value at high depth', ...
     'control shannon mean( sd) at full depth', 'case shannon mean( sd) at full depth', 'Count of sig p_value at full depth', 'Counts of nonsig p_value at full depth'});

%% ===================== EXPORT TABLE =====================
summaryOut = fullfile(outDir, "Shannon_KW_summary_" + safeName + ".csv");
writetable(SummaryTable, summaryOut);

fprintf('Done. Final summary table saved in %s\n', outDir);

%% =================== Local functions ===================

function Xp = to_composition_for_shannon(X, isCLR)
    if isCLR
        Xp = exp(X);
        Xp(~isfinite(Xp)) = 0;
        Xp(Xp < 0) = 0;
    else
        Xp = X;
        Xp(Xp < 0) = 0;
        Xp(~isfinite(Xp)) = 0;
    end

    rs = sum(Xp, 2);
    rs(rs == 0) = 1;
    Xp = Xp ./ rs;
end

function H = shannon_index_rows(Xp)
    H = nan(size(Xp,1),1);
    for i = 1:size(Xp,1)
        p = Xp(i,:);
        p = p(p > 0);
        if isempty(p)
            H(i) = 0;
        else
            H(i) = -sum(p .* log(p));
        end
    end
end