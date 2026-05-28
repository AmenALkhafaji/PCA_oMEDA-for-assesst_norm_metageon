% clc; clear; close all;
% 
% load('genus.mat');   % must contain lab

dataPath = "Genus/True";

lab = {
    'Bifidobacterium'
    'Butyricimonas'
    'Odoribacter'
    'Paraprevotella'
    'Bacteroides'
    'Parabacteroides'
    'Prevotella'
    'unknown1'
    'unknown2'
    'unknown3'
    'Clostridium'
    'unknown4'
    'unknown5'
    'Ruminococcus'
    'Anaerostipes'
    'Blautia'
    'Coprococcus'
    'Dorea'
    'Lachnospira'
    'Roseburia'
    'unknown6'
    'unknown7'
    'unknown8'
    'unknown9'
    'Faecalibacterium'
    'Oscillospira'
    'Ruminococcus1'
    'Dialister'
    'unknown10'
    'unknown11'
    'unknown12'
    'Sutterella'
    'Escherichia'
    'Klebsiella'
    'Haemophilus'
    'Akkermansia'
};
lab = cellstr(string(lab(:)));
outDirPDF = fullfile(dataPath, 'Figures_PDF');
outDirPNG = fullfile(dataPath, 'Figures_PNG');

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

    if contains(lower(methodName), 'edger')
        pcsToUse = [1 3];
    else
        pcsToUse = 1:2;
    end

    csvFiles = dir(fullfile(methodFullPath, '*.csv'));

    if isempty(csvFiles)
        fprintf('No CSV files found in: %s\n', methodFullPath);
        continue;
    end

    fprintf('Processing method: %s\n', methodName);

    for r = 1:min(R, length(csvFiles))

        currentFile = fullfile(methodFullPath, csvFiles(r).name);
        X = readmatrix(currentFile);

        % Group labels
        y = repmat([ones(50,1)*-1; ones(50,1)], 3, 1);

        % Remove last column ONLY if it exactly matches group coding
        if size(X,1) == numel(y) && size(X,2) > 1
            lastCol = X(:,end);
            if isequal(lastCol, y) || isequal(lastCol, -y)
                X(:,end) = [];
            end
        end

        % Extra safeguard
        if size(X,2) ~= numel(lab)
            if size(X,2) == numel(lab) + 1
                X(:,end) = [];
            else
                error('Mismatch: X has %d columns but lab has %d labels.', ...
                    size(X,2), numel(lab));
            end
        end

        if any(isnan(X(:))) || any(isinf(X(:)))
            error('X contains NaN or Inf.');
        end

        if rank(X) < 2
            error('X has rank < 2, PCA score plot will be degenerate.');
        end

        s = cellstr(categorical(y, [-1, 1], ["Control", "Case"]));
        mylab = repmat({''}, size(X,1), 1);

        % -----------------------------
        % Scores by class
        % -----------------------------
        figBefore = findall(0, 'Type', 'figure');

        scoresPca(X, ...
            'PCs', pcsToUse, ...
            'Preprocessing', 1, ...
            'ObsClass', s, ...
            'ObsLabel', mylab);

        drawnow;

        figAfter = findall(0, 'Type', 'figure');
        fig1 = setdiff(figAfter, figBefore);

        if isempty(fig1)
            error('scoresPca did not create a new figure for class plot.');
        end

        fig1 = fig1(1);
        ax1 = findMainAxes(fig1);

        if isempty(ax1)
            error('No plotted axes found in class figure.');
        end

        addClassLegend(ax1, {'Control','Case'});

        % -----------------------------
        % Scores by mean
        % -----------------------------
        figBefore = findall(0, 'Type', 'figure');

        scoresPca(X, ...
            'PCs', pcsToUse, ...
            'Preprocessing', 1, ...
            'ObsClass', mean(X,2), ...
            'ObsLabel', mylab);

        drawnow;

        figAfter = findall(0, 'Type', 'figure');
        fig2 = setdiff(figAfter, figBefore);

        if isempty(fig2)
            error('scoresPca did not create a new figure for mean plot.');
        end

        fig2 = fig2(1);
        ax2 = findMainAxes(fig2);

        if isempty(ax2)
            error('No plotted axes found in mean figure.');
        end

        % -----------------------------
        % Loadings
        % -----------------------------
        figBefore = findall(0, 'Type', 'figure');

        loadingsPca(X, ...
            'PCs', pcsToUse, ...
            'VarsLabel', lab, ...
            'Preprocessing', 1, ...
            'VarsClass', lab);

        drawnow;

        figAfter = findall(0, 'Type', 'figure');
        fig3 = setdiff(figAfter, figBefore);

        if isempty(fig3)
            error('loadingsPca did not create a new figure.');
        end

        fig3 = fig3(1);
        ax3 = findMainAxes(fig3);

        if isempty(ax3)
            error('No plotted axes found in loadings figure.');
        end

        % -----------------------------
        % Combine into one final figure
        % Panel A = scores by mean
        % Panel B = scores by class
        % Panel C = loadings
        % -----------------------------
        comboFig = figure('Color','w', ...
            'Units','pixels', ...
            'Position',[100 100 1800 550]);

        tl = tiledlayout(comboFig, 1, 3, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        % Panel A
        dst1 = nexttile(tl, 1);
        copyAxesContent(ax2, dst1);
        addPanelLabel(dst1, 'A');

        % Panel B
        dst2 = nexttile(tl, 2);
        copyAxesContent(ax1, dst2);
        addPanelLabel(dst2, 'B');
        addClassLegend(dst2, {'Control','Case'});

        % Panel C
        dst3 = nexttile(tl, 3);
        copyAxesContent(ax3, dst3);
        addPanelLabel(dst3, 'C');

        % -----------------------------
        % Export final figure as high-quality PDF (vector)
        % -----------------------------
        pdfFile = fullfile(outDirPDF, ...
            sprintf('%s_rep%d_scores_loadings_combined.pdf', methodSaveName, r));

        exportgraphics(comboFig, pdfFile, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');

        % -----------------------------
        % Optional: export high-resolution PNG
        % -----------------------------
        pngFile = fullfile(outDirPNG, ...
            sprintf('%s_rep%d_scores_loadings_combined.png', methodSaveName, r));

        exportgraphics(comboFig, pngFile, ...
            'Resolution', 600, ...
            'BackgroundColor', 'white');

        close(fig1);
        close(fig2);
        close(fig3);
        close(comboFig);

    end
end

fprintf('Done.\n');

% =========================================================
% Helper functions
% =========================================================

function ax = findMainAxes(fig)

    axs = findall(fig, 'Type', 'axes');
    ax = [];

    maxChildren = -1;

    for i = 1:numel(axs)
        n = numel(allchild(axs(i)));
        if n > maxChildren
            maxChildren = n;
            ax = axs(i);
        end
    end

end

function copyAxesContent(srcAx, dstAx)

    kids = allchild(srcAx);
    copyobj(kids, dstAx);
    hold(dstAx, 'on');

    dstAx.XLim = srcAx.XLim;
    dstAx.YLim = srcAx.YLim;
    dstAx.XTick = srcAx.XTick;
    dstAx.YTick = srcAx.YTick;
    dstAx.XScale = srcAx.XScale;
    dstAx.YScale = srcAx.YScale;
    dstAx.XDir = srcAx.XDir;
    dstAx.YDir = srcAx.YDir;
    dstAx.Box = srcAx.Box;
    dstAx.LineWidth = srcAx.LineWidth;
    dstAx.FontSize = srcAx.FontSize;

    xlabel(dstAx, srcAx.XLabel.String, ...
        'Interpreter', srcAx.XLabel.Interpreter);

    ylabel(dstAx, srcAx.YLabel.String, ...
        'Interpreter', srcAx.YLabel.Interpreter);

    title(dstAx, '');

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

function addClassLegend(ax, labels)

    lgd = legend(ax);
    if ~isempty(lgd) && isvalid(lgd)
        delete(lgd);
    end

    h = findobj(ax, 'Type', 'Scatter');
    if isempty(h)
        h = findobj(ax, 'Type', 'line');
    end

    h = flipud(h);

    if numel(h) >= 2
        legend(ax, h(1:2), labels, 'Location', 'best');
    else
        warning('Could not create class legend: fewer than two plot objects found.');
    end

end