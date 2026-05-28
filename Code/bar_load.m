%% Ground truth vs oMEDA (RawData) in one panel
% Left: Ground truth
% Right: oMEDA for RawData
% Compact layout with small top white space so upper labels are not cut.
% Exported as high-quality vector PDF.

close all;

%% ---------------- INPUTS ----------------
load Meta_genus.mat;   % must contain GGT

% This assumes raw data matrix X is already in the workspace.
% This assumes omedaPca(...) is available in the MATLAB path.

y = repmat([ones(50,1)*-1; ones(50,1)], 3, 1);

lab = [
    "Bifidobacterium"
    "Butyricimonas"
    "Odoribacter"
    "Paraprevotella"
    "Bacteroides"
    "Parabacteroides"
    "Prevotella"
    "Unknown1"
    "Unknown2"
    "Unknown3"
    "Clostridium"
    "Unknown4"
    "Unknown5"
    "Ruminococcus"
    "Anaerostipes"
    "Blautia"
    "Coprococcus"
    "Dorea"
    "Lachnospira"
    "Roseburia"
    "Unknown6"
    "Unknown7"
    "Unknown8"
    "Unknown9"
    "Faecalibacterium"
    "Oscillospira"
    "Ruminococcus_dup"
    "Dialister"
    "Unknown10"
    "Unknown11"
    "Unknown12"
    "Sutterella"
    "Escherichia"
    "Klebsiella"
    "Haemophilus"
    "Akkermansia"
];

lab = string(lab(:));

%% ---------------- GROUND TRUTH ----------------
GTP1 = GGT;
GTP1 = GTP1 - mean(GTP1(:,1:2), 2) * ones(1, 3);

ground_truth = sign(GTP1(:,2)) .* (GTP1(:,2).^2) ...
             - sign(GTP1(:,1)) .* (GTP1(:,1).^2);

ground_truth = ground_truth / norm(ground_truth);
GT = ground_truth(:);

%% ---------------- oMEDA FOR RAW DATA ----------------
rj = omedaPca(X, [], X, y, 'VarsLabel', lab, 'Preprocessing', 1);
close all;

if ~isvector(rj)
    rj = rj(:,1);   % change component if needed
end
rj = rj(:);

%% ---------------- SAFE ALIGNMENT ----------------
n = min([numel(lab), numel(GT), numel(rj)]);

lab = lab(1:n);
GT  = GT(1:n);
rj  = rj(1:n);

%% ---------------- BUILD TABLES ----------------
Tgt = table(lab, double(GT), 'VariableNames', {'Taxon','Diff'});
Tom = table(lab, double(rj), 'VariableNames', {'Taxon','Diff'});

% Flip so top label is Akkermansia and bottom is Bifidobacterium
Tgt = flipud(Tgt);
Tom = flipud(Tom);

% Display labels to match attached figure style
Tgt.TaxonPlot = Tgt.Taxon;
Tom.TaxonPlot = Tom.Taxon;

Tgt.TaxonPlot = replace(Tgt.TaxonPlot, "Unknown", "unknown");
Tom.TaxonPlot = replace(Tom.TaxonPlot, "Unknown", "unknown");

Tgt.TaxonPlot = replace(Tgt.TaxonPlot, "Ruminococcus_dup", "Ruminococcus");
Tom.TaxonPlot = replace(Tom.TaxonPlot, "Ruminococcus_dup", "Ruminococcus");

% Scale each panel independently to [-1, 1]
if max(abs(Tgt.Diff)) > 0
    Tgt.DiffPlot = Tgt.Diff ./ max(abs(Tgt.Diff));
else
    Tgt.DiffPlot = Tgt.Diff;
end

if max(abs(Tom.Diff)) > 0
    Tom.DiffPlot = Tom.Diff ./ max(abs(Tom.Diff));
else
    Tom.DiffPlot = Tom.Diff;
end

%% ---------------- PRINT RESULTS ----------------
disp('--- GROUND TRUTH: HEALTHY / CONTROL (Diff < 0) ---');
disp(Tgt(Tgt.Diff < 0, {'Taxon','Diff'}));

disp('--- GROUND TRUTH: SICK / CASE (Diff > 0) ---');
disp(Tgt(Tgt.Diff > 0, {'Taxon','Diff'}));

disp('--- oMEDA RAW DATA: HEALTHY / CONTROL (Diff < 0) ---');
disp(Tom(Tom.Diff < 0, {'Taxon','Diff'}));

disp('--- oMEDA RAW DATA: SICK / CASE (Diff > 0) ---');
disp(Tom(Tom.Diff > 0, {'Taxon','Diff'}));

%% ---------------- FIGURE SETTINGS ----------------
figW_px = 800;
figH_px = 820;

f = figure('Color','w', ...
           'Units','pixels', ...
           'Position',[100 100 figW_px figH_px], ...
           'Renderer','painters');

% Slightly lower and shorter axes to create top white space.
% This prevents the top y-label from being cut.
ax1 = axes('Parent',f, ...
           'Units','normalized', ...
           'Position',[0.18 0.10 0.26 0.78]);   % left panel

ax2 = axes('Parent',f, ...
           'Units','normalized', ...
           'Position',[0.53 0.10 0.26 0.78]);   % right panel

%% ---------------- COLORS ----------------
negColor = [11/255, 250/255, 253/255];   % healthy / control
posColor = [1, 0, 0];                    % sick / case

%% ---------------- LEFT PANEL: GROUND TRUTH ----------------
hold(ax1,'on');

for i = 1:height(Tgt)
    if Tgt.DiffPlot(i) >= 0
        barh(ax1, i, Tgt.DiffPlot(i), ...
            'FaceColor', posColor, ...
            'EdgeColor', 'none', ...
            'BarWidth', 0.80);
    else
        barh(ax1, i, Tgt.DiffPlot(i), ...
            'FaceColor', negColor, ...
            'EdgeColor', 'none', ...
            'BarWidth', 0.80);
    end
end

set(ax1, ...
    'YTick', 1:height(Tgt), ...
    'YTickLabel', Tgt.TaxonPlot, ...
    'YDir', 'reverse', ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'FontAngle', 'italic', ...
    'LineWidth', 0.75, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'TickDir', 'out', ...
    'TickLength', [0 0]);

xlim(ax1, [-1 1]);
ylim(ax1, [0.2 height(Tgt)+0.8]);  % small top/bottom white space

xticks(ax1, [-1 -0.5 0 0.5 1]);
xticklabels(ax1, {'-1','-0.5','0','0.5','1'});

xline(ax1, 0, 'k-', 'LineWidth', 0.75);

xlabel(ax1, '');
ylabel(ax1, '');


box(ax1, 'on');
grid(ax1, 'off');

%% ---------------- RIGHT PANEL: oMEDA RAW DATA ----------------
hold(ax2,'on');

for i = 1:height(Tom)
    if Tom.DiffPlot(i) >= 0
        barh(ax2, i, Tom.DiffPlot(i), ...
            'FaceColor', posColor, ...
            'EdgeColor', 'none', ...
            'BarWidth', 0.80);
    else
        barh(ax2, i, Tom.DiffPlot(i), ...
            'FaceColor', negColor, ...
            'EdgeColor', 'none', ...
            'BarWidth', 0.80);
    end
end

set(ax2, ...
    'YTick', 1:height(Tom), ...
    'YTickLabel', [], ...
    'YDir', 'reverse', ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'FontAngle', 'italic', ...
    'LineWidth', 0.75, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'TickDir', 'out', ...
    'TickLength', [0 0]);

xlim(ax2, [-1 1]);
ylim(ax2, [0.2 height(Tom)+0.8]);  % same white space as left panel

xticks(ax2, [-1 -0.5 0 0.5 1]);
xticklabels(ax2, {'-1','-0.5','0','0.5','1'});

xline(ax2, 0, 'k-', 'LineWidth', 0.75);

xlabel(ax2, '');
ylabel(ax2, '');

box(ax2, 'on');
grid(ax2, 'off');

%% ---------------- PANEL LETTERS ----------------
annotation(f,'textbox',[0.06 0.92 0.03 0.03], ...
    'String','A', ...
    'LineStyle','none', ...
    'FontName','Arial', ...
    'FontSize',13, ...
    'FontWeight','bold');

annotation(f,'textbox',[0.48 0.92 0.03 0.03], ...
    'String','B', ...
    'LineStyle','none', ...
    'FontName','Arial', ...
    'FontSize',13, ...
    'FontWeight','bold');

%% ---------------- EXPORT ----------------
exportgraphics(f, 'GroundTruth_vs_oMEDA_RawData.pdf', ...
    'ContentType', 'vector');

exportgraphics(f, 'GroundTruth_vs_oMEDA_RawData.png', ...
    'Resolution', 600);