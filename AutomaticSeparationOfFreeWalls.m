%% AutoSegment_LV_RV_Septum_AxisBandSimple.m
% Automatic LV / RV / Septum segmentation using:
%  - LV and RV endocardial facet tags (auto-detected from STLs)
%  - An LV<->RV axis
%  - A mid-band along that axis for the septum
%
% Every element becomes either LV, RV, or Septum.

clc; clear; close all;

%% ====================== USER INPUT ======================================
ratname = 'W283W8';
rootDir = 'C:\Users\k1garcia\Desktop\PAH_FEA\Volume_Meshes';
baseDir = fullfile(rootDir, ratname);
h5file  = fullfile(baseDir, [ratname '_With_Fibers.h5']);

% Septum thickness control (fractions of LV-RV separation along the axis)
% bandFrac_LV  = how far the band extends toward the LV side
% bandFrac_RV  = how far the band extends toward the RV side
bandFrac_LV = 0.30;   % thinner on LV side  (less LV labeled as S)
bandFrac_RV = 0.50;   % thicker on RV side  (keep RV septal curvature)

%% ====================== LOAD MESH ======================================
assert(exist(h5file,'file')==2, 'Missing H5 file: %s', h5file);

V = h5read(h5file,'/Mesh/Mesh/geometry');          % [3 x N] or [N x 3]
if size(V,1)==3, V = V.'; end                      % [N x 3]

E = h5read(h5file,'/Mesh/Mesh/topology');          % [4 x Ne] or [Ne x 4]
if size(E,1)==4, E = E.'; end                      % [Ne x 4]
E = double(E);
if min(E(:))==0, E = E + 1; end                    % 1-based

[Nn, ~] = size(V);
Ne      = size(E,1);
fprintf('Loaded mesh for %s: %d nodes, %d tets\n', ratname, Nn, Ne);

% Outer surface for visualization
TR       = triangulation(E,V);
Fb_outer = freeBoundary(TR);

%% ================= FACET TAGS FROM H5 ==================================
Fb = h5read(h5file,'/MeshTags/Facet tags/topology');  % [3 x Nf] or [Nf x 3]
if size(Fb,1)==3, Fb = Fb.'; end
Fb = double(Fb);
if min(Fb(:))==0, Fb = Fb + 1; end

Cb = h5read(h5file,'/MeshTags/Facet tags/Values');    % [Nf x 1]
Cb = double(Cb(:));

fprintf('--- FACET TAG SUMMARY ---\n');
uTags = unique(Cb);
for k = 1:numel(uTags)
    fprintf('Tag %d  →  %d faces\n', uTags(k), sum(Cb==uTags(k)));
end

%% ---------- AUTO-DETECT WHICH FACET TAG IS LV vs RV --------------------
% Root folder where your LV/RV STLs live:
stlRoot = 'C:\Users\k1garcia\Dropbox\Kristen Garcia Working Documents\Cardiac MRI\PAH_Imaging_FEA\Segmented stl files';

% Try a couple of common patterns; pick the first that exists
cand1  = fullfile(stlRoot, ratname);
cand2  = fullfile(stlRoot, [ratname ' Segmentations']);
if exist(cand1,'dir')
    stlDir = cand1;
elseif exist(cand2,'dir')
    stlDir = cand2;
else
    error('Could not find STL folder for %s. Checked:\n  %s\n  %s', ...
          ratname, cand1, cand2);
end
fprintf('\nUsing STL directory: %s\n', stlDir);

% Find LV/RV STL files (case-insensitive)
stlFiles = dir(fullfile(stlDir,'*.stl'));
if isempty(stlFiles)
    error('No STL files found in %s. Update stlRoot or naming.', stlDir);
end

LVfile = [];
RVfile = [];
for k = 1:numel(stlFiles)
    nm = lower(stlFiles(k).name);
    if contains(nm,'_lv.stl')
        LVfile = fullfile(stlFiles(k).folder, stlFiles(k).name);
    elseif contains(nm,'_rv.stl')
        RVfile = fullfile(stlFiles(k).folder, stlFiles(k).name);
    end
end
if isempty(LVfile) || isempty(RVfile)
    error('Could not find *_LV.stl and *_RV.stl in %s', stlDir);
end

% Load LV/RV STL and compute centers
TR_LV   = stlread(LVfile);
TR_RV   = stlread(RVfile);
cLV_stl = mean(TR_LV.Points,1);
cRV_stl = mean(TR_RV.Points,1);

fprintf('STL centers: LV=[%.2f %.2f %.2f], RV=[%.2f %.2f %.2f]\n', ...
        cLV_stl, cRV_stl);

% Compute centers for *every* facet tag in this H5
uTags = unique(Cb);
tagCenters = nan(numel(uTags),3);
for i = 1:numel(uTags)
    t = uTags(i);
    F_t = Fb(Cb==t,:);          % faces with this tag
    V_t = V(unique(F_t(:)),:);  % corresponding nodes
    tagCenters(i,:) = mean(V_t,1);
end

fprintf('\nFacet tag centers:\n');
for i = 1:numel(uTags)
    t = uTags(i);
    c = tagCenters(i,:);
    fprintf('  tag %d : [%.2f %.2f %.2f]\n', t, c(1), c(2), c(3));
end

% Distances from each tag center to LV and RV STL centers
nT = numel(uTags);
dLV = zeros(nT,1);
dRV = zeros(nT,1);
for i = 1:nT
    c = tagCenters(i,:);
    dLV(i) = norm(c - cLV_stl);
    dRV(i) = norm(c - cRV_stl);
end

fprintf('\nDistances tag→STL (mm):\n');
for i = 1:nT
    fprintf('  tag %d : to LV = %.2f, to RV = %.2f\n', ...
            uTags(i), dLV(i), dRV(i));
end

% Choose the best pair of distinct tags for LV and RV
bestScore = inf;
bestLV = NaN;
bestRV = NaN;

for i = 1:nT
    for j = 1:nT
        if i==j, continue; end
        % require that each tag is closer to its assigned cavity
        if dLV(i) < dRV(i) && dRV(j) < dLV(j)
            score = dLV(i) + dRV(j);
            if score < bestScore
                bestScore = score;
                bestLV   = uTags(i);
                bestRV   = uTags(j);
            end
        end
    end
end

if isnan(bestLV) || isnan(bestRV)
    warning('Could not find a clean LV/RV tag pair; falling back to tags 3/4.');
    tag_LV = 3;
    tag_RV = 4;
else
    tag_LV = bestLV;
    tag_RV = bestRV;
end

fprintf('\n>>> Auto-detected facet tags: LV = %d, RV = %d\n\n', tag_LV, tag_RV);

%% ================= FACET TAGS → LV/RV ENDO NODES =======================
F_LV_endo = Fb(Cb==tag_LV,:);
F_RV_endo = Fb(Cb==tag_RV,:);
assert(~isempty(F_LV_endo) && ~isempty(F_RV_endo), ...
    'LV or RV endocardial faces (tags %d/%d) are empty.', tag_LV, tag_RV);

LV_nodes = unique(F_LV_endo(:));
RV_nodes = unique(F_RV_endo(:));

V_LV = V(LV_nodes,:);
V_RV = V(RV_nodes,:);

fprintf('LV endo nodes: %d   RV endo nodes: %d\n', size(V_LV,1), size(V_RV,1));

%% ================= LV<->RV AXIS & ELEMENT PROJECTIONS ==================
% Approximate LV and RV cavity centers (from H5 endocardial nodes)
cLV = mean(V_LV,1);   % [1 x 3]
cRV = mean(V_RV,1);   % [1 x 3]

fprintf('\nLV center (H5): [%.2f  %.2f  %.2f]\n', cLV);
fprintf('RV center (H5): [%.2f  %.2f  %.2f]\n', cRV);
fprintf('Center differences (LV - RV): dX=%.2f, dY=%.2f, dZ=%.2f\n', ...
        cLV(1)-cRV(1), cLV(2)-cRV(2), cLV(3)-cRV(3));

axisVec = cLV - cRV;          % from RV -> LV
axisLen = norm(axisVec);
assert(axisLen > 0, 'LV/RV centers are identical; cannot define axis.');
axisVec = axisVec / axisLen;  % unit vector

% Element centroids
C_elem = ( V(E(:,1),:) + V(E(:,2),:) + V(E(:,3),:) + V(E(:,4),:) ) / 4;  % [Ne x 3]

% Project LV and RV endo nodes on axis
proj_LV_nodes = V_LV * axisVec.';    % [nLV x 1]
proj_RV_nodes = V_RV * axisVec.';    % [nRV x 1]

meanLV = mean(proj_LV_nodes);
meanRV = mean(proj_RV_nodes);

% Midplane along the axis
midval = 0.5*(meanLV + meanRV);

% Project all element centroids on axis
proj_elem = C_elem * axisVec.';      % [Ne x 1]

% Axis coordinate and span
dAxis = proj_elem - midval;       % negative ~ RV side, positive ~ LV side
dSpan = abs(meanLV - meanRV);     % typical LV-RV separation

fprintf('\nAxis separation |meanLV - meanRV| = %.3f (mesh units)\n', dSpan);

% Asymmetric septum band widths (in coordinate units)
bandWidth_LV = bandFrac_LV * dSpan;   % toward LV
bandWidth_RV = bandFrac_RV * dSpan;   % toward RV

fprintf('Using bandWidth_LV = %.3f (frac=%.2f), bandWidth_RV = %.3f (frac=%.2f)\n', ...
        bandWidth_LV, bandFrac_LV, bandWidth_RV, bandFrac_RV);

%% ================= CLASSIFY ELEMENTS: LV / RV / SEPTUM =================
% Septum: between -bandWidth_RV and +bandWidth_LV
isS  = (dAxis >= -bandWidth_RV) & (dAxis <= bandWidth_LV);
isLV =  dAxis >  bandWidth_LV;
isRV =  dAxis < -bandWidth_RV;

% Safety: if anything unclassified (should be rare), assign by sign
unclassified = ~(isLV | isRV | isS);
if any(unclassified)
    fprintf('Warning: %d unclassified elems; assigning by sign of dAxis.\n', ...
        sum(unclassified));
    isLV(unclassified & (dAxis >= 0)) = true;
    isRV(unclassified & (dAxis <  0)) = true;
end

LV_ids = find(isLV);
RV_ids = find(isRV);
S_ids  = find(isS);

fprintf('\nFinal classification (bandFrac_LV=%.2f, bandFrac_RV=%.2f):\n', ...
        bandFrac_LV, bandFrac_RV);

fprintf('  LV: %d elems (%.1f %%)\n', numel(LV_ids), 100*numel(LV_ids)/Ne);
fprintf('  RV: %d elems (%.1f %%)\n', numel(RV_ids), 100*numel(RV_ids)/Ne);
fprintf('  S : %d elems (%.1f %%)\n', numel(S_ids),  100*numel(S_ids)/Ne);

%% ==================== SAVE TO CSV ======================================
writematrix(LV_ids, fullfile(baseDir,'LV_freewall.csv'));
writematrix(RV_ids, fullfile(baseDir,'RV_freewall.csv'));
writematrix(S_ids,  fullfile(baseDir,'S_freewall.csv'));

fprintf('\nSaved element sets:\n  %s\n  %s\n  %s\n', ...
    fullfile(baseDir,'LV_freewall.csv'), ...
    fullfile(baseDir,'RV_freewall.csv'), ...
    fullfile(baseDir,'S_freewall.csv'));

%% ==================== VISUALIZATION: OUTER SURFACE =====================
% Color LV = blue, RV = green, S = red, on outer surface.
nodeRGB   = zeros(Nn,3);
nodeCount = zeros(Nn,1);

colLV = [0 0.3 1];
colRV = [0 0.7 0];
colS  = [1 0 0];

for e = LV_ids.'
    nd = E(e,:);
    nodeRGB(nd,:) = nodeRGB(nd,:) + colLV;
    nodeCount(nd) = nodeCount(nd) + 1;
end
for e = RV_ids.'
    nd = E(e,:);
    nodeRGB(nd,:) = nodeRGB(nd,:) + colRV;
    nodeCount(nd) = nodeCount(nd) + 1;
end
for e = S_ids.'
    nd = E(e,:);
    nodeRGB(nd,:) = nodeRGB(nd,:) + colS;
    nodeCount(nd) = nodeCount(nd) + 1;
end

nodeCount(nodeCount==0) = 1;
C_rgb = nodeRGB ./ nodeCount;

figure('Color','w');
patch('Faces',Fb_outer,'Vertices',V, ...
      'FaceVertexCData',C_rgb, ...
      'FaceColor','interp', ...
      'EdgeColor','k','EdgeAlpha',0.05);
axis equal vis3d; view(3);
camlight headlight; lighting gouraud;
title(sprintf('%s: LV (blue) / RV (green) / Septum (red) – axis band', ratname), ...
      'FontSize',14);
xlabel('X'); ylabel('Y'); zlabel('Z');
set(gca,'DataAspectRatio',[1 1 1]);

%% ==================== VISUALIZATION: XY VIEW + LEGEND ==================
figure('Color','w','Position',[100 100 1600 600]);

patch('Faces',Fb_outer,'Vertices',V, ...
      'FaceVertexCData',C_rgb, ...
      'FaceColor','interp', ...
      'EdgeColor','k','EdgeAlpha',0.15);
view(2);                          % XY plane
axis equal tight;
xlabel('X (mm)');
ylabel('Y (mm)');
title('XY view','FontSize',20,'FontWeight','bold');

hold on;
pLV = patch(nan,nan,colLV,'FaceColor',colLV,'EdgeColor','none');
pRV = patch(nan,nan,colRV,'FaceColor',colRV,'EdgeColor','none');
pS  = patch(nan,nan,colS, 'FaceColor',colS, 'EdgeColor','none');

legend([pLV pRV pS], {'LV (blue)','RV (green)','Septum (red)'}, ...
       'Location','bestoutside');

ax = gca;
ax.LabelFontSizeMultiplier = 1.3;
ax.FontSize = 14;
set(ax,'YAxisLocation','left','XAxisLocation','bottom');

%% ==================== DIAGNOSTIC LV/RV ENDO SCATTER ====================
figure('Color','w','Position',[100 100 1400 450]);

% 3D view
subplot(1,3,1); hold on;
scatter3(V_LV(:,1), V_LV(:,2), V_LV(:,3), 10, 'r', 'filled');
scatter3(V_RV(:,1), V_RV(:,2), V_RV(:,3), 10, 'b', 'filled');
axis equal; grid on;
view(120,20);
title('LV/RV endo – 3D');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend({'LV','RV'});

% XY (top) view
subplot(1,3,2); hold on;
scatter(V_LV(:,1), V_LV(:,2), 10, 'r', 'filled');
scatter(V_RV(:,1), V_RV(:,2), 10, 'b', 'filled');
axis equal; grid on;
title('LV/RV endo – XY');
xlabel('X'); ylabel('Y');

% YZ (side) view
subplot(1,3,3); hold on;
scatter(V_LV(:,2), V_LV(:,3), 10, 'r', 'filled');
scatter(V_RV(:,2), V_RV(:,3), 10, 'b', 'filled');
axis equal; grid on;
title('LV/RV endo – YZ');
xlabel('Y'); ylabel('Z');
