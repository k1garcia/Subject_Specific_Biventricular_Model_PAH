%% Step1_NewMethodComputeCoordLongCircRad_OPTION1_LOCKED.m
% OPTION 1 (fixed): eL = base-plane normal (sign-corrected), then LOCK eL.
% Only apexCtr/baseCtr are refined along eL (no axis drifting).
%
% Outputs:
%  - FigA: chosen base slab highlighted
%  - Fig1: surface colored by s with eL arrow + apex/base markers
%  - Fig2: histogram of s
%  - Fig3: local axes quiver with legend (eL/eR/eC)
%
% Saves: outMat with X,E,Fs,surfNodeIds,eL,eR,eC,apexCtr,baseCtr,s,best,n_base

clear; clc; close all;

%% ================= USER EDIT =================
root    = 'C:\Users\k1garcia\Desktop\PAH_FEA\ZeroLoadState_InverseFEA\FeBio';
febFile = fullfile(root,'W283W0.feb');   % <-- change
outMat  = fullfile(root,'W283W0_geomCoords.mat');

% Base-plane search settings
nCandPerSide   = 6;
slabFrac       = 0.020;   % slab half-thickness = slabFrac*bboxDiag
minSlabNodes   = 100;
planarTolFrac  = 0.010;

% Candidate scoring weights
w_rms  = 0.65;
w_area = 0.25;
w_n    = 0.10;

% Force base slab to be in TOP of heart in GLOBAL Z
useGlobalZFilter = true;
zQuantileKeep    = 0.80;  % keep candidates with ctr.z >= top 20%

% Option 1 refinement (centers only) ALONG FIXED eL
qBase = 0.90;
qApex = 0.10;
nIter = 5;

% Plot settings
nQuiver = 450;
saveFigs = true;
figDir   = fullfile(root,'FIGS_COORDSYS_CHECKS');
%% ============================================

assert(exist(febFile,'file')==2, 'Missing .feb file: %s', febFile);
if saveFigs && ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('Reading FEB mesh:\n  %s\n', febFile);
[X, E] = read_feb_nodes_elements(febFile);
N = size(X,1);
M = size(E,1);
fprintf('  Nodes: %d | Elements: %d (tet4)\n', N, M);

%% ---- Surface extraction ----
[Fs, ~] = boundary_faces_from_tets(E);
surfNodeIds = unique(Fs(:));
Xs = X(surfNodeIds,:);
fprintf('Surface triangles: %d | Surface nodes: %d\n', size(Fs,1), numel(surfNodeIds));

%% ---- Bounding box ----
bb = max(X) - min(X);
bboxDiag = norm(bb);
slabHalfThickness = slabFrac * bboxDiag;
planarTol = planarTolFrac * bboxDiag;

%% =========================================================
%  GUARDED base-plane detection
%% =========================================================
fprintf('\n==== Base-plane detection (GUARDED search) ====\n');

axVecs = eye(3);
cand = struct('id',{},'axis',{},'side',{},'u',{},'c0',{},'ids',{}, ...
              'ctr',{},'n',{},'rms',{},'maxd',{},'area',{},'nPts',{});
candID = 0;

for ax = 1:3
    u = axVecs(:,ax)';       % 1x3
    proj = Xs*u';            % Ns x 1

    qList = linspace(0.85, 0.99, nCandPerSide);
    for side = 1:2
        if side == 1, q = 1 - qList; else, q = qList; end

        for k = 1:numel(q)
            c0 = quantile(proj, q(k));
            idsLocal = find(abs(proj - c0) <= slabHalfThickness);
            if numel(idsLocal) < minSlabNodes, continue; end

            P = Xs(idsLocal,:);
            Pc = P - mean(P,1);
            [~,~,V] = svd(Pc,'econ');
            n = V(:,3); n = n / norm(n); % plane normal (sign ambiguous)

            d = (P - mean(P,1)) * n;
            rmsd = sqrt(mean(d.^2));
            maxd = max(abs(d));

            [t1,t2] = plane_basis(n);
            uv = [ (P-mean(P,1))*t1, (P-mean(P,1))*t2 ];
            area2D = 0;
            if size(uv,1) >= 10
                try
                    K = convhull(uv(:,1), uv(:,2));
                    area2D = polyarea(uv(K,1), uv(K,2));
                catch
                    area2D = 0;
                end
            end

            candID = candID + 1;
            cand(candID).id   = candID;
            cand(candID).axis = ax;
            cand(candID).side = side;
            cand(candID).u    = u;
            cand(candID).c0   = c0;
            cand(candID).ids  = surfNodeIds(idsLocal); % global IDs
            cand(candID).ctr  = mean(P,1);
            cand(candID).n    = n(:)';                 % 1x3
            cand(candID).rms  = rmsd;
            cand(candID).maxd = maxd;
            cand(candID).area = area2D;
            cand(candID).nPts = numel(idsLocal);

            fprintf('  cand %3d axis %d side %d: n=%4d | rms=%.3f | area=%.2f\n', ...
                candID, ax, side, numel(idsLocal), rmsd, area2D);
        end
    end
end

assert(~isempty(cand), 'No candidates found. Try increasing nCandPerSide/slabFrac or lowering minSlabNodes.');

%% ---- Global-Z filter (keep top-of-heart slabs) ----
if useGlobalZFilter
    zThresh = quantile(Xs(:,3), zQuantileKeep);
    keep = arrayfun(@(c) c.ctr(3) >= zThresh, cand);
    candF = cand(keep);

    fprintf('\n==== Global-Z filter ====\n');
    fprintf('  Keeping candidates with ctr.z >= quantile(z, %.2f) = %.3f\n', zQuantileKeep, zThresh);
    fprintf('  Before: %d candidates | After: %d candidates\n', numel(cand), numel(candF));

    if ~isempty(candF), cand = candF;
    else
        warning('Z-filter removed all candidates. Relax zQuantileKeep or disable filter. Using unfiltered candidates.');
    end
end

%% ---- Score candidates ----
rmsAll  = reshape([cand.rms],  [], 1);
areaAll = reshape([cand.area], [], 1);
nAll    = reshape([cand.nPts], [], 1);

rmsN  = (rmsAll  - min(rmsAll )) ./ max(eps, (max(rmsAll )-min(rmsAll )));
areaN = (areaAll - min(areaAll)) ./ max(eps, (max(areaAll)-min(areaAll)));
nN    = (nAll    - min(nAll   )) ./ max(eps, (max(nAll   )-min(nAll   )));

score = -w_rms*rmsN + w_area*areaN + w_n*nN;
[~, iBest] = max(score);
best = cand(iBest);

fprintf('\n==== Base-plane result (GUARDED) ====\n');
fprintf('  bboxDiag=%.2f | slabHalfThickness=%.3f\n', bboxDiag, slabHalfThickness);
fprintf('  BEST cand %d: axis %d side %d\n', best.id, best.axis, best.side);
fprintf('  slab nodes=%d | rms=%.4f (tol=%.3f) | area=%.2f\n', best.nPts, best.rms, planarTol, best.area);
fprintf('  base plane center = [%.3f %.3f %.3f]\n', best.ctr(1), best.ctr(2), best.ctr(3));
fprintf('  base plane normal (raw) = [%.3f %.3f %.3f]\n', best.n(1), best.n(2), best.n(3));

%% ---- PROOF FIG A: highlight chosen base slab ----
figA = figure('Color','w','Name','PROOF: chosen base slab (highlight)');
trisurf(Fs, X(:,1), X(:,2), X(:,3), 'FaceAlpha',0.10, 'EdgeColor','none');
axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
hold on;
plot3(X(best.ids,1), X(best.ids,2), X(best.ids,3), 'r.', 'MarkerSize',12);
plot3(best.ctr(1), best.ctr(2), best.ctr(3), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
title(sprintf('Chosen base slab (red): cand %d | axis %d side %d | ctr z=%.3f | rms=%.3f', ...
    best.id, best.axis, best.side, best.ctr(3), best.rms));
hold off;

fprintf('  Base slab z-range: [%.3f, %.3f] | mean=%.3f\n', ...
    min(X(best.ids,3)), max(X(best.ids,3)), mean(X(best.ids,3)));
fprintf('  Global z-range    : [%.3f, %.3f]\n', min(X(:,3)), max(X(:,3)));

if saveFigs
    exportgraphics(figA, fullfile(figDir,'FigA_baseSlab_highlight.png'), 'Resolution', 300);
end

%% =========================================================
%  OPTION 1 (FIXED): sign-correct normal + LOCK eL
%% =========================================================
baseCtr = best.ctr(:);                 % 3x1
n_base  = best.n(:); n_base = n_base / norm(n_base);  % 3x1, sign ambiguous

% ---- SIGN FIX: make normal point from heart center toward the baseCtr ----
Xmean = mean(X,1)';                    % 3x1
toBase = baseCtr - Xmean;              % vector from heart center to base centroid
if dot(n_base, toBase) < 0
    n_base = -n_base;                  % flip to point toward base region
end

% Also (optional) ensure it roughly points to +Z (since base is top)
if dot(n_base, [0;0;1]) < 0
    n_base = -n_base;
end

eL = n_base;  % LOCKED longitudinal axis (apex -> base)

fprintf('\n==== Option 1 initialization (LOCKED eL) ====\n');
fprintf('  baseCtr = [%.3f %.3f %.3f]\n', baseCtr(1), baseCtr(2), baseCtr(3));
fprintf('  n_base (sign-fixed) = [%.3f %.3f %.3f]\n', n_base(1), n_base(2), n_base(3));
fprintf('  dot(n_base, +Z) = %.3f\n', dot(n_base,[0;0;1]));

% Initial apex guess: go opposite direction from base along eL
apexCtr = baseCtr - eL*(0.60*bboxDiag);

% ---- Refinement: update centers ONLY; keep eL fixed ----
for it = 1:nIter
    s = (X - Xmean') * eL;            % Nx1 along fixed eL
    sBase = quantile(s, qBase);
    sApex = quantile(s, qApex);

    baseIds = find(s >= sBase);
    apexIds = find(s <= sApex);

    baseCtr = mean(X(baseIds,:),1)';  % update centers
    apexCtr = mean(X(apexIds,:),1)';
end

% Ensure orientation is apex->base along fixed eL:
% (base should have larger s than apex)
s_ap = (apexCtr' - Xmean') * eL;
s_ba = (baseCtr' - Xmean') * eL;
if s_ba < s_ap
    eL = -eL;
    s_ap = (apexCtr' - Xmean') * eL;
    s_ba = (baseCtr' - Xmean') * eL;
end

fprintf('\n==== After refinement (centers only; eL locked) ====\n');
fprintf('  apexCtr = [%.3f %.3f %.3f]\n', apexCtr(1), apexCtr(2), apexCtr(3));
fprintf('  baseCtr = [%.3f %.3f %.3f]\n', baseCtr(1), baseCtr(2), baseCtr(3));
fprintf('  eL      = [%.3f %.3f %.3f]\n', eL(1), eL(2), eL(3));
fprintf('  Alignment |n_base · eL| = %.4f (angle %.2f deg)\n', ...
    abs(dot(n_base,eL)), acosd(min(1,max(-1,abs(dot(n_base,eL))))));

%% ---- Build eR and eC using centerline through apexCtr along eL ----
Xa = X - apexCtr';
t = Xa * eL;                      % Nx1
Xline = apexCtr' + t.*eL';         % Nx3
eR = X - Xline;                    % Nx3
eR = normalize_rows(eR);

eL_node = repmat(eL', N, 1);
eC = cross(eL_node, eR, 2);
eC = normalize_rows(eC);

% Re-orthogonalize eR to avoid drift
eR = cross(eC, eL_node, 2);
eR = normalize_rows(eR);

%% ---- PROOF checks ----
s = (X - Xmean') * eL;

Pbest = X(best.ids,:);
dBest = (Pbest - mean(Pbest,1)) * n_base;

fprintf('\n==== Coordinate-system confirmations (PROOF) ====\n');
fprintf('  mean |dot(eL,eC)| = %.3e\n', mean(abs(sum(eL_node.*eC,2))));
fprintf('  mean |dot(eL,eR)| = %.3e\n', mean(abs(sum(eL_node.*eR,2))));
fprintf('  mean |dot(eC,eR)| = %.3e\n', mean(abs(sum(eC.*eR,2))));
fprintf('  s-range: min(s)=%.3f  max(s)=%.3f\n', min(s), max(s));
fprintf('  s(apexCtr)=%.3f  s(baseCtr)=%.3f  (base should be larger)\n', s_ap, s_ba);
fprintf('  Base planarity: RMS dist to plane = %.4f | max dist = %.4f\n', sqrt(mean(dBest.^2)), max(abs(dBest)));

%% =========================================================
%  FIGURES
%% =========================================================

% FIG 1: Surface colored by s + axis arrow
fig1 = figure('Name','Surface + longitudinal axis','Color','w');
trisurf(Fs, X(:,1), X(:,2), X(:,3), s, 'EdgeColor','none');
axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Surface colored by s, with e_L (apex \rightarrow base)');
colormap turbo; cb = colorbar; ylabel(cb,'s (mm along e_L)');

hold on;
plot3(apexCtr(1), apexCtr(2), apexCtr(3), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
plot3(baseCtr(1), baseCtr(2), baseCtr(3), 'ks', 'MarkerFaceColor','k', 'MarkerSize',7);

Lscale = 0.35 * bboxDiag;
quiver3(apexCtr(1), apexCtr(2), apexCtr(3), eL(1)*Lscale, eL(2)*Lscale, eL(3)*Lscale, ...
    0, 'LineWidth',3, 'Color',[0.95 0.65 0.10]);
text(apexCtr(1), apexCtr(2), apexCtr(3), '  apex', 'FontSize',12, 'Color','k');
text(baseCtr(1), baseCtr(2), baseCtr(3), '  base', 'FontSize',12, 'Color','k');
hold off;

% FIG 2: Histogram of s
fig2 = figure('Name','Histogram of s','Color','w');
histogram(s, 60);
xlabel('s'); ylabel('count');
title('Distribution of longitudinal coordinate s');

% FIG 3: Local axes sampled + legend
rng(0);
idx = surfNodeIds;
if numel(idx) > nQuiver
    idx = idx(randperm(numel(idx), nQuiver));
end

fig3 = figure('Name','Local axes (sampled)','Color','w');
trisurf(Fs, X(:,1), X(:,2), X(:,3), s, 'FaceAlpha',0.12, 'EdgeColor','none');
axis equal; grid on; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Sampled local axes: e_L (long), e_R (rad), e_C (circ)');
%colormap turbo; cb3 = colorbar; ylabel(cb3,'s (mm along e_L)');
hold on;

Lq = 0.12 * mean(bb);

hL = quiver3(X(idx,1), X(idx,2), X(idx,3), ...
    eL_node(idx,1)*Lq, eL_node(idx,2)*Lq, eL_node(idx,3)*Lq, ...
    0, 'LineWidth',1.2, 'Color',[0.00 0.45 0.74]);

hR = quiver3(X(idx,1), X(idx,2), X(idx,3), ...
    eR(idx,1)*Lq, eR(idx,2)*Lq, eR(idx,3)*Lq, ...
    0, 'LineWidth',1.2, 'Color',[0.93 0.69 0.13]);

hC = quiver3(X(idx,1), X(idx,2), X(idx,3), ...
    eC(idx,1)*Lq, eC(idx,2)*Lq, eC(idx,3)*Lq, ...
    0, 'LineWidth',1.2, 'Color',[0.85 0.33 0.10]);

plot3(apexCtr(1), apexCtr(2), apexCtr(3), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
plot3(baseCtr(1), baseCtr(2), baseCtr(3), 'ks', 'MarkerFaceColor','k', 'MarkerSize',7);

legend([hL hR hC], ...
    {'e_L : longitudinal (apex \rightarrow base)', ...
     'e_R : radial (centerline \rightarrow wall)', ...
     'e_C : circumferential (e_L \times e_R)'}, ...
    'Location','southeastoutside');

hold off;

%% ---- Save outputs ----
save(outMat, 'X','E','Fs','surfNodeIds','best','n_base','eL','apexCtr','baseCtr','s','eR','eC');
fprintf('\n✅ Saved coordinate system to:\n  %s\n', outMat);

if saveFigs
    exportgraphics(fig1, fullfile(figDir,'Fig1_surface_s_eL.png'), 'Resolution', 300);
    exportgraphics(fig2, fullfile(figDir,'Fig2_hist_s.png'), 'Resolution', 300);
    exportgraphics(fig3, fullfile(figDir,'Fig3_localAxes_legend.png'), 'Resolution', 300);
    fprintf('✅ Saved figures to:\n  %s\n', figDir);
end

%% ========================= Helper functions =========================
function [X, E] = read_feb_nodes_elements(febFile)
    txt = fileread(febFile);

    nBlock = regexp(txt, '<Nodes[^>]*>(.*?)</Nodes>', 'tokens', 'once');
    if isempty(nBlock), error('Could not find <Nodes> block in %s', febFile); end
    nLines = regexp(nBlock{1}, '<node[^>]*>.*?</node>', 'match');

    X = zeros(numel(nLines),3);
    for i=1:numel(nLines)
        xyz = regexp(nLines{i}, '>([^<]+)<', 'tokens', 'once');
        vals = sscanf(strrep(xyz{1},',',' '), '%f %f %f');
        X(i,:) = vals(:)';
    end

    eBlock = regexp(txt, '<Elements[^>]*type="tet4"[^>]*>(.*?)</Elements>', 'tokens', 'once');
    if isempty(eBlock)
        error('Could not find <Elements type="tet4"> block. If your mesh is not tet4, tell me the type.');
    end
    eLines = regexp(eBlock{1}, '<elem[^>]*>.*?</elem>', 'match');

    E = zeros(numel(eLines),4);
    for i=1:numel(eLines)
        conn = regexp(eLines{i}, '>([^<]+)<', 'tokens', 'once');
        vals = sscanf(strrep(conn{1},',',' '), '%d %d %d %d');
        E(i,:) = vals(:)';
    end
end

function [F, surfNodes] = boundary_faces_from_tets(E)
    f1 = E(:,[1 2 3]);
    f2 = E(:,[1 2 4]);
    f3 = E(:,[1 3 4]);
    f4 = E(:,[2 3 4]);
    Fall = [f1; f2; f3; f4];

    Fsorted = sort(Fall,2);
    [uF, ~, ic] = unique(Fsorted, 'rows');
    counts = accumarray(ic, 1);

    isBoundary = counts == 1;
    F = uF(isBoundary,:);
    surfNodes = unique(F(:));
end

function A = normalize_rows(A)
    n = vecnorm(A,2,2);
    n(n==0) = 1;
    A = A ./ n;
end

function [t1,t2] = plane_basis(n)
    n = n(:) / norm(n);
    if abs(dot(n,[1;0;0])) < 0.9
        a = [1;0;0];
    else
        a = [0;1;0];
    end
    t1 = cross(n,a); t1 = t1 / norm(t1);
    t2 = cross(n,t1); t2 = t2 / norm(t2);
end
