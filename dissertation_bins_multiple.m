%% A.6 Circumferential Stress Binning (core script)
% Computes bin-averaged circumferential and longitudinal deviatoric stresses
% within LV free wall, RV free wall, and septum using 20 circumferential bins
% (implemented as quantile bins along the short-axis X coordinate within each region).
%
% REQUIRED INPUT FILES (per ratID):
%   1) <geomDir>\<ratID>_geomStates.mat  containing:
%        - V_ED_final : [nNodes x 3] end-diastolic nodal coordinates
%        - E          : [nElem  x 4] tet4 connectivity (1-based)
%        - (optional) eL : [3 x 1] or [1 x 3] longitudinal axis unit vector
%   2) <geomDir>\<ratID>_stress_out.txt  FEBio element stress logfile (Cauchy stress)
%   3) <volMeshParent>\<meshID>\LV_freewall.csv / RV_freewall.csv / S_freewall.csv
%        Each CSV contains element IDs for that region (1-based or 0-based; 0-based auto-fixed).
%
% OUTPUTS:
%   - Results struct array "R" containing element stresses + bin IDs + bin means (per region)
%   - (optional) writes per-rat CSV summaries to outDir

clear; clc;

%% ---------------- USER SETTINGS ----------------
ratIDs    = {'W282W0','W282W8_ControlParams','W282W0_SuHxParams','W282W8'};
ratLabels = {'Female Baseline','Geometry Only','Parameters Only','Full Remodeling'};

% IMPORTANT: meshID selects the folder that contains the region CSVs
meshIDs   = {'W282W0','W282W8','W282W0','W282W8'};

geomDir       = 'C:\Users\k1garcia\Desktop\PAH_FEA\ZeroLoadState_InverseFEA\FeBio';
volMeshParent = 'C:\Users\k1garcia\Desktop\PAH_FEA\Volume_Meshes';
outDir        = fullfile(geomDir,'APPENDIX_A6_BINNING_CORE');
if ~exist(outDir,'dir'); mkdir(outDir); end

% Binning settings
nBins     = 20;
binMethod = 'quantile';  % 'quantile' recommended to avoid empty bins

% Regions and CSV names
regions      = {'LV','RV','S'};
regionFiles  = struct('LV','LV_freewall.csv','RV','RV_freewall.csv','S','S_freewall.csv');

% Metrics to bin-average (element-based)
metricsToAverage = {'sigLp','sigCp','diffLC'};  % σ'_L, σ'_C, σ_L - σ_C

% GIBBON logfile importer must be on path
% addpath(genpath('C:\Users\k1garcia\Desktop\GIBBON-master\GIBBON-master'));

assert(numel(ratIDs)==numel(meshIDs) && numel(ratIDs)==numel(ratLabels), ...
    'ratIDs, meshIDs, ratLabels must have the same length.');

%% ---------------- PROCESS ALL RATS ----------------
nR = numel(ratIDs);
R  = repmat(struct(), nR, 1);
valid = false(nR,1);

for r = 1:nR
    ratID  = ratIDs{r};
    label  = ratLabels{r};
    meshID = meshIDs{r};

    fprintf('\n--- %s (%s) ---\n', ratID, label);

    geomFile   = fullfile(geomDir, [ratID '_geomStates.mat']);
    stressFile = fullfile(geomDir, [ratID '_stress_out.txt']);

    if ~isfile(geomFile) || ~isfile(stressFile)
        warning('Missing geom or stress file for %s. Skipping.', ratID);
        continue;
    end

    % ---- Load geometry ----
    Sg = load(geomFile);
    if ~isfield(Sg,'V_ED_final') || ~isfield(Sg,'E')
        warning('%s geomStates missing V_ED_final or E. Skipping.', ratID);
        continue;
    end

    V = double(Sg.V_ED_final);
    E = double(Sg.E);

    if ~isempty(V) && size(V,2)~=3 && size(V,1)==3, V = V.'; end
    if isempty(V) || size(V,2)~=3 || any(~isfinite(V(:)))
        warning('%s invalid V_ED_final. Skipping.', ratID);
        continue;
    end

    nElem  = size(E,1);
    nNodes = size(V,1);

    % ---- Load stresses (last time point) ----
    D = importFEBio_logfile(stressFile, 0, 1);
    if ~isfield(D,'data') || isempty(D.data)
        warning('%s stress log empty. Skipping.', ratID);
        continue;
    end
    Sall = D.data(:,:,end);

    if size(Sall,1) ~= nElem || size(Sall,2) < 6
        warning('%s stress size mismatch (need nElem x >=6). Skipping.', ratID);
        continue;
    end

    sx  = Sall(:,1); sy  = Sall(:,2); sz  = Sall(:,3);
    sxy = Sall(:,4); syz = Sall(:,5); sxz = Sall(:,6);

    % ---- Element centroids ----
    Xe = (V(E(:,1),:) + V(E(:,2),:) + V(E(:,3),:) + V(E(:,4),:)) / 4;

    % ---- Hydrostatic pressure p = tr(sigma)/3 ----
    p_elem = (sx + sy + sz) / 3;

    % ---- Longitudinal axis eL (saved or PCA fallback) ----
    if isfield(Sg,'eL') && ~isempty(Sg.eL) && numel(Sg.eL)>=3 && all(isfinite(Sg.eL(:)))
        eL = double(Sg.eL(:)); eL = eL(1:3);
    else
        Xc = V - mean(V,1);
        [~,~,Vsvd] = svd(Xc,'econ');
        eL = Vsvd(:,1);
    end
    eL = eL / max(norm(eL), eps);
    if eL(3) < 0, eL = -eL; end

    % ---- Build local circumferential direction eC per element ----
    % Project element centroid to long-axis line, then define radial eR and circumferential eC = eL x eR
    ctr   = mean(V,1);
    d     = Xe - ctr;
    t     = d*eL;               % scalar projection onto eL
    xline = ctr + t.*eL';       % closest point on the long-axis line
    eR    = Xe - xline;
    eR    = eR ./ max(vecnorm(eR,2,2), eps);
    eC    = cross(repmat(eL',nElem,1), eR, 2);
    eC    = eC ./ max(vecnorm(eC,2,2), eps);

    % ---- Directional normal stresses ----
    ex = eL(1); ey = eL(2); ez = eL(3);
    sigmaL = sx.*ex.^2 + sy.*ey.^2 + sz.*ez.^2 + 2*(sxy.*ex.*ey + sxz.*ex.*ez + syz.*ey.*ez);

    cx = eC(:,1); cy = eC(:,2); cz = eC(:,3);
    sigmaC = sx.*cx.^2 + sy.*cy.^2 + sz.*cz.^2 + 2*(sxy.*cx.*cy + sxz.*cx.*cz + syz.*cy.*cz);

    % Deviatoric (remove hydrostatic) + anisotropy metrics
    sigLp  = sigmaL - p_elem;     % σ'_L
    sigCp  = sigmaC - p_elem;     % σ'_C
    diffLC = sigmaL - sigmaC;     % σ_L - σ_C

    % ---- Load region element IDs from CSVs ----
    regionDir = fullfile(volMeshParent, meshID);
    regElems = struct('LV',[],'RV',[],'S',[]);

    for rr = 1:numel(regions)
        reg = regions{rr};
        f   = fullfile(regionDir, regionFiles.(reg));
        if ~isfile(f), continue; end

        L = readmatrix(f); L = L(:); L = L(isfinite(L));
        if isempty(L), continue; end

        % allow 0-based IDs
        if any(L==0), L = L + 1; end

        L = unique(round(L));
        L = L(L>=1 & L<=nElem);
        regElems.(reg) = L;
    end

    % ---- Bin assignment within each region along X (quantiles recommended) ----
    binID_elem = struct();
    for rr = 1:numel(regions)
        reg = regions{rr};
        binID_elem.(reg) = zeros(nElem,1);

        L = regElems.(reg);
        if isempty(L), continue; end

        xR = Xe(L,1);

        if strcmpi(binMethod,'quantile')
            edges = quantile(xR, linspace(0,1,nBins+1));
        else
            edges = linspace(min(xR), max(xR), nBins+1);
        end
        edges = make_edges_strict(edges);

        b = discretize(xR, edges);
        b(isnan(b)) = 0;

        tmp = zeros(nElem,1);
        tmp(L) = b;
        binID_elem.(reg) = tmp;
    end

    % ---- Bin-averaged metrics per region ----
    binStats = struct();
    for rr = 1:numel(regions)
        reg = regions{rr};
        L   = regElems.(reg);
        if isempty(L), continue; end

        bAll = binID_elem.(reg);

        for m = 1:numel(metricsToAverage)
            met = metricsToAverage{m};
            binStats.(reg).(met).mean = nan(nBins,1);
            binStats.(reg).(met).n    = zeros(nBins,1);

            switch lower(met)
                case 'siglp',  v = sigLp;
                case 'sigcp',  v = sigCp;
                case 'difflc', v = diffLC;
                otherwise,     continue;
            end

            for b = 1:nBins
                idx = L(bAll(L) == b);
                vv  = v(idx);
                vv  = vv(isfinite(vv));
                if isempty(vv), continue; end

                binStats.(reg).(met).mean(b) = mean(vv);
                binStats.(reg).(met).n(b)    = numel(vv);
            end
        end
    end

    % ---- Store results ----
    R(r).ratID       = ratID;
    R(r).label       = label;
    R(r).meshID      = meshID;
    R(r).V           = V;
    R(r).E           = E;
    R(r).Xe          = Xe;
    R(r).eL          = eL;

    R(r).p_elem      = p_elem;
    R(r).sigmaL      = sigmaL;
    R(r).sigmaC      = sigmaC;
    R(r).sigLp       = sigLp;
    R(r).sigCp       = sigCp;
    R(r).diffLC      = diffLC;

    R(r).regElems    = regElems;
    R(r).binID_elem  = binID_elem;
    R(r).binStats    = binStats;

    valid(r) = true;

    % ---- Optional: write a compact CSV summary per rat ----
    outCSV = fullfile(outDir, sprintf('%s_A6_binMeans.csv', ratID));
    writeBinMeansCSV(outCSV, binStats, regions, metricsToAverage, nBins);
    fprintf('Saved: %s\n', outCSV);
end

fprintf('\nDone. Valid rats processed: %d/%d\n', sum(valid), nR);

%% ========================= LOCAL FUNCTIONS =========================
function edges = make_edges_strict(edges)
% Ensures strictly increasing bin edges (avoids discretize() failures)
edges = double(edges(:))';
for i = 2:numel(edges)
    if edges(i) <= edges(i-1)
        edges(i) = edges(i-1) + 1e-9;
    end
end
end

function writeBinMeansCSV(outCSV, binStats, regions, metrics, nBins)
% Writes one table with rows = (region, bin) and columns = mean metrics + n
regCol  = strings(nBins*numel(regions),1);
binCol  = nan(nBins*numel(regions),1);

T = table(regCol, binCol, 'VariableNames', {'Region','Bin'});

row = 0;
for rr = 1:numel(regions)
    reg = regions{rr};
    for b = 1:nBins
        row = row + 1;
        regCol(row) = string(reg);
        binCol(row) = b;
    end
end
T.Region = regCol;
T.Bin    = binCol;

for m = 1:numel(metrics)
    met = metrics{m};
    meanCol = nan(height(T),1);
    nCol    = nan(height(T),1);

    row = 0;
    for rr = 1:numel(regions)
        reg = regions{rr};
        if ~isfield(binStats, reg) || ~isfield(binStats.(reg), met)
            row = row + nBins; %#ok<AGROW>
            continue;
        end
        mu = binStats.(reg).(met).mean(:);
        nn = binStats.(reg).(met).n(:);

        meanCol(row+1:row+nBins) = mu;
        nCol(row+1:row+nBins)    = nn;
        row = row + nBins;
    end

    T.([met '_mean']) = meanCol;
    T.([met '_n'])    = nCol;
end

writetable(T, outCSV);
end