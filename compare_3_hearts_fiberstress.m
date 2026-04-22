%% compare_sigmaF_three_females_pretty.m
% Defense-style layout:
%   Row 1: Female CONTROL   (3D, short-axis slice, histogram)
%   Row 2: Female SuHx W4   (3D, short-axis slice, histogram)
%   Row 3: Female SuHx W8   (3D, short-axis slice, histogram)
% Left 3D panels are recentered for nicer figure layout (not anatomical).
% Here we plot fiber-direction Cauchy stress σ_f instead of s_1.

clear; clc; close all;

%% --- USER SETTINGS ------------------------------------------------------
ratIDs    = {'W282W0','W282W4','W282W8'};   % control, W4, W8
ratLabels = {'Female control', 'Female SuHx W4', 'Female SuHx W8'};

geomDir    = fullfile('C:\Users\k1garcia\Desktop\PAH_FEA\ZeroLoadState_InverseFEA','FeBio');
gibbonRoot = 'C:\Users\k1garcia\Desktop\GIBBON-master\GIBBON-master';
addpath(genpath(gibbonRoot));

fontSize = 14;
slab_mm  = 0.8; %0.8  % slice half-thickness in mm

%% --- Load all rats & compute sigma_f + slices ---------------------------
nRats = numel(ratIDs);
R = cell(nRats,1);

for i = 1:nRats
    R{i} = get_sigmaF_and_slice_three(ratIDs{i}, geomDir, slab_mm);
end

% pooled color scale for σ_f (5–95% to avoid extreme outliers), symmetric
all_sf   = [];
allV_plot = [];
for i = 1:nRats
    all_sf    = [all_sf;    R{i}.sigma_f_nodes(:)];
    allV_plot = [allV_plot; R{i}.V_plot];     % centered coords for plotting
end

nBins_hist = 30;
xlo_hist   = -3.0;   % or min(all_sf)
xhi_hist   =  6.0;   % or max(all_sf)
histEdges  = linspace(xlo_hist, xhi_hist, nBins_hist+1);

% Color scale: 0 to global max positive sigma_f across all animals
pos_sf = all_sf(all_sf > 0)  % only positive values
if isempty(pos_sf)
    % fallback: if everything is <= 0, just use full range
    pos_sf = all_sf;
end

cminF = 0;
cmaxF = 4;

if cmaxF == 0
    cmaxF = 1e-6;   % avoid degenerate color axis
end

% shared symmetric axis limits in centered coordinates (for nice layout)
xmax = max(abs(allV_plot(:,1)));
ymax = max(abs(allV_plot(:,2)));
zmax = max(abs(allV_plot(:,3)));
pad = 0.05 * max([xmax ymax zmax]); % small padding
Lx = xmax + pad; Ly = ymax + pad; Lz = zmax + pad;
axL = [-Lx Lx -Ly Ly -Lz Lz];


%% --- Big 3x3 layout -----------------------------------------------------
figure('Color','w','Position',[100 100 1400 1000]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for i = 1:nRats
    Ri = R{i};

    % ---------- (1) 3D panel (using centered coords) ----------
    nexttile; hold on;
    %title(sprintf('%s (%s): \\sigma_f 3D', ratLabels{i}, Ri.ratname), ...
     %     'FontSize',fontSize);

    hp = gpatch(Ri.Fb, Ri.V_plot, Ri.sigma_f_nodes, 'k', 0.3);  % V_plot here
    hp.FaceColor  = 'interp';
    hp.EdgeColor  = 'none';
    hp.Marker     = '.';
    hp.MarkerSize = 4;

    axisGeom(gca,fontSize);
    axis(axL);                       % same centered box for all rows
    view(135,25);
    colormap(gjet(250));
    caxis([cminF cmaxF]);
    camlight headlight;
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % ---------- (2) short-axis slice (physical coords, unchanged) ----------
    nexttile; hold on;
    %title(sprintf('%s short-axis slice', ratLabels{i}), 'FontSize',fontSize);
    scatter(Ri.Xs, Ri.Ys, 12, Ri.sigma_f_slice, 'filled');
    axis equal tight;
    set(gca,'YDir','normal');
    xlabel('X (mm)'); ylabel('Y (mm)');
    colormap(gjet(250));
    caxis([cminF cmaxF]);
    grid on;

    % ---------- (3) histogram ----------
    nexttile; hold on;
    %title(sprintf('%s \\sigma_f distribution', ratLabels{i}), 'FontSize',fontSize);
    histogram(Ri.sigma_f_nodes, histEdges, ...
        'Normalization','probability', ...
        'FaceAlpha',0.7);
    xlabel('\sigma_f (kPa)');
    ylabel('Probability');
    % xlim([xlo xhi]);
    % ylim([0 ymax]);
    % Force axis limits big enough to show desired ticks
    xlim([-3 6]);      % this allows -2, -1, 0, 1, 2, 3 to appear
    ylim([0 0.15]);    % this allows 0.15 probability to appear
    yticks(0:0.05:0.15);
    xticks(-2:1:3);

    grid on;
end

% Global colorbar along the bottom
cb = colorbar('Location','southoutside');
cb.Layout.Tile = 'south';
ylabel(cb,'Fiber-direction Cauchy stress \sigma_f (kPa)');

fprintf('\n[3-female \\sigma_f comparison – pretty layout] complete.\n');

%% ===================== LOCAL FUNCTIONS ==================================
function R = get_sigmaF_and_slice_three(ratname, geomDir, slab_mm)
    % Load geometry + connectivity
    geomFile = fullfile(geomDir, [ratname '_geomStates.mat']);
    if ~isfile(geomFile)
        error('geomStates missing: %s', geomFile);
    end
    S = load(geomFile);

    V_ED = double(S.V_ED_final);
    E    = double(S.E);
    Fb   = double(S.Fb);
    f0   = double(S.f0);

    Fb = patchNormalFix(Fb);   % ensure consistent outward normals

    nNodes = size(V_ED,1);
    nElem  = size(E,1);

    % ---- full stress tensor log (for σ_f) ----
    log_stress = fullfile(geomDir, [ratname '_stress_out.txt']);
    if ~isfile(log_stress)
        error('Stress log missing: %s', log_stress);
    end
    dataS  = importFEBio_logfile(log_stress,0,1);
    S_all  = dataS.data;                    % [nElem x nComp x nSteps]
    [nElemS,nCompS,~] = size(S_all);
    if nElemS~=nElem
        error('Stress log nElem mismatch for %s', ratname);
    end
    S_last = S_all(:,:,end);

    sx = S_last(:,1);
    sy = S_last(:,2);
    sz = S_last(:,3);
    if nCompS >= 4
        sxy = S_last(:,4);
        syz = S_last(:,5);
        sxz = S_last(:,6);
    else
        sxy = zeros(nElem,1);
        syz = zeros(nElem,1);
        sxz = zeros(nElem,1);
    end

    % ---- fiber directions at elements ----
    if size(f0,1)==nElem
        f_elem = f0;
    elseif size(f0,1)==nNodes
        f_elem = (f0(E(:,1),:)+f0(E(:,2),:)+f0(E(:,3),:)+f0(E(:,4),:))/4;
    else
        error('f0 has incompatible dimensions for %s', ratname);
    end

    f_norm = max(eps, vecnorm(f_elem,2,2));
    f_elem = f_elem ./ f_norm;
    fx = f_elem(:,1); fy = f_elem(:,2); fz = f_elem(:,3);

    % σ_f = fᵀ σ f at elements
    sigma_f_elem = sx.*fx.^2 + sy.*fy.^2 + sz.*fz.^2 + ...
                   2*(sxy.*fx.*fy + sxz.*fx.*fz + syz.*fy.*fz);

    % element -> node + light smoothing
    sigma_f_nodes = elemToNodeMean_local(E, sigma_f_elem, nNodes);
    sigma_f_nodes = lapSmoothScalarField_local(E, sigma_f_nodes, 5, 0.5);

    % ---- recenter geometry for pretty plotting (does NOT affect stresses) ----
    ctr    = mean(V_ED,1);           % centroid
    V_plot = V_ED - ctr;             % centered coordinates for 3D panels

    % % % ---- short-axis slice (still using original coords for z) ----
    % % z_all  = V_ED(:,3);
    % % z0     = mean(z_all);
    % % inSlab = abs(z_all - z0) <= slab_mm;
    % % 
    % % Xs = V_ED(inSlab,1);
    % % Ys = V_ED(inSlab,2);
    % % sigma_f_slice = sigma_f_nodes(inSlab);

        % ---- short-axis slice (choose more basal slice) ----
    z_all = V_ED(:,3);

    % alphaBase: 0 = apex, 1 = base
    alphaBase = 0.85;   % try 0.7–0.85 for a clearly basal slice

    z_min = min(z_all);
    z_max = max(z_all);
    z0    = z_min + alphaBase*(z_max - z_min);

    inSlab = abs(z_all - z0) <= slab_mm;

    Xs = V_ED(inSlab,1);
    Ys = V_ED(inSlab,2);
    sigma_f_slice = sigma_f_nodes(inSlab);


    fprintf('%s: sigma_f range [%.3f , %.3f] kPa | slice nodes = %d\n', ...
        ratname, min(sigma_f_nodes), max(sigma_f_nodes), nnz(inSlab));

    % Pack up
    R.ratname        = ratname;
    R.V_ED           = V_ED;
    R.V_plot         = V_plot;   % centered version for plotting
    R.E              = E;
    R.Fb             = Fb;
    R.sigma_f_nodes  = sigma_f_nodes;
    R.Xs             = Xs;
    R.Ys             = Ys;
    R.sigma_f_slice  = sigma_f_slice;
end

function Vnod = elemToNodeMean_local(E, elemVals, nNodes)
    E        = double(E);
    elemVals = double(elemVals(:));
    Vsum   = accumarray(E(:), repmat(elemVals,4,1), [nNodes,1], @sum, 0);
    Vcount = accumarray(E(:), 1,                  [nNodes,1], @sum, 0);
    Vcount(Vcount==0) = 1;
    Vnod = Vsum ./ Vcount;
end

function field_s = lapSmoothScalarField_local(E, field, nIter, alpha)
    nNodes  = numel(field);
    field_s = field(:);

    E = double(E);
    pairs = [E(:,[1 2]); E(:,[1 3]); E(:,[1 4]); ...
             E(:,[2 3]); E(:,[2 4]); E(:,[3 4])];
    A = sparse(pairs(:,1),pairs(:,2),1,nNodes,nNodes);
    A = A + A.';
    A = A - diag(diag(A));

    deg   = sum(A,2);
    DinvA = spdiags(1./max(deg,1),0,nNodes,nNodes)*A;

    for k = 1:nIter
        field_s = (1-alpha)*field_s + alpha*(DinvA*field_s);
    end
end
