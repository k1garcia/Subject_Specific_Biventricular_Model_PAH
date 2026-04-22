%% compare_s1_three_females_pretty.m
% Defense-style layout:
%   Row 1: Female CONTROL   (3D, short-axis slice, histogram)
%   Row 2: Female SuHx W4   (3D, short-axis slice, histogram)
%   Row 3: Female SuHx W8   (3D, short-axis slice, histogram)
% Left 3D panels are recentered for nicer figure layout (not anatomical).

clear; clc; close all;

%% --- USER SETTINGS ------------------------------------------------------
ratIDs    = {'W282W0','W282W4','W282W8'};   % control, W4, W8
ratLabels = {'Female control', 'Female SuHx W4', 'Female SuHx W8'};

geomDir    = fullfile('C:\Users\k1garcia\Desktop\PAH_FEA\ZeroLoadState_InverseFEA','FeBio');
gibbonRoot = 'C:\Users\k1garcia\Desktop\GIBBON-master\GIBBON-master';
addpath(genpath(gibbonRoot));

fontSize = 14;
slab_mm  = 0.8;   % slice half-thickness in mm

%% --- Load all rats & compute s1 + slices -------------------------------
nRats = numel(ratIDs);
R = cell(nRats,1);

for i = 1:nRats
    R{i} = get_s1_and_slice_three(ratIDs{i}, geomDir, slab_mm);
end

% pooled color scale for s1 (5–95% to avoid outliers)
all_s1 = [];
allV_plot = [];
for i = 1:nRats
    all_s1    = [all_s1;    R{i}.s1_nodes(:)];
    allV_plot = [allV_plot; R{i}.V_plot];     % centered coords for plotting
end
p_lo   = prctile(all_s1,5);
p_hi   = prctile(all_s1,95);
cminS1 = p_lo;
cmaxS1 = p_hi;
if cminS1==cmaxS1, cmaxS1=cminS1+1e-6; end

% shared symmetric axis limits in centered coordinates (for nice layout)
xmax = max(abs(allV_plot(:,1)));
ymax = max(abs(allV_plot(:,2)));
zmax = max(abs(allV_plot(:,3)));
pad = 0.05 * max([xmax ymax zmax]); % small padding
Lx = xmax + pad;
Ly = ymax + pad;
Lz = zmax + pad;
axL = [-Lx Lx -Ly Ly -Lz Lz];

%% --- Big 3x3 layout -----------------------------------------------------
figure('Color','w','Position',[100 100 1400 1000]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for i = 1:nRats
    Ri = R{i};

    % ---------- (1) 3D panel (using centered coords) ----------
    nexttile; hold on;
    %title(sprintf('%s (%s): s_1 3D', ratLabels{i}, Ri.ratname), 'FontSize',fontSize);
    hp = gpatch(Ri.Fb, Ri.V_plot, Ri.s1_nodes, 'k', 0.3);  % V_plot here
    hp.FaceColor  = 'interp';
    hp.EdgeColor  = 'none';
    hp.Marker     = '.';
    hp.MarkerSize = 4;

    axisGeom(gca,fontSize);
    axis(axL);                       % same centered box for all rows
    view(135,25);
    colormap(gjet(250));
    caxis([cminS1 cmaxS1]);
    camlight headlight;
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % ---------- (2) short-axis slice (physical coords, unchanged) ----------
    nexttile; hold on;
    %title(sprintf('%s short-axis slice', ratLabels{i}), 'FontSize',fontSize);
    scatter(Ri.Xs, Ri.Ys, 12, Ri.s1_slice, 'filled');
    axis equal tight;
    set(gca,'YDir','normal');
    xlabel('X (mm)'); ylabel('Y (mm)');
    colormap(gjet(250));
    caxis([cminS1 cmaxS1]);
    grid on;

    % ---------- (3) histogram ----------
    nexttile; hold on;
    %title(sprintf('%s s_1 distribution', ratLabels{i}), 'FontSize',fontSize);
    histogram(Ri.s1_nodes, 40, ...
        'Normalization','probability', ...
        'FaceAlpha',0.7);
    xlabel('s_1 (kPa)');
    ylabel('Probability');
    xlim([-2 3.5]);      % this allows -2, -1, 0, 1, 2, 3 to appear
    ylim([0 0.2]);    % this allows 0.15 probability to appear
    yticks(0:0.05:0.2);
    xticks(-2:1:3.5);
    grid on;
end

% Global colorbar along the bottom
cb = colorbar('Location','southoutside');
cb.Layout.Tile = 'south';
ylabel(cb,'Max principal stress s_1 (kPa)');

fprintf('\n[3-female s1 comparison – pretty layout] complete.\n');

%% ===================== LOCAL FUNCTIONS ==================================
function R = get_s1_and_slice_three(ratname, geomDir, slab_mm)
    % Load geometry + connectivity
    geomFile = fullfile(geomDir, [ratname '_geomStates.mat']);
    if ~isfile(geomFile)
        error('geomStates missing: %s', geomFile);
    end
    S = load(geomFile);

    V_ED = double(S.V_ED_final);
    E    = double(S.E);
    Fb   = double(S.Fb);

    Fb = patchNormalFix(Fb);   % ensure consistent outward normals

    nNodes = size(V_ED,1);
    nElem  = size(E,1);

    % ---- principal stresses log ----
    log_princ = fullfile(geomDir, [ratname '_prinstress_out.txt']);
    dataP  = importFEBio_logfile(log_princ,0,1);
    P_elem = dataP.data;                    % [nElem x 3 x nSteps]
    [nElemP,nCompP,~] = size(P_elem);
    if nElemP~=nElem || nCompP~=3
        error('Principal stress log has unexpected size for %s', ratname);
    end
    s1_elem  = P_elem(:,1,end);             % max principal at final step

    % element -> node + light smoothing
    s1_nodes = elemToNodeMean_local(E, s1_elem, nNodes);
    s1_nodes = lapSmoothScalarField_local(E, s1_nodes, 5, 0.5);

    % ---- recenter geometry for pretty plotting (does NOT affect stresses) ----
    ctr   = mean(V_ED,1);           % centroid
    V_plot = V_ED - ctr;            % centered coordinates for 3D panels

    % ---- short-axis slice (still using original coords for z) ----
    z_all = V_ED(:,3);
    z0    = mean(z_all);
    inSlab = abs(z_all - z0) <= slab_mm;

    Xs = V_ED(inSlab,1);
    Ys = V_ED(inSlab,2);
    s1_slice = s1_nodes(inSlab);

    fprintf('%s: s1 range [%.3f , %.3f] kPa | slice nodes = %d\n', ...
        ratname, min(s1_nodes), max(s1_nodes), nnz(inSlab));

    % Pack up
    R.ratname   = ratname;
    R.V_ED      = V_ED;
    R.V_plot    = V_plot;   % centered version for plotting
    R.E         = E;
    R.Fb        = Fb;
    R.s1_nodes  = s1_nodes;
    R.Xs        = Xs;
    R.Ys        = Ys;
    R.s1_slice  = s1_slice;
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
