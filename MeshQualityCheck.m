%% check_all_meshes.m
% Batch quality check for all Gmsh .msh files in a folder
% Skips any meshes that fail to load or parse.

clear; clc;

meshDir = '';  % leave '' to choose interactively

if isempty(meshDir)
    meshDir = uigetdir(pwd, 'Select folder containing .msh files');
    if meshDir == 0
        error('No folder selected. Exiting.');
    end
end

mshFiles = dir(fullfile(meshDir, '*.msh'));
if isempty(mshFiles)
    error('No .msh files found in folder: %s', meshDir);
end

fprintf('Found %d .msh files in %s\n\n', numel(mshFiles), meshDir);

for k = 1:numel(mshFiles)
    fname = fullfile(meshDir, mshFiles(k).name);
    fprintf('\n=============================================\n');
    fprintf('Checking mesh %d of %d: %s\n', k, numel(mshFiles), mshFiles(k).name);
    fprintf('=============================================\n');

    try
        meshQualityReport(fname);
    catch ME
        fprintf('\n⚠️  Skipping file due to error: %s\n', ME.message);
        fprintf('   File skipped: %s\n\n', mshFiles(k).name);
        continue; % move to next mesh
    end

    fprintf('\nPress any key for next mesh, or Ctrl+C to stop.\n');
    pause;
end

fprintf('\nDone processing all meshes.\n');


function meshQualityReport(mshFile)
% meshQualityReport
%   Reads a Gmsh v4 ASCII .msh file, extracts nodes + tets, computes
%   element quality and aspect ratio, and generates a basic visualization.
%
%   Usage:
%       meshQualityReport('myMesh.msh')

    if ~isfile(mshFile)
        error('File not found: %s', mshFile);
    end

    %% --- Read mesh from Gmsh v4 ASCII ---
    fprintf('Reading mesh from %s ...\n', mshFile);
    [V, E] = readMshGmsh_ascii(mshFile);    % NEW (handles v2 or v4)

    nNodes = size(V,1);
    nTets  = size(E,1);
    fprintf('  Nodes: %d\n', nNodes);
    fprintf('  Tets : %d\n', nTets);

    if nTets == 0
        warning('No tetrahedral elements found. Skipping quality computation.');
        return;
    end

    %% --- Compute basic quality metrics ---
    fprintf('Computing quality metrics ...\n');

    [q, aspectRatio, vol] = tetQualitySimple(V,E);

    % Summary stats
    q_min   = min(q);
    q_5     = prctile(q,5);
    q_med   = median(q);
    q_mean  = mean(q);

    ar_max  = max(aspectRatio);
    ar_95   = prctile(aspectRatio,95);

    vol_min = min(vol);
    vol_med = median(vol);
    vol_max = max(vol);

    %% --- Print summary to command window ---
    fprintf('\n--- Element Quality (radius-ratio style, 0–1, higher = better) ---\n');
    fprintf('  min      = %.3f\n', q_min);
    fprintf('  5th pct  = %.3f\n', q_5);
    fprintf('  median   = %.3f\n', q_med);
    fprintf('  mean     = %.3f\n', q_mean);

    % Threshold warnings
    frac_bad_0p1 = sum(q < 0.10) / numel(q);
    frac_bad_0p2 = sum(q < 0.20) / numel(q);
    fprintf('  %% elements q < 0.10: %.2f %%\n', 100*frac_bad_0p1);
    fprintf('  %% elements q < 0.20: %.2f %%\n', 100*frac_bad_0p2);

    fprintf('\n--- Aspect Ratio (max edge / min edge, lower = better) ---\n');
    fprintf('  max      = %.2f\n', ar_max);
    fprintf('  95th pct = %.2f\n', ar_95);

    fprintf('\n--- Volume (in mesh units^3) ---\n');
    fprintf('  min      = %.3e\n', vol_min);
    fprintf('  median   = %.3e\n', vol_med);
    fprintf('  max      = %.3e\n', vol_max);

    %% --- Quick plots ---

    % Histogram of quality
    figure('Name', ['Quality histogram: ' mshFile], 'Color', 'w');
    histogram(q, 30);
    xlabel('Element quality (0–1, higher = better)');
    ylabel('Count');
    title(sprintf('Quality histogram: %s', mshFile), 'Interpreter', 'none');

    % 3D plot colored by nodal-averaged quality on outer surface
    fprintf('Computing surface and nodal-averaged quality for plotting ...\n');

    % 1) Compute nodal quality as average of connected element qualities
    nNodes = size(V,1);
    nTets  = size(E,1);
    q_node = zeros(nNodes,1);
    counts = zeros(nNodes,1);

    for ei = 1:nTets
        nodes = E(ei,:);
        q_e   = q(ei);
        q_node(nodes) = q_node(nodes) + q_e;
        counts(nodes) = counts(nodes) + 1;
    end

    nonzero = counts > 0;
    q_node(nonzero) = q_node(nonzero) ./ counts(nonzero);

    % 2) Get outer surface triangles using boundary on node coordinates
    shrinkFactor = 0.9;  % must be in (0,1], lower -> tighter wrap
    K = boundary(V(:,1), V(:,2), V(:,3), shrinkFactor);  % K: [nSurfTri x 3]


    % 3D plot
    figure('Name', ['Quality map: ' mshFile], 'Color', 'w');
    hold on;
    title(sprintf('Surface quality map (nodal avg): %s', mshFile), ...
          'Interpreter', 'none');
    axis equal; axis off;

    patch('Faces', K, ...
          'Vertices', V, ...
          'FaceVertexCData', q_node, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    colorbar;
    colormap(parula);
    caxis([0 1]);
    camlight; lighting gouraud;

end
function [V,E] = readMshGmsh_ascii(filename)
% readMshGmsh_ascii
%   Minimal reader for Gmsh ASCII .msh files.
%   Automatically detects version (2.x or 4.x) and calls
%   the appropriate parser.
%
%   Returns:
%       V: [nNodes x 3] node coordinates
%       E: [nTets x 4]  tetrahedral elements (node indices, type 4)

    fid = fopen(filename,'r');
    if fid < 0
        error('Could not open file: %s', filename);
    end

    gmshVersion = [];
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if strcmp(line,'$MeshFormat')
            fmtLine = strtrim(fgetl(fid));  % e.g. "4.1 0 8"
            vals    = sscanf(fmtLine,'%f');
            gmshVersion = vals(1);
            break;
        end
    end
    fclose(fid);

    if isempty(gmshVersion)
        error('Could not find $MeshFormat in file: %s', filename);
    end

    fprintf('  Detected Gmsh version: %.2f\n', gmshVersion);

    if gmshVersion < 3
        [V,E] = readGmsh2_core(filename);
    else
        [V,E] = readGmsh4_core(filename);
    end
end

%% ================== Gmsh 4.x core parser ==================
function [V,E] = readGmsh4_core(filename)

    fid = fopen(filename,'r');
    if fid < 0
        error('Could not open file: %s', filename);
    end

    V = [];
    E = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if strcmp(line,'$Nodes')
            V = parseNodesV4(fid);
        elseif strcmp(line,'$Elements')
            E = parseElementsV4(fid);
        end
    end

    fclose(fid);

    if isempty(V)
        error('No nodes read from file: %s', filename);
    end
    if isempty(E)
        warning('No tetrahedral elements (type 4) found in file: %s', filename);
        E = zeros(0,4);
    end
end

function V = parseNodesV4(fid)
    % First line after $Nodes:
    %   numEntityBlocks numNodes minNodeTag maxNodeTag
    hdr = fscanf(fid,'%d %d %d %d',4);
    numEntityBlocks = hdr(1);
    numNodesTotal   = hdr(2);

    V = zeros(numNodesTotal,3);

    for b = 1:numEntityBlocks
        % entityDim entityTag parametric numNodesInBlock
        blk = fscanf(fid,'%d %d %d %d',4);
        numNodesInBlock = blk(4);

        % Node tags (ints, one per line but fscanf doesn't care)
        nodeTags = fscanf(fid,'%d',numNodesInBlock);

        % Coordinates: x y z for each node
        coords = fscanf(fid,'%f',3*numNodesInBlock);
        coords = reshape(coords,3,numNodesInBlock).';

        V(nodeTags,:) = coords;
    end

    % consume "$EndNodes"
    endTag = strtrim(fgetl(fid)); %#ok<NASGU>
end

function E = parseElementsV4(fid)
    % First line after $Elements:
    %   numEntityBlocks numElements minElementTag maxElementTag
    hdr = fscanf(fid,'%d %d %d %d',4);
    numEntityBlocks = hdr(1);
    numElemsTotal   = hdr(2); %#ok<NASGU>

    tetList = [];

    for b = 1:numEntityBlocks
        % entityDim entityTag elementType numElementsInBlock
        blk = fscanf(fid,'%d %d %d %d',4);
        elemType        = blk(3);
        numElemsInBlock = blk(4);

        % number of nodes per element for this type
        nPerElem = getNumNodesV4(elemType);

        % read all ints for this block:
        %   elemTag node1 ... nodeN  (repeated)
        data = fscanf(fid,'%d',(1+nPerElem)*numElemsInBlock);
        data = reshape(data,1+nPerElem,numElemsInBlock).';

        if elemType == 4
            % keep tets only: columns 2–5
            tetList = [tetList; data(:,2:end)]; %#ok<AGROW>
        end
        % other element types are skipped but still consumed
    end

    E = tetList;

    % consume "$EndElements"
    endTag = strtrim(fgetl(fid)); %#ok<NASGU>
end

function n = getNumNodesV4(elemType)
    % minimal mapping for types seen in your meshes
    switch elemType
        case 1      % 2-node line
            n = 2;
        case 2      % 3-node triangle
            n = 3;
        case 3      % 4-node quad
            n = 4;
        case 4      % 4-node tetra
            n = 4;
        case 5      % 8-node hex (not used here, but safe)
            n = 8;
        otherwise
            % You can extend this if you ever use higher-order elements.
            error('Unhandled element type in v4 parser: %d', elemType);
    end
end

%% ================== Gmsh 2.x core parser ==================
function [V,E] = readGmsh2_core(filename)

    fid = fopen(filename,'r');
    if fid < 0
        error('Could not open file: %s', filename);
    end

    V = [];
    E = [];

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if strcmp(line,'$Nodes')
            V = parseNodesV2(fid);
        elseif strcmp(line,'$Elements')
            E = parseElementsV2(fid);
        end
    end

    fclose(fid);

    if isempty(V)
        error('No nodes read from file: %s', filename);
    end
    if isempty(E)
        warning('No tetrahedral elements (type 4) found in file: %s', filename);
        E = zeros(0,4);
    end
end

function V = parseNodesV2(fid)
    % First line after $Nodes: numNodes
    numNodes = fscanf(fid,'%d',1);
    V = zeros(numNodes,3);

    for i = 1:numNodes
        % nodeTag x y z
        line = fgetl(fid);
        if ~ischar(line), break; end
        nums = sscanf(line,'%f');
        if numel(nums) < 4
            continue;
        end
        tag    = round(nums(1));
        coords = nums(2:4).';
        if tag >= 1 && tag <= numNodes
            V(tag,:) = coords;
        end
    end

    % consume "$EndNodes"
    endTag = strtrim(fgetl(fid)); %#ok<NASGU>
end

function E = parseElementsV2(fid)
    % First line after $Elements: numElements
    numElems = fscanf(fid,'%d',1);

    tetList = [];

    for i = 1:numElems
        % elm-number elm-type number-of-tags <tags> node1 ... nodeN
        line = fgetl(fid);
        if ~ischar(line), break; end
        nums = sscanf(line,'%d');
        if numel(nums) < 4
            continue;
        end

        elemType = nums(2);
        nTags    = nums(3);

        if elemType ~= 4
            % skip non-tets
            continue;
        end

        % last 4 entries are node indices for tet4
        nodeIds = nums(end-3:end).';
        tetList = [tetList; nodeIds]; %#ok<AGROW>
    end

    E = tetList;

    % consume "$EndElements"
    endTag = strtrim(fgetl(fid)); %#ok<NASGU>
end

function [q, aspectRatio, Vtet] = tetQualitySimple(V,E)
% tetQualitySimple
%   Compute simple quality metrics for 4-node tetrahedra:
%       q           - radius-ratio-style quality (0–1, higher = better)
%       aspectRatio - max edge length / min edge length
%       Vtet        - element volumes
%
%   Good for flagging slivers or near-degenerate elements.

    % Node coordinates per tet
    V1 = V(E(:,1),:);
    V2 = V(E(:,2),:);
    V3 = V(E(:,3),:);
    V4 = V(E(:,4),:);

    % Edge vectors
    e12 = V2 - V1;
    e13 = V3 - V1;
    e14 = V4 - V1;
    e23 = V3 - V2;
    e24 = V4 - V2;
    e34 = V4 - V3;

    % Edge lengths
    L12 = sqrt(sum(e12.^2,2));
    L13 = sqrt(sum(e13.^2,2));
    L14 = sqrt(sum(e14.^2,2));
    L23 = sqrt(sum(e23.^2,2));
    L24 = sqrt(sum(e24.^2,2));
    L34 = sqrt(sum(e34.^2,2));

    allL = [L12, L13, L14, L23, L24, L34];
    Lmin = min(allL,[],2);
    Lmax = max(allL,[],2);

    % Tet volume
    Vtet = abs(dot(e12, cross(e13,e14,2), 2)) / 6;

    % Aspect ratio
    aspectRatio = Lmax ./ Lmin;

    % Simple radius-ratio-like quality
    q = (6*sqrt(2) * Vtet) ./ (Lmax.^3);
    q(Vtet <= 0) = 0;  % inverted or zero-volume

    q(q < 0) = 0;
    q(q > 1) = 1;
end
