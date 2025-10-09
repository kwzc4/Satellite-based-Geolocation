%% Geometry and visibilty calculation
%% Find distance
% Ground vantage (positive elevations): convenient for RSS/TDoA/Doppler masks
[AZg, ELg, RHOg, ~] = aer(GSTrue, LeoSats); % fix
if size(ELg,1)==walker.LeoNum && size(ELg,2)~=walker.LeoNum
    AZg = AZg.'; ELg = ELg.'; RHOg = RHOg.';     % [Nt x LeoNum]
end
VisMaskG = (ELg > minElDeg);                      % [Nt x LeoNum]

% Satellite vantage (negative elevations): needed for AoA-at-satellite
[AZs, ELs, ~, ~] = aer(LeoSats, GSTrue);
if size(ELs,1)==walker.LeoNum && size(ELs,2)~=walker.LeoNum
    AZs = AZs.'; ELs = ELs.';                    % [Nt x LeoNum]
end
VisMaskS = (ELs < -minElDeg);                    % [Nt x LeoNum]

% Debug
nVisG = sum(VisMaskG,2);
fprintf('Ground-vantage visible sats (min/med/max): %d / %d / %d (Total time steps=%d, Total number of satellites=%d)\n', ...
    min(nVisG), median(nVisG), max(nVisG), size(ELg,1), walker.LeoNum);
%% For TDoA — per-epoch visible set and reference (highest elevation)
visIdx      = cell(Nt,1);           % indices of visible sats per epoch
refPerEpoch = nan(Nt,1);            % reference sat per epoch
for k = 1:Nt
    idx_k = find(VisMaskG(k,:));
    if numel(idx_k) >= 2
        [~, imax] = max(ELg(k, idx_k));  % <-- use ELg here (not EL_all)
        refPerEpoch(k) = idx_k(imax);
        visIdx{k}      = idx_k;
    end
end
%% For Doppler — pick the satellite with the most visible samples
nVisPerSat = sum(VisMaskG,1);
[bestCount, satIdx] = max(nVisPerSat);
fprintf('Chosen sat %d with %d visible epochs.\n', satIdx, bestCount);
if bestCount==0
    error('No visible satellite in window. Adjust time or elevation mask.');
end
VisCol = VisMaskG(:,satIdx);        % [Nt x 1] visibility of the chosen sat
%% Now get position/velocity for the most visible satellite
[pos, vel] = states(LeoSats(satIdx), 'CoordinateFrame','ecef');  % pos/vel: 3 x Nt
SatPos = pos.';    % Nt x 3
SatVel = vel.';    % Nt x 3