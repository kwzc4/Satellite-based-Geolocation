%% S201_Measurements
function [MeasAll] = S201_Measurements( walker, Param, Device, Const, Nt, ELg, RHOg , VisMaskG, ...
    DevicePos, SatPos, SatVel, VisCol, AZs, ELs, VisMaskS)
    
    % ---------- RSS (no reference; per-epoch offset profiled in estimator) ----------
    Lfs_true   = 20*log10(RHOg) + 20*log10(Param.f) - 147.55;   % [Nt x leoNum] dB
    EIRP_true  = Device.muEIRP + Device.sigmaEIRP*randn();               % one constant (will be cancelled by K^k)
    RSS_ideal  = EIRP_true - Lfs_true + Param.GRx;                 % [Nt x leoNum] dBm
    RSS_GT   = NaN(Nt, walker.LeoNum);
    RSS_GT(VisMaskG) = RSS_ideal(VisMaskG) + Device.sigmaRSS*randn(nnz(VisMaskG),1);
    
    % ---------- TDoA (per-epoch reference = highest elevation) ----------
    TOA_true = RHOg ./ Const.c;                                    % s
    TOA_meas = NaN(Nt, walker.LeoNum);
    TOA_meas(VisMaskG) = TOA_true(VisMaskG) + Device.sigmaT*randn(nnz(VisMaskG),1);
    
    refPerEpoch = nan(Nt,1);
    for k=1:Nt
        idx = find(VisMaskG(k,:));
        if numel(idx)>=2
            [~,imax] = max(ELg(k,idx));
            refPerEpoch(k) = idx(imax);
        end
    end
    
    TDOA = NaN(Nt, walker.LeoNum);
    for k=1:Nt
        r = refPerEpoch(k);
        if ~isnan(r)
            visk = VisMaskG(k,:);
            TDOA(k,visk) = TOA_meas(k,visk) - TOA_meas(k,r);
            % optional: include 0 at ref
            TDOA(k,r) = 0;
        end
    end
    TDoA_GT = struct('TDOA',TDOA,'refPerEpoch',refPerEpoch,'VisMaskG',VisMaskG);
    
    % ---------- Doppler (multi- and single-satellite) ----------
    DevicePos_true = lla2ecef(DevicePos);                     % 1 x 3 (ensure row)
    DevicePos_true = DevicePos_true(:).';                         % force 1x3
    
    LOS   = SatPos - DevicePos_true;                          % Nt x 3 (implicit expansion)
    rng   = sqrt(sum(LOS.^2,2));
    uLOS  = LOS ./ rng;
    vr    = sum(SatVel .* uLOS, 2);                       % Nt x 1  (m/s)
    fD_true = -(Param.f/Const.c).*vr + Device.DoppBias;   % Hz
    
    % Simulated measurements (only at visible epochs)
    fD_meas = NaN(Nt,1);
    fD_meas(VisCol) = fD_true(VisCol) + Device.sigmaf*randn(nnz(VisCol),1);

    Doppler_GT = struct( ...
        'fD',     fD_meas, ...               % [Nt x 1] Hz (NaN when not visible)
        'VisCol', VisCol(:), ...             % [Nt x 1] logical
        'posSat', SatPos, ...                % [Nt x 3] ECEF
        'velSat', SatVel );                  % [Nt x 3] ECEF

    % ---------- AoA at satellite (angles noisy, wrapped/clamped) ----------
    MeasAz = NaN(Nt, walker.LeoNum);
    MeasEl = NaN(Nt, walker.LeoNum);
    
    % Add noise then wrap/clamp
    Az_noisy = AZs + Device.sigmaAz*randn(Nt, walker.LeoNum);
    El_noisy = ELs + Device.sigmaEl*randn(Nt, walker.LeoNum);
    
    % Wrap azimuth to (-180,180] for consistency
    Az_noisy = wrapTo180(Az_noisy);
    % Clamp elevation to [-90,90]
    El_noisy = max(min(El_noisy, 90), -90);
    
    MeasAz(VisMaskS) = Az_noisy(VisMaskS);
    MeasEl(VisMaskS) = El_noisy(VisMaskS);
    AoA_GT = struct('Az',MeasAz,'El',MeasEl,'VisMaskS',VisMaskS);

    % Build the big measurement struct once (from your S201_Measurements)
    MeasAll.RSS     = RSS_GT;            % [Nt x Ns] dBm (NaNs outside vis)
    MeasAll.TDoA    = TDoA_GT;           % struct with .TDOA, .TOA, .refPerEpoch, .VisMaskG
    MeasAll.Doppler = Doppler_GT;        % struct with .fD_multi, .fD_single, .satBest, .VisMaskG
    MeasAll.AoA     = AoA_GT;            % struct with .Az_sat, .El_sat, .VisMaskS
    MeasAll.VisMaskG = TDoA_GT.VisMaskG; % reuse mask
    MeasAll.VisMaskS = AoA_GT.VisMaskS;

end




