%% S201_Measurements
function [MeasAll] = S201_Measurements( LeoSats, walker, Param, Device, Const, Nt, ELg, RHOg , VisMaskG, ...
    satIdx, DevicePos, VisCol, AZs, ELs, VisMaskS)
    
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
    [POS, VEL] = states(LeoSats, 'CoordinateFrame','ecef');  % 3 x Nt x Ns
    posSat = permute(POS, [2 3 1]);                           % [Nt x Ns x 3]
    velSat = permute(VEL, [2 3 1]);                           % [Nt x Ns x 3]

    % Device ECEF repeated to match dims
    r_dev = reshape(DevicePos_true, 1,1,3);
    r_dev = repmat(r_dev, Nt, walker.LeoNum, 1);
    LOS   = posSat - r_dev;                                   % [Nt x Ns x 3]
    rngNm = sqrt(sum(LOS.^2, 3));                             % [Nt x Ns]
    uLOS  = LOS ./ max(rngNm, eps);                           % unit LOS
    vr    = sum(velSat .* uLOS, 3);                           % [Nt x Ns] m/s
    bias  = 0; if isfield(Device,'DoppBias'), bias = Device.DoppBias; end
    fD_true  = -(Param.f/Const.c).*vr + bias;                 % Hz
    
    % Multi-Sat
    fD_multi = NaN(Nt, walker.LeoNum);
    fD_multi(VisMaskG) = fD_true(VisMaskG) + Device.sigmaf*randn(nnz(VisMaskG),1);

    % Choose single satellite with most visibility
    fD_single = NaN(Nt,1);
    fD_single(VisCol) = fD_multi(VisCol, satIdx)+ Device.sigmaf*randn(nnz(VisCol),1);


    Doppler_GT = struct('fD_multi',fD_multi,'fD_single',fD_single,'Nt',Nt, ...
        'satBest',satIdx, ...
        'VisMaskG',VisMaskG,...
        'VisCol', VisCol(:), ...             % [Nt x 1] logical
        'posSat', posSat, ...                % [Nt x 3] ECEF
        'velSat', velSat );                  % [Nt x 3] ECEF

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