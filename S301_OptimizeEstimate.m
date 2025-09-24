%% S301_Estimation
function [EstimPos] = S301_OptimizeEstimate(Selector, sc, LeoSats, InitPos, ...
                            Meas, Param, Const, Device, varargin)
% Unified estimator for RSS / TDoA / Doppler / AoA with joint weighted cost.

% Selector = [useRSS useTDoA useDopp useAoA] (0/1 each)
% Meas.RSS            : [Nt x Ns] dBm (NaNs outside vis)
% Meas.VisMaskG       : [Nt x Ns] logical for ground-vantage visibility
% Meas.TDoA           : struct with .TDOA [Nt x Ns] s (NaN else), .refPerEpoch [Nt x 1], .VisMaskG
% Meas.Doppler        : struct with .fD_multi [Nt x Ns], .fD_single [Nt x 1], .satBest, .VisMaskG
% Meas.AoA            : struct with .Az_sat, .El_sat, .VisMaskS (satellite vantage)
% Param               : struct with f (Hz), GRx (dBi)
% Const               : struct with c (m/s)
% Device              : sigmas (sigmaRSS, sigmaT, sigmaf, sigmaAz, sigmaEl), optional DoppBias
%
% Name-Value:
%   'Weights'            [wRSS wTDoA wDopp wAoA], default = ones(1,4)
%   'Display'             'iter'|'off' for fminsearch
%
% Output:
%   EstimPos : [lat lon alt]

% ---------- Parse options
p = inputParser;
p.addParameter('Weights', [1 1 1 1], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('Display','iter', @ischar);
p.parse(varargin{:});
W = p.Results.Weights(:).';
dispMode      = p.Results.Display;

useRSS   = logical(Selector(1));
useTDoA  = logical(Selector(2));
useDOP   = logical(Selector(3));
useAoA   = logical(Selector(4));


% ---------- Optimizer
opts = optimset('Display', dispMode);

% Wrap the joint cost
x0 = [InitPos(1), InitPos(2)];
costFun = @(xy) joint_cost([xy(1) xy(2) InitPos(3)]);

[sol] = fminsearch(costFun, x0, opts);
EstimPos = [sol(1) sol(2) InitPos(3)];

% ==================== Nested: joint cost ====================
    function J = joint_cost(LLA)
        [Jrss, Jtdoa, Jdopp, Jaoa] = modality_costs(LLA);
        parts = [Jrss, Jtdoa, Jdopp, Jaoa];
        sel   = [useRSS, useTDoA, useDOP, useAoA] & (W>0);
        if ~any(sel)
            J = inf; return;
        end
        w = W(sel);
        q = parts(sel);
        % Weighted RMS of selected parts
        J = sqrt( sum( (w .* (q.^2)) ) / sum(w) );
    end

% ==================== Nested: per-modality costs ====================
    function [Jrss, Jtdoa, Jdopp, Jaoa] = modality_costs(LLA)
        Jrss  = NaN;  Jtdoa = NaN;  Jdopp = NaN;  Jaoa  = NaN;
        if useRSS
            Jrss = cost_rss(LLA);
        end
        if useTDoA
            Jtdoa = cost_tdoa(LLA);
        end
        if useDOP
            Jdopp = cost_doppler(LLA);
        end
        if useAoA
            Jaoa = cost_aoa(LLA);
        end
        % Replace NaNs with +inf to avoid contaminating joint cost
        if ~isfinite(Jrss),  Jrss  = inf; end
        if ~isfinite(Jtdoa), Jtdoa = inf; end
        if ~isfinite(Jdopp), Jdopp = inf; end
        if ~isfinite(Jaoa),  Jaoa  = inf; end
    end

% ==================== RSS: per-epoch K^k ====================
    function J = cost_rss(LLA)
        if ~isfield(Meas,'RSS') || ~isfield(Meas,'VisMaskG')
            J = inf; return;
        end
        RSS_meas = Meas.RSS;
        VisMask  = Meas.VisMaskG;

        gs_tmp = groundStation(sc,"Latitude",LLA(1),"Longitude",LLA(2),"Altitude",LLA(3));
        [~, ~, rho_mod, ~] = aer(gs_tmp, LeoSats);
        delete(gs_tmp);

        Ns = numel(LeoSats);
        if size(rho_mod,1)==Ns && size(rho_mod,2)~=Ns, rho_mod = rho_mod.'; end

        Lfs = 20*log10(rho_mod) + 20*log10(Param.f) - 147.55;
        S   = -Lfs + Param.GRx;   % model up to K^k

        res_all = [];
        Nt = size(S,1);
        for k=1:Nt
            idx = VisMask(k,:);
            if nnz(idx) < 2, continue; end
            m = RSS_meas(k,idx); s = S(k,idx);
            good = isfinite(m) & isfinite(s);
            if any(good)
                Khat = mean(m(good) - s(good));
                res_all = [res_all, (s(good)+Khat) - m(good)]; %#ok<AGROW>
            end
        end
        if isempty(res_all), J = inf; else, J = sqrt(mean(res_all.^2)); end % dB
    end

% ==================== TDoA: per-epoch reference ====================
    function J = cost_tdoa(LLA)
        if ~isfield(Meas,'TDoA') || ~isfield(Meas.TDoA,'TDOA')
            J = inf; return;
        end
        TDOA  = Meas.TDoA.TDOA;          % [Nt x Ns] NaN except vis + ref
        refPE = Meas.TDoA.refPerEpoch;   % [Nt x 1]
        VisG  = Meas.TDoA.VisMaskG;

        gs_tmp = groundStation(sc,"Latitude",LLA(1),"Longitude",LLA(2),"Altitude",LLA(3));
        [~, ~, rho_mod, ~] = aer(gs_tmp, LeoSats);
        delete(gs_tmp);

        Ns = numel(LeoSats);
        if size(rho_mod,1)==Ns && size(rho_mod,2)~=Ns, rho_mod = rho_mod.'; end
        TOA_mod = rho_mod ./ Const.c;

        res_all = [];
        Nt = size(TDOA,1);
        for k=1:Nt
            r = refPE(k);
            if isnan(r), continue; end
            idx = VisG(k,:);
            good = idx & isfinite(TDOA(k,:));
            if ~any(good), continue; end
            y_meas = TDOA(k, good);
            y_mod  = TOA_mod(k, good) - TOA_mod(k, r);
            res_all = [res_all, (y_mod - y_meas)]; %#ok<AGROW>
        end

        if isempty(res_all), J = inf; else, J = sqrt(mean(res_all.^2)); end % seconds
    end

% ==================== Doppler: single or multi ====================
    function J = cost_doppler(LLA)
        if ~isfield(Meas,'Doppler'), J = inf; return; end

        % Build LOS & radial velocity for all sats/epochs
        Nt = size(posSatAll,1); Ns = size(posSatAll,2);
        r_dev = reshape(lla2ecef(LLA), 1,1,3);
        r_dev = repmat(r_dev, Nt, Ns, 1);
        LOS   = posSatAll - r_dev;
        rngNm = sqrt(sum(LOS.^2,3));
        uLOS  = LOS ./ max(rngNm, eps);
        vr    = sum(velSatAll .* uLOS, 3);         % [Nt x Ns] m/s
        fDmod = -(Param.f/Const.c) .* vr;          % Hz

        % single-sat with global bias
        satBest = Meas.Doppler.satBest;
        mask    = Meas.Doppler.VisMaskG(:,satBest);
        y       = Meas.Doppler.fD_single(mask);
        s       = fDmod(mask, satBest);
        if isempty(y), J = inf; return; end
        b0 = mean(y - s);
        r  = (s + b0) - y;
        J  = sqrt(mean(r.^2));                 % Hz
    end

% ==================== AoA at satellite ====================
    function J = cost_aoa(LLA)
        if ~isfield(Meas,'AoA'), J = inf; return; end

        Az_meas = Meas.AoA.Az_sat;
        El_meas = Meas.AoA.El_sat;
        VisS    = Meas.AoA.VisMaskS;

        gs_tmp = groundStation(sc,"Latitude",LLA(1),"Longitude",LLA(2),"Altitude",LLA(3));
        [AZ_mod, EL_mod, ~, ~] = aer(LeoSats, gs_tmp);  % satellite vantage
        delete(gs_tmp);

        Ns = numel(LeoSats);
        if size(AZ_mod,1)==Ns && size(AZ_mod,2)~=Ns
            AZ_mod = AZ_mod.'; EL_mod = EL_mod.'; 
        end

        M = VisS & isfinite(Az_meas) & isfinite(El_meas);
        if ~any(M,'all'), J = inf; return; end

        az_res = angdiff_deg(AZ_mod(M), Az_meas(M));  % wrap diff
        el_res = EL_mod(M) - El_meas(M);

        % Optional noise-normalized residuals
        if isfield(Device,'sigmaAz') && isfield(Device,'sigmaEl')
            r = [az_res(:)/max(Device.sigmaAz,eps); el_res(:)/max(Device.sigmaEl,eps)];
        else
            r = [az_res(:); el_res(:)];
        end
        J = sqrt(mean(r.^2));                         % degrees (normalized if sigmas given)
    end

end 


% ===== helper: minimal signed angular difference in degrees =====
function d = angdiff_deg(a, b)
    d = atan2d(sind(a - b), cosd(a - b));  % wrap to (-180,180]
end
