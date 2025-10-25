%% S301_OptimizationEstimation
function [EstimPos, OptPath] = S301_OptimizeEstimate(Selector, sc, LeoSats, InitPos, walker, ...
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
%   'UseSingleSatDoppler' logical, default true  (false => multi-sat Doppler)
%   'FixAltitude'         logical, default true  (optimize lat/lon only)
%   'Display'             'iter'|'off' for fminsearch
%
% Output:
%   EstimPos : [lat lon alt]

% ---------- Parse options
p = inputParser;
p.addParameter('Weights', [1 1 1 1], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('UseSingleSatDoppler', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('FixAltitude', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('Display','iter', @ischar);
p.parse(varargin{:});
W = p.Results.Weights(:).';
useDoppSingle = p.Results.UseSingleSatDoppler;
fixAlt        = p.Results.FixAltitude;
dispMode      = p.Results.Display;

useRSS   = logical(Selector(1));
useTDoA  = logical(Selector(2));
useDOP   = logical(Selector(3));
useAoA   = logical(Selector(4));

% ---------- Path capture container ----------
OptPath.lat = []; OptPath.lon = []; OptPath.alt = []; OptPath.J = [];

% ---------- Optimizer
opts = optimset('Display', dispMode, 'OutputFcn', @outfun);

% Wrap the joint cost
if fixAlt
    x0 = [InitPos(1), InitPos(2)];
    costFun = @(xy) joint_cost([xy(1) xy(2) InitPos(3)]);
    [sol, ~] = fminsearch(costFun, x0, opts);
    EstimPos = [sol(1) sol(2) InitPos(3)];
else
    x0 = InitPos(:).';
    costFun = @(lla) joint_cost(lla);
    [sol, ~] = fminsearch(costFun, x0, opts);
    EstimPos = sol;
end

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
        RSS_meas = Meas.RSS;            % Measured P_r (dBm) incl. noise; size [Nt x Ns]
        VisMask  = Meas.VisMaskG;       % Visibility mask (ground vantage); logical [Nt x Ns]
        % ---- Forward model: compute ranges to all sats from the trial LLA ----
        gs_tmp = groundStation(sc,"Latitude",LLA(1),"Longitude",LLA(2),"Altitude",LLA(3));  % Create a temporary GS at trial position
        [~, ~, rho_mod, ~] = aer(gs_tmp, LeoSats);  % Model slant ranges ρ(t,s) from trial GS to satellites
        delete(gs_tmp);

        Ns = numel(LeoSats);
        if size(rho_mod,1)==Ns && size(rho_mod,2)~=Ns, rho_mod = rho_mod.'; end

        % ---- Free-space path loss (FSPL) and geometry-only model term ----
        Lfs = 20*log10(rho_mod) + 20*log10(Param.f) - 147.55;   % FSPL in dB for each (t,s)
        S   = -Lfs + Param.GRx;                                 % Geometry term up to a per-epoch offset K^k

        % ---- Accumulate residuals across all epochs/sats ----
        res_all = [];
        NtLocal = size(S,1);
        for k=1:NtLocal
            idx = VisMask(k,:);         % Visible satellites at epoch k
            if nnz(idx) < 2, continue; end  % Need ≥2 links to separate geometry from offset K^k
            m = RSS_meas(k,idx);        % Measured P_r for visible sats at epoch k
            s = S(k,idx);               % Modelled geometry term at epoch k
            good = isfinite(m) & isfinite(s);
            if any(good)
                % --- Profile out per-epoch offset K^k via least squares (closed form = mean residual) ---
                % Khat minimizes sum_s (m_s - (s_s + K))^2 over visible sats at epoch k
                Khat = mean(m(good) - s(good));
                res_all = [res_all, (s(good)+Khat) - m(good)]; %#ok<AGROW>  % stack residuals across all epochs and satellites and compute
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
        for k=1:size(TDOA,1)
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
        % Guard: no Doppler block present
        if ~isfield(Meas,'Doppler'), J = inf; return; end
    
        % Build LOS & radial velocity for all sats/epochs
        r_dev = reshape(lla2ecef(LLA), 1,1,3);
        r_dev = repmat(r_dev, Meas.Doppler.Nt, walker.LeoNum, 1);
    
        % LOS and radial velocity (single satellite over time)
        % posSat, velSat: [Nt x 3]
        LOS   = Meas.Doppler.posSat - r_dev;              % [Nt x 3]
        rngNm = sqrt(sum(LOS.^2, 3));                     % [Nt x 1]
        uLOS  = LOS ./ max(rngNm, eps);                   % [Nt x 3] (implicit expansion)
        vr    = sum(Meas.Doppler.velSat .* uLOS, 3);      % [Nt x LeoNum] m/s
        fDmod = -(Param.f/Const.c) .* vr;                 % [Nt x 1] Hz
        
        if useDoppSingle
            % single-sat with global bias
            satBest = Meas.Doppler.satBest;
            mask    = Meas.Doppler.VisMaskG(:,satBest);
            y       = Meas.Doppler.fD_single(mask);
            s       = fDmod(mask, satBest);
            if isempty(y), J = inf; return; end
            b0 = mean(y - s);
            r  = (s + b0) - y;
            J  = sqrt(mean(r.^2));                          % RMSE in Hz
        else
            % multi-sat with per-epoch bias
            fDmeas = Meas.Doppler.fD_multi;
            Vis    = Meas.Doppler.VisMaskG;
            res_all = [];
            for k=1:size(fDmeas,1)
                idx = Vis(k,:) & isfinite(fDmeas(k,:));
                if nnz(idx) < 2, continue; end
                y = fDmeas(k, idx);
                s = fDmod(k, idx);
                b = mean(y - s);                    % per-epoch bias
                res_all = [res_all, (s + b) - y]; %#ok<AGROW>
            end
            if isempty(res_all), J = inf; else, J = sqrt(mean(res_all.^2)); end % Hz
        end
    end


% ==================== AoA at satellite ====================
    function J = cost_aoa(LLA)
        if ~isfield(Meas,'AoA'), J = inf; return; end

        Az_meas = Meas.AoA.Az;
        El_meas = Meas.AoA.El;
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

% ---------- OutputFcn to capture optimizer path ----------
    function stop = outfun(x, optimValues, state)
        stop = false;  % never stop early here
        if strcmp(state,'iter')
            if fixAlt
                lat = x(1); lon = x(2); alt = InitPos(3);
            else
                lat = x(1); lon = x(2); alt = x(3);
            end
            OptPath.lat(end+1) = lat;
            OptPath.lon(end+1) = lon;
            OptPath.alt(end+1) = alt;
            if isfield(optimValues,'fval')
                OptPath.J(end+1)   = optimValues.fval;
            else
                % As a fallback, evaluate once (rarely needed)
                OptPath.J(end+1) = joint_cost([lat,lon,alt]);
            end
        end
    end

end 


% ===== helper: minimal signed angular difference in degrees =====
function d = angdiff_deg(a, b)
    d = atan2d(sind(a - b), cosd(a - b));  % wrap to (-180,180]
end