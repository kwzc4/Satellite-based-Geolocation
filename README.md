# Satellite-based-Geolocation
End-to-end scripts for estimating a ground device's position from a LEO Walker-star constellation using four classical methods:

- **RSS** — Received Signal Strength (with per-epoch offset profiling, no fixed reference)
- **TDoA** — Time Difference of Arrival (per-epoch reference satellite)
- **Doppler** — Frequency shift; single-satellite over time (global bias)
- **AoA** — Angle of Arrival from azimuth/elevation geometry at the satellites

The pipeline:
1) build/load the constellation and times,
2) compute **per-epoch visibility**,
3) simulate **noisy measurements** for *all* modalities in one pass, and
4) run a **batch** optimizer that can use **any combination** of methods to
   return a single best-fit LLA point. Plots show **Truth vs Estimate vs Initial**.
5) Plots show Truth vs Estimate vs Initial.

---

## Repository layout (flat)

Main.m                          % Orchestration: runs the full pipeline
S101_Parameters.m               % Global constants, time window, noise, priors/search window
S102_CreateSatGeometry.m        % satelliteScenario + Walker-star constellation + true GS
S103_GeometryAndVisibility.m    % Per-epoch visibility (ground & satellite vantage) + Doppler sat pos/vel
S201_Measurements.m             % Simulate RSS, TDoA, Doppler, AoA (truth -> noisy)
S301_OptimizeEstimate.m         % Joint optimizer: choose methods via Selector, optional weights
S401_Plot.m                     % Geoplot: truth / estimate / initial + summary text
README.md



## Requirements

- **MATLAB** R2022b+ (R2025a recommended)
- **Toolboxes**
  - Satellite Communications Toolbox (`satelliteScenario`, `walkerStar`, `access`, `aer`, `states`)
  - Aerospace Toolbox (`lla2ecef`, etc.)
  - Mapping Toolbox (for `geoplot` basemaps)

---

## Quick start

1. Open **`S00_Parameters.m`** and review:
   - Simulation window: startTime, duration_min, sampleTime
   - RF & noise: Param.f, IoT/Device noise sigmas (sigmaRSS, sigmaT, sigmaf, sigmaAz, sigmaEl)
   - Search box / priors: IoTLat, IoTLon, IoTAlt, latMin/Max, lonMin/Max
   - Walker config: walker.NPlanes, walker.SatsPerPlane, walker.Inc, etc.
   - Elevation mask: minElDeg (default ~10°)
  
2. Run Main.m. By default it:
   - builds geometry,
   - computes visibility,
   - simulates measurements,
   - solves for the best-fit position,
   - plots results.
   
## How to choose methods

In Main.m set the Selector (on/off per method) and optional weights:

% [ RSS  TDoA  Doppler  AoA ]
Selector = [ 1    1      0     0 ];                 % enable RSS + TDoA
EstimPos = S301_OptimizeEstimate(Selector, sc, LeoSats, InitPos, ...
                                 MeasAll, Param, Const, Device, ...
                                 'Weights', [0.5 1 1 1], ...     % W_RSS, W_TDoA, W_Dop, W_AoA
                                 'FixAltitude', true, ...        % 2D solve at IoTAlt
                                 'Display','iter');              % or 'off'

Single method: set one element to 1 (e.g., [1 0 0 0] for RSS only).

Fusion: turn on multiple methods and optionally adjust 'Weights'.

## What each script does

#### S102_CreateSatGeometry.m

Creates a Walker-star constellation (satelliteScenario, walkerStar) and a true ground station at IoTPos.

#### S103_GeometryAndVisibility.m

Computes az/el/range histories two ways:

Ground vantage (aer(GS, sats)): for RSS, TDoA, Doppler. Visibility mask VisMaskG = ELg > minElDeg.

Satellite vantage (aer(sats, GS)): for AoA at satellites. Visibility mask VisMaskS = ELs < -minElDeg.
Also selects a Doppler satellite with the most visible epochs, and extracts its ECEF pos/vel time series via states.

#### S201_Measurements.m
Creates noisy measurements only where visible:

RSS: RSS = EIRP - FSPL + GRx + noise. Per-epoch constant later profiled out.

TDoA: build TOA from ranges; per-epoch reference = highest elevation sat; TDoA = TOA - TOA_ref.

Doppler (single sat): fD = -(f/c)·v_r + bias + noise, with v_r = dot(v_sat, u_LOS).

AoA: az/el at satellites with noise and angle wrapping.

#### S301_OptimizeEstimate.m
Joint cost = weighted sum of modality-specific RMSE terms. Options:

RSS: profiles one offset per epoch (cancels unknown EIRP & constants).

TDoA: RMSE of modeled vs measured TDoA using the same per-epoch visible set.

Doppler (single sat): profiles one global bias (e.g., LO offset).

AoA: bearing-only residuals with azimuth wrapping to (−180°, 180°].
Supports 2D solve ('FixAltitude', true) or full 3D.

#### S401_Plot.m
geoplot of Truth (green), Estimate (red X), Initial (blue) + error text.

## Model snippets (for reference)

Free-space path loss (dB):

FSPL = 20·log10(R[m]) + 20·log10(f[Hz]) − 147.55

RSS model: RSS ≈ −FSPL + GRx + K_k (K_k profiled out per epoch)

Doppler: f_D ≈ −(f/c)·v_r, v_r = dot(v_sat, u_LOS), u_LOS = (r_sat − r_dev)/||·||

TDoA: TDOA = TOA − TOA_ref, TOA = range/c

## Tips & troubleshooting

Visibility: If no satellites are visible, increase duration_min, shift startTime, or lower minElDeg (for debugging only).

Initial guess: InitPos is jittered around (IoTLat, IoTLon) with ~1σ and clamped to the search box. Tighter boxes → easier convergence.

Altitude: Fixing altitude ('FixAltitude', true) often stabilizes RSS/TDoA/Doppler.

Units: Param.f in Hz; ranges in meters; velocities in m/s; Doppler in Hz; angles in degrees.

Noise realism: Start with small sigmas to validate geometry; then scale to your scenario.

Weights: When fusing methods, begin with equal weights, then tune (e.g., increase TDoA weight if clocks are precise).

## License

MIT

## Acknowledgements

Built with MATLAB Satellite Communications, Aerospace, and Mapping toolboxes. Contributions & pull requests welcome!
