%% S101_Parameters
%% Physical Constants
Const.E = wgs84Ellipsoid('meters');     % Reference ellipsoid (WGS-84)
Const.c = physconst('LightSpeed');
Const.kb = physconst('Boltzmann');      % Boltzmann constant [J/K]
Const.TempK = 293;                      % System noise temperature [K]
%% General Simulation Parameters
fprintf('Initializing simulation parameters ...\n');
startTime = datetime(2025, 7, 12, 0, 0, 0); % Simulation start
duration = 10;                              % 10 min simulation in seconds
sampleTime = 1;                             % Time step in seconds
stopTime = startTime + minutes(duration);
ts = startTime:seconds(sampleTime):stopTime;
Nt = numel(ts);                             % Time vector
%% LEO Walker-Star Constellation Parameters
walker.a = 500e3 + earthRadius;    % Semi-major axis
walker.alfa = earthRadius / walker.a;
walker.Inc = 87.9;                 % Inclination in degrees (typical for OneWeb)
walker.NPlanes = 12;               % Number of orbital planes 
walker.SatsPerPlane = 12;          % Number of satellites per plane 
walker.PhaseOffset = 8;            % Phase offset for phasing between planes
walker.LeoNum = walker.NPlanes * walker.SatsPerPlane;
Ref = 1;                           % Reference satellite
minElDeg = 10;                     % small cutoff (tweak if needed)
% add maximum nuber of satellites ==> 1 to inf (use all)
% epoch length ==> see the satellite which has visibility every  for
% example - divide the duration into number of epochs as an option in the
% optimization
% TDoA, 
%% Radar parameters
Param.f = 2e9;                     % 2 GHz uplink
Param.muLoS = 0;
Param.sigmaLoS = 2;
Param.muNLoS = 10;
Param.sigmaNLoS = 8;
Param.beta = 0.3;                  % Urban environment constant, can vary 0.0-0.57
Param.sigmaD = 10;                 % Doppler noise std deviation [Hz]
Param.GRx = 6;                     % dBi
%% Device Power Parameters
Device.Number = 1;                    % single unknown device
Device.muEIRP = 14;                    % Mean of RF Intensity
Device.sigmaEIRP = 2;                  % Standard deviation of RF Intensity
Device.sigmaT = 100e-9;                % Standard deviation of time
Device.sigmaf = 2;                     % Standard deviation of Doppler frequency
Device.sigmaRSS =2;                    % Standard deviation of RSS
Device.sigmaAz = 1.0;                  % Standard deviation of azimuth error (degrees)
Device.sigmaEl = 1.0;                  % Standard deviationf elevation error (degrees)
Device.DoppBias = 20;                  % Hz, constant bias (for simulation)
%% Device Position
DeviceLat = -37.8136;
DeviceLon = 144.9631;
DeviceAlt = 0;
DevicePos = [DeviceLat, DeviceLon, DeviceAlt];
latMin = DeviceLat - 2;
latMax = DeviceLat + 2;
lonMin = DeviceLon - 2;
lonMax = DeviceLon + 2;
% Use the 3σ rule: ±3σ spans the box  ⇒  σ = width/6
SigmaLat = (latMax - latMin)/6;
SigmaLon = (lonMax - lonMin)/6;
% Initial location guess (1σ jitter, clamped)
lat0 = min(max(DeviceLat + SigmaLat*randn(), latMin), latMax);
lon0 = min(max(DeviceLon + SigmaLon*randn(), lonMin), lonMax);
InitPos = [lat0, lon0, DeviceAlt];