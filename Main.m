%% Main
clc; clear; close all;
%% STEP 1 — Load parameters, create geometry, and Ground truth
S101_Parameters
[sc,LeoSats,GSTrue] = S102_CreateSatGeometry(walker,startTime,stopTime,sampleTime,DevicePos);
S103_GeometryAndVisibility % 
%% STEP 2 — Measurement Generation
MeasAll =S201_Measurements( LeoSats, walker, Param, Device, Const, Nt, ELg, RHOg , VisMaskG, ...
    satIdx, DevicePos, VisCol, AZs, ELs, VisMaskS);
%% STEP 3 — Optimization and estimation
Selector = [0 1 0 0]; % [RSS / TDoA / Doppler / AoA]
EstimPos = S301_OptimizeEstimate(Selector, sc, LeoSats, InitPos, walker,...
                       MeasAll, Param, Const, Device,'Weights', [1 1 1 1], ...
                       'UseSingleSatDoppler', true, 'FixAltitude', true);
%% STEP 4 — Visualization and Evaluation
S401_VisualizationEvaluation;
%% Save Figure
% saveLocation = 'C:\Users\nermi\OneDrive - RMIT University\28. Bisma Project\6. Paper\1. Figures';
% exportgraphics(Fig, fullfile(saveLocation,'Figure15.png'), 'Resolution',600);