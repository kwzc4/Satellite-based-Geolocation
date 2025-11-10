%% S401_VisualizationEvaluation
co = colororder;
Scale = 0.7;
figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
geoplot(DevicePos(1),DevicePos(2),  'go', 'color',co(4,:),'MarkerSize',8,'LineWidth',2);hold on;   % Truth
geoplot(EstimPos(1),   EstimPos(2),    'x', 'color',co(2,:),'MarkerSize',10,'LineWidth',2);  % Estimate
geoplot(InitPos(1),  InitPos(2),  'o', 'color',co(5,:),'MarkerSize',8,'LineWidth',2);   % Initial
geoplot(Path.lat, Path.lon, 'k^:','LineWidth',1.5,'MarkerSize',0.5);
legend({'Truth','Estimate','Initial','Optimizer path'}, 'Location','northwest');
geobasemap streets;

ax = gca; 
ax.Scalebar.Visible = 'off';
ax.LatitudeAxis.TickValues = []; 
ax.LongitudeAxis.TickValues = [];
ax.LatitudeAxis.Label.String  = 'Latitude ^\circ';
ax.LongitudeAxis.Label.String = 'Longitude ^\circ';
ax.Box = 'on';
ax.Grid="on";
drawnow;                          % let geoaxes compute its limits
latlim = ax.LatitudeLimits;       % freeze them
lonlim = ax.LongitudeLimits;
geolimits(ax, latlim, lonlim);
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);

% --- after you plot Truth/Estimate/Initial ---
ell    = wgs84Ellipsoid('kilometer');              % distance in km
err_km = distance(DevicePos(1),DevicePos(2), ...
                  EstimPos(1),EstimPos(2), ell);
% Put a label near the estimated point
latOff = 0.2;  lonOff = 0.4;                     % small offset in degrees
text(DevicePos(1)+latOff, DevicePos(2)+lonOff, ...
     sprintf('Estimation error: %.4f km', err_km), ...
     'FontSize',10,'FontWeight','bold', ...
     'BackgroundColor','w','Margin',2);
%% Plot Optimization path
% Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
% geobasemap streets; hold on
% geoplot(DevicePos(1), DevicePos(2), 'go','MarkerSize',8,'LineWidth',1.2);
% geoplot(EstimPos(1), EstimPos(2), 'rx','MarkerSize',10,'LineWidth',1.5);
% geoplot(Path.lat, Path.lon, 'k^:','LineWidth',2.0,'MarkerSize',2);
% legend('Truth','Estimate','Optimizer path');
%% Evaluation
fprintf('Estimated LLA: [%.6f, %.6f, %.1f m]\n', EstimPos(1), EstimPos(2), EstimPos(3));
fprintf('True LLA:      [%.6f, %.6f, %.1f m]\n', DevicePos(1), DevicePos(2), DevicePos(3));
% True vs estimated (lat,lon in deg; alt in m)
latT = DevicePos(1); lonT = DevicePos(2); altT = DevicePos(3);
latE = EstimPos(1);  lonE = EstimPos(2);  altE = EstimPos(3);

% Horizontal error (great-circle) in Km (Mapping Toolbox)
hErr_km = deg2km(distance(latT, lonT, latE, lonE)) ;

% Vertical error (Km)
vErr_Km = (altE - altT)/1000;

% 3D error (Km): geodesy-accurate using ECEF
pT = lla2ecef([latT, lonT, altT]);   % [x y z] (m)
pE = lla2ecef([latE, lonE, altE]);
err3D_Km = norm(pE - pT)/1000;

fprintf('Horizontal error: %.4f Km | Vertical error: %.4f Km | 3D error: %.4f Km\n', ...
        hErr_km, vErr_Km, err3D_Km);




