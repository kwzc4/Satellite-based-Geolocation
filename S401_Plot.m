%% S401_Plot
figure; 
geoplot(DevicePos(1),DevicePos(2),  'go', 'MarkerSize',8,'LineWidth',2);hold on;   % Truth
geoplot(EstimPos(1),   EstimPos(2),    'rx', 'MarkerSize',10,'LineWidth',2);  % Estimate
geoplot(InitPos(1),  InitPos(2),  'bo', 'MarkerSize',8,'LineWidth',2);   % Initial
legend({'Truth','Estimate','Initial'}, 'Location','best');
geobasemap streets;
