%% S401_Plot
co = colororder;
Scale = 0.7;
Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
geoplot(DevicePos(1),DevicePos(2),  'go', 'color',co(4,:),'MarkerSize',8,'LineWidth',2);hold on;   % Truth
geoplot(EstimPos(1),   EstimPos(2),    'x', 'color',co(2,:),'MarkerSize',10,'LineWidth',2);  % Estimate
geoplot(InitPos(1),  InitPos(2),  'o', 'color',co(5,:),'MarkerSize',8,'LineWidth',2);   % Initial
legend({'Truth','Estimate','Initial'}, 'Location','northwest');
geobasemap streets;

ax = gca; 
ax.Scalebar.Visible = 'off';
ax.LatitudeAxis.TickValues = []; 
ax.LongitudeAxis.TickValues = [];
ax.LatitudeAxis.Label.String  = 'Latitude ^\circ';
ax.LongitudeAxis.Label.String = 'Longitude ^\circ';
ax.Box = 'on';

drawnow;                          % let geoaxes compute its limits
latlim = ax.LatitudeLimits;       % freeze them
lonlim = ax.LongitudeLimits;
geolimits(ax, latlim, lonlim);
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);

