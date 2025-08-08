function GP_CBF_Delay_Postprocessing_DisturbanceMap

%% Arrow Map
rx = linspace(-1.5,1.5,20);
ry = linspace(-1.5,1.5,20);
[rx_grid,ry_grid] = meshgrid(rx,ry);
d_1_grid = nan(size(rx_grid));
d_2_grid = nan(size(rx_grid));
for i = 1:numel(rx_grid)
	x = zeros(4,1);
	x(1) = rx_grid(i);
	x(2) = ry_grid(i);
	d = ISSF_CBF_disturbance(x);
	d_1_grid(i) = d(1);
	d_2_grid(i) = d(2);
end
%
Disturbance_FigureObj = figure('Name','Disturbance');
%
Disturbance_AxesObj = subplot(2,2,[1;3],'Parent',Disturbance_FigureObj);
quiver(Disturbance_AxesObj,rx_grid,ry_grid,d_1_grid,d_2_grid);
axis(Disturbance_AxesObj,2 * [-1;1;-1;1]);
xlabel(Disturbance_AxesObj,'x');ylabel(Disturbance_AxesObj,'y');
%%
rx = linspace(-1.5,1.5,100);
ry = linspace(-1.5,1.5,100);
[rx_grid,ry_grid] = meshgrid(rx,ry);
d_1_grid = nan(size(rx_grid));
d_2_grid = nan(size(rx_grid));
for i = 1:numel(rx_grid)
	x = zeros(4,1);
	x(1) = rx_grid(i);
	x(2) = ry_grid(i);
	d = ISSF_CBF_disturbance(x);
	d_1_grid(i) = d(1);
	d_2_grid(i) = d(2);
end
%
d1_AxesObj = subplot(2,2,2,'Parent',Disturbance_FigureObj);
surf(d1_AxesObj,rx_grid,ry_grid,d_1_grid, ...
	'EdgeAlpha',0,'FaceColor','interp');
xlabel(d1_AxesObj,'x');ylabel(d1_AxesObj,'y');zlabel(d1_AxesObj,'d_1');
axis(d1_AxesObj,1.5 * [-1;1;-1;1]);
view(d1_AxesObj,-30,75)
%
d2_AxesObj = subplot(2,2,4,'Parent',Disturbance_FigureObj);
surf(d2_AxesObj,rx_grid,ry_grid,d_2_grid, ...
	'EdgeAlpha',0,'FaceColor','interp');
xlabel(d2_AxesObj,'x');ylabel(d2_AxesObj,'y');zlabel(d2_AxesObj,'d_2');
axis(d2_AxesObj,1.5 * [-1;1;-1;1]);
view(d2_AxesObj,-30,75)
% quiver(d1_AxesObj,rx_grid,ry_grid,d_1_grid,d_2_grid);
% axis(d1_AxesObj,2 * [-1;1;-1;1]);
% xlabel(d1_AxesObj,'x');ylabel(d1_AxesObj,'y');
end