%% GP_CBF_Delay_BatchTest_OneComparison
clc; clear all; close all;
rng(0);
%%
SaveFolderName = 'results';
%% Configuration
SimulationConfigurationSet = cell(4,1);
%
SimulationConfigurationSet{1}.do_Consider_GPError = false;
SimulationConfigurationSet{1}.do_Consider_Delay = false;
SimulationConfigurationSet{1}.Name = 'NoGP_NoDelay';
SimulationConfigurationSet{1}.Color = [1 0 1];
%
SimulationConfigurationSet{2}.do_Consider_GPError = true;
SimulationConfigurationSet{2}.do_Consider_Delay = false;
SimulationConfigurationSet{2}.Name = 'WithGP_NoDelay';
SimulationConfigurationSet{2}.Color = [0 0 1];
%
SimulationConfigurationSet{3}.do_Consider_GPError = false;
SimulationConfigurationSet{3}.do_Consider_Delay = true;
SimulationConfigurationSet{3}.Name = 'NoGP_WithDelay';
SimulationConfigurationSet{3}.Color = [0 1 0];
%
SimulationConfigurationSet{4}.do_Consider_GPError = true;
SimulationConfigurationSet{4}.do_Consider_Delay = true;
SimulationConfigurationSet{4}.Name = 'WithGP_WithDelay';
SimulationConfigurationSet{4}.Color = [1 0 0];
%% Simulation
do_simulation = true;
if do_simulation
	for SimulationConfigurationNr = 1:numel(SimulationConfigurationSet)
		SimulationConfiguration = SimulationConfigurationSet{SimulationConfigurationNr};
		do_Consider_GPError = SimulationConfiguration.do_Consider_GPError;
		do_Consider_Delay = SimulationConfiguration.do_Consider_Delay;
		SaveFileName = [SimulationConfiguration.Name];
		%
		GP_CBF_Delay_TestFunc(do_Consider_GPError,do_Consider_Delay, ...
			SaveFolderName,SaveFileName);
	end
end
%% Plotting
RawResult = load([SaveFolderName,'/',SimulationConfigurationSet{1}.Name,'.mat']);
t_set = RawResult.t_set;
x_ref_set = GP_CBF_Delay_Reference(t_set);
%
x_set_all = cell(numel(SimulationConfigurationSet),1);
h_set_all = nan(numel(SimulationConfigurationSet),numel(t_set));
h_star_set_all = nan(numel(SimulationConfigurationSet),numel(t_set));
delta_h_all = nan(numel(SimulationConfigurationSet),1);
PredictionError_set_all = nan(numel(SimulationConfigurationSet),numel(t_set));
eta_GP_set_all = nan(numel(SimulationConfigurationSet),numel(t_set));
eta_Delay_set_all = nan(numel(SimulationConfigurationSet),numel(t_set));
for SimulationConfigurationNr = 1:numel(SimulationConfigurationSet)
	SimulationConfiguration = SimulationConfigurationSet{SimulationConfigurationNr};
	SaveFileName = SimulationConfiguration.Name;
	RawResult = load([SaveFolderName,'/',SaveFileName,'.mat']);
	x_set = RawResult.x_set(:,1:numel(t_set));
	delta_h = RawResult.delta_h;
	mu_set = RawResult.mu_set;
	eta_GP_set = RawResult.eta_GP_set;
	eta_Delay_set = RawResult.eta_Delay_set;
	%
	x_set_all{SimulationConfigurationNr} = x_set;
	h_set_all(SimulationConfigurationNr,:) = GP_CBF_Delay_HOCBF_h(x_set);
	h_star_set_all(SimulationConfigurationNr,:) = GP_CBF_Delay_HOCBF_h(x_set) - delta_h;
	delta_h_all(SimulationConfigurationNr) = delta_h;
	d_set = ISSF_CBF_disturbance(x_set);
	PredictionError_set = sqrt(sum((d_set - mu_set) .^ 2));
	PredictionError_set_all(SimulationConfigurationNr,:) = PredictionError_set;
	eta_GP_set_all(SimulationConfigurationNr,:) = eta_GP_set;
	eta_Delay_set_all(SimulationConfigurationNr,:) = eta_Delay_set;
end
eta_GP_set_all(eta_GP_set_all < 1e-10) = 1e-10;
eta_Delay_set_all(eta_Delay_set_all < 1e-10) = 1e-10;
% x_ref_set = GP_CBF_Delay_Reference(t_set);
% h_set_dist = GP_CBF_Delay_HOCBF_h(x_set(:,1:numel(t_set)));
% h_star_set = h_set_dist - delta_h;
%%
figure; 
hold on;
plot(x_ref_set(1,:), x_ref_set(2,:), 'c--', 'LineWidth', 1.2);
theta = linspace(0,2*pi,100);
h_bound = [cos(theta); sin(theta)];
plot(h_bound(1,:), h_bound(2,:), 'k--', 'LineWidth', 3);
%plot(x_set_nodist(1,:), x_set_nodist(2,:), 'b-', 'LineWidth', 1.5);
for SimulationConfigurationNr = 1:numel(SimulationConfigurationSet)
	SimulationConfiguration = SimulationConfigurationSet{SimulationConfigurationNr};
	x_set = x_set_all{SimulationConfigurationNr};
	plot(x_set(1,:), x_set(2,:), '-', ...
		'LineWidth', 1.5,'Color',SimulationConfiguration.Color);
	h_star_bound = h_bound * sqrt(1 - delta_h_all(SimulationConfigurationNr));
	plot(h_star_bound(1,:), h_star_bound(2,:), '--', ...
		'LineWidth', 1.2,'Color',SimulationConfiguration.Color);
end

% theta = linspace(0,2*pi,100);
% h_bound = [cos(theta); sin(theta)];
% plot(h_bound(1,:), h_bound(2,:), 'k--', 'LineWidth', 1.2);
% 
% circle_star = sqrt(1 - delta_h) * [cos(theta); sin(theta)];
% plot(circle_star(1,:), circle_star(2,:), 'm--', 'LineWidth', 1.2);

% legend('reference', 'with disturbance', 'Safe set','Conservative zone');
title('Trajectory x Without Delay: With vs. Without Disturbance');
xlabel('q_1'); ylabel('q_2'); axis equal; grid on;
%%
figure;
hold on;
for SimulationConfigurationNr = 1:numel(SimulationConfigurationSet)
	SimulationConfiguration = SimulationConfigurationSet{SimulationConfigurationNr};
	h_set = h_set_all(SimulationConfigurationNr,:);
	plot(t_set, h_set, '-', ...
		'LineWidth', 1.5,'Color',SimulationConfiguration.Color);
	h_star_set = h_star_set_all(SimulationConfigurationNr,:);
	plot(t_set, h_star_set, '--', ...
		'LineWidth', 1.5,'Color',SimulationConfiguration.Color);
end
legend('h(x) Disturbed', 'h0*(x)');
xlabel('Time (s)'); ylabel('h(q)');
title('Safety constraint h(q) evolution'); grid on;

%%

figure;

for SimulationConfigurationNr = 1:numel(SimulationConfigurationSet)
	SimulationConfiguration = SimulationConfigurationSet{SimulationConfigurationNr};
	subplot(2,2,SimulationConfigurationNr);

	PredictionError_set = PredictionError_set_all(SimulationConfigurationNr,:);
	semilogy(t_set, PredictionError_set, '-', ...
		'LineWidth', 1.2,'Color',SimulationConfiguration.Color);
	hold on;
	%
	eta_GP_set = eta_GP_set_all(SimulationConfigurationNr,:);
	semilogy(t_set, eta_GP_set, '-.', ...
		'LineWidth', 1,'Color',SimulationConfiguration.Color);
	%
	eta_Delay_set = eta_Delay_set_all(SimulationConfigurationNr,:);
	semilogy(t_set, eta_Delay_set, ':', ...
		'LineWidth', 1,'Color',SimulationConfiguration.Color);
	%
	eta_set = eta_GP_set + eta_Delay_set;
	semilogy(t_set, eta_set, '--', ...
		'LineWidth', 2,'Color',[SimulationConfiguration.Color,0.5]);
	%
	ylim([1e-2;1e1]);
	legend({'$\| e \|$', '$\bar{e}_{GP}$', '$\bar{e}_{Delay}$', '$\bar{e}$'}, ...
	'Interpreter', 'latex', 'NumColumns',2);
end
% plot(t_set, d_set(1,:), 'r-', 'LineWidth', 1.2);
% hold on;
% plot(t_set, mu_set(1,:), 'b--', 'LineWidth', 1.2);hold on;
% plot(t_set, d_set(2,:), 'g-', 'LineWidth', 1.2);
% plot(t_set, mu_set(2,:), 'm--', 'LineWidth', 1.2);

xlabel('Time (s)'); ylabel('Disturbance');
title('True disturbance vs GP prediction');
grid on;
%%
GP_CBF_Delay_Postprocessing_DisturbanceMap;
