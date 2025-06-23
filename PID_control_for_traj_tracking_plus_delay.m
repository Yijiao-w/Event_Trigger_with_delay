%% Robot PID control for trajectory tracking within circular safe set
clc; clear all; close all; 
rng(0);
%% system setting
dt = 0.01;  % time step
T = 20;     % time
t_set = 0:dt:(T-dt);

%% eta set
eta_samples = 1;
beta = 5;
SigmaN = 0.02;
eta_star_max = beta * SigmaN;

%% Simulation setting
sub_mode = 2;  % change sub_mode for delay
% 1 = No Delay, 2 = With Delay

MonteCarlo_runs = 1;  % 10-times monte carlo simulation

eta_star_list = linspace(1e-5, 1, eta_samples) * eta_star_max;  % (0, eta_star_max]

result_list = cell(eta_samples, MonteCarlo_runs);

h_bar = nan(1, eta_samples);    % 平均值
h_std = nan(1, eta_samples);    %varianz
h_list = nan(1,eta_samples);

h_bar_nodelay = nan(1, eta_samples);

eta_actual_mean = nan(eta_samples, MonteCarlo_runs);

%% eta list loop
for p = 1:eta_samples
	eta_p = eta_star_list(p);
	h_samples = nan(1, MonteCarlo_runs);
	%% simulation for monte carlo
	for mc = 1:MonteCarlo_runs
		[min_h0_star, x_set_nodist, x_set_dist, delta_h,...
			h_set_nodist, h_set_dist, u_set, d_set, mu_d_set,ET_count,GP,eta_mean_used]...
			= GP_CBF_Delay_TestFunc(eta_p, sub_mode,mc);

		h_samples(mc) = min_h0_star;
		fprintf('Computed: eta %d, run %d\n', p, mc);
	end

	h_bar(p) = mean(h_samples,'omitnan');
	h_std(p) = std(h_samples,'omitnan');
	h_sort = sort(h_bar, 'descend');

	h_bar_smooth = smooth(h_bar, 5);  % 滑动窗口为5
	eta_actual_mean(p, mc) = eta_mean_used;
	eta_actual_mean_bar = mean(eta_actual_mean, 2, 'omitnan');  % size = [eta_samples x 1]
end
if ~exist('results', 'dir')
	mkdir('results');
end
save('results/all_simulation_results.mat', 'result_list');
%% Plotting
x_set_dist = x_set_nodist;
x_ref_set = GP_CBF_Delay_Reference(t_set);
h_set_dist = GP_CBF_Delay_HOCBF_h(x_set_dist(:,1:numel(t_set)));
h_star_set = h_set_dist - delta_h;
%% Monte Carlo Plot for h0_star vs eta_star

figure; hold on;
%plot(eta_star_list, h0_star_theory, 'b-o', 'LineWidth', 2);
%plot(eta_star_list, h_list, 'b-o', 'LineWidth', 1.5);  % 平均值线

plot(eta_star_list, h_bar,'b-o','LineWidth', 1.5);  % 平均值线

%yline(mean(h_bar_sorted), 'r--', 'LineWidth', 1.2);  % 总体平均线

xlabel('\eta^*');
ylabel('h_0^* (mean)');
title('Monte Carlo: Mean of h_0^* under Different \eta^*');
legend('±1 std deviation', 'Mean of h_0^*', 'Global Mean');
grid on;

%%
figure; hold on;
plot(x_ref_set(1,:), x_ref_set(2,:), 'g--', 'LineWidth', 1.2);
%plot(x_set_nodist(1,:), x_set_nodist(2,:), 'b-', 'LineWidth', 1.5);
plot(x_set_dist(1,:), x_set_dist(2,:), 'r-', 'LineWidth', 1.5);

theta = linspace(0,2*pi,100);
circle = [cos(theta); sin(theta)];
plot(circle(1,:), circle(2,:), 'k--', 'LineWidth', 1.2);

circle_star = sqrt(1 - delta_h) * [cos(theta); sin(theta)];
plot(circle_star(1,:), circle_star(2,:), 'm--', 'LineWidth', 1.2);
% valid_h0_star = delta_h(~isnan(delta_h) & delta_h >= 0 & delta_h <= 1);
% if ~isempty(valid_h0_star)
% 	h0_star_mean = mean(valid_h0_star);
% 	r_h0_star_mean = sqrt(1 - h0_star_mean);
% 	h0_star_circle = r_h0_star_mean * [cos(theta); sin(theta)];
% 	plot(h0_star_circle(1,:), h0_star_circle(2,:), 'm--', 'LineWidth', 1.2);
% end

legend('reference', 'with disturbance', 'Safe set','Conservative zone');
title('Trajectory x Without Delay: With vs. Without Disturbance');
xlabel('q_1'); ylabel('q_2'); axis equal; grid on;
%%
figure;
plot(t_set, h_set_dist, 'r', t_set, h_star_set, 'g', 'LineWidth', 1.5);
legend('h(x) Disturbed', 'h0*(x)');
xlabel('Time (s)'); ylabel('h(q)');
title('Safety constraint h(q) evolution'); grid on;

%%
figure;
plot(t_set, d_set(1,:), 'r-', 'LineWidth', 1.2);
hold on;
plot(t_set, mu_d_set(1,:), 'b--', 'LineWidth', 1.2);hold on;
plot(t_set, d_set(2,:), 'g-', 'LineWidth', 1.2);
plot(t_set, mu_d_set(2,:), 'm--', 'LineWidth', 1.2);
legend('d_1 (true)', 'μ_1 (predicted)', 'd_2 (true)', 'μ_2 (predicted)');
xlabel('Time (s)'); ylabel('Disturbance');
title('True disturbance vs GP prediction');
grid on;

