%% Robot PID control for trajectory tracking within circular safe set
clear all; close all; clc;
rng(0);
%% system setting
dt = 0.01;  % time step
T = 20;     % time
t = 0:dt:T-dt;
steps = length(t);
n = 2;

%for delay
Ld = 3;
F = 1.5;
delta = 10 * dt;
eta_delay = Ld * sqrt(2 * F * delta);

% system initialization
q = [0; 0];
q_dot = [1; 0];
x = [q; q_dot]; % x=[q;q_dot] q=[-1;1]; q_dot=[1;0]

% Euler-Lagrange form: M*q_ddot + C*q_dot + G = u
M = eye(n);
C = zeros(n);
G = zeros(n, 1);

%% PID controller gains
Kp = 10 * eye(n);
Kd = 10 * eye(n);

%% 控制参数
alpha_0 = 10;
alpha_1 = 13;
beta = 5;
gamma = 0;
epsilon = 2;
%% Safe set: h(q) = 1 - q1^2 - q2^2 > 0
h = @(q) 1 - q(1)^2 - q(2)^2;
dh_dq = @(q) [-2*q(1), -2*q(2)];

% x_dot=f(x)+g(x)(u+d(x))
f = @(x) [x(3:4); -M \ (C * x(3:4) + G) ];
g = @(x) [zeros(n); M \ eye(n)];

%% x reference
a = 1;   % x 方向半轴
b = 0.6; % y 方向半轴
c = 0.3;

x_ref_set = @(t) [a * cos(c*t)-0.3; b * sin(c*t)];
dx_ref_set = @(t) [-c*a*sin(c*t); c*b*cos(c*t)];
ddx_ref_set = @(t) [-c^2*a*cos(c*t); -c^2*b*sin(c*t)];

%% GP
x_dim = 4;
y_dim = 2;
SigmaF = 1;
SigmaL = 0.6;
SigmaN = 0.02;
MaxDataQuantity = 200;

GP = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, SigmaN, SigmaF, SigmaL);
for j = 1:18
    x_sample = randn(4,1);  % random x in 4-dim
    d_sample = ISSF_CBF_disturbance(x_sample);
    y_sample = d_sample + SigmaN * randn(y_dim,1);
    GP.addPoint(x_sample, d_sample);
end
%% eta set
eta_samples = 1;
eta_star_max = beta * SigmaN;

%% Computation Time Coefficient
ComputationPowerCoefficient = 0.25e-6; % 5e-7 % 1e-8
ComputationPowerCoefficient_update = 0/2 * 0.25e-6;
ComputationPowerCoefficient_delete = 0/2 * 0.25e-6;
ComputationPowerBias = 0.000;

DataQuantity = GP.DataQuantity;
t_step = ComputationPowerCoefficient * DataQuantity ^ 2 + ...
    ComputationPowerBias;
t_step = 1;
%% 初始化 step_counter
step_counter = 0;  % 计数器从0开始

%% Data storage
x_set_nodist = nan(n*2, steps);
x_set_dist = nan(n*2, steps);
h_set_nodist = nan(1, steps);
h_set_dist = nan(1, steps);
h0_star_set = nan(1, steps);
h1_star_set = nan(1, steps);
r_h0_star_set = nan(1, steps); % h0*(x) 对应的半径
u_set = nan(n, steps);

mu_d_pre = zeros(y_dim, 1);  % 初始的 mu_d（可以是零向量或其他合适的值）
sigma_d_pre = zeros(y_dim, 1);  % 初始的 sigma_d（可以是零向量或其他合适的值）
eta_pre = zeros(y_dim, 1);  % 初始的 eta（可以是零向量或其他合适的值）

mu_d_set = nan(y_dim, steps);

d_set = nan(y_dim, steps);       % 实际扰动值

%% Simulation setting
ET_count = 0;
delay_steps = round(delta / dt);
trigger_queue = false(1, delay_steps);



sub_mode = 2;  % change sub_mode for delay
% 1 = No Delay, 2 = With Delay

MonteCarlo_runs = 1;  % 10-times monte carlo simulation

eta_star_list = linspace(1e-5, 1, eta_samples) * eta_star_max;  % (0, eta_star_max]
%eta_star_list = linspace(0.001249, 0.001255, eta_samples);  % (0.001247, 0.001255]


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
    h_samples_nodelay = nan(1, MonteCarlo_runs);


    %% simulation for monte carlo
    for mc = 1:MonteCarlo_runs
      

        % simulation for 2 mode

       file_name = sprintf('results/simulation_eta_%.4f_run%d.mat', eta_p, mc);
       file_name = sprintf('results/simulation_eta_%.4f_submode%d_run%d.mat', eta_p, sub_mode, mc);

        % if exist(file_name, 'file')
        %     % 加载已有结果
        %     data = load(file_name);
        %     result_list{p, mc}.eta_p = eta_p;
        %     result_list{p, mc}.h0_star_set = data.h0_star_set;
        %     h_samples(mc) = mean(data.h0_star_set);
        %     fprintf('Skipped: eta %d, run %d (already done)\n', p, mc);
        % else

[min_h0_star, x_ref_traj, x_set_nodist, x_set_dist, h0_star_set,...
    h_set_nodist, h_set_dist, u_set, d_set, mu_d_set,ET_count,GP,eta_mean_used]...
= simulation_result(eta_p, sub_mode,GP,ET_count,mc);  

h_samples(mc) = min_h0_star;

            fprintf('Computed: eta %d, run %d\n', p, mc);
    end
   
h_bar(p) = mean(h_samples,'omitnan');
h_std(p) = std(h_samples,'omitnan');
h_sort = sort(h_bar, 'descend');

%     valid_samples = h_samples(~isnan(h_samples));
% if ~isempty(valid_samples)
%     h_bar(p) = sum(valid_samples) / length(valid_samples);
% else
%     h_bar(p) = NaN;  % 或设置为0/默认值
% end
    h_bar_smooth = smooth(h_bar, 5);  % 滑动窗口为5
    eta_actual_mean(p, mc) = eta_mean_used;
    eta_actual_mean_bar = mean(eta_actual_mean, 2, 'omitnan');  % size = [eta_samples x 1]
    %h0_star_theory = -epsilon^2 * (eta_star_list + eta_delay).^2 / (alpha_0 * alpha_1);
end
if ~exist('results', 'dir')
    mkdir('results');
end
save('results/all_simulation_results.mat', 'result_list');
%% Plotting
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
plot(x_ref_traj(1,:), x_ref_traj(2,:), 'g--', 'LineWidth', 1.2);
%plot(x_set_nodist(1,:), x_set_nodist(2,:), 'b-', 'LineWidth', 1.5);
plot(x_set_dist(1,:), x_set_dist(2,:), 'r-', 'LineWidth', 1.5);

 theta = linspace(0,2*pi,100);
circle = [cos(theta); sin(theta)];
plot(circle(1,:), circle(2,:), 'k--', 'LineWidth', 1.2);

valid_h0_star = h0_star_set(~isnan(h0_star_set) & h0_star_set >= 0 & h0_star_set <= 1);
if ~isempty(valid_h0_star)
    h0_star_mean = mean(valid_h0_star);
    r_h0_star_mean = sqrt(1 - h0_star_mean);
    h0_star_circle = r_h0_star_mean * [cos(theta); sin(theta)];
    plot(h0_star_circle(1,:), h0_star_circle(2,:), 'm--', 'LineWidth', 1.2);
end

legend('reference', 'with disturbance', 'Safe set','Conservative zone');
title('Trajectory x Without Delay: With vs. Without Disturbance');
xlabel('q_1'); ylabel('q_2'); axis equal; grid on;
%%
figure;
plot(t, h_set_dist, 'r', t, h0_star_set, 'g', 'LineWidth', 1.5);
legend('h(x) Disturbed', 'h0*(x)');
xlabel('Time (s)'); ylabel('h(q)');
title('Safety constraint h(q) evolution'); grid on;

%%
figure;
plot(t, d_set(1,:), 'r-', 'LineWidth', 1.2);
hold on;
plot(t, mu_d_set(1,:), 'b--', 'LineWidth', 1.2);hold on;
plot(t, d_set(2,:), 'g-', 'LineWidth', 1.2);
plot(t, mu_d_set(2,:), 'm--', 'LineWidth', 1.2);
legend('d_1 (true)', 'μ_1 (predicted)', 'd_2 (true)', 'μ_2 (predicted)');
xlabel('Time (s)'); ylabel('Disturbance');
title('True disturbance vs GP prediction');
grid on;



% figure;
% plot(t, u_set, 'r-', 'LineWidth', 1.2);

% figure; hold on;
% plot(eta_star_list, h_bar, 'b-o', 'LineWidth', 1.5);  % 原有: 有 delay
% plot(eta_star_list, h_bar_nodelay, 'r--^', 'LineWidth', 1.5);  % 新增: 无 delay
% xlabel('\eta^*');
% ylabel('h_0^* (mean)');
% legend('With delay', 'No delay');
% title('Monte Carlo: Mean of h_0^* under Different \eta^*');
% grid on;
% 
