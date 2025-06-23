function [min_h0_star, x_ref_traj, x_set_nodist, x_set_dist, ...
     h0_star_set, h_set_nodist, h_set_dist, u_set, d_set, ...
     mu_d_set,ET_count,GP,eta_mean_used] = simulation_result(eta_p, sub_mode,GP,ET_count,run_id);

dt = 0.01;  % time step
T = 20;     % time
t = 0:dt:T-dt;
steps = length(t);
n = 2;

%for delay
Ld = 3;
F = 2.5;
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
Kp = 100 * eye(n);
Kd = 150 * eye(n);

%% 控制参数
alpha_0 = 20;
alpha_1 = 25;
beta = 4.4;
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

% x_ref_traj = arrayfun(x_ref_set, t, 'UniformOutput', false); % reference x_ref(t) in cell array
% x_ref_traj = cell2mat(x_ref_traj);

%% GP
x_dim = 4;
y_dim = 2;
SigmaF = 1.5;
SigmaL = 0.35;
SigmaN = 0.015;
MaxDataQuantity = 200;

GP = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, SigmaN, SigmaF, SigmaL);
% for j = 1:20
%     x_sample = randn(4,1);  % random x in 4-dim
%     d_sample = ISSF_CBF_disturbance(x_sample);
%     y_sample = d_sample + SigmaN * randn(y_dim,1);
%     GP.addPoint(x_sample, d_sample);
% end
for j = 1:0
    t_j = (j-1) * T/20;
    x_ref_j = [x_ref_set(t_j); dx_ref_set(t_j)];
    d_sample = ISSF_CBF_disturbance(x_ref_j);
    y_sample = d_sample + SigmaN * randn(y_dim,1);
    GP.addPoint(x_ref_j, y_sample);
end
%% eta
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
eta_used_record = zeros(y_dim, steps);

ET_count = 0;
delay_steps = round(delta / dt);
trigger_queue = false(1, delay_steps);
ET_index_list = [];
%% simulation for 2 mode
for mode = 1:2  % 1 = no disturbance, 2 = with disturbance

    % 每个 mode 都从相同初始条件开始
    q = [0; 0];
    q_dot = [1; 0];
    x = [q; q_dot];
    step_counter = 0;
    mu_d_pre = zeros(y_dim, 1);
    sigma_d_pre = zeros(y_dim, 1);
    eta_pre = zeros(y_dim, 1);
    for i = 1:steps
        q = x(1:2);
        q_dot = x(3:4);
        q_ref = x_ref_set(t(i));
        dq_ref = dx_ref_set(t(i));
        ddq_ref = ddx_ref_set(t(i));

        e = q - q_ref;
        de = q_dot - dq_ref;
        PID = ddq_ref - Kp * e - Kd * de;

        % PID Control law (computed torque style)
        u_nom = M * PID + C*q_dot + G;

        % HOCBF derivation
        h0 = h(q);
        dh0_star_dx = dh_dq(q) * q_dot;
        Lf_h0 = dh_dq(q) * q_dot;
        h1 = Lf_h0 + alpha_0 * h0;
        dh1_dx = [zeros(1,n), dh_dq(q)];
        Lf_h1 = dh1_dx * f(x);
        Lg_h1 = dh1_dx * g(x);


        if mode == 1

            %QP problem
            Q = 2*eye(n);
            f_qp = -2 * u_nom;
            b_qp = Lf_h1 + alpha_1 * h1;
            A_qp = -Lg_h1;

            % Solve QP
            u = quadprog(Q, f_qp, A_qp, b_qp, [], [], [], [], [], optimoptions('quadprog','Display','off'));

            % Integrate system
            odefun = @(t, x) f(x) + g(x) * u;
            [t_ode, x_ode] = ode45(odefun, [t(i), t(i) + dt], x);
            x = x_ode(end, :)';

            x_set_nodist(:,i) = x;
            h_set_nodist(i) = h(q);

        end

        % Disturbance (if mode == 2)
        if mode == 2
            step_counter = step_counter + 1;
            d = ISSF_CBF_disturbance(x);
            d_set(:, i) = d;  % 存储实际扰动

            delete_point_step = 200;  % 在第200步时删点
            num_delete = 20;          % 删掉最旧的20个点
            deleted_flag = false;     % 标志位，防止多次执行

            %% delay setting
            if step_counter < 10  % check for GP update condition

                mu_d = mu_d_pre;
                sigma_d = sigma_d_pre;
                eta = eta_pre;                
             else
                      [mu_d, sigma_d, ~, ~, ~, ~] = GP.predict(x);
                 eta = beta * sigma_d;
                 %eta = min(beta*sigma_d,eta_p);

                psi_r1_star = h1 - (epsilon^2 * norm(eta)^2) / alpha_1;
                alpha_val = alpha_1 * psi_r1_star;
                eta_bound = norm(eta) + eta_delay;  % 延迟时的 eta 上界
                e_mu = eta_p + eta_delay;  % GP 的均值估计误差上界

                %% event trigger setting
                % ξ_ET condition
                rho = 0.7;
                L_delta_alpha = 0.1;

                xi_ET = - rho * alpha_val - epsilon * e_mu^2 ...
                    + 2 * rho * L_delta_alpha * delta ...
                    + (1 - rho) * epsilon * eta_bound^2;

                trigger_flag = (psi_r1_star < 0) || (xi_ET >= 0);

                % 将当前 trigger_flag 存入循环队列
                idx_now = mod(step_counter - 1, delay_steps) + 1;
                trigger_queue(idx_now) = trigger_flag;

                if mod(step_counter, delay_steps) == 0
                    idx_check = mod(step_counter - delay_steps, delay_steps) + 1;
                    % 组合 trigger 判据：只要进入保守区 或 原始判据成立，就触发

                    fprintf('[step %d] xi_ET = %.5f, alpha_val = %.5f, e_mu^2 = %.5f, eta_bound^2 = %.5f\n', ...
                        i, xi_ET, alpha_val, e_mu^2, eta_bound^2);

                    % [mu_d, sigma_d, ~, ~, ~, ~] = GP.predict(x);

                    %if xi_ET >= 0
                    if trigger_queue(idx_check)
                        ET_count = ET_count + 1;
                        ET_index_list(ET_count) = i;
                        if GP.check_Saturation
                            GP.downdateParam(1);
                            t_step = t_step + ComputationPowerCoefficient_delete * DataQuantity^2;
                        end
                        y = d + SigmaN * randn(y_dim,1);
                        GP.addPoint(x, y);

                        % === 这里插入删除旧点的逻辑（仅执行一次）===
                        if step_counter == delete_point_step && ~deleted_flag
                            for k = 1:num_delete
                                if GP.DataQuantity > 1
                                    GP.downdateParam(1);  % 始终删除最旧的点
                                end
        end
        deleted_flag = true;
        fprintf('[step %d] Deleted first %d GP points.\n', step_counter, num_delete);
    end
                    end
                    %GP
                    [mu_d, sigma_d, ~, ~, ~, ~] = GP.predict(x);
                    eta = min(beta * sigma_d, eta_p);  % 控制 GP η 不超过 η*
                else
                    % 非 trigger 步，使用上一步 GP
                    mu_d = mu_d_pre;
                    sigma_d = sigma_d_pre;
                    eta = eta_pre;
                end
            end

                eta_used_record(:, i) = eta;  % i 是时间步
                

                mu_d_pre = mu_d;
                sigma_d_pre = sigma_d;
                eta_pre = eta;

                

                %  dh1_star/dx part
                % 从 GP 中提取训练数据
                X_train = GP.X(:, 1:GP.DataQuantity)';
                N = size(X_train, 1);
                w = x_dim;  % input dimension
                sigma_f = SigmaF;
                Lambda = diag(SigmaL.^2);

                % 核矩阵
                K = GP.K(1:N, 1:N) + SigmaN^2 * eye(N);
                invK = inv(K);

                % 核函数计算
                K_x = zeros(N, 1);
                grad_K_x = zeros(N, w);

                for j = 1:N
                    xi = X_train(j, :)';
                    diff = x - xi;
                    r2 = diff' * (Lambda \ diff);

                    K_x(j) = sigma_f^2 * exp(-0.5 * r2);

                    % dKx/dx
                    grad_K_x(j, :) = K_x(j) * (-Lambda \ diff)';
                end

                % dsigma^2/dx
                grad_sigma2 = -2 * (grad_K_x' * invK * K_x);

                %% h0_star & h1_star
                dh1_star_dx = dh1_dx - (epsilon^2 * beta^2/ alpha_1) * grad_sigma2;
                Lf_h1_star = dh1_star_dx * f(x);
                Lg_h1_star = dh1_star_dx * g(x);


                mu_d_set(:,i) = mu_d;
                a1_x = Lf_h1_star + Lg_h1_star * mu_d;
                b_x = dh1_star_dx;
        

                if sub_mode == 1
                    h0_star = h0 - epsilon^2 *  eta_p^2 / (alpha_0 * alpha_1);
                    h1_star = dh_dq(q) * q_dot + alpha_0 * h0_star;

                elseif sub_mode == 2
                    h0_star = h0 - epsilon^2 * (eta_p + eta_delay)^2 / (alpha_0 * alpha_1); % for delay
                    h1_star = dh_dq(q) * q_dot + alpha_0 * h0_star;

                end
                %% QP problem
                Q = 2*eye(n);
                f_qp = -2 * u_nom;

                %  b_qp for 2 sub mode
                if sub_mode == 1  % No Delay
                    b_qp = a1_x - norm(b_x)^2 / (4 * epsilon^2)  + alpha_1 * h1_star;
                elseif sub_mode == 2  % With Delay
                    b_qp = a1_x - norm(b_x)^2 / (4 * epsilon^2) + alpha_1 * h1_star - epsilon^2 * eta_delay^2;
                end

                A_qp = -Lg_h1_star;

                u = quadprog(Q, f_qp, A_qp, b_qp, [], [], [], [], [], optimoptions('quadprog','Display','off'));

                u(u >= 80) = 80;
                u(u <= -80) = -80;
                u_set(:,i) = u;


                if h0_star >= 0 && h0_star <= 1
                    r_h0_star_set(i) = sqrt(1 - h0_star);  % h0* radius
                else
                    r_h0_star_set(i) = NaN;  % 记录为 NaN，表示不在有效范围内
                end

                %% system setting
                odefun = @(t, x) f(x) + g(x) * (u + d);
                [t_ode, x_ode] = ode45(odefun, [t(i), t(i) + dt], x);
                x = x_ode(end, :)';

                h_set_dist(i) = h(q);
                x_set_dist(:,i) = x;
                h0_star_set(i) = h0_star;
                h1_star_set(i) = h1_star;

                x_ref_traj = arrayfun(x_ref_set, t, 'UniformOutput', false); % reference x_ref(t) in cell array
                x_ref_traj = cell2mat(x_ref_traj);
            end
        end
        fprintf('t = %6.4f\n',t(i));

    end
% end
min_h0_star = mean(h0_star_set,'omitnan');
%min_h0_star = min(h0_star_set(end-500:end), [], 'omitnan');
% valid_h0_star = h0_star_set(h0_star_set > 0);  % 仅保留正的值
% min_h0_star = mean(valid_h0_star, 'omitnan');

% valid_h0_star = h0_star_set(h0_star_set > 0 & h0_star_set < 1);
% min_h0_star = mean(valid_h0_star, 'omitnan');  % 或用 median()
 %min_h0_star = min(h0_star_set, [], 'omitnan');

 % 剔除无效值
% valid_h0 = h0_star_set(h0_star_set > 0 & h0_star_set < 1);
% 
% if numel(valid_h0) >= 500
%     % 末段均值（稳定）
%     min_h0_star = mean(valid_h0(end-499:end), 'omitnan');
% elseif numel(valid_h0) >= 10
%     % 中位数作为次优替代
%     min_h0_star = median(valid_h0, 'omitnan');
% else
%     % 极端情况：全局均值
%     min_h0_star = mean(valid_h0, 'omitnan');
% end


eta_mean_used = mean(vecnorm(eta_used_record(:, 1:steps), 2, 1), 'omitnan');

file_name = sprintf('results/simulation_eta_%.4f_run%d.mat', eta_p, run_id);  
if ~exist('results', 'dir')
    mkdir('results');
end
save(file_name, 'h0_star_set', 'eta_p');

end