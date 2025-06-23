function [x_set_nodist, x_set_dist, ...
	delta_h, h_set_nodist, h_set_dist, u_set, d_set, ...
	mu_set,ET_count,GP,eta_mean_used] = GP_CBF_Delay_TestFunc( ...
	eta_p, sub_mode,run_id)
%%
t_start = 0;
t_end = 20;     % time
t_step = 0.01;  % time step
t_set = t_start:t_step:(t_end-t_step);
steps = length(t_set);
q_dim = 2;
%% system initialization
q0 = [0; 0];
q0_dot = [1; 0];
x0 = [q0; q0_dot]; % x=[q;q_dot] q=[-1;1]; q_dot=[1;0]

% Euler-Lagrange form: M*q_ddot + C*q_dot + G = u
M = eye(q_dim);
C = zeros(q_dim);
G = zeros(q_dim, 1);
Model_Parameter.M = M;
Model_Parameter.C = C;
Model_Parameter.G = G;
%% PID controller gains
Kp = 1 * eye(q_dim);
Kd = 1.5 * eye(q_dim);
PD_Controller_Parameter.Kp = Kp;
PD_Controller_Parameter.Kd = Kd;
%% 控制参数
alpha_0 = 2.0;
alpha_1 = 2.5;
HOCBF_Parameter.alpha_0 = alpha_0;
HOCBF_Parameter.alpha_1 = alpha_1;

epsilon = 20;
%% GP
x_dim = 2 * q_dim;
y_dim = q_dim;
SigmaF = 1.5;
SigmaL = 0.35;
SigmaN = 0.2;
MaxDataQuantity = 200;

GP = LocalGP_MultiOutput(x_dim, y_dim, MaxDataQuantity, SigmaN, SigmaF, SigmaL);

beta = 0.2;
eta_desired = 1 * beta * SigmaN;
%% delay
Ld = 0.1;
F = 0.1;
Max_GP_Delay_Step = 10;
	Max_GP_Delay_Step = max(1,Max_GP_Delay_Step); % if Max_GP_Delay_Step = 1, then almost no delay
delta = Max_GP_Delay_Step * t_step;
eta_delay = Ld * sqrt(2 * F * delta);

delta_h = epsilon^2 * (eta_desired + 1 * eta_delay)^2 / (alpha_0 * alpha_1);
%% Data storage
x_set_nodist = nan(q_dim*2, steps);
x_set_nodist(:,1) = x0;
x_set_dist = nan(q_dim*2, steps);
h_set_nodist = nan(1, steps);
h_set_dist = nan(1, steps);
h0_star_set = nan(1, steps);
h1_star_set = nan(1, steps);
r_h0_star_set = nan(1, steps); % h0*(x) 对应的半径
u_set = nan(q_dim, steps);

mu_set = nan(y_dim, steps);

d_set = nan(y_dim, steps);       % 实际扰动值
eta_used_record = zeros(y_dim, steps);

ET_count = 0;
delay_steps = round(delta / t_step);
trigger_queue = false(1, delay_steps);
ET_index_list = [];
%% simulation for 2 mode
for mode = 1:1  % 1 = no disturbance, 2 = with disturbance
	step_counter = 0;
	mu_d_pre = zeros(y_dim, 1);
	sigma_d_pre = zeros(y_dim, 1);
	eta_pre = zeros(y_dim, 1);

	GP_CountDown = 0;
	

	mu = zeros(y_dim,1);
	for t_Nr = 1:steps
		x = x_set_nodist(:,t_Nr);
		t = t_set(t_Nr);
		u_nom = GP_CBF_Delay_NominalController( ...
			t,x,Model_Parameter,PD_Controller_Parameter);

		if mode == 1
			% Solve QP
			% 			mu = ISSF_CBF_disturbance(x);
			if GP_CountDown == 0 % GP is free now
				% Output previous prediction
				if t_Nr > 1
					[mu,GP] = GP_CBF_Delay_Previous_Prediction_Output( ...
						x_GP,GP,do_GP_update);
				end
				% Start new computation
				x_GP = x;
				do_GP_update = true; % Determine whether Gp should be updated

				GP_CountDown = Max_GP_Delay_Step - 1;
			else
				GP_CountDown = GP_CountDown - 1;
			end
			
			d_set(:,t_Nr) = ISSF_CBF_disturbance(x);
			mu_set(:,t_Nr) = mu;

			h = GP_CBF_Delay_HOCBF_h(x);
			h_star = h - delta_h;

			SolverType = 'analytic';
			u = GP_CBF_Delay_ISSf_HOCBF_QP_Solver(SolverType, ...
				x,u_nom,mu,h_star,epsilon,Model_Parameter,HOCBF_Parameter);
			% Integrate system
			[~, x_temp_set] = ode45(@(t,x)GP_CBF_Delay_Dynamics(t,x,u,Model_Parameter), ...
				[t, t + t_step], x);

			x_set_nodist(:,t_Nr + 1) = x_temp_set(end, :)';
			h_set_nodist(t_Nr) = GP_CBF_Delay_HOCBF_h(x);

		end

		% Disturbance (if mode == 2)
% 		if mode == 2
% 			step_counter = step_counter + 1;
% 			d = ISSF_CBF_disturbance(x);
% 			d_set(:, t_Nr) = d;  % 存储实际扰动
% 
% 			delete_point_step = 200;  % 在第200步时删点
% 			num_delete = 20;          % 删掉最旧的20个点
% 			deleted_flag = false;     % 标志位，防止多次执行
% 
% 			%% delay setting
% 			if step_counter < 10  % check for GP update condition
% 
% 				mu_d = mu_d_pre;
% 				sigma_d = sigma_d_pre;
% 				eta = eta_pre;
% 			else
% 				[~, sigma_d] = GP.predict(x);
% 				eta = beta * sigma_d;
% 				%eta = min(beta*sigma_d,eta_p);
% 
% 				psi_r1_star = h1 - (epsilon^2 * norm(eta)^2) / alpha_1;
% 				alpha_val = alpha_1 * psi_r1_star;
% 				eta_bound = norm(eta) + eta_delay;  % 延迟时的 eta 上界
% 				e_mu = eta_p + eta_delay;  % GP 的均值估计误差上界
% 
% 				%% event trigger setting
% 				% ξ_ET condition
% 				rho = 0.7;
% 				L_delta_alpha = 0.1;
% 
% 				xi_ET = - rho * alpha_val - epsilon * e_mu^2 ...
% 					+ 2 * rho * L_delta_alpha * delta ...
% 					+ (1 - rho) * epsilon * eta_bound^2;
% 
% 				trigger_flag = (psi_r1_star < 0) || (xi_ET >= 0);
% 
% 				% 将当前 trigger_flag 存入循环队列
% 				idx_now = mod(step_counter - 1, delay_steps) + 1;
% 				trigger_queue(idx_now) = trigger_flag;
% 
% 				if mod(step_counter, delay_steps) == 0
% 					idx_check = mod(step_counter - delay_steps, delay_steps) + 1;
% 					% 组合 trigger 判据：只要进入保守区 或 原始判据成立，就触发
% 
% 					fprintf('[step %d] xi_ET = %.5f, alpha_val = %.5f, e_mu^2 = %.5f, eta_bound^2 = %.5f\n', ...
% 						t_Nr, xi_ET, alpha_val, e_mu^2, eta_bound^2);
% 
% 					% [mu_d, sigma_d, ~, ~, ~, ~] = GP.predict(x);
% 
% 					%if xi_ET >= 0
% 					if trigger_queue(idx_check)
% 						ET_count = ET_count + 1;
% 						ET_index_list(ET_count) = t_Nr;
% 						if GP.check_Saturation
% 							GP.downdateParam(1);
% 							t_step = t_step + ComputationPowerCoefficient_delete * DataQuantity^2;
% 						end
% 						y = d + SigmaN * randn(y_dim,1);
% 						GP.addPoint(x, y);
% 
% 						% === 这里插入删除旧点的逻辑（仅执行一次）===
% 						if step_counter == delete_point_step && ~deleted_flag
% 							for k = 1:num_delete
% 								if GP.DataQuantity > 1
% 									GP.downdateParam(1);  % 始终删除最旧的点
% 								end
% 							end
% 							deleted_flag = true;
% 							fprintf('[step %d] Deleted first %d GP points.\n', step_counter, num_delete);
% 						end
% 					end
% 					%GP
% 					[mu_d, sigma_d] = GP.predict(x);
% 					eta = min(beta * sigma_d, eta_p);  % 控制 GP η 不超过 η*
% 				else
% 					% 非 trigger 步，使用上一步 GP
% 					mu_d = mu_d_pre;
% 					sigma_d = sigma_d_pre;
% 					eta = eta_pre;
% 				end
% 			end
% 
% 			eta_used_record(:, t_Nr) = eta;  % i 是时间步
% 
% 
% 			mu_d_pre = mu_d;
% 			sigma_d_pre = sigma_d;
% 			eta_pre = eta;
% 
% 
% 
% 			%  dh1_star/dx part
% 			% 从 GP 中提取训练数据
% 			X_train = GP.X(:, 1:GP.DataQuantity)';
% 			N = size(X_train, 1);
% 			w = x_dim;  % input dimension
% 			sigma_f = SigmaF;
% 			Lambda = diag(SigmaL.^2);
% 
% 			% 核矩阵
% 			K = GP.K(1:N, 1:N) + SigmaN^2 * eye(N);
% 			invK = inv(K);
% 
% 			% 核函数计算
% 			K_x = zeros(N, 1);
% 			grad_K_x = zeros(N, w);
% 
% 			for j = 1:N
% 				xi = X_train(j, :)';
% 				diff = x - xi;
% 				r2 = diff' * (Lambda \ diff);
% 
% 				K_x(j) = sigma_f^2 * exp(-0.5 * r2);
% 
% 				% dKx/dx
% 				grad_K_x(j, :) = K_x(j) * (-Lambda \ diff)';
% 			end
% 
% 			% dsigma^2/dx
% 			grad_sigma2 = -2 * (grad_K_x' * invK * K_x);
% 
% 			%% h0_star & h1_star
% 			dh1_star_dx = dh1_dx - (epsilon^2 * beta^2/ alpha_1) * grad_sigma2;
% 			Lf_h1_star = dh1_star_dx * f(x);
% 			Lg_h1_star = dh1_star_dx * g(x);
% 
% 
% 			mu_d_set(:,t_Nr) = mu_d;
% 			a1_x = Lf_h1_star + Lg_h1_star * mu_d;
% 			b_x = dh1_star_dx;
% 
% 
% 			if sub_mode == 1
% 				h0_star = h0 - epsilon^2 *  eta_p^2 / (alpha_0 * alpha_1);
% 				h1_star = dh_dq(q) * q_dot + alpha_0 * h0_star;
% 
% 			elseif sub_mode == 2
% 				h0_star = h0 - epsilon^2 * (eta_p + eta_delay)^2 / (alpha_0 * alpha_1); % for delay
% 				h1_star = dh_dq(q) * q_dot + alpha_0 * h0_star;
% 
% 			end
% 			%% QP problem
% 			Q = 2*eye(n);
% 			f_qp = -2 * u_nom;
% 
% 			%  b_qp for 2 sub mode
% 			if sub_mode == 1  % No Delay
% 				b_qp = a1_x - norm(b_x)^2 / (4 * epsilon^2)  + alpha_1 * h1_star;
% 			elseif sub_mode == 2  % With Delay
% 				b_qp = a1_x - norm(b_x)^2 / (4 * epsilon^2) + alpha_1 * h1_star - epsilon^2 * eta_delay^2;
% 			end
% 
% 			A_qp = -Lg_h1_star;
% 
% 			u = quadprog(Q, f_qp, A_qp, b_qp, [], [], [], [], [], optimoptions('quadprog','Display','off'));
% 
% 			u(u >= 80) = 80;
% 			u(u <= -80) = -80;
% 			u_set(:,t_Nr) = u;
% 
% 
% 			if h0_star >= 0 && h0_star <= 1
% 				r_h0_star_set(t_Nr) = sqrt(1 - h0_star);  % h0* radius
% 			else
% 				r_h0_star_set(t_Nr) = NaN;  % 记录为 NaN，表示不在有效范围内
% 			end
% 
% 			%% system setting
% 			odefun = @(t, x) f(x) + g(x) * (u_nom + d);
% 			[~, x_temp_set] = ode45(@(t,x)GP_CBF_Delay_Dynamics(t,x,u), ...
% 				[t, t + t_step], x);
% 			x = x_temp_set(end, :)';
% 
% 			h_set_dist(t_Nr) = h(q);
% 			x_set_dist(:,t_Nr) = x;
% 			h0_star_set(t_Nr) = h0_star;
% 			h1_star_set(t_Nr) = h1_star;
% 
% 			x_ref_set = arrayfun(q_ref_set, t_set, 'UniformOutput', false); % reference x_ref(t) in cell array
% 			x_ref_set = cell2mat(x_ref_set);
% 		end
	end
	fprintf('t = %6.4f\n',t_set(t_Nr));

end
% end

eta_mean_used = mean(vecnorm(eta_used_record(:, 1:steps), 2, 1), 'omitnan');

file_name = sprintf('results/simulation_eta_%.4f_run%d.mat', eta_p, run_id);
if ~exist('results', 'dir')
	mkdir('results');
end
save(file_name, 'h0_star_set', 'eta_p');

end