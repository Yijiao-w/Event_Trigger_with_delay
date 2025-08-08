% function GP_CBF_Delay_Manipulator2D_TestFunc(do_Consider_GPError,do_Consider_Delay, ...
% 	SaveFolderName,SaveFileName)
%% GP_CBF_Delay_Manipulator2D_Test
do_Consider_GPError = true;
do_Consider_Delay = true;
%%
t_start = 0;
t_end = 20;     % time
t_step = 0.01;  % time step
t_set = t_start:t_step:(t_end-t_step);
q_dim = 2;
%% system initialization
q0 = [0; 0];
q0_dot = [0; 0];
x0 = [q0; q0_dot]; % x=[q;q_dot] q=[-1;1]; q_dot=[1;0]
%% Euler-Lagrange form: M*q_ddot + C*q_dot + G = u
% M = 1 * eye(q_dim);
% C = zeros(q_dim);
% G = zeros(q_dim, 1);
L1 = 0.75;
L2 = 0.75;
m1 = 1;
m2 = 1;
% Model_Parameter = GP_CBF_Delay_Manipulator2D_get_Model_Parameter( ...
% 	x0,L1,L2,m1,m2);
% Model_Parameter.M = M;
% Model_Parameter.C = C;
% Model_Parameter.G = G;
%% PID controller gains
Kp = 20 * eye(q_dim);
Kd = 15 * eye(q_dim);
PD_Controller_Parameter.Kp = Kp;
PD_Controller_Parameter.Kd = Kd;
%% 控制参数
alpha_0 = 2.0;
alpha_1 = 2.5;
HOCBF_Parameter.alpha_0 = alpha_0;
HOCBF_Parameter.alpha_1 = alpha_1;

epsilon = 2;
%% GP
x_dim = 2 * q_dim;
y_dim = q_dim;
SigmaF = 1.5;
SigmaL = 0.35;
SigmaN = 0.2;
MaxDataQuantity = 200;

GP = LocalGP_MultiOutput(x_dim, y_dim, ...
	MaxDataQuantity, SigmaN, SigmaF, SigmaL);

beta = 1.2;
eta_desired = beta * sqrt(y_dim) * SigmaN;
%% delay
Ld = 1.5;
F = 0.1;
Max_GP_Delay_Step = 20;
Max_GP_Delay_Step = max(1,Max_GP_Delay_Step); % if Max_GP_Delay_Step = 1, then almost no delay
delta = Max_GP_Delay_Step * t_step;
eta_delay = Ld * sqrt(2 * F * delta);
%% delta_h
if do_Consider_GPError
	GPError_Coefficient = 1;
else
	GPError_Coefficient = 0;
end
if do_Consider_Delay
	Delay_Coefficient = 1;
else
	Delay_Coefficient = 0;
end
delta_h = epsilon^2 * (GPError_Coefficient * eta_desired + ...
	Delay_Coefficient * eta_delay)^2 / (alpha_0 * alpha_1);
%% Data storage
x_set = nan(q_dim*2, numel(t_set));
x_set(:,1) = x0;
x_GP = x0;
u_set = nan(q_dim, numel(t_set));
mu_set = nan(y_dim, numel(t_set));
mu = zeros(y_dim,1);
eta_GP_set = nan(1, numel(t_set));
eta_Delay_set = nan(1, numel(t_set));
sigma = sqrt(y_dim) * SigmaF;
GP_CountDown = 0;
PrintCounter = 0;
%% simulation
fprintf([SaveFileName,': \t']);
for t_Nr = 1:numel(t_set)
	x = x_set(:,t_Nr);
	t = t_set(t_Nr);
	Model_Parameter = GP_CBF_Delay_Manipulator2D_get_Model_Parameter( ...
		x,L1,L2,m1,m2);
	u_nom = GP_CBF_Delay_NominalController( ...
		t,x,Model_Parameter,PD_Controller_Parameter);

	% Solve QP
	if GP_CountDown == 0 % GP is free now
		% Output previous prediction
		if t_Nr > 1
			[mu,sigma,GP] = GP_CBF_Delay_Previous_Prediction_Output( ...
				x_GP,GP,do_GP_update);
		end
		% Start new computation
		x_GP = x;
		do_GP_update = true; % Determine whether GP should be updated
		GP_CountDown = Max_GP_Delay_Step - 1;
	else
		GP_CountDown = GP_CountDown - 1;
	end
	mu_set(:,t_Nr) = mu;
	eta_GP_set(t_Nr) = beta * norm(sigma);
	eta_Delay_set(t_Nr) = Ld * sqrt(norm(x - x_GP));

	h = GP_CBF_Delay_HOCBF_h(x);
	h_star = h - delta_h;
	SolverType = 'analytic';
	u = GP_CBF_Delay_ISSf_HOCBF_QP_Solver(SolverType, ...
		x,u_nom,mu,h_star,epsilon,Model_Parameter,HOCBF_Parameter);
	u_set(:,t_Nr) = u;
	% Integrate system
	[~, x_temp_set] = ode45(@(t,x)GP_CBF_Delay_Dynamics(t,x,u,Model_Parameter), ...
		[t, t + t_step], x);
	x_set(:,t_Nr + 1) = x_temp_set(end, :)';
	%
	if t >= PrintCounter
		fprintf('*');
		PrintCounter = PrintCounter + 1;
	end
end
fprintf('#\n');
%%
% if ~exist(SaveFolderName, 'dir')
% 	mkdir(SaveFolderName);
% end
% save([SaveFolderName,'/',SaveFileName,'.mat'], ...
% 	't_set','x_set','u_set','mu_set','delta_h','eta_GP_set','eta_Delay_set');
% end