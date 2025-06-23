function [H,f,A,b] = GP_CBF_Delay_HOCBF_QP_parameter(x,u_nom, ...
	mu,h,Model_Parameter,HOCBF_Parameter)
q = x(1:(numel(x) / 2));
dq = x((numel(x) / 2 + 1):end);
%%
M = Model_Parameter.M;
C = Model_Parameter.C;
G = Model_Parameter.G;
%%
alpha_0 = HOCBF_Parameter.alpha_0;
alpha_1 = HOCBF_Parameter.alpha_1;
%% Constraint
inv_M = inv(M);
h0 = h;
dh0dt = - 2 * q' * dq;
h1 = dh0dt + alpha_0 * h0;
b = - 2 * dq' * (alpha_0 * q + dq) + 2 * q' * inv_M * (C * dq + G - mu) + ...
	alpha_1 * h1;
A = 2 * q' * inv_M;
%% Cost Function
H = 2 * eye(numel(u_nom));
f = - 2 * u_nom;
end