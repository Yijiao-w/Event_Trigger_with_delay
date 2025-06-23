function u_nom = GP_CBF_Delay_NominalController( ...
	t,x,Model_Parameter,PD_Controller_Parameter)
%%
[x_ref,ddq_ref] = GP_CBF_Delay_Reference(t);
q_ref = x_ref(1:(numel(x_ref) / 2));
dq_ref = x_ref((numel(x_ref) / 2 + 1):end);
%%
q = x(1:(numel(x) / 2));
dq = x((numel(x) / 2 + 1):end);
%%
Kp = PD_Controller_Parameter.Kp;
Kd = PD_Controller_Parameter.Kd;
%%
e = q - q_ref;
de = dq - dq_ref;
PID = ddq_ref - Kp * e - Kd * de;
%%
M = Model_Parameter.M;
C = Model_Parameter.C;
G = Model_Parameter.G;
u_nom = M * PID + C * dq + G - ISSF_CBF_disturbance(x);
end