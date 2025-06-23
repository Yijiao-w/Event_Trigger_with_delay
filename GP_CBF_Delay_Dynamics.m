function dx = GP_CBF_Delay_Dynamics(t,x,u,Model_Parameter)
dq = x((numel(x) / 2 + 1):end);
%%
d = ISSF_CBF_disturbance(x);
M = Model_Parameter.M;
C = Model_Parameter.C;
G = Model_Parameter.G;
%%
ddq = M \ (- C * dq - G + u + d);
dx = [dq; ddq];
end