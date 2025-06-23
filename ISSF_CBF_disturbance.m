function d = ISSF_CBF_disturbance(x)
q1 = x(1);
q2 = x(2);
%%
% d = 10 * [2 * sin(2 * q1); 
% 	3 * sin(1*q2)] + q1 + q2;
d1 = 10 * 2 * sin(2 * q1) + q1;
d2 = 10 * 3 * sin(1 * q2) + q2;
d = [d1; d2];
end