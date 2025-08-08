function d = GP_CBF_Delay_Disturbance(x)
q1 = x(1,:);
q2 = x(2,:);
%%
% d = 10 * [2 * sin(2 * q1); 
% 	3 * sin(1*q2)] + q1 + q2;
d1 = 1 * 5 * cos(5 * q1) + 5 * q1;
d2 = 1 * 10 * cos(8 * q2) + 10 * q2;
d = [d1; d2];
end