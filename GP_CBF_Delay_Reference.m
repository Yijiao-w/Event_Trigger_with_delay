function [x_ref,ddq_ref] = GP_CBF_Delay_Reference(t)
a = 1;   % x 方向半轴
b = 0.6; % y 方向半轴
c = 0.3;

q_ref = [a * cos(c * t) - 0.3; b * sin(c * t)];
dq_ref = [-c * a * sin(c * t); c * b * cos(c * t)];
ddq_ref = [-c^2 * a * cos(c * t); -c^2 * b * sin(c * t)];

x_ref = [q_ref; dq_ref];
end