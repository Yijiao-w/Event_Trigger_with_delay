function Model_Parameter = GP_CBF_Delay_Manipulator2D_get_Model_Parameter( ...
	x,L1,L2,m1,m2)
r = x(1:2);
dr = x(3:4);
%% Convert to joint coordinate
q = Manipulator_2D_2DoF_get_q_from_r(r, L1, L2);
q1 = q(1);
q2 = q(2);

invT = Manipulator_2D_2DoF_get_invT(q, L1, L2);
dq = invT * dr;
dq1 = dq(1);
dq2 = dq(2);
%% Mass and Force in joint coordinate
[Mass_q,DampingMatrix_q,Gravity_q] = Manipulator_2D_2DoF_get_MassForce_q( ...
	q,dq,L1,L2,m1,m2);

Force_Coordinate_Convertion_1 = L1 * cos(q1) * dq1^2 + L2 * cos(q1 + q2) * (dq1 + dq2)^2;
Force_Coordinate_Convertion_2 = L1 * sin(q1) * dq1^2 + L2 * sin(q1 + q2) * (dq1 + dq2)^2;
Force_Coordinate_Convertion = [Force_Coordinate_Convertion_1; Force_Coordinate_Convertion_2];
%%
Mass_r = invT' * Mass_q * invT;

%%
Model_Parameter.M = Mass_r;
Model_Parameter.C = invT' * DampingMatrix_q;
Model_Parameter.G = Gravity_q + invT' * Mass_q * invT * Force_Coordinate_Convertion;

end