function u = GP_CBF_Delay_HOCBF_QP_Solver(SolverType, ...
	x,u_nom,mu,h,Model_Parameter,HOCBF_Parameter)

switch SolverType
	case 'numeric'
		[Q,f_qp,A_qp,b_qp] = GP_CBF_Delay_HOCBF_QP_parameter( ...
			x,u_nom,mu,h,Model_Parameter,HOCBF_Parameter);
		u = quadprog(Q, f_qp, A_qp, b_qp, [], [], ...
			[], [], u_nom, optimoptions('quadprog','Display','off'));
	case 'analytic'
		[~,~,A_qp,b_qp] = GP_CBF_Delay_HOCBF_QP_parameter( ...
			x,u_nom,mu,h,Model_Parameter,HOCBF_Parameter);
		if all(A_qp == 0)
			u = u_nom;
		else
			u = u_nom - ...
				max(0, - b_qp + A_qp * u_nom) * A_qp' / (A_qp * A_qp');
		end
end

end