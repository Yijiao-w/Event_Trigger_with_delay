function [mu,GP] = GP_CBF_Delay_Previous_Prediction_Output( ...
	x_GP,GP,do_GP_update)
SigmaN = GP.SigmaN;
y_dim = GP.y_dim;
%%
if do_GP_update == true
	y_GP = ISSF_CBF_disturbance(x_GP) + SigmaN * randn(y_dim,1);
	if GP.DataQuantity >= GP.MaxDataQuantity
		GP.downdateParam(1);
	end
	GP.addPoint(x_GP,y_GP);
end
mu = GP.predict(x_GP);

end