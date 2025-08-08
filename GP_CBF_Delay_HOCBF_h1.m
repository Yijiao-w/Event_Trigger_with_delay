function h1_set = GP_CBF_Delay_HOCBF_h1(h0_set,x_set,HOCBF_Parameter)

alpha_0 = HOCBF_Parameter.alpha_0;
SampleQuantity = size(x_set,2);
h1_set = nan(SampleQuantity,1);
for SampleNr = 1:SampleQuantity
	x = x_set(:,SampleNr);
	q = x(1:(numel(x) / 2));
	dq = x((numel(x) / 2 + 1):end);

	dh0dt = - 2 * q' * dq;
	h0 = h0_set(SampleNr);
	h1_set(SampleNr) = dh0dt + alpha_0 * h0;
end

end