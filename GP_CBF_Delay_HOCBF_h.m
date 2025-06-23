function h_set = GP_CBF_Delay_HOCBF_h(x_set)
SampleQuantity = size(x_set,2);
h_set = nan(SampleQuantity,1);
for SampleNr = 1:SampleQuantity
	x = x_set(:,SampleNr);
	q = x(1:(numel(x) / 2));
	h_set(SampleNr) = 1 - q' * q;
end
end