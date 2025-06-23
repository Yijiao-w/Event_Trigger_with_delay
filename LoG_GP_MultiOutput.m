classdef LoG_GP_MultiOutput < handle
	properties (Access = public)
		Max_LocalGP_DataQuantity;
		Max_LocalGP_Quantity;
		x_dim;
		y_dim;
		SigmaF;
		SigmaL;
		SigmaN;
		o_ratio = 1/10;

		tau = 1e-6;
		delta = 0.05;
		% 		Lf = 10;

		DataQuantity = 0;
		MaxDataQuantity;
	end
	properties (Access = public)
		LocalGP_set;
		ActivatedGPNr;
		UnactivatedGPNr;
		ActivatedNodeNr;
		UnactivateNodeNr;

		ActivatedNodeQuantity
		ActivatedGPQuantity;
		Node_GP_Map;

		RootModel;
		children;
		% 		LocalGPNr_set;

		HyperplaneMean;
		HyperplaneDimension;
		HyperplaneOverlap;

		AgeOfLocalGP;

		AggregationMethod = 'MOE';
		DataSaturation = false;
	end
	methods
		function obj = LoG_GP_MultiOutput( ...
				Max_LocalGP_DataQuantity,Max_LocalGP_Quantity, ...
				x_dim_in,y_dim_in, ...
				SigmaN,SigmaF,SigmaL)
			obj.Max_LocalGP_DataQuantity = Max_LocalGP_DataQuantity;
			obj.Max_LocalGP_Quantity = Max_LocalGP_Quantity;
			obj.x_dim =	x_dim_in;
			obj.y_dim = y_dim_in;
			obj.DataQuantity = 0;

			obj.MaxDataQuantity = obj.Max_LocalGP_Quantity * obj.Max_LocalGP_DataQuantity;

			obj.SigmaN = SigmaN;
			obj.SigmaF = SigmaF;
			obj.SigmaL = SigmaL;

			obj.RootModel = 1;
			obj.LocalGP_set = cell(obj.Max_LocalGP_Quantity,1);
			for i = 1:obj.Max_LocalGP_Quantity
				obj.LocalGP_set{i} = LocalGP_MultiOutput( ...
					obj.x_dim,obj.y_dim,obj.Max_LocalGP_DataQuantity, ...
					obj.SigmaN,obj.SigmaF,obj.SigmaL);
				obj.LocalGP_set{i}.tau = obj.tau;
			end
			obj.LocalGP_set{1}.ActivateState = true;

			obj.children = nan(2 * obj.Max_LocalGP_Quantity - 1,2); % first column: Left Node Nr; second column: Right Node Nr
			obj.children(1,:) = -1;

			obj.ActivatedGPNr = [1;nan(obj.Max_LocalGP_Quantity-1,1)];
			obj.UnactivatedGPNr = [(2:obj.Max_LocalGP_Quantity)';nan];
			obj.ActivatedGPQuantity = 1;

			obj.ActivatedNodeNr = nan(2 * obj.Max_LocalGP_Quantity - 1,1);
			obj.ActivatedNodeNr(1) = 1;
			obj.UnactivateNodeNr = nan(2 * obj.Max_LocalGP_Quantity - 1,1);
			obj.UnactivateNodeNr(1:end-1) = 2:(2 * obj.Max_LocalGP_Quantity - 1);
			obj.ActivatedNodeQuantity = 1;
			obj.Node_GP_Map = nan(obj.Max_LocalGP_Quantity,1);
			obj.Node_GP_Map(1) = 1;

			obj.HyperplaneDimension = nan(2 * obj.Max_LocalGP_Quantity - 1,1);
			obj.HyperplaneMean = nan(2 * obj.Max_LocalGP_Quantity - 1,1);
			obj.HyperplaneOverlap = nan(2 * obj.Max_LocalGP_Quantity - 1,1);

			obj.AgeOfLocalGP = nan(obj.Max_LocalGP_Quantity,1);
			obj.AgeOfLocalGP(1) = 0;
		end
		%% Set Global X range
		function set_gloabl_x_range(obj,xMin,xMax)
			for i = 1:obj.Max_LocalGP_Quantity
				obj.LocalGP_set{i}.xMax = xMax;
				obj.LocalGP_set{i}.xMin = xMin;
			end
		end
		%% Update
		function flag = update(obj,x,y,Kickout)
			if nargin < 4
				Kickout = true;
			end
			flag = 0;
			addpoint_flag = inf;
			deledt_flag = inf;
			saturation_flag = inf;
			devide_flag = inf;
			%
			if obj.DataQuantity >= obj.MaxDataQuantity
				obj.DataSaturation = true;
				saturation_flag = -3;
				flag = min([addpoint_flag,deledt_flag,saturation_flag,devide_flag]);
				return;
			end
			%
			obj.DataQuantity = obj.DataQuantity + 1;
			NodeNr = obj.RootModel;
			while obj.children(NodeNr,1) ~= -1 %if model is a parent
				%search for the leaf to asign the point
				[pL, ~] = obj.activation(x, NodeNr);
				if pL >= rand()
					NodeNr = obj.children(NodeNr,1);%left child
				else
					NodeNr = obj.children(NodeNr,2);
				end
			end
			%add the model to the randomly selected model
			LocalGPNr = find(obj.Node_GP_Map == NodeNr);
			if obj.LocalGP_set{LocalGPNr}.DataQuantity < obj.Max_LocalGP_DataQuantity
				obj.LocalGP_set{LocalGPNr}.addPoint(x,y);
				addpoint_flag = 1; % Successfully add data into the local GP
			else
				% warning('wrong!')
			end
			if obj.LocalGP_set{LocalGPNr}.DataQuantity == obj.Max_LocalGP_DataQuantity 
				if (obj.Max_LocalGP_Quantity - obj.ActivatedGPQuantity) <= 0 || ...
						(2 * obj.Max_LocalGP_Quantity - 1 - obj.ActivatedNodeQuantity) < 2
					% all GP nodes are full and occupied
					if Kickout
						obj.delete_Node_GP(LocalGPNr);
						deledt_flag = -4;
					else
						obj.DataSaturation = true;
						saturation_flag = -3;
						flag = min([addpoint_flag,deledt_flag,saturation_flag,devide_flag]);
						return;
					end
				end
				% divide the current GP
				if (obj.Max_LocalGP_Quantity - obj.ActivatedGPQuantity) > 0 && ...
						(2 * obj.Max_LocalGP_Quantity - 1 - obj.ActivatedNodeQuantity) >= 2
					% there is at least one GP not activated
					% there are at least 2 free nodes, which can be used as children
					% determine the new activated GP1;
					obj.devide(NodeNr,LocalGPNr);
					devide_flag = -2;
				end
			end
			flag = min([addpoint_flag,deledt_flag,saturation_flag,devide_flag]);
		end
		%% Delete Node and GP
		function delete_Node_GP(obj,LocalGPNr,Unused_LocalGPNr)
			if nargin < 3 % obj,LocalGPNr
				All_AgeOfLocalGP = obj.AgeOfLocalGP;
				All_AgeOfLocalGP(LocalGPNr) = -1;
				Unused_LocalGPNr = find(All_AgeOfLocalGP == max(All_AgeOfLocalGP));
				Unused_LocalGPNr = Unused_LocalGPNr(1);
			end
			obj.AgeOfLocalGP(Unused_LocalGPNr) = nan;
			Unused_NodeNr = obj.Node_GP_Map(Unused_LocalGPNr);
			% delete_LocalGP
			DeleteDataQuantity = obj.LocalGP_set{Unused_LocalGPNr}.DataQuantity;
			obj.DataQuantity = obj.DataQuantity - DeleteDataQuantity;

			obj.Node_GP_Map(Unused_LocalGPNr) = nan;
			obj.LocalGP_set{Unused_LocalGPNr}.resetGP;
			new_ActivatedGPNr = rmmissing(obj.ActivatedGPNr);
			new_ActivatedGPNr(new_ActivatedGPNr == Unused_LocalGPNr) = nan;
			new_UnactivatedGPNr = rmmissing(obj.UnactivatedGPNr);
			new_UnactivatedGPNr = [new_UnactivatedGPNr;Unused_LocalGPNr];

			obj.ActivatedGPQuantity = obj.ActivatedGPQuantity - 1;
			obj.ActivatedGPNr(1:obj.ActivatedGPQuantity) = rmmissing(new_ActivatedGPNr);
			obj.ActivatedGPNr(obj.ActivatedGPQuantity+1:end) = nan;
			obj.UnactivatedGPNr(1:(obj.Max_LocalGP_Quantity - obj.ActivatedGPQuantity)) = new_UnactivatedGPNr;
			obj.UnactivatedGPNr((obj.Max_LocalGP_Quantity - obj.ActivatedGPQuantity+1):end) = nan;
			% change node relationship
			% find cousin node
			[ParentNodeNr,ChildrenPos] = find(obj.children == Unused_NodeNr);
			if ChildrenPos == 1
				CousinNodeNr = obj.children(ParentNodeNr,2);
			elseif ChildrenPos == 2
				CousinNodeNr = obj.children(ParentNodeNr,1);
			end
			% Cousin Node becomes Parent Node
			[GrandParentNode,ParentPos] = find(obj.children == ParentNodeNr);
			if isempty(GrandParentNode) % Parent Node is root node
				obj.RootModel = CousinNodeNr;
			else % Parent Node isn't root node, Grand parent node exists
				% Grand Parent Node connect directly to cousin node
				if ParentPos == 1
					obj.children(GrandParentNode,1) = CousinNodeNr;
				elseif ParentPos == 2
					obj.children(GrandParentNode,2) = CousinNodeNr;
				end
			end
			obj.children(ParentNodeNr,:) = nan;
			obj.children(Unused_NodeNr,:) = nan;

			new_ActivatedNodeNr = rmmissing(obj.ActivatedNodeNr);
			new_ActivatedNodeNr(new_ActivatedNodeNr == ParentNodeNr) = nan;
			new_ActivatedNodeNr(new_ActivatedNodeNr == Unused_NodeNr) = nan;
			new_UnactivatedNodeNr = [...
				obj.UnactivateNodeNr(1:(2 * obj.Max_LocalGP_Quantity - 1 -obj.ActivatedNodeQuantity)); ...
				ParentNodeNr;Unused_NodeNr];
			obj.ActivatedNodeQuantity = obj.ActivatedNodeQuantity - 2;
			obj.ActivatedNodeNr(1:obj.ActivatedNodeQuantity) = rmmissing(new_ActivatedNodeNr);
			obj.ActivatedNodeNr(obj.ActivatedNodeQuantity+1:end) = nan;
			obj.UnactivateNodeNr(1:(2 * obj.Max_LocalGP_Quantity - 1 -obj.ActivatedNodeQuantity)) = new_UnactivatedNodeNr;
			obj.UnactivateNodeNr((2 * obj.Max_LocalGP_Quantity - 1 -obj.ActivatedNodeQuantity + 1):end) = nan;
		end
		%% Delete one data in one local GP
		function delete_Point_GP(obj,LocalGPNr,PointNr)
			obj.LocalGP_set{LocalGPNr}.downdateParam(PointNr);
			obj.DataQuantity = obj.DataQuantity - 1;
		end
		%% Devide
		function devide(obj,NodeNr,LocalGPNr)
			Extend_LocalGPNr = rmmissing(obj.UnactivatedGPNr);
			Extend_LocalGPNr = Extend_LocalGPNr(1);
			% divide data to new activated GP
			obj.divideGP(NodeNr,LocalGPNr,Extend_LocalGPNr);
			% update activated and unactivated node information
			obj.ActivatedGPQuantity = obj.ActivatedGPQuantity + 1;
			obj.ActivatedGPNr(obj.ActivatedGPQuantity) = Extend_LocalGPNr;
			obj.UnactivatedGPNr(obj.UnactivatedGPNr == Extend_LocalGPNr) = nan;
			% determine the new activated node
			new_ActivadteNodeNr = rmmissing(obj.UnactivateNodeNr);
			new_ActivadteNodeNr = new_ActivadteNodeNr(1:2);
			% update activated and unactivated node information
			UnactivateNodeNr_temp = obj.UnactivateNodeNr;
			UnactivateNodeNr_temp(UnactivateNodeNr_temp == new_ActivadteNodeNr(1)) = nan;
			UnactivateNodeNr_temp(UnactivateNodeNr_temp == new_ActivadteNodeNr(2)) = nan;
			UnactivateNodeNr_temp = rmmissing(UnactivateNodeNr_temp);
			obj.UnactivateNodeNr(1:numel(UnactivateNodeNr_temp)) = UnactivateNodeNr_temp;

% 			obj.UnactivateNodeNr(obj.UnactivateNodeNr == new_ActivadteNodeNr(1)) = nan;
% 			obj.UnactivateNodeNr(obj.UnactivateNodeNr == new_ActivadteNodeNr(2)) = nan;
			obj.ActivatedNodeNr(obj.ActivatedNodeQuantity + (1:2)) = new_ActivadteNodeNr;
			obj.ActivatedNodeQuantity = obj.ActivatedNodeQuantity + 2;
			% change the parent-children relationship
			obj.children(NodeNr,:) = new_ActivadteNodeNr;
			obj.children(new_ActivadteNodeNr,:) = -1;
			%
			obj.Node_GP_Map(LocalGPNr) = new_ActivadteNodeNr(1);
			obj.Node_GP_Map(Extend_LocalGPNr) = new_ActivadteNodeNr(2);
			%
			obj.children(NodeNr,1) = new_ActivadteNodeNr(1);
			obj.children(NodeNr,2) = new_ActivadteNodeNr(2);
			obj.children(new_ActivadteNodeNr(1),:) = -1;
			obj.children(new_ActivadteNodeNr(2),:) = -1;
		end
		%% Divide Gaussian Process
		function divideGP(obj,ModelNr,LocalGPNr,Extend_LocalGPNr)
			All_X = obj.LocalGP_set{LocalGPNr}.X;
			All_y = obj.LocalGP_set{LocalGPNr}.Y;
			% Find the Hyperplane in one dimension
			[~,cutD]=max((max(All_X,[],2) - min(All_X,[],2)));
			mP = mean(All_X(cutD,:));
			maxV = max(All_X(cutD,:));
			minV = min(All_X(cutD,:));
			o  = (maxV-minV) * obj.o_ratio;
			% Distribute the data
			lcount = 0;
			rcount = 0;
			iL = nan(1,obj.Max_LocalGP_DataQuantity); %left index order vector
			iR = nan(1,obj.Max_LocalGP_DataQuantity); %right index order vector
			while true
				for i = 1:obj.Max_LocalGP_DataQuantity
					xD = All_X(cutD,i);
					if xD < mP - o/2 %if in left set
						lcount = lcount+1;
						iL(lcount) = i;
					elseif xD >= mP - o/2 && xD <= mP + o/2 %if in overlapping
						pL = 0.5 + (xD - mP) / o;
						if pL>=rand() %left side
							lcount = lcount+1;
							iL(lcount) = i;
						else
							rcount = rcount + 1;
							iR(rcount) = i;
						end
					elseif xD > mP + o/2 %if in right
						rcount = rcount + 1;
						iR(rcount) = i;
					end
					iL = rmmissing(iL);
					iR = rmmissing(iR);

				end
				if numel(iL) > 0 && numel(iR) > 0
					break;
				end
			end
			% Put data in new Local GP Object
			obj.LocalGP_set{Extend_LocalGPNr}.ActivateState = true;
			obj.LocalGP_set{Extend_LocalGPNr}.DataQuantity = numel(iR);
			obj.LocalGP_set{Extend_LocalGPNr}.X(:,1:numel(iR)) = All_X(:,iR);
			obj.LocalGP_set{Extend_LocalGPNr}.Y(1:numel(iR),:) = All_y(iR,:);
			obj.LocalGP_set{Extend_LocalGPNr}.updateParam_Alldata;
			obj.AgeOfLocalGP(Extend_LocalGPNr) = 0;
			% Change data in old Local GP Object
			obj.LocalGP_set{LocalGPNr}.ActivateState = true;
			obj.LocalGP_set{LocalGPNr}.DataQuantity = numel(iL);
			obj.LocalGP_set{LocalGPNr}.X(:,1:numel(iL)) = All_X(:,iL);
			obj.LocalGP_set{LocalGPNr}.Y(1:numel(iL),:) = All_y(iL,:);
			obj.LocalGP_set{LocalGPNr}.updateParam_Alldata;
			% Record Hyperplanr Information
			obj.HyperplaneDimension(ModelNr) = cutD;
			obj.HyperplaneMean(ModelNr) = mP;
			obj.HyperplaneOverlap(ModelNr) = o;
		end
		%% Determine which children x belongs
		function [pL,pR] = activation(obj, x, model)
			mP = obj.HyperplaneMean(model);
			cutD = obj.HyperplaneDimension(model);
			xD = x(cutD); %x value in cut dimension
			o = obj.HyperplaneOverlap(model); %half of the overlapping region
			if xD < mP - o/2
				pL = 1;
			elseif  xD >= mP - o/2 && xD <= mP + o/2 %if in overlapping
				pL = 0.5 + (xD - mP) / o;
				if(pL <= 1e-12)%avoid numerical errors
					pL = 0;
				elseif(pL >= 1 - 1e-12)
					pL = 1;
				end
			else
				pL = 0;
			end
			pR = 1 - pL;
		end
		%% prediction
		function [mu,var,likelyhood,eta,eta_max,Na, ...
				eta_pre,beta_m_set,gamma_set,p_set] = predict(obj,x,y)
			if obj.DataQuantity == 0
				mu = zeros(obj.y_dim,1);
				var = obj.SigmaF ^ 2 * ones(obj.y_dim,1);
				likelyhood = nan;
				eta = nan;
				eta_max = nan;
				Na = 0;
				eta_pre = nan;
				beta_m_set = nan;
				gamma_set = nan;
				return;
			end
			% calculate the probability
			moP = nan(2 * obj.Max_LocalGP_Quantity - 1,2);
			mCount = 1;% start from the root
			moP(1,1) = obj.RootModel;
			moP(1,2) = 1;
			while ~isequal( obj.children(moP(1:mCount,1),:) , -1*ones(mCount,2) )
				for j = 1:mCount
					if ~isequal( obj.children(moP(j,1),:) , -1*ones(1,2) )
						[pL, pR] = obj.activation(x,moP(j,1));
						if pL > 0 && pR == 0
							moP(j,1) = obj.children(moP(j,1),1);
							moP(j,2) = moP(j,2)*pL;
						elseif pR > 0 && pL == 0
							moP(j,1) = obj.children(moP(j,1),2);
							moP(j,2) = moP(j,2) * pR;
						elseif pL > 0 && pR > 0
							mCount = mCount + 1;
							moP(mCount,1) = obj.children(moP(j,1),2);
							moP(mCount,2) = moP(j,2) * pR;

							moP(j,1) = obj.children(moP(j,1),1);
							moP(j,2) = moP(j,2) * pL;
						end
					end
				end
			end
			% update the age of GP
			obj.AgeOfLocalGP = obj.AgeOfLocalGP + 1;


			% aggregation
			mu = zeros(obj.y_dim,1);
			var = zeros(obj.y_dim,1);
			eta = 0;
			eta_max = 0;
			eta_pre = inf;
			beta_m_set = nan(mCount,1);
			gamma_set = nan(mCount,1);
			p_set = moP(:,2);
			for i=1:mCount
				NodeNr = moP(i,1);
				LocalGPNr = obj.Node_GP_Map == NodeNr;
				% Note before prediction, update the delta of activated local GP
				obj.LocalGP_set{LocalGPNr}.delta = obj.delta / mCount;
				[mu_m,var_m,eta_m,beta_m,gamma_m,eta_max_m,eta_pre_m] = obj.LocalGP_set{LocalGPNr}.predict(x);
				switch obj.AggregationMethod
					case 'MOE'
						% MOE
						mu = mu + moP(i,2) * mu_m;
						var = var + (var_m + mu_m.^ 2) * moP(i,2);
					case 'GPOE'
						mu = mu + moP(i,2) * mu_m ./ var_m;
						var = var + moP(i,2) ./ var_m;
				end
				eta = eta + moP(i,2) * eta_m;
				eta_max = eta_max + moP(i,2) * eta_max_m;
				eta_pre = min(eta_pre,eta_pre_m);
				beta_m_set(i) = beta_m;
				gamma_set(i) = gamma_m;
				% update age of activated GP
				obj.AgeOfLocalGP(obj.Node_GP_Map ==	moP(i,1)) = 0;
			end
			switch obj.AggregationMethod
				case 'MOE'
					% MOE
					var = var - mu .^ 2;
				case 'GPOE'
					var = 1 ./ var;
					mu = mu .* var;
			end
			% likelyhood
			likelyhood = log(max(1e-300,normpdf(y,mu,sqrt(var+obj.SigmaN^2))));
			% Activated Number
			Na = mCount;
		end
	end
end