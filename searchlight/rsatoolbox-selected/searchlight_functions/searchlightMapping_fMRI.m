function [smm_dists, smm_corrs, n, searchlightRDMs] = searchlightMapping_fMRI(fullBrainVolumes, models, mask, userOptions, localOptions)

	% ARGUMENTS
	% fullBrainVolumes	A voxel x condition x session matrix of activity
	% 				patterns.
	%
	% models		A struct of model RDMs.
	%
	% mask     		A 3d or 4d mask to perform the searchlight in.
	%
	% userOptions and localOptions
	%
	% RETURN VALUES
	%               ToDo
    %
	% n             an array of the same dimensions as the volume, which
	%               indicates for each position how many voxels contributed
	%               data to the corresponding values of the infomaps.
	%               this is the number of searchlight voxels, except at the
	%               fringes, where the searchlight may illuminate voxels
	%               outside the input-data mask or voxel with all-zero
	%               time-courses (as can arise from head-motion correction).
	%
	% mappingMask_actual
	%               3D mask indicating locations for which valid searchlight
	%               statistics have been computed.
	%
	% Based on Niko Kriegeskorte's searchlightMapping_RDMs.m
	%
	% Additions by Cai Wingfield 2-2010:
	% 	- Now skips points in the searchlight where there's only one voxel inside.
	% 	- Now takes a userOptions struct for the input parameters.

	localOptions = setIfUnset(localOptions, 'averageSessions', true);

	%% Figure out whether to average over sessions or not
	if localOptions.averageSessions
		for sessionNumber = 1:size(fullBrainVolumes,3)
			thisSessionId = ['s' num2str(sessionNumber)];
			t_patsPerSession.(thisSessionId) = fullBrainVolumes(:,:,sessionNumber)';
		end%for:sessionNumber
	else
		justThisSession = 1;
		t_pats = fullBrainVolumes(:,:,justThisSession)';
		
		fprintf(['\nYou have selected not to average over sessions.\n         Only session number ' num2str(justThisSession) ' will be used.\n']);
		
	end%if

	%% Get parameters
	voxSize_mm = userOptions.voxelSize;
	searchlightRad_mm = userOptions.searchlightRadius;
	monitor = localOptions.monitor;
	nConditions = size(fullBrainVolumes, 2);
	
	clear fullBrainVolumes;

	% Prepare models
	modelRDMs_ltv = permute(unwrapRDMs(vectorizeRDMs(models)), [3 2 1]);

	% Prepare masks
	mask(isnan(mask)) = 0; % Just in case!
	if ndims(mask)==3
		inputDataMask=logical(mask);
		mappingMask_request=logical(mask);
	else
		inputDataMask=logical(mask(:,:,:,1));
		mappingMask_request=logical(mask(:,:,:,2));
	end

	% Check to see if there's more data than mask...
	if localOptions.averageSessions
		for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
			thisSessionId = ['s' num2str(sessionNumber)];
			t_patsPerSession.(thisSessionId) = t_patsPerSession.(thisSessionId)(:, inputDataMask(:));
		end%for:sessionNumber
	else
		if (size(t_pats,2)>sum(inputDataMask(:)))
			t_pats=t_pats(:,inputDataMask(:));
		end%if
	end%if

	% Other data
	volSize_vox=size(inputDataMask);
	nModelRDMs=size(modelRDMs_ltv,1);
	rad_vox=searchlightRad_mm./voxSize_mm;
	minMargin_vox=floor(rad_vox);


	%% create spherical multivariate searchlight
	[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
	sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
	sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )

	if monitor, figure(50); clf; showVoxObj(sphere); end % show searchlight in 3D

	% compute center-relative sphere SUBindices
	[sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
	sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
	ctrSUB=sphereSize_vox/2+[.5 .5 .5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
	ctrRelSphereSUBs=sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices

	nSearchlightVox=size(sphereSUBs,1);


	%% define masks
	validInputDataMask=inputDataMask;

	if localOptions.averageSessions
		for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
			thisSessionId = ['s' num2str(sessionNumber)];
			sumAbsY=sum(abs(t_patsPerSession.(thisSessionId)),1);
		end%for:sessionNumber
	else
		sumAbsY=sum(abs(t_pats),1);
	end%if

	validYspace_logical= (sumAbsY~=0) & ~isnan(sumAbsY); clear sumAbsY;
	validInputDataMask(inputDataMask)=validYspace_logical; % define valid-input-data brain mask

	if localOptions.averageSessions
		for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
			thisSessionId = ['s' num2str(sessionNumber)];
			t_patsPerSession.(thisSessionId) = t_patsPerSession.(thisSessionId)(:,validYspace_logical);
			nVox_validInputData=size(t_patsPerSession.(thisSessionId),2);
		end%for:sessionNumber
	else
		t_pats=t_pats(:,validYspace_logical); % reduce t_pats to the valid-input-data brain mask
		nVox_validInputData=size(t_pats,2);
	end%if

	mappingMask_request_INDs=find(mappingMask_request);
	nVox_mappingMask_request=length(mappingMask_request_INDs);

	if monitor
		disp([num2str(round(nVox_mappingMask_request/prod(volSize_vox)*10000)/100),'% of the cuboid volume requested to be mapped.']);
		disp([num2str(round(nVox_validInputData/prod(volSize_vox)*10000)/100),'% of the cuboid volume to be used as input data.']);
		disp([num2str(nVox_validInputData),' of ',num2str(sum(inputDataMask(:))),' declared input-data voxels included in the analysis.']);
	end

	volIND2YspaceIND=nan(volSize_vox);
	volIND2YspaceIND(validInputDataMask)=1:nVox_validInputData;

	% n voxels contributing to infobased t at each location
	n=nan(volSize_vox);

	%% similarity-graph-map the volume with the searchlight
	smm_bestModel=nan(volSize_vox);
    
	smm_dists.Riem_bc=nan([volSize_vox, nModelRDMs]);
    smm_dists.Frob_bc=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Spea_bc=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Pear_bc=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Kend_bc=nan([volSize_vox, nModelRDMs]);
    
    smm_dists.Riem=nan([volSize_vox, nModelRDMs]);
    smm_dists.Frob=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Spea=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Pear=nan([volSize_vox, nModelRDMs]);
    smm_corrs.Kend=nan([volSize_vox, nModelRDMs]);
    
	searchlightRDMs = nan([nConditions, nConditions, volSize_vox]);

	if monitor
		h_progressMonitor=progressMonitor(1, nVox_mappingMask_request,  'Similarity-graph-mapping...');
	end

	%% THE BIG LOOP! %%
    
	for cMappingVoxI=1:nVox_mappingMask_request
		
		if mod(cMappingVoxI,10000)==0
			if monitor
				progressMonitor(cMappingVoxI, nVox_mappingMask_request, 'Searchlight mapping Mahalanobis distance...', h_progressMonitor);
				%                 cMappingVoxI/nVox_mappingMask_request
			else
				fprintf('.');
			end%if
		end%if

		[x, y, z]=ind2sub(volSize_vox,mappingMask_request_INDs(cMappingVoxI));

		% compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;

		% exclude out-of-volume voxels
		outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
					cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
					cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));

		cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);

		% list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));

		% restrict searchlight to voxels inside validDataBrainMask
		cIllValidVox_volINDs=cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
		cIllValidVox_YspaceINDs=volIND2YspaceIND(cIllValidVox_volINDs);

		% note how many voxels contributed to this locally multivariate stat
		n(x,y,z)=length(cIllValidVox_YspaceINDs);
		
		if n(x,y,z) < 2, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
		
        % ---------------- TBD ----------------
		if localOptions.averageSessions
			searchlightRDM = zeros(localOptions.nConditions, localOptions.nConditions);
			for session = 1:localOptions.nSessions
				sessionId = ['s' num2str(session)];
				searchlightRDM = searchlightRDM + squareform(pdist(t_patsPerSession.(sessionId)(:,cIllValidVox_YspaceINDs),'correlation'));
			end%for:sessions
			searchlightRDM = searchlightRDM / localOptions.nSessions;
		else
			searchlightRDM = squareform(pdist(t_pats(:,cIllValidVox_YspaceINDs), 'correlation'));
		end%if
		
		searchlightRDM = vectorizeRDM(searchlightRDM);
		
		% Locally store the full brain's worth of indexed RDMs.
		searchlightRDMs(:, :, x, y, z) = squareform(searchlightRDM);
		
        switch userOptions.RDMcorrelationType
            case 'all'
                dists.Riem_bc = nan(1, nModelRDMs);
                dists.Frob_bc = nan(1, nModelRDMs);
                corrs.Spea_bc = nan(1, nModelRDMs);
                corrs.Pear_bc = nan(1, nModelRDMs);
                corrs.Kend_bc = nan(1, nModelRDMs);
                
                dists.Riem = nan(1, nModelRDMs);
                dists.Frob = nan(1, nModelRDMs);
                corrs.Spea = nan(1, nModelRDMs);
                corrs.Pear = nan(1, nModelRDMs);
                corrs.Kend = nan(1, nModelRDMs);
                
                for model_i = 1 : nModelRDMs
                    sl_RSM = 1 - squareform(searchlightRDM);
                    mdl_RSM = 1 - squareform(modelRDMs_ltv(model_i, :));
                    [~, dists.Riem(model_i), ~, dists.Riem_bc(model_i)] = pval_perm_Riemann(sl_RSM, mdl_RSM, true, 20);
                    [~, dists.Frob(model_i), ~, dists.Frob_bc(model_i)] = pval_perm_Frobenius(sl_RSM, mdl_RSM, true, 20);
                    [~, corrs.Spea(model_i), ~, corrs.Spea_bc(model_i)] = pval_perm_Spearman(sl_RSM, mdl_RSM, true, 20);
                    [~, corrs.Pear(model_i), ~, corrs.Pear_bc(model_i)] = pval_perm_Pearson(sl_RSM, mdl_RSM, true, 20);
                    [~, corrs.Kend(model_i), ~, corrs.Kend_bc(model_i)] = pval_perm_Kendall(sl_RSM, mdl_RSM, true, 20);                
                end
                
                if localOptions.fisher
                    %ToDo
                end%if
                
            otherwise
                warning('RDMCorrelationType is not valid!')
        end
		
        % ToDo: if we had multiple models
		smm_dists.Riem_bc(x,y,z,:) = dists.Riem_bc;
        smm_dists.Frob_bc(x,y,z,:) = dists.Frob_bc;
        smm_corrs.Spea_bc(x,y,z,:) = corrs.Spea_bc;
        smm_corrs.Pear_bc(x,y,z,:) = corrs.Pear_bc;
        smm_corrs.Kend_bc(x,y,z,:) = corrs.Kend_bc;
        
        smm_dists.Riem(x,y,z,:) = dists.Riem;
        smm_dists.Frob(x,y,z,:) = dists.Frob;
        smm_corrs.Spea(x,y,z,:) = corrs.Spea;
        smm_corrs.Pear(x,y,z,:) = corrs.Pear;
        smm_corrs.Kend(x,y,z,:) = corrs.Kend;
		
	end%for:cMappingVoxI

	%% END OF THE BIG LOOP! %%
    if monitor
        fprintf('\n');
		close(h_progressMonitor);
    end
    % ---------- Commented by Ali ------------
	%mappingMask_actual=mappingMask_request;
	%mappingMask_actual(isnan(sum(smm_rs,4)))=0;

end%function
