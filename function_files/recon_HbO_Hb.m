function [HbO_c Hb_c HbT_c] = recon_HbO_Hb(PD_data,weekN, goodch_idx)

    if weekN == 30;
        load('Jacobian_infant30week_850nm.mat')
        load('mesh_infant30week_850nm.mat')
    end
    J_amp = J.complete; %amp - We just need this part of the Jacobian as it describes the intensity
    %EACH ROW IS A MEASUREMENT CORROSPONDING TO THE MESH.LINK
    %EACH COLUMN IS A NODE CORROSPONDING TO THE MESH.NODE

    %need to change remCh, to index the row number. e.g 1 2 3 4 7 11 ... not
    %just 1 1 1 0 0 0 1 1 etc. 
    Ylhs_wv1 = [PD_data.dod(:,goodch_idx(1:end/2,1))]'; %780
    Ylhs_wv2 = [PD_data.dod(:,goodch_idx((end/2)+1:end,1))]'; %850

    % % mua Amplitude only
    J = [J_amp(goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.

    %do inversion 
    % Obtain specific absorption coefficients
    wavelengths = PD_data.SD.Lambda; % wavelengths of the system
    nWavs = length(PD_data.SD.Lambda); % n of wavelengths
    Eall = [];
    for i = 1:nWavs
            Etmp = GetExtinctions(wavelengths(i));
            Etmp = Etmp(1:2); %HbO and HbR only
            Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
    end
    E=Eall;

    %For re running this section, set J. mua Amplitude only
    J = [J_amp(goodch_idx(1:end/2,1),:)]; %amp, mua only %goodch idx is symmetrical.
    % normalise for voxels
    JTJ = sum(J(:,1:end/2).^2,1);
    L1 = sqrt(JTJ + (1e-2*max(JTJ)));
    JTJ = sum(J(:,end/2+1:end).^2,1);
    L2 = sqrt(JTJ + (1e-2*max(JTJ)));
    %Can set to 1 or normalize
    %L = [L1 L2]; 
    L(:) = 1; 

    J = bsxfun(@rdivide,J,L);
    % normalise for data magnitude
    JJT = sum(J.^2,2);
    %Can set to 1 or normalize
    %M = sqrt(JJT + (1e-2*max(JJT)));
    M(:) = 1; 
    J = bsxfun(@rdivide,J,M);

    lambda = 1E-2;
    % get update
    Hess = J*J';
    reg = eye(size(Hess)).*lambda.*sqrt(norm(Hess)); %can set to 1 or use regularisation
    %higher reg = more blur , lower reg = high contrast but may overfit to
    %noise
    %
    tmp = J'*((Hess+reg)\(diag(1./M)*Ylhs_wv1));
    tmp = tmp./L';

    %mus1_lhs = tmp(1:end/2,:); %for using mua and mus' in Jac
    %mua1_lhs = tmp(end/2+1:end,:); %for using mua and mus' in Jac

    mua1_lhs = tmp(:,:); %for just using mu a in Jac

    tmp = J'*((Hess+reg)\(diag(1./M)*Ylhs_wv2));
    tmp = tmp./L'; %Not included before, just added 08:26 30 06 22

    %mus2_lhs = tmp(1:end/2,:);%for using mua and mus' in Jac
    %mua2_lhs = tmp(end/2+1:end,:);%for using mua and mus' in Jac

    mua2_lhs = tmp(:,:); %for just using mu a in Jac
    %
    % unmix chromophores
    for i = 1 : size(mua1_lhs,2)
        tmp = inv(E)*[mua1_lhs(:,i)'; mua2_lhs(:,i)']; %just using mu a

    %tmp = inv(E)*[mus1_lhs(:,i)'; mus2_lhs(:,i)']; %just using kappa

    %tmp = inv(E)*[mua1_lhs(:,i)' mus1_lhs(:,i)'; mua2_lhs(:,i)' mus2_lhs(:,i)']; %Using mua and mus'
    HbO_c(:,i) = tmp(1,:)';
    Hb_c(:,i) = tmp(2,:)';
    HbT_c(:,i) = sum(tmp);
    end
    
end