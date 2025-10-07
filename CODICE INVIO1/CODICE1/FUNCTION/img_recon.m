function [img] = img_recon(lambda,JCropped,SD,dodFilt,headVolumeMesh,gmSurfaceMesh,vol2gm,save_r)
% Perform the images reconstraction for both HbO and HbR

% INPUT  = lambda         --> Tikhonov regularization 
%          JCropped       --> J matrices for both wavelength only good channels
%          SD             --> NIRS data info structure
%          dodFilt        --> ALL channels definitive dod signal after BP filtering
%          headVolumeMesh --> Age matched volumetric head mesh
%          gmSurfaceMesh  --> Age matched GM surface mesh
%          vol2gm         --> Age matched transfer matrix (volumetric to GM)
%          save_r         --> 1 for saving img results

% OUTPUT = img --> img.hbo and img.hbr images reconstruct over GM

%%  INVERSE OF JACOBIAN

% Initialization of usefull variable
invJ = cell(length(SD.Lambda),1);

% For each Jacobian
for i = 1:length(SD.Lambda)
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    % Single Value Decomposition
    S=svd(JJT);
    % Theoretical formula fo J inverted
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda*max(S)));
end

%% OPTICAL PROPERTIES

% Initialize matrices and load useful data
nNodeGM  = size(gmSurfaceMesh.node,1);  % The node count of the GM mesh
nNodeVol = size(headVolumeMesh.node,1); % The node count of the volume mesh
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(SD.Lambda); % n of wavelengths

% Obtain specific absorption coefficients serve in MBLL
Eall = [];
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    % La formula estrae come primi HbO e HbR poi estrae i valori dei coef.
    % di assorbimento specifici di altri tessuti a noi inutili. Prendiamo
    % solo HbO e HbR
    % 'GetExtinctions' extracts optical property for each tissue of the
    % brain. We are interested only in HbO and HbR.
    Etmp = Etmp(1:2); %HbO and HbR only
    % Conversion of units of measurement 
    % [Mol*L/cm] --> [uMol*L/mm] 
    % Just dived for 1e7
    Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
end

%% RECONSTRUCTION FOR EACH FRAME

% Obtain changes in HbO and HbR starting from changes of optical density

% Change of variable name
datarecon = dodFilt;
% Number of samples to reconstruct 
nFrames  = size(datarecon,1); 

% Initialize final results matrices (HbO and HbR for GM surface mesh)
hbo.gm = zeros(nFrames,nNodeGM);
hbr.gm = zeros(nFrames,nNodeGM);

% For each frame
disp('START RECONSTRUCION')
disp(' ')
disp('Progression:')
for frame = 1:nFrames
    if mod(frame,100) == 0
        disp(num2str(frame))
    end

    % Reconstruct absorption changes
    muaImageAll = zeros(nWavs,nNodeVol);
    for wav = 1:nWavs
        dataTmp = datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1);
        invJtmp = invJ{wav};
        tmp = invJtmp * dataTmp';
        muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
    end

    % Convert to concentration changes (Inverted MBLL)
    hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
    hbr_tmpVol = (muaImageAll(2,:)-Eall(2,1)*hbo_tmpVol)/Eall(2,2);

    % Map to GM surface mesh.
    hbo_tmpGM = (vol2gm*hbo_tmpVol');                                           
    hbr_tmpGM = (vol2gm*hbr_tmpVol');                                           

    % Book-keeping and saving
    hbo.gm(frame,:) = hbo_tmpGM;
    hbr.gm(frame,:) = hbr_tmpGM;

end

disp('RECONSTRUCTION DONE')

%% SAVE OF THE RESULTS

% Function output
img.hbo = hbo;
img.hbr = hbr;

% Save of the results
if save_r == 1
    save img_recon.mat img
    disp('IMAGES SAVED')
end