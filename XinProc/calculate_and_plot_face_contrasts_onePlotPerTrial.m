%%

saving = false;

% monkey = 'Cadbury';
% load('D:\XINTRINSIC\Cadbury_20220405d\20220406d004446t_Recording_JoinedCollection_300x480@5fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d';
% load('D:\XINTRINSIC\Cadbury_20220405d_LessBin\20230126d195611t_Recording_JoinedCollection_300x480@10fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d_LessBin';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230129d\20230129d231636t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Crumpet_20230129d';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230828d\20230828d170205t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Crumpet_20230828d';

% monkey = 'Dali';
% load('D:\XINTRINSIC\Dali_20230505d_proc\20230505d134732t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Dali_20230505d_proc';

% monkey = 'Dali';
% load('D:\XINTRINSIC\Dali_20220523d\20220523d132338t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Dali_20220523d';

% monkey = 'Hershey';
%load('D:\XINTRINSIC\Hershey_20220405d\20220406d003223t_Recording_JoinedCollection_300x480@5fps_14repeats.mat')
% save_path = 'D:\XINTRINSIC\Hershey_20220405d';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230811d\20230812d120014t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Larry_combined_proc';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images\20230912d081854t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images\results';
% 
% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images_dup\setA\20230912d201726t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images_dup\setA\results';
% 
% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images_dup\setB\20230912d202539t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images_dup\setB\results';

monkey = 'Curly';
load('D:\XINTRINSIC\Curly_20231004d\20231004d151631t_Recording_JoinedCollection_300x480@5fps')
save_path = 'D:\XINTRINSIC\Curly_20231004d\results';

%%
if saving
    if isfolder(save_path)
        cd(save_path);
    else
        mkdir(save_path);
        cd(save_path);
    end
end

cycle_num = num2str(size(P.ProcDataMat, 1));

% P.ProcDataMat [cycles x trials x pixels(H) x pixels(W) x frames]

[A.NumC, A.NumT, A.NumH, A.NumW, A.NumF] = size(P.ProcDataMat);

% as in Song et al Wang 2022 NatCommun Fig6d...
% blank pre-stim period is 2 sec (0-2 sec)
% stim period is 10 sec (2-12 sec)
% post period is 8 sec (12-20 sec)
% P.ProcFrameRate is generally 5 fps, so...
%   IdxPre = 0-2 sec (matching stm.Vis.TrlDurPreStim = 2)
%   IdxRes = 6-12 sec (probably 6 sec to factor in 4 sec hemodynamic delay)
A.IdxPre = 1:10;
A.IdxRes = 31:60;

% R.snapshot [cycles x trials x height x width]
R.snapshot = zeros(A.NumC, A.NumT, A.NumH, A.NumW);

for i = 1:A.NumT
    clear Rtrltemp
    Rtrltemp = zeros(A.NumC, A.NumH, A.NumW, A.NumF);
    disp([ 'Getting trial# ', num2str(i), ' done']);
    Rtrltemp(:, :, :, :) = -(squeeze(...
        P.ProcDataMat(:,i,:,:,:) ./ ...
        reshape(repmat(mean(P.ProcDataMat(:,i,:,:,A.IdxPre),5), ...
        [1 1 1 A.NumF]), ...
        A.NumC, 1, A.NumH, A.NumW, A.NumF) ) - 1);
end

%% 
%trial types 'FBAVOUPS';
    % 'F';  Faces
    % 'B';  Body parts
    % 'A';  Animals               
    % 'V';  Fruits & Vegetables	
    % 'O';  Familiar Objects 
    % 'U';  Unfamiliar Objects 
    % 'P';  Phase scrambled Faces
    % 'S';  Spatially scrambled Faces
    
% P.ProcDataMat [cycles x trials x pixels(H) x pixels(W) x frames]
R.windowavg = squeeze(mean(mean(mean(P.ProcDataMat(:,:,:,:,:), 5), 2), 1));
figure;
f1h = imagesc(R.windowavg); 
colormap gray; axis equal; colorbar; axis off;
title('window average');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((R.windowavg/max(R.windowavg(:))),...
        [save_path filesep datestrf '_' monkey '_Window_' cycle_num 'reps_AvgAll.png'])
    saveas(f1h,[save_path filesep datestrf '_' monkey '_Window_' cycle_num 'reps_AvgAll_FigSaveAs.png'])
end

% R.snapshot [cycles x trials x height x width]
R.snapshot(:,:,:,:) = -(squeeze(...
    mean(P.ProcDataMat(:,:,:,:,A.IdxRes), 5) ./ ...
    mean(P.ProcDataMat(:,:,:,:,A.IdxPre), 5)) - 1);

[D.R.NumC, ~, D.R.NumH, D.R.NumW] = size(R.snapshot);

R.AmpF = squeeze(R.snapshot(:,1,:,:));
R.AmpB = squeeze(R.snapshot(:,2,:,:));
R.AmpA = squeeze(R.snapshot(:,3,:,:));
R.AmpV = squeeze(R.snapshot(:,4,:,:));
R.AmpO = squeeze(R.snapshot(:,5,:,:));
R.AmpU = squeeze(R.snapshot(:,6,:,:));
R.AmpP = squeeze(R.snapshot(:,7,:,:));
R.AmpS = squeeze(R.snapshot(:,8,:,:));
R.AmpOU = reshape(R.snapshot(:,5:6,:,:), D.R.NumC*2, D.R.NumH, D.R.NumW);
R.AmpPS = reshape(R.snapshot(:,7:8,:,:), D.R.NumC*2, D.R.NumH, D.R.NumW);
R.AmpBAVOU = reshape(R.snapshot(:,2:6,:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpBAVOUPS = reshape(R.snapshot(:,2:8,:,:), D.R.NumC*7, D.R.NumH, D.R.NumW);
R.AmpFAVOU = reshape(R.snapshot(:,[1 3:6],:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpFAVOUPS = reshape(R.snapshot(:,[1 3:8],:,:), D.R.NumC*7, D.R.NumH, D.R.NumW);
R.AmpFBAOU = reshape(R.snapshot(:,[1:3 5:6],:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpFBAOUPS = reshape(R.snapshot(:,[1:3 5:8],:,:), D.R.NumC*7, D.R.NumH, D.R.NumW);
R.AmpFBAVU = reshape(R.snapshot(:,[1:4 6],:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpFBAVUPS = reshape(R.snapshot(:,[1:4 6:8],:,:), D.R.NumC*7, D.R.NumH, D.R.NumW);
R.AmpFBAVO = reshape(R.snapshot(:,1:5,:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpFBAVOPS = reshape(R.snapshot(:,[1:5 7:8],:,:), D.R.NumC*7, D.R.NumH, D.R.NumW);
R.AmpFBAV = reshape(R.snapshot(:,1:4,:,:), D.R.NumC*4, D.R.NumH, D.R.NumW);
R.AmpFBAVPS = reshape(R.snapshot(:,[1:4 7:8],:,:), D.R.NumC*6, D.R.NumH, D.R.NumW);

D.R.MeanF = squeeze(mean(R.AmpF, 1));
D.R.MeanB = squeeze(mean(R.AmpB, 1));
D.R.MeanA = squeeze(mean(R.AmpA, 1));
D.R.MeanV = squeeze(mean(R.AmpV, 1));
D.R.MeanO = squeeze(mean(R.AmpO, 1));
D.R.MeanU = squeeze(mean(R.AmpU, 1));
D.R.MeanP = squeeze(mean(R.AmpP, 1));
D.R.MeanS = squeeze(mean(R.AmpS, 1));
D.R.MeanOU = squeeze(mean(R.AmpOU, 1));
D.R.MeanPS = squeeze(mean(R.AmpPS, 1));
D.R.MeanBAVOU = squeeze(mean(R.AmpBAVOU, 1));
D.R.MeanBAVOUPS = squeeze(mean(R.AmpBAVOUPS, 1));
D.R.MeanFAVOU = squeeze(mean(R.AmpFAVOU, 1));
D.R.MeanFAVOUPS = squeeze(mean(R.AmpFAVOUPS, 1));
D.R.MeanFBAOU = squeeze(mean(R.AmpFBAOU, 1));
D.R.MeanFBAOUPS = squeeze(mean(R.AmpFBAOUPS, 1));
D.R.MeanFBAVU = squeeze(mean(R.AmpFBAVU, 1));
D.R.MeanFBAVUPS = squeeze(mean(R.AmpFBAVUPS, 1));
D.R.MeanFBAVO = squeeze(mean(R.AmpFBAVO, 1));
D.R.MeanFBAVOPS = squeeze(mean(R.AmpFBAVOPS, 1));
D.R.MeanFBAV = squeeze(mean(R.AmpFBAV, 1));
D.R.MeanFBAVPS = squeeze(mean(R.AmpFBAVPS, 1));

D.R.StdF = squeeze(std(R.AmpF, 0, 1));
D.R.StdB = squeeze(std(R.AmpB, 0, 1));
D.R.StdA = squeeze(std(R.AmpA, 0, 1));
D.R.StdV = squeeze(std(R.AmpV, 0, 1));
D.R.StdO = squeeze(std(R.AmpO, 0, 1));
D.R.StdU = squeeze(std(R.AmpU, 0, 1));
D.R.StdP = squeeze(std(R.AmpP, 0, 1));
D.R.StdS = squeeze(std(R.AmpS, 0, 1));
D.R.StdOU = squeeze(std(R.AmpOU, 0, 1));
D.R.StdPS = squeeze(std(R.AmpPS, 0, 1));
D.R.StdBAVOU = squeeze(std(R.AmpBAVOU, 0, 1));
D.R.StdBAVOUPS = squeeze(std(R.AmpBAVOUPS, 0, 1));
D.R.StdBAVOU = squeeze(std(R.AmpBAVOU, 0, 1));
D.R.StdBAVOUPS = squeeze(std(R.AmpBAVOUPS, 0, 1));
D.R.StdFAVOU = squeeze(std(R.AmpFAVOU, 0, 1));
D.R.StdFAVOUPS = squeeze(std(R.AmpFAVOUPS, 0, 1));
D.R.StdFBAOU = squeeze(std(R.AmpFBAOU, 0, 1));
D.R.StdFBAOUPS = squeeze(std(R.AmpFBAOUPS, 0, 1));
D.R.StdFBAVU = squeeze(std(R.AmpFBAVU, 0, 1));
D.R.StdFBAVUPS = squeeze(std(R.AmpFBAVUPS, 0, 1));
D.R.StdFBAVO = squeeze(std(R.AmpFBAVO, 0, 1));
D.R.StdFBAVOPS = squeeze(std(R.AmpFBAVOPS, 0, 1));
D.R.StdFBAV = squeeze(std(R.AmpFBAV, 0, 1));
D.R.StdFBAVPS = squeeze(std(R.AmpFBAVPS, 0, 1));

D.R.AmpMeanMaxObjs = squeeze(max(mean(R.snapshot(:,5:6,:,:),1), [], 2)); % familiar and unfamiliar objects
D.R.AmpMeanDiff = D.R.MeanF - D.R.MeanBAVOUPS;
D.R.AmpMeanDiffObjs = D.R.MeanF - D.R.AmpMeanMaxObjs;
D.R.AmpMeanS = squeeze(mean(R.snapshot(:,8,:,:),1));
D.R.AmpMeanMaxNonS = squeeze(max(mean(R.snapshot(:,1:7,:,:),1), [], 2)); % all others
D.R.AmpMeanDiffS = D.R.AmpMeanS - D.R.AmpMeanMaxNonS;


%%
close all
saving = false;

% Calculate and plot each category vs the maximum of all other categories
category_list = 1:8;
category_string = 'FBAVOUPS';

[trial_n, ~, ~, ~] = size(R.snapshot);
for trial_i = 1:trial_n
    figure(100+trial_i)
    for category_i = 1:8 
        D.R.MeanF_tmp = squeeze(R.snapshot(trial_i,category_i,:,:));
        D.R.MeanBAVOUPS_tmp = squeeze(max(R.snapshot(trial_i,category_list (category_list ~= category_i),:,:), [], 2)); % all others
        D.R.AmpMeanDiff_tmp = D.R.MeanF_tmp - D.R.MeanBAVOUPS_tmp;
        subplot(2,4,category_i)

        %f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.050,0.02]); % through-window
        f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.0015,0.0015]); %through-skull
        title(strcat('category = ', category_string(category_i)))
        axis equal; colorbar;
        axis off;
    end
end
%If you want to save the figure as pdf without weird-Matlab-cropping:
if saving
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters AllVsMax(All)
end

figure(101);

for category_i = 1:8 
    D.R.MeanF_tmp = imgaussfilt(squeeze(mean(R.snapshot(:,category_i,:,:),1)),2);
    D.R.MeanBAVOUPS_tmp = imgaussfilt(squeeze(max(mean(R.snapshot(:,category_list(category_list ~= category_i),:,:),1), [], 2)),2); % all others
    D.R.AmpMeanDiff_tmp = D.R.MeanF_tmp - D.R.MeanBAVOUPS_tmp;
    subplot(2,4,category_i)
    
    %f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.008,0.008]); % through-window
    f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.0015,0.0015]); %through-skull
    title(strcat('(blurred) category = ', category_string(category_i)))
    axis equal; colorbar;
    axis off;
end

%If you want to save the figure as pdf without weird-Matlab-cropping:
if saving
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters AllVsMax(All)__Blurred
end

figure(102);

for category_i = 1 
    D.R.MeanF_tmp = imgaussfilt(squeeze(mean(R.snapshot(:,category_i,:,:),1)),2);
    D.R.MeanBAVOUPS_tmp = imgaussfilt(squeeze(max(mean(R.snapshot(:,category_list (category_list ~= category_i),:,:),1), [], 2)),2); % all others
    D.R.AmpMeanDiff_tmp = D.R.MeanF_tmp - D.R.MeanBAVOUPS_tmp;
    
    %f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.008,0.008]); % through-window
    f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.0015,0.0015]); %through-skull
    title('FacesVsMax(All)__Blurred')
    axis equal; colorbar;
    axis off;
end

%If you want to save the figure as pdf without weird-Matlab-cropping:
if saving
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters FacesVsMax(All)__Blurred
end
% figure(103);
% fused_average_and_FacesVsAll = imfuse(R.windowavg*0.000002,squeeze(D.R.AmpMeanDiff_tmp)*800,'falsecolor','Scaling','none','ColorChannels',[1 2 0]);
% f1h1 = imagesc(fused_average_and_FacesVsAll); 
% axis equal; axis off;
% title('Window average fused with FacesVsMax(All)');
% 
% %If you want to save the figure as pdf without weird-Matlab-cropping:
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters FacePatchesOnAnatomical_Blended


figure(103);
bgim = gray2rgb(R.windowavg*0.000002);
fgim = zeros(size(bgim));
fgim(:,:,3) = squeeze(D.R.AmpMeanDiff_tmp)*800;
% fused_average_and_FacesVsAll = imfuse(R.windowavg*0.000002,squeeze(D.R.AmpMeanDiff_tmp)*800,'falsecolor','Scaling','none','ColorChannels',[1 2 0]);
% blended = imblend(squeeze(D.R.AmpMeanDiff_tmp)*800, bgim, 0.5, 'normal');
blended = imblend(fgim, bgim, 0.6, 'normal');
%f1h1 = imagesc(fused_average_and_FacesVsAll); 
f1h1 = imshow(blended);
axis equal; axis off;
title('Window average blended with FacesVsMax(All)');

%If you want to save the figure as pdf without weird-Matlab-cropping:
if saving
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print -dpdf -painters FacePatchesOnAnatomical_Blended
end