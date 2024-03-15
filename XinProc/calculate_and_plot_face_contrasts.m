%%

saving = true;

% monkey = 'Cadbury';
% load('D:\XINTRINSIC\Cadbury_20220405d\20220406d004446t_Recording_JoinedCollection_300x480@5fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d';
% load('D:\XINTRINSIC\Cadbury_20220405d_LessBin\20230126d195611t_Recording_JoinedCollection_300x480@10fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d_LessBin';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230129d\20230129d231636t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Crumpet_20230129d';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230828d_Intrinsic\20230828d170205t_Recording_JoinedCollection_300x480@5fps.mat')
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
% load('D:\XINTRINSIC\20230812d120014t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Larry_combined_proc';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images\20230912d081854t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images\results';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images_dup\setA\20230912d201726t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images_dup\setA\results';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230911d_org\images_dup\setB\20230912d202539t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Larry_20230911d_org\images_dup\setB\results';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20230913d_org\images\20230913d215830t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Larry_20230913d_org\images\results';

% monkey = 'Larry';
% load('D:\XINTRINSIC\Larry_20231011d_PDproctest\20231011d172353t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Larry_20231011d_PDproctest\results';

% monkey = 'Curly';
% load('D:\XINTRINSIC\Curly_20231004d_Intrinsic\20231004d151631t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Curly_20231004d\results';

% monkey = 'Curly';
% load('D:\XINTRINSIC\Curly_20231016d\20231016d215420t_Recording_JoinedCollection_300x480@5fps')
% save_path = 'D:\XINTRINSIC\Curly_20231016d\results';

monkey = 'Coconut';
% load('D:\XINTRINSIC\Coconut_20240112d\images_pxBin1_frBin5_gaussFalse\20240112d191105t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Coconut_20240112d\results';
load('D:\XINTRINSIC\Coconut_20240205d_org\images_pxBin1_frBin5_gaussFalse\20240205d172049t_Recording_JoinedCollection_300x480@5fps.mat')
save_path = 'D:\XINTRINSIC\Coconut_20240205d_org\results';
% load('D:\XINTRINSIC\Coconut_20240207d_org\images_pxBin1_frBin5_gaussFalse\20240207d153353t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Coconut_20240207d_org\results';


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

%% plot response/amplitude difference maps
figure;
f2h = imagesc(squeeze(D.R.AmpMeanDiff), [-0.000,0.007]);
axis equal; colorbar; axis off;
title('response difference F against max of all others');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps_FigSaveAs.png'])
end

figure;
f2h = imagesc(squeeze(D.R.AmpMeanDiffObjs), [-0.00,0.015]);
axis equal; colorbar; axis off;
title('response difference F against max of objects only');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiffObjs / max(D.R.AmpMeanDiffObjs(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiffObjs_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiffObjs_' cycle_num 'reps_FigSaveAs.png'])
end

figure;
f2h = imagesc(squeeze(D.R.AmpMeanDiffS), [-0.02,0.010]);
axis equal; colorbar; axis off;
title('response difference S against max of all others');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiffS_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiffS_' cycle_num 'reps_FigSaveAs.png'])
end


%% calculate and plot Tvalue maps contrasting F with other categories
D.R.TvalueFvsBAVOUPS = (D.R.MeanF - D.R.MeanBAVOUPS) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdBAVOUPS.^2 / (D.R.NumC * 7));
figure;
% can also try plotting something like (D.R.TvalueFvsBAVOUPS>=1.5)
fh = imagesc(squeeze(D.R.TvalueFvsBAVOUPS),[0,3]);
axis equal; colorbar; axis off;
title('T contrast FvsBAVOUPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_FvsBAVOUPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsBAVOUPS / max(D.R.TvalueFvsBAVOUPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueFvsBAVOU = (D.R.MeanF - D.R.MeanBAVOU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdBAVOU.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueFvsBAVOU),[-1,2]);
axis equal; colorbar; axis off;
title('T contrast FvsBAVOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_FvsBAVOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsBAVOU / max(D.R.TvalueFvsBAVOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueFvsOU = (D.R.MeanF - D.R.MeanOU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdOU.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueFvsOU),[-2,1.5]);
axis equal; colorbar; axis off;
title('T contrast FvsOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_FvsOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsOU / max(D.R.TvalueFvsOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueFvsO = (D.R.MeanF - D.R.MeanO) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdO.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueFvsO));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast FvsO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_FvsO_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsO / max(D.R.TvalueFvsO(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueFvsU = (D.R.MeanF - D.R.MeanU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdU.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueFvsU));%[-2,2]);
axis equal; colorbar; axis off;
title('T contrast FvsU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_FvsU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsU / max(D.R.TvalueFvsU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueFvsPS = (D.R.MeanF - D.R.MeanPS) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueFvsPS));
axis equal; colorbar; axis off;
title('T contrast FvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_FvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueFvsPS / max(D.R.TvalueFvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end


%% calculate and plot Tvalue maps contrasting B with other categories
D.R.TvalueBvsFAVOUPS = (D.R.MeanB - D.R.MeanFAVOUPS) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdFAVOUPS.^2 / (D.R.NumC * 7));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsFAVOUPS));%,[-2,4]);
axis equal; colorbar; axis off;
title('T contrast BvsFAVOUPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_BvsFAVOUPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsFAVOUPS / max(D.R.TvalueBvsFAVOUPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueBvsFAVOU = (D.R.MeanB - D.R.MeanFAVOU) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdFAVOU.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsFAVOU));%,[-2,2]);
axis equal; colorbar; axis off;
title('T contrast BvsFAVOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_BvsFAVOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsFAVOU / max(D.R.TvalueBvsFAVOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueBvsOU = (D.R.MeanB - D.R.MeanOU) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdOU.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsOU));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast BvsOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_BvsOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsOU / max(D.R.TvalueBvsOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueBvsO = (D.R.MeanB - D.R.MeanO) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdO.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsO));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast BvsO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_BvsO_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsO / max(D.R.TvalueBvsO(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueBvsU = (D.R.MeanB - D.R.MeanU) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdU.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsU));%[-2,2]);
axis equal; colorbar; axis off;
title('T contrast BvsU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_BvsU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsU / max(D.R.TvalueBvsU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueBvsPS = (D.R.MeanB - D.R.MeanPS) ./ sqrt(D.R.StdB.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueBvsPS));
axis equal; colorbar; axis off;
title('T contrast BvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_BvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueBvsPS / max(D.R.TvalueBvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

%% calculate and plot Tvalue maps contrasting V with other categories
D.R.TvalueVvsFBAOUPS = (D.R.MeanV - D.R.MeanFBAOUPS) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdFBAOUPS.^2 / (D.R.NumC * 7));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsFBAOUPS),[-1,2]);
axis equal; colorbar; axis off;
title('T contrast VvsFBAOUPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_VvsFBAOUPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsFBAOUPS / max(D.R.TvalueVvsFBAOUPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueVvsFBAOU = (D.R.MeanV - D.R.MeanFBAOU) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdFBAOU.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsFBAOU));%,[-1,1]);
axis equal; colorbar; axis off;
title('T contrast VvsFBAOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_VvsFBAOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsFBAOU / max(D.R.TvalueVvsFBAOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueVvsOU = (D.R.MeanV - D.R.MeanOU) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdOU.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsOU));%,[-3,1]);
axis equal; colorbar; axis off;
title('T contrast VvsOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_VvsOU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsOU / max(D.R.TvalueVvsOU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueVvsO = (D.R.MeanV - D.R.MeanO) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdO.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsO));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast VvsO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_VvsO_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsO / max(D.R.TvalueVvsO(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueVvsU = (D.R.MeanV - D.R.MeanU) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdU.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsU));%[-2,2]);
axis equal; colorbar; axis off;
title('T contrast VvsU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_VvsU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsU / max(D.R.TvalueVvsU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueVvsPS = (D.R.MeanV - D.R.MeanPS) ./ sqrt(D.R.StdV.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueVvsPS));
axis equal; colorbar; axis off;
title('T contrast VvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_VvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueVvsPS / max(D.R.TvalueVvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end


%% calculate and plot Tvalue maps contrasting O with other categories
D.R.TvalueOvsFBAVUPS = (D.R.MeanO - D.R.MeanFBAVUPS) ./ sqrt(D.R.StdO.^2 / D.R.NumC + D.R.StdFBAVUPS.^2 / (D.R.NumC * 7));
figure;
fh = imagesc(squeeze(D.R.TvalueOvsFBAVUPS));%,[-2,4]);
axis equal; colorbar; axis off;
title('T contrast OvsFBAVUPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_OvsFBAVUPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOvsFBAVUPS / max(D.R.TvalueOvsFBAVUPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueOvsFBAVU = (D.R.MeanO - D.R.MeanFBAVU) ./ sqrt(D.R.StdO.^2 / D.R.NumC + D.R.StdFBAVU.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueOvsFBAVU));%,[-2,2]);
axis equal; colorbar; axis off;
title('T contrast OvsFBAVU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_OvsFBAVU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOvsFBAVU / max(D.R.TvalueOvsFBAVU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueOvsU = (D.R.MeanO - D.R.MeanU) ./ sqrt(D.R.StdO.^2 / D.R.NumC + D.R.StdU.^2 / D.R.NumC);
figure;
fh = imagesc(squeeze(D.R.TvalueOvsU));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast OvsU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_OvsU_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOvsU / max(D.R.TvalueOvsU(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueOvsPS = (D.R.MeanO - D.R.MeanPS) ./ sqrt(D.R.StdO.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueOvsPS));
axis equal; colorbar; axis off;
title('T contrast OvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_OvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOvsPS / max(D.R.TvalueOvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end


%% calculate and plot Tvalue maps contrasting U with other categories
D.R.TvalueUvsFBAVOPS = (D.R.MeanU - D.R.MeanFBAVOPS) ./ sqrt(D.R.StdU.^2 / D.R.NumC + D.R.StdFBAVOPS.^2 / (D.R.NumC * 7));
figure;
fh = imagesc(squeeze(D.R.TvalueUvsFBAVOPS));%,[-2,4]);
axis equal; colorbar; axis off;
title('T contrast UvsFBAVOPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_UvsFBAVOPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueUvsFBAVOPS / max(D.R.TvalueUvsFBAVOPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueUvsFBAVO = (D.R.MeanU - D.R.MeanFBAVO) ./ sqrt(D.R.StdU.^2 / D.R.NumC + D.R.StdFBAVO.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueUvsFBAVO));%,[-2,2]);
axis equal; colorbar; axis off;
title('T contrast UvsFBAVO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_UvsFBAVO_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueUvsFBAVO / max(D.R.TvalueUvsFBAVO(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueUvsO = (D.R.MeanU - D.R.MeanO) ./ sqrt(D.R.StdU.^2 / D.R.NumC + D.R.StdO.^2 / D.R.NumC);
figure;
fh = imagesc(squeeze(D.R.TvalueUvsO));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast UvsO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_UvsO_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueUvsO / max(D.R.TvalueUvsO(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueUvsPS = (D.R.MeanU - D.R.MeanPS) ./ sqrt(D.R.StdU.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueUvsPS));
axis equal; colorbar; axis off;
title('T contrast UvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_UvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueUvsPS / max(D.R.TvalueUvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end


%% calculate and plot Tvalue maps contrasting OU with other categories
D.R.TvalueOUvsFBAVPS = (D.R.MeanOU - D.R.MeanFBAVPS) ./ sqrt(D.R.StdOU.^2 / (D.R.NumC * 2) + D.R.StdFBAVUPS.^2 / (D.R.NumC * 6));
figure;
fh = imagesc(squeeze(D.R.TvalueOUvsFBAVPS));%,[-2,4]);
axis equal; colorbar; axis off;
title('T contrast OUvsFBAVPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_OUvsFBAVPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOUvsFBAVPS / max(D.R.TvalueOUvsFBAVPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueOUvsFBAV = (D.R.MeanOU - D.R.MeanFBAV) ./ sqrt(D.R.StdOU.^2 / (D.R.NumC * 2) + D.R.StdFBAV.^2 / (D.R.NumC * 4));
figure;
fh = imagesc(squeeze(D.R.TvalueOUvsFBAV));%,[-2,2]);
axis equal; colorbar; axis off;
title('T contrast OUvsFBAV');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_OUvsFBAV_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOUvsFBAV / max(D.R.TvalueOUvsFBAV(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end

D.R.TvalueOUvsPS = (D.R.MeanOU - D.R.MeanPS) ./ sqrt(D.R.StdOU.^2 / (D.R.NumC * 2) + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueOUvsPS));
axis equal; colorbar; axis off;
title('T contrast OUvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    fn = [datestrf '_' monkey '_Tmap_OUvsPS_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueOUvsPS / max(D.R.TvalueOUvsPS(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
end


%%

saving = false;

% Calculate and plot each category vs the maximum of all other categories
category_list = 1:8;
category_string = 'FBAVOUPS';
figure(100);

for category_i = 1:8 
    D.R.MeanF_tmp = squeeze(mean(R.snapshot(:,category_i,:,:),1));
    D.R.MeanBAVOUPS_tmp = squeeze(max(mean(R.snapshot(:,category_list (category_list ~= category_i),:,:),1), [], 2)); % all others
    D.R.AmpMeanDiff_tmp = D.R.MeanF_tmp - D.R.MeanBAVOUPS_tmp;
    subplot(2,4,category_i)
    
    f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.020,0.006]); % through-window
    %f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.0015,0.0015]); %through-skull
    title(strcat('category = ', category_string(category_i)))
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
% IM = R.windowavg*0.000002;
% red_IM = cast(cat(3, IM, zeros(size(IM)), zeros(size(IM))), class(IM));
% green_IM = cast(cat(3, zeros(size(IM)), IM, zeros(size(IM))), class(IM));
% blue_IM = cast(cat(3, zeros(size(IM)), zeros(size(IM)), IM), class(IM));
% bgim = uint8(IM * 256);
bgim = gray2rgb(R.windowavg*0.000002);
fgim = zeros(size(bgim));
fgim(:,:,3) = squeeze(D.R.AmpMeanDiff_tmp)*150;
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



figure(104);
bgim = gray2rgb(double(R.windowavg/max(R.windowavg(:))));
% bgim = gray2rgb(R.windowavg*0.000002);
fgim = zeros(size(bgim));
AmpMeanThresh = max(D.R.AmpMeanDiff(:)) - 0.65*(max(D.R.AmpMeanDiff(:)) - min(D.R.AmpMeanDiff(:)));
% AmpMeanDifftmp = D.R.AmpMeanDiff(D.R.AmpMeanDiff >= AmpMeanThresh);
%AmpMeanObjThresh = max(D.R.AmpMeanDiffObjs(:)) - 0.40*(max(D.R.AmpMeanDiffObjs(:)) - min(D.R.AmpMeanDiffObjs(:)));
fgim(:,:,2) = D.R.AmpMeanDiff >= AmpMeanThresh;
% fgim(:,:,2) = D.R.AmpMeanDiffObjs >= AmpMeanObjThresh;
DotAvgThresh = max(average_map(:)) - 0.55*(max(average_map(:)) - min(average_map(:)));
fgim(:,:,3) = average_map >= DotAvgThresh;
% fgim(:,:,3) = D.R.AmpMeanDiff/max(D.R.AmpMeanDiff(:));
% fgim(:,:,3) = squeeze(D.R.AmpMeanDiff)*800;
% fused_average_and_FacesVsAll = imfuse(R.windowavg*0.000002,squeeze(D.R.AmpMeanDiff_tmp)*800,'falsecolor','Scaling','none','ColorChannels',[1 2 0]);
% blended = imblend(squeeze(D.R.AmpMeanDiff_tmp)*800, bgim, 0.5, 'normal');
blended = imblend(0.6*fgim, 1.5*bgim, 0.6, 'normal');
%f1h1 = imagesc(fused_average_and_FacesVsAll); 
f1h1 = imshow(blended);
axis equal; axis off;
title('blended FacesVsMaxAll (b) and DotsMovVsNotResp (g)');

f2h = imagesc(squeeze(D.R.AmpMeanDiff), [-0.0,0.015]);
axis equal; colorbar; axis off;
title('response difference F against max of all others');
% if saving
%     datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
%     imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
%         [save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps.png'])
%     saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps_FigSaveAs.png'])
% end

figure;
f2h = imagesc(squeeze(D.R.AmpMeanDiffObjs), [-0.01,0.015]);
axis equal; colorbar; axis off;
title('response difference F against max of objects only');


