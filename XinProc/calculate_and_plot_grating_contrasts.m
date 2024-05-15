%%

saving = false;

monkey = 'Louwho';
load('D:\XINTRINSIC\Louwho_20240513d_proc\20240513d140841t_Recording_JoinedCollection_300x480@5fps');
% load('D:\XINTRINSIC\Louwho_20240513d_proc\20240513d143017t_Recording_JoinedCollection_300x480@5fps');
save_path = 'D:\XINTRINSIC\Louwho_20240513d_proc\results';



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
%trial types 'LURD';
    % 'U';  Up
    % 'R';  Right
    % 'D';  Down               
    % 'L';  Left
    
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

    % 'U';  Up
    % 'R';  Right
    % 'D';  Down               
    % 'L';  Left

UDidx = [1 3];
LRidx = [2 4];
R.AmpU = squeeze(R.snapshot(:,1,:,:));
R.AmpR = squeeze(R.snapshot(:,2,:,:));
R.AmpD = squeeze(R.snapshot(:,3,:,:));
R.AmpL = squeeze(R.snapshot(:,4,:,:));
R.AmpUD = reshape(R.snapshot(:,UDidx,:,:), D.R.NumC*2, D.R.NumH, D.R.NumW);
R.AmpLR = reshape(R.snapshot(:,LRidx,:,:), D.R.NumC*2, D.R.NumH, D.R.NumW);

D.R.MeanU = squeeze(mean(R.AmpU, 1));
D.R.MeanR = squeeze(mean(R.AmpR, 1));
D.R.MeanD = squeeze(mean(R.AmpD, 1));
D.R.MeanL = squeeze(mean(R.AmpL, 1));
D.R.MeanUD = squeeze(mean(R.AmpUD, 1));
D.R.MeanLR = squeeze(mean(R.AmpLR, 1));

D.R.StdU = squeeze(std(R.AmpU, 0, 1));
D.R.StdR = squeeze(std(R.AmpR, 0, 1));
D.R.StdD = squeeze(std(R.AmpD, 0, 1));
D.R.StdL = squeeze(std(R.AmpL, 0, 1));
D.R.StdUD = squeeze(std(R.AmpUD, 0, 1));
D.R.StdLR = squeeze(std(R.AmpLR, 0, 1));

D.R.AmpMeanMaxUD = squeeze(max(mean(R.snapshot(:,UDidx,:,:),1), [], 2));
D.R.AmpMeanMaxLR = squeeze(max(mean(R.snapshot(:,LRidx,:,:),1), [], 2));

D.R.AmpMeanDiff = D.R.MeanUD - D.R.MeanLR;
D.R.AmpMeanMaxDiff = D.R.AmpMeanMaxUD - D.R.AmpMeanMaxLR;
% D.R.AmpMeanDiffObjs = D.R.MeanF - D.R.AmpMeanMaxObjs;
% D.R.AmpMeanS = squeeze(mean(R.snapshot(:,8,:,:),1));
% D.R.AmpMeanMaxNonS = squeeze(max(mean(R.snapshot(:,1:7,:,:),1), [], 2)); % all others
% D.R.AmpMeanDiffS = D.R.AmpMeanS - D.R.AmpMeanMaxNonS;

%% plot response/amplitude difference maps
figure;
f2h = imagesc(squeeze(D.R.AmpMeanDiff), [-0.02,0.04]);
axis equal; colorbar; axis off;
title('response mean difference UD against LR');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps_FigSaveAs.png'])
end

figure;
f3h = imagesc(squeeze(D.R.AmpMeanMaxDiff)); %, [-0.000,0.08]);
axis equal; colorbar; axis off;
title('response mean max difference UD against LR');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanMaxDiff / max(D.R.AmpMeanMaxDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanMaxDiff_' cycle_num 'reps.png'])
    saveas(f3h,[save_path filesep datestrf '_' monkey '_AmpMeanMaxDiff_' cycle_num 'reps_FigSaveAs.png'])
end


%% calculate and plot Tvalue maps contrasting UD with LR
D.R.TvalueUDvsLR = (D.R.MeanUD - D.R.MeanLR) ./ sqrt(D.R.StdUD.^2 / (D.R.NumC * 2) + D.R.StdLR.^2 / (D.R.NumC * 2));
figure;
% can also try plotting something like (D.R.TvalueFvsBAVOUPS>=1.5)
fh = imagesc(squeeze(D.R.TvalueUDvsLR),[0,3]);
axis equal; colorbar; axis off;
title('T contrast UDvsLR');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        fn = [datestrf '_' monkey '_Tmap_UDvsLR_' cycle_num 'reps'];
    vals = squeeze(D.R.TvalueUDvsLR / max(D.R.TvalueUDvsLR(:)));
    imwrite(vals, [save_path filesep fn '.png'])
    saveas(fh, [save_path filesep fn '_FigSaveAs.png'])
    % saveas(fh, [save_path filesep fn '_FigSaveAs.svg'])
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
    
    f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.010,0.004]); % through-window
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
% if saving
%     set(gcf,'Units','inches');
%     screenposition = get(gcf,'Position');
%     set(gcf,...
%         'PaperPosition',[0 0 screenposition(3:4)],...
%         'PaperSize',[screenposition(3:4)]);
%     print -dpdf -painters FacePatchesOnAnatomical_Blended
% end



% figure(104);
% bgim = gray2rgb(double(R.windowavg/max(R.windowavg(:))));
% % bgim = gray2rgb(R.windowavg*0.000002);
% fgim = zeros(size(bgim));
% AmpMeanThresh = max(D.R.AmpMeanDiff(:)) - 0.65*(max(D.R.AmpMeanDiff(:)) - min(D.R.AmpMeanDiff(:)));
% % AmpMeanDifftmp = D.R.AmpMeanDiff(D.R.AmpMeanDiff >= AmpMeanThresh);
% %AmpMeanObjThresh = max(D.R.AmpMeanDiffObjs(:)) - 0.40*(max(D.R.AmpMeanDiffObjs(:)) - min(D.R.AmpMeanDiffObjs(:)));
% fgim(:,:,2) = D.R.AmpMeanDiff >= AmpMeanThresh;
% % fgim(:,:,2) = D.R.AmpMeanDiffObjs >= AmpMeanObjThresh;
% DotAvgThresh = max(average_map(:)) - 0.55*(max(average_map(:)) - min(average_map(:)));
% fgim(:,:,3) = average_map >= DotAvgThresh;
% % fgim(:,:,3) = D.R.AmpMeanDiff/max(D.R.AmpMeanDiff(:));
% % fgim(:,:,3) = squeeze(D.R.AmpMeanDiff)*800;
% % fused_average_and_FacesVsAll = imfuse(R.windowavg*0.000002,squeeze(D.R.AmpMeanDiff_tmp)*800,'falsecolor','Scaling','none','ColorChannels',[1 2 0]);
% % blended = imblend(squeeze(D.R.AmpMeanDiff_tmp)*800, bgim, 0.5, 'normal');
% blended = imblend(0.6*fgim, 1.5*bgim, 0.6, 'normal');
% %f1h1 = imagesc(fused_average_and_FacesVsAll); 
% f1h1 = imshow(blended);
% axis equal; axis off;
% title('blended FacesVsMaxAll (b) and DotsMovVsNotResp (g)');
