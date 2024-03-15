%%

saving = true;

% monkey = 'Dali';
% load('D:\XINTRINSIC\Dali_20230505d_proc\20230505d134732t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Dali_20230505d_proc';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230129d\20230129d231636t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Crumpet_20230129d';

% monkey = 'Crumpet';
% load('D:\XINTRINSIC\Crumpet_20230127d\20230127d180107t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Crumpet_20230127d';

% monkey = 'Cadbury';
% load('D:\XINTRINSIC\Cadbury_20220405d\20220406d004446t_Recording_JoinedCollection_300x480@5fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d';
% load('D:\XINTRINSIC\Cadbury_20220405d_LessBin\20230126d195611t_Recording_JoinedCollection_300x480@10fps_24repeats.mat');
% save_path = 'D:\XINTRINSIC\Cadbury_20220405d_LessBin';

monkey = 'Cashew';
load('D:\XINTRINSIC\Cashew_combined\20220406d093729t_Recording_JoinedCollection_300x480@5fps_18repeats.mat')
save_path = 'D:\XINTRINSIC\Cashew_combined_proc';

% monkey = 'Dali';
% load('D:\XINTRINSIC\Dali_20220523d\20220523d132338t_Recording_JoinedCollection_300x480@5fps.mat')
% save_path = 'D:\XINTRINSIC\Dali_20220523d';

% monkey = 'Hershey';
%load('D:\XINTRINSIC\Hershey_20220405d\20220406d003223t_Recording_JoinedCollection_300x480@5fps_14repeats.mat')
% save_path = 'D:\XINTRINSIC\Hershey_20220405d';

%%
if saving
    if isfolder(save_path)
        cd(save_path);
    else
        mkdir(save_path);
        cd(save_path);
    end
end

cycle_num = num2str(size(P.ProcDataMat,1));

% P.ProcDataMat [cycles x trials x pixels(H) x pixels(W) x frames]

% A.NumCses = 2;
% A.NumCall = A.NumCses * size(A.FileName,2);
[A.NumC, A.NumT, A.NumH, A.NumW, A.NumF] = size(P.ProcDataMat);

% as in Song et al Wang 2022 NatCommun Fig6d...
% blank pre-stim period is 2 sec (0-2 sec)
% stim period is 10 sec (2-12 sec)
% post period is 8 sec (12-20 sec)
% P.ProcFrameRate is generally 5 fps, so...
%   IdxPre = 0-2 sec (matching stm.Vis.TrlDurPreStim = 2)
%   IdxRes = 6-12 sec (probably 6 sec to factor in 4 sec hemodynamic delay)
A.IdxPre =  1:10;
A.IdxRes =  31:60;

% R.snapshot [cycles x trials x height x width]
R.snapshot =    zeros(A.NumC,    A.NumT,        A.NumH, A.NumW);
R.trlmean =     zeros(           A.NumT,        A.NumH, A.NumW, A.NumF);
R.trlstd =      zeros(           A.NumT,        A.NumH, A.NumW, A.NumF);

for i = 1:A.NumT
    clear Rtrltemp
    Rtrltemp = zeros(A.NumC, A.NumH, A.NumW, A.NumF);
        disp([  'Getting trial# ', num2str(i),  ' done']);
    %for j = 1:length(A.FileName)
        %disp([  '  Reading file #', num2str(j), ': "', A.FileName{j},  '"']);
        %Pcurses = load([A.PathName, A.FileName{j}]);
        %Pcurses = Pcurses.P;
        Rtrltemp(:, :, :, :) = -(squeeze(...
            P.ProcDataMat(:,i,:,:,:) ./ ...
            reshape(repmat(mean(P.ProcDataMat(:,i,:,:,A.IdxPre),5), [1 1 1 A.NumF]), ...
                A.NumC, 1, A.NumH, A.NumW, A.NumF) ) - 1);
        %if i == 1
        %    R.snapshot(:, :, :, :) = -(squeeze(...
        %        mean(P.ProcDataMat(:,:,:,:,A.IdxRes),5) ./ ...
        %        mean(P.ProcDataMat(:,:,:,:,A.IdxPre),5)) - 1);
        %end
    %end
    %save(   [   A.PathName, datestr(now, 'yymmddTHHMMSS'),...
    %            A.FileName{1}(14:end-7), '_Trial#', num2str(i),'_R.mat'],...
    %            'Rtrltemp', '-v7.3'); 
    R.trlmean(i,:,:,:) = mean(Rtrltemp, 1);
    R.trlstd(i,:,:,:) =  std(Rtrltemp, 0, 1);
end
%save(   [   A.PathName, datestr(now, 'yymmddTHHMMSS'),...
%            A.FileName{1}(14:end-7), '_all_R.mat'],...
%            'R', '-v7.3');

%R.snapshot(:, :, :, :) = -(squeeze(...
%    mean(P.ProcDataMat(:,:,:,:,A.IdxRes),5) ./ ...
%    mean(P.ProcDataMat(:,:,:,:,A.IdxPre),5)) - 1);
%figure
%for i = 1:8
%    subplot(1,8,i);
%    imagesc(squeeze(mean(R.snapshot(:,i,:,:),1)));
%    caxis(0.5e-2*[-1 1]);
%end

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
%[A.NumC, A.NumT, A.NumH, A.NumW, A.NumF] =	size(P.ProcDataMat);
R.windowavg = squeeze(mean(mean(mean(P.ProcDataMat(:,:,:,:,:),5),2),1));
figure(1);
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
R.snapshot(:, :, :, :) = -(squeeze(...
    mean(P.ProcDataMat(:,:,:,:,A.IdxRes),5) ./ ...
    mean(P.ProcDataMat(:,:,:,:,A.IdxPre),5)) - 1);
D.R.AmpMeanFace = squeeze(mean(R.snapshot(:,1,:,:),1));
D.R.AmpMeanMaxOthers = squeeze(max(mean(R.snapshot(:,2:8,:,:),1), [], 2)); % all others
D.R.AmpMeanMaxObjs = squeeze(max(mean(R.snapshot(:,5:6,:,:),1), [], 2)); % familiar and unfamiliar objects
D.R.AmpMeanDiff = D.R.AmpMeanFace - D.R.AmpMeanMaxOthers;
D.R.AmpMeanDiffObjs = D.R.AmpMeanFace - D.R.AmpMeanMaxObjs;
D.R.AmpMeanS = squeeze(mean(R.snapshot(:,8,:,:),1));
D.R.AmpMeanMaxNonS = squeeze(max(mean(R.snapshot(:,1:7,:,:),1), [], 2)); % all others
D.R.AmpMeanDiffS = D.R.AmpMeanS - D.R.AmpMeanMaxNonS;

figure();
f2h = imagesc(squeeze(D.R.AmpMeanDiff), [-0.02,0.012]);
axis equal; colorbar; axis off;
title('response difference F against max of all others');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiff_' cycle_num 'reps_FigSaveAs.png'])
end

figure();
f2h = imagesc(squeeze(D.R.AmpMeanDiffObjs), [-0.02,0.01]);
axis equal; colorbar; axis off;
title('response difference F against max of objects only');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiffObjs / max(D.R.AmpMeanDiffObjs(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiffObjs_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiffObjs_' cycle_num 'reps_FigSaveAs.png'])
end

figure();
f2h = imagesc(squeeze(D.R.AmpMeanDiffS), [-0.02,0.010]);
axis equal; colorbar; axis off;
title('response difference S against max of all others');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite((D.R.AmpMeanDiff / max(D.R.AmpMeanDiff(:))),...
        [save_path filesep datestrf '_' monkey '_AmpMeanDiffS_' cycle_num 'reps.png'])
    saveas(f2h,[save_path filesep datestrf '_' monkey '_AmpMeanDiffS_' cycle_num 'reps_FigSaveAs.png'])
end

% D.R.AmpMeanMaxOthers2 = squeeze(max(mean(R.snapshot(:,[2:3 5:8],:,:),1), [], 2)); % except fruits & veggies
% D.R.AmpMeanDiff2 = D.R.AmpMeanFace - D.R.AmpMeanMaxOthers2;
[D.R.NumC, ~, D.R.NumH, D.R.NumW] = size(R.snapshot); % Cycles, Height, Width
R.trlsem = R.trlstd / sqrt(D.R.NumC);
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

D.R.TvalueFvsBAVOUPS = (D.R.MeanF - D.R.MeanBAVOUPS) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdBAVOUPS.^2 / (D.R.NumC * 7));
figure(4);
% (D.R.TvalueFvsBAVOUPS>=1.5)
fh = imagesc(squeeze(D.R.TvalueFvsBAVOUPS));%,[-2,4]);
axis equal; colorbar; axis off;
title('T contrast FvsBAVOUPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvsBAVOUPS / max(D.R.TvalueFvsBAVOUPS(:))),...
        [save_path filesep datestrf '_' monkey '_FvsBAVOUPS_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsBAVOUPS_' cycle_num 'reps_FigSaveAs.png'])
end

D.R.TvalueFvBAVOU = (D.R.MeanF - D.R.MeanBAVOU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdBAVOU.^2 / (D.R.NumC * 5));
figure;
fh = imagesc(squeeze(D.R.TvalueFvBAVOU));%,[-2,2]);
axis equal; colorbar; axis off;
title('T contrast FvsBAVOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvBAVOU/max(D.R.TvalueFvBAVOU(:))),...
        [save_path filesep datestrf '_' monkey '_FvsBAVOU_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsBAVOU_' cycle_num 'reps_FigSaveAs.png'])
end

D.R.TvalueFvOU = (D.R.MeanF - D.R.MeanOU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdOU.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueFvOU));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast FvsOU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvOU/max(D.R.TvalueFvOU(:))),...
        [save_path filesep datestrf '_' monkey '_FvsOU_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsOU_' cycle_num 'reps_FigSaveAs.png'])
end

D.R.TvalueFvU = (D.R.MeanF - D.R.MeanU) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdU.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueFvU)),;%[-2,2]);
axis equal; colorbar; axis off;
title('T contrast FvsU');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvU/max(D.R.TvalueFvU(:))),...
        [save_path filesep datestrf '_' monkey '_FvsU_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsU_' cycle_num 'reps_FigSaveAs.png'])
end

D.R.TvalueFvO = (D.R.MeanF - D.R.MeanO) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdO.^2 / (D.R.NumC));
figure;
fh = imagesc(squeeze(D.R.TvalueFvO));%,[-2,3]);
axis equal; colorbar; axis off;
title('T contrast FvsO');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvO/max(D.R.TvalueFvO(:))),...
        [save_path filesep datestrf '_' monkey '_FvsO_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsO_' cycle_num 'reps_FigSaveAs.png'])
end

D.R.TvalueFvPS = (D.R.MeanF - D.R.MeanPS) ./ sqrt(D.R.StdF.^2 / D.R.NumC + D.R.StdPS.^2 / (D.R.NumC * 2));
figure;
fh = imagesc(squeeze(D.R.TvalueFvPS));
axis equal; colorbar; axis off;
title('T contrast FvsPS');
if saving
    datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
    imwrite(squeeze(D.R.TvalueFvPS/max(D.R.TvalueFvPS(:))),...
        [save_path filesep datestrf '_' monkey '_FvsPS_' cycle_num 'reps.png'])
    saveas(fh,[save_path filesep datestrf '_' monkey '_FvsPS_' cycle_num 'reps_FigSaveAs.png'])
end
%D.R.PixelTrlMean =  squeeze(R.trlmean(:,D.R.PxlIdx(1),D.R.PxlIdx(2),:));
%D.R.PixelTrlSem =   squeeze(R.trlsem( :,D.R.PxlIdx(1),D.R.PxlIdx(2),:));

% *** Does not yet respect 'saving' boolean.
% % Calculate and plot each category vs the maximum of all other categories
% category_list = 1:8;
% category_string = 'FBAVOUPS';
% figure(27);
% 
% for category_i = 1:8 
%     D.R.AmpMeanFace_tmp = squeeze(mean(R.snapshot(:,category_i,:,:),1));
%     D.R.AmpMeanMaxOthers_tmp = squeeze(max(mean(R.snapshot(:,category_list (category_list ~= category_i),:,:),1), [], 2)); % all others
%     D.R.AmpMeanDiff_tmp = D.R.AmpMeanFace_tmp - D.R.AmpMeanMaxOthers_tmp;
%     subplot(2,4,category_i)
%     
%     f2h = imagesc(squeeze(D.R.AmpMeanDiff_tmp),[-0.0015,0.0015]);
%     title(strcat('category = ', category_string(category_i)))
%     axis equal; colorbar;
%     axis off;
% end
% 
% %If you want to save the figure as pdf without weird-Matlab-cropping:
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters epsFig
