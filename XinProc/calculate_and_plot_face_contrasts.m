%load('D:\XINTRINSIC\Cashew_combined\20220406d093729t_Recording_JoinedCollection_300x480@5fps_18repeats.mat')
%load('D:\XINTRINSIC\Hershey_20220405d\20220406d003223t_Recording_JoinedCollection_300x480@5fps_14repeats.mat')
load('D:\XINTRINSIC\Cadbury_20220405d\20220406d004446t_Recording_JoinedCollection_300x480@5fps_24repeats.mat')

% P.ProcDataMat [cycles x trials x pixels(H) x pixels(W) x frames]

% A.NumCses =     2;
% A.NumCall =     A.NumCses*size(A.FileName,2);
[A.NumC, A.NumT, A.NumH, A.NumW, A.NumF] =	size(P.ProcDataMat);
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

%trial types 'FBAVOUPS';
    % 'F';  Faces
    % 'B';  Body parts
    % 'A';  Animals               
    % 'V';  Fruits & Vegetables	
    % 'O';  Familiar Objects 
    % 'U';  Unfamiliar Objects 
    % 'P';  Phase scrambled Faces
    % 'S';  Spatially scrambled Faces
    
    % current comparisons include
    % FvsBAVOU and FvsPS
    
% P.ProcDataMat [cycles x trials x pixels(H) x pixels(W) x frames]
%[A.NumC, A.NumT, A.NumH, A.NumW, A.NumF] =	size(P.ProcDataMat);
R.windowavg = squeeze(mean(mean(mean(P.ProcDataMat(:,:,:,:,:),5),2),1));
figure(1);
imagesc(R.windowavg); colormap gray;
axis equal; colorbar; axis off;
%imwrite((R.windowavg/max(R.windowavg(:))),...
%   'D:\XINTRINSIC\Cashew_combined\20220406d095500t_Cashew_Window_AvgAll.png')
%imwrite((R.windowavg/max(R.windowavg(:))),...
%    'D:\XINTRINSIC\Hershey_20220405d\20220406d101000t_Hershey_Window_AvgAll.png')
imwrite((R.windowavg/max(R.windowavg(:))),...
  'D:\XINTRINSIC\Cadbury_20220405d\20220406d105000t_Cadbury_Window_AvgAll.png')

% R.snapshot [cycles x trials x height x width]
R.snapshot(:, :, :, :) = -(squeeze(...
    mean(P.ProcDataMat(:,:,:,:,A.IdxRes),5) ./ ...
    mean(P.ProcDataMat(:,:,:,:,A.IdxPre),5)) - 1);
D.R.AmpMeanFace = squeeze(mean(R.snapshot(:,1,:,:),1));
D.R.AmpMeanMaxOthers = squeeze(max(mean(R.snapshot(:,2:end,:,:),1), [], 2)); % all others
D.R.AmpMeanDiff = D.R.AmpMeanFace - D.R.AmpMeanMaxOthers;
figure(2);
imagesc(squeeze(D.R.AmpMeanDiff));
axis equal; axis off;
% D.R.AmpMeanMaxOthers2 = squeeze(max( mean(R.snapshot(:,[2:3 5:8],:,:),1), [], 2));  % except fruits & veggies
% D.R.AmpMeanDiff2 =      D.R.AmpMeanFace - D.R.AmpMeanMaxOthers2;
[D.R.NumC, ~, D.R.NumH, D.R.NumW] = size(R.snapshot);   % Cycles, Height, Width
R.trlsem = R.trlstd / sqrt(D.R.NumC);
R.AmpFace = squeeze(R.snapshot(:,1,:,:));
R.AmpObjects = reshape(R.snapshot(:,2:6,:,:), D.R.NumC*5, D.R.NumH, D.R.NumW);
R.AmpScmbld = reshape(R.snapshot(:,7:8,:,:), D.R.NumC*2, D.R.NumH, D.R.NumW);
D.R.MeanFace = squeeze(mean(R.AmpFace, 1));
D.R.MeanObjects = squeeze(mean(R.AmpObjects, 1));
D.R.MeanScmbld = squeeze(mean(R.AmpScmbld, 1));
D.R.StdFace = squeeze(std(R.AmpFace, 0, 1));
D.R.StdObjects = squeeze(std(R.AmpObjects, 0, 1));
D.R.StdScmbld = squeeze(std(R.AmpScmbld, 0, 1));
D.R.TvalueFO = (D.R.MeanFace - D.R.MeanObjects) ./ sqrt(D.R.StdFace.^2 / D.R.NumC + D.R.StdObjects.^2 / (D.R.NumC * 5));
figure(3);
imagesc(squeeze(D.R.TvalueFO));
axis equal; colorbar; axis off;
D.R.TvalueFS = (D.R.MeanFace - D.R.MeanScmbld) ./ sqrt(D.R.StdFace.^2 / D.R.NumC + D.R.StdScmbld.^2 / (D.R.NumC * 2));
figure(4);
imagesc(squeeze(D.R.TvalueFS));
axis equal; colorbar; axis off;
%D.R.PixelTrlMean =  squeeze(R.trlmean(:,D.R.PxlIdx(1),D.R.PxlIdx(2),:));
%D.R.PixelTrlSem =   squeeze(R.trlsem( :,D.R.PxlIdx(1),D.R.PxlIdx(2),:));