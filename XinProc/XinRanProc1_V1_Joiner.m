function XinRanProc1_FaceTrial_Joiner(varargin)
% Xintrinsic preProcessing 1 
% DATA BINNING

clear global
global A P S

%% Get preprocessed ('*.rec') file
[~, A.Sys.pcname] = system('hostname');
if strcmp(A.Sys.pcname(1:end-1), 'Intrinsic')
        A.Sys.folder = 'D:\XINTRINSIC\';    
else
        A.Sys.folder = 'D:\XINTRINSIC\';     
end

if nargin ==0
    % Calling from direct running of the function
    A.RunningSource =   'D';
    [A.FileName, A.PathName, A.FilterIndex] = uigetfile(...
        [A.Sys.folder '*_Processed*.mat'],...
        'Select processed recording files to join',...
        'MultiSelect',              'On');
    if A.FilterIndex == 0
        clear A;                    % nothing selected
        return
    end
    if iscell(A.FileName) == 0      % single file selected
        A.FileName = {A.FileName};
    end
else
    A.RunningSource =   'S';
    % Calling from another script
    [A.PathName, A.FileName, FileExt] = fileparts(varargin{1});
    A.PathName =        [A.PathName, filesep];
    A.FileName =        {[A.FileName, FileExt]};
end

disp(['Xintrinsic Processing Stage 1 (spatiotemporal binning) is about to start on ' ...
    num2str(length(A.FileName)) ' files']);

%% DATA JOINING
for i = 1:length(A.FileName)
   
    %% Load 'S'
    A.curfilename = [A.PathName, A.FileName{i}];
    [fp,fn,fe] = fileparts(A.curfilename);
    ofn = regexprep(fn,'_Processed_.*','');
    S = load(fullfile(fp, [ofn, '_StimulusInformation', fe]));  
    P = load(A.curfilename);
    S = S.S;
    P = P.P;
    if i==1
        Sall = S;
        Pall = P;
    else
        Sall.TrlDurTotal =      Sall.TrlDurTotal +      S.TrlDurTotal;
        Sall.SesCycleNumTotal =	Sall.SesCycleNumTotal + S.SesCycleNumTotal;
        Sall.SesDurTotal =      Sall.SesDurTotal +      S.SesDurTotal;  
        Sall.SesTrlOrderMat = [ Sall.SesTrlOrderMat;    S.SesTrlOrderMat];
        Sall.SesTrlOrderVec = [ Sall.SesTrlOrderVec,    S.SesTrlOrderVec];
        Sall.SesTrlOrderSoundVec = [	Sall.SesTrlOrderSoundVec, S.SesTrlOrderSoundVec];
        Pall.ProcFrameNumTotal = 	Pall.ProcFrameNumTotal +    P.ProcFrameNumTotal;
        Pall.RawMeanPixel = [       Pall.RawMeanPixel,          P.RawMeanPixel];
        Pall.RawMeanPower = [       Pall.RawMeanPower,          P.RawMeanPower];
        Pall.ProcMeanPixel = [      Pall.ProcMeanPixel,         P.ProcMeanPixel];
        Pall.ProcMeanPower = [      Pall.ProcMeanPower,         P.ProcMeanPower];
        Pall.ProcDataMat = [        Pall.ProcDataMat;           P.ProcDataMat];
    end
    disp([  'Reading: "', A.FileName{i},  '"']);
end
    P = Pall;
    S = Sall;
    nowthen = now;
    [~,fn,~] = fileparts(A.FileName{1});
    extrastr = regexprep(fn,'.*_','');
    A.combinedname = [A.PathName, datestr(nowthen, 'yyyymmdd') 'd' datestr(nowthen, 'HHMMSS'), 't_Recording_JoinedCollection_' extrastr];
    save(A.combinedname, 'P', 'S', '-v7.3');  
    %save([A.combinedname, A.FileName{1}(14:36), '.mat'], 'S', '-v7.3');   

disp('All files are processed');
return;
