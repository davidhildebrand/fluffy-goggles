
clearvars stm sys
global stm sys

ScanImagePath = 'E:\FreiwaldSync\MarmoScope\ScanImage\SI-Basic_2021.1.0_(2021-12-22)_2af5d7cfec'
addpath(genpath(ScanImagePath))

sys.TempStimFile = ['D:\XINTRINSIC\', 'TempStimData.mat'];
if isfile(sys.TempStimFile)
    load(sys.TempStimFile, 'ExpDataDir', 'StimDir', 'CycleNum', 'CycleDur');
    stm.DataDir = ExpDataDir;
    stm.StimDir = StimDir;
    stm.CycleNum = CycleNum;
    stm.CycleDur = CycleDur;
    clear ExpDataDir StimDir CycleNum CycleDur;
else
    stm.DataDir = 'D:\XINTRINSIC\TESTING_20220402d\';
    stm.StimDir = 'E:\FreiwaldSync\XINTRINSIC\Stimuli\';
end


%% Switch multi-display mode
if max(Screen('Screens')) ~= 2
    opts = struct(  'WindowStyle',  'modal',... 
                    'Interpreter',  'tex');
    errordlg(...
        [   '\fontsize{20} The monitors are not in extended mode \newline ',...
            'Close this MATLAB instance \newline',...
            'Extend all screens in windows \newline' ,...
            'And re-open MATLAB'],...
        '', opts)
    return
end
dos('C:\Windows\System32\DisplaySwitch.exe /extend');
sca;       
pause(0.5);

%% Specify session parameters
stm.Vis.SesTime = now;
stm.Vis.SesTimeStr = [datestr(stm.Vis.SesTime, 'yyyymmdd'), ...
    'd', datestr(stm.Vis.SesTime, 'HHMMSS'), 't'];
stm.Vis.Name = 'FacePatchCycle';
% locate the screen number or default
sys.screenNumber = max(Screen('Screens'));
for i = 1: max(Screen('Screens'))
    info = Screen('Resolution', i);
    if info.hz == 144
        sys.screenNumber = i;
        break
    end
end

stm.SR =                100e3;
stm.SesDurTotal =       stm.CycleNum * stm.CycleDur; %60;  
stm.Vis.CycleDurTotal = stm.CycleDur; %20; % in second
stm.Vis.PicSource = '\Last\shifthalf';
stm.Vis.PicBackground = 'gray';
stm.Vis.PicDur =        0.5;

%   F:Face;	B:Body;     O:Object;
%   P:Phase scambled;   S:Spatial scambled

% stm.Vis.SesOptionContrast = 'AvF';  % Animals               vs Faces
% stm.Vis.SesOptionContrast = 'VvF';  % Fruits & Vegetables	vs Faces
% stm.Vis.SesOptionContrast = 'OvF';  % Familiar Objects  	vs Faces
 stm.Vis.SesOptionContrast = 'UvF';  % Unfamiliar Objects  	vs Faces
% stm.Vis.SesOptionContrast = 'BvF';	% Body Parts            vs Faces
% stm.Vis.SesOptionContrast = 'PvF';  % Phase SCRBD Faces     vs Faces 
% stm.Vis.SesOptionContrast = 'SvF';  % Spatial SCRBD Faces   vs Faces 

% stm.Vis.SesOptionContrast = 'BvF';	% Body Parts            vs Faces
% stm.Vis.SesOptionContrast = 'OvF';  % Familiar Objects      vs Faces
% stm.Vis.SesOptionContrast = 'PvF';  % Phase SCRBD Faces     vs Faces 
% stm.Vis.SesOptionContrast = 'SvF';  % Spatial SCRBD Faces   vs Faces 
% stm.Vis.SesOptionContrast = 'OvB';  % Familiar Object       vs Body

%% Visual CO Parameters
% Session Timer 

stm.Vis.CycleNumTotal =         stm.SesDurTotal / stm.Vis.CycleDurTotal;     % rep # total
stm.Vis.CycleNumCurrent =       0;
stm.Vis.CycleDurCurrentTimer =	tic;
stm.Vis.CyclePicNum =           round((stm.Vis.CycleDurTotal / 2) / stm.Vis.PicDur);
stm.Vis.CyclePicNumHf =         round(stm.Vis.CyclePicNum / 2);

%% Prepare the Psychtoolbox window
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% Define shades
white = WhiteIndex(sys.screenNumber);
black = BlackIndex(sys.screenNumber);   
gray =  GrayIndex(sys.screenNumber, 0.5);
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
% Open an on screen window
switch stm.Vis.PicBackground
    case 'white'
        bgcolor = white;
    case 'black'
        bgcolor = black;
    case 'gray'
        bgcolor = gray;
end
[stm.Vis.windowPtr, windowRect] = PsychImaging('OpenWindow', sys.screenNumber, bgcolor);

% Check frame duration
stm.Vis.TrialIFI = Screen('GetFlipInterval', stm.Vis.windowPtr);
if (abs(stm.Vis.TrialIFI - 1/144) / (1/144)) > 0.05 % 1/59)/(1/59) > 0.05 
    errordlg('screen is not at right refresh rate!');
    return;
end
vbl = Screen('Flip', stm.Vis.windowPtr);

%% Compare stimulus size in Wang and Freiwald labs and adjust accordingly

%Wang
stm.Vis.MonitorDistance =   75;             % in cm
stm.Vis.MonitorHeight =     0.02724 * 1440;	% in cm
stm.Vis.MonitorWidth =      0.02724 * 2560;	% in cm
stm.Vis.MonitorPixelNumX =  2560;
stm.Vis.MonitorPixelNumY =  1440;

stm.Vis.MonitorAngleX =         2 * atan((stm.Vis.MonitorWidth / 2) / stm.Vis.MonitorDistance)/pi * 180;  
stm.Vis.MonitorAngleY =         2 * atan((stm.Vis.MonitorHeight / 2) / stm.Vis.MonitorDistance)/pi * 180;
stm.Vis.MonitorPixelAngleX =    stm.Vis.MonitorAngleX / stm.Vis.MonitorPixelNumX;
stm.Vis.MonitorPixelAngleY =    stm.Vis.MonitorAngleY / stm.Vis.MonitorPixelNumY;

stm.Vis.MonitorPixelAngle_Wang = mean([stm.Vis.MonitorPixelAngleX stm.Vis.MonitorPixelAngleY]);

%Freiwald
stm.Vis.MonitorDistance =   20;         % in cm
stm.Vis.MonitorHeight =     193.5 / 10; % in cm  % or 0.026875 cm/px
stm.Vis.MonitorWidth =      344 / 10;   % in cm  % or 0.026875 cm/px
[stm.Vis.MonitorPixelNumX, stm.Vis.MonitorPixelNumY] = ...
                                                Screen('WindowSize', stm.Vis.windowPtr);
                                           
stm.Vis.MonitorAngleX =         2 * atan((stm.Vis.MonitorWidth/2) / stm.Vis.MonitorDistance)/pi * 180;  
stm.Vis.MonitorAngleY =         2 * atan((stm.Vis.MonitorHeight/2) / stm.Vis.MonitorDistance)/pi * 180;
stm.Vis.MonitorPixelAngleX =    stm.Vis.MonitorAngleX / stm.Vis.MonitorPixelNumX;
stm.Vis.MonitorPixelAngleY =    stm.Vis.MonitorAngleY / stm.Vis.MonitorPixelNumY;

stm.Vis.MonitorPixelAngle_Freiwald = mean([stm.Vis.MonitorPixelAngleX stm.Vis.MonitorPixelAngleY]);


%% Initialize parameters

stm.Vis.MonitorPixelAngle =     stm.Vis.MonitorPixelAngle_Freiwald;
stm.Vis.MonitorCenter =         [stm.Vis.MonitorPixelNumX/2 stm.Vis.MonitorPixelNumY/2];

% texture patch display size
stm.Vis.MonitorTexPatchPos = round([1280 -720 -1280 720] * ...
    stm.Vis.MonitorPixelAngle_Wang / stm.Vis.MonitorPixelAngle_Freiwald); %SOC.%[1280 -720 -1280 720]; 
                                        
% read in texture patches
stm.Vis.TexImDir = [stm.StimDir, filesep, 'Visual', filesep, stm.Vis.PicSource, filesep];
for i = 1:20
    if strcmp(stm.Vis.PicSource(1:5), '\P002') || strcmp(stm.Vis.PicSource(1:5), '\Last')   
        stm.Vis.TexImFaceOri{i} =	imread([stm.Vis.TexImDir 'SHINEd_m'  num2str(i) '.tif']);
        stm.Vis.TexImAnimals{i} =	imread([stm.Vis.TexImDir 'SHINEd_a'  num2str(i) '.tif']);
        stm.Vis.TexImFruitVe{i} =	imread([stm.Vis.TexImDir 'SHINEd_f'  num2str(i) '.tif']);
        stm.Vis.TexImObjFami{i} =	imread([stm.Vis.TexImDir 'SHINEd_o'  num2str(i) '.tif']);
        stm.Vis.TexImObjUnfm{i} =	imread([stm.Vis.TexImDir 'SHINEd_u'  num2str(i) '.tif']);
        stm.Vis.TexImBodyPrt{i} =	imread([stm.Vis.TexImDir 'SHINEd_b'  num2str(i) '.tif']);
        stm.Vis.TexImFacePhs{i} =	imread([stm.Vis.TexImDir 'SHINEd_p'  num2str(i) '.tif']);
        stm.Vis.TexImFaceSpt{i} =	imread([stm.Vis.TexImDir 'SHINEd_s'  num2str(i) '.tif']);
        stm.Vis.TexInd(i,1) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFaceOri{i}); % 1: Face, Frontal
        stm.Vis.TexInd(i,2) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImAnimals{i}); % 2: Animals
        stm.Vis.TexInd(i,3) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFruitVe{i}); % 3: Fruits & Vegetables
        stm.Vis.TexInd(i,4) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImObjFami{i}); % 4: Objects, Familiar
        stm.Vis.TexInd(i,5) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImObjUnfm{i}); % 5: Objects, Unfamiliar 
        stm.Vis.TexInd(i,6) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImBodyPrt{i}); % 6: Body Parts  
        stm.Vis.TexInd(i,7) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFacePhs{i}); % Phase Scrambled faces
        stm.Vis.TexInd(i,8) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFaceSpt{i}); % Space Scrabled faces
%     else
%         stm.Vis.TexImFaceOri{i} =	imread([stm.Vis.TexImDir 'm'  num2str(i) '.png']);
%         stm.Vis.TexImFacePhs{i} =	imread([stm.Vis.TexImDir 'mp' num2str(i) '.png']);
%         stm.Vis.TexImFaceSpt{i} =	imread([stm.Vis.TexImDir 'ms' num2str(i) '.png']);
%         stm.Vis.TexImBodyPrt{i} =	imread([stm.Vis.TexImDir 'b'  num2str(i) '.png']);
%         stm.Vis.TexImObjects{i} =	imread([stm.Vis.TexImDir 'h'  num2str(i) '.png']);
%         stm.Vis.TexInd(i,1) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFaceOri{i}); % F:1
%         stm.Vis.TexInd(i,2) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFacePhs{i}); % P:2
%         stm.Vis.TexInd(i,3) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImFaceSpt{i}); % S:3
%         stm.Vis.TexInd(i,4) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImBodyPrt{i}); % B:4
%         stm.Vis.TexInd(i,5) =	Screen('MakeTexture', stm.Vis.windowPtr, stm.Vis.TexImObjects{i}); % O:5
    end
end
    % pre-arrange image sequence
    stm.Vis.TexIdxCurrent =     0;  
	stm.Vis.TexIdxAll =         [];
    for i = 1:stm.Vis.CycleNumTotal
        for j = 1:size(stm.Vis.TexInd, 2)
            stm.Vis.TexIdxAll(i, j, :) =    randperm(20, stm.Vis.CyclePicNum);
        end
    end
    if strcmp(stm.Vis.PicSource(1:5), '\P002') || strcmp(stm.Vis.PicSource(1:5), '\Last')          
        switch stm.Vis.SesOptionContrast
            case 'AvF';	stm.Vis.TexSeq = [ 2*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 2*ones(1,stm.Vis.CyclePicNumHf)];
            case 'VvF';	stm.Vis.TexSeq = [ 3*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 3*ones(1,stm.Vis.CyclePicNumHf)];
            case 'OvF';	stm.Vis.TexSeq = [ 4*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 4*ones(1,stm.Vis.CyclePicNumHf)];
            case 'UvF';	stm.Vis.TexSeq = [ 5*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 5*ones(1,stm.Vis.CyclePicNumHf)];
            case 'BvF';	stm.Vis.TexSeq = [ 6*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 6*ones(1,stm.Vis.CyclePicNumHf)];
            case 'PvF';	stm.Vis.TexSeq = [ 7*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 7*ones(1,stm.Vis.CyclePicNumHf)];
            case 'SvF'; stm.Vis.TexSeq = [ 8*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 8*ones(1,stm.Vis.CyclePicNumHf)];
     
            otherwise
        end       
%     else
%         stm.Vis.TexIdxAll =     stm.Vis.TexIdxAll(:,[1 1 1 2 3],:);
%         switch stm.Vis.SesOptionContrast
%             case 'PvF';	stm.Vis.TexSeq = [ 2*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 2*ones(1,stm.Vis.CyclePicNumHf)];
%             case 'SvF';	disp('hi'); stm.Vis.TexSeq = [ 3*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 3*ones(1,stm.Vis.CyclePicNumHf)];
%             case 'BvF';	stm.Vis.TexSeq = [ 4*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 4*ones(1,stm.Vis.CyclePicNumHf)];
%             case 'OvF';	stm.Vis.TexSeq = [ 5*ones(1,stm.Vis.CyclePicNumHf) 1*ones(1,stm.Vis.CyclePicNum) 5*ones(1,stm.Vis.CyclePicNumHf)];
%             case 'OvB';	stm.Vis.TexSeq = [ 5*ones(1,stm.Vis.CyclePicNumHf) 4*ones(1,stm.Vis.CyclePicNum) 5*ones(1,stm.Vis.CyclePicNumHf)];
%             otherwise
%         end
    end
    stm.Vis.TexSeq = [  stm.Vis.TexSeq; 
                        [1:stm.Vis.CyclePicNumHf    1:stm.Vis.CyclePicNum    1:stm.Vis.CyclePicNumHf]   ];

%% NI-DAQ

sys.NIDAQ.D.InTimebaseRate =    stm.SR;
sys.NIDAQ.D.CycleSmplHigh =     2;
sys.NIDAQ.D.CycleSmplLow =      sys.NIDAQ.D.InTimebaseRate * stm.Vis.CycleDurTotal - ...
                                sys.NIDAQ.D.CycleSmplHigh;
import dabs.ni.daqmx.*
sys.NIDAQ.TaskCO = Task('Recording Session Cycle Switcher');
sys.NIDAQ.TaskCO.createCOPulseChanTicks(...
    'Intrinsic_PCIe6323', 1, 'Cycle Counter', '100kHzTimebase', ...
    sys.NIDAQ.D.CycleSmplLow, sys.NIDAQ.D.CycleSmplHigh,...
    0, 'DAQmx_Val_Low');
sys.NIDAQ.TaskCO.cfgImplicitTiming(...
    'DAQmx_Val_FiniteSamps',	stm.Vis.CycleNumTotal+1);
sys.NIDAQ.TaskCO.cfgDigEdgeStartTrig(...
    'RTSI6',            'DAQmx_Val_Rising');
sys.NIDAQ.TaskCO.registerSignalEvent(...
    @XinStimEx_VisSom_Localizer_Callback, 'DAQmx_Val_CounterOutputEvent');
sys.NIDAQ.TaskCO.start()
stm.Vis.Running =               1;     
    opts = struct(  'WindowStyle',  'non-modal',... 
                    'Interpreter',  'tex');
sys.MsgBox =	msgbox('\fontsize{20} Click to terminate the session after current visual cycle','', opts);

%% Play 
while stm.Vis.Running
    
    % Session timing
    stm.Vis.SesOn =           (   stm.Vis.CycleNumCurrent>0 && ...
                                    stm.Vis.CycleNumCurrent<=stm.Vis.CycleNumTotal);
    stm.Vis.CycleDurTotalCurrent =    toc(stm.Vis.CycleDurCurrentTimer); 
    
    %% Display 
    if stm.Vis.SesOn
        a = stm.Vis.TexIdxCurrent;
        stm.Vis.TexIdxCurrent =	...
            min(ceil(stm.Vis.CycleDurTotalCurrent/stm.Vis.PicDur), length(stm.Vis.TexSeq));
        if stm.Vis.TexIdxCurrent ~= a
            fprintf('.'); %disp([sprintf('frame time:   '), datestr(now, 'HH:MM:SS.FFF')]);
            Screen('DrawTextures', stm.Vis.windowPtr,...
                stm.Vis.TexInd(...
                    stm.Vis.TexIdxAll(...
                        stm.Vis.CycleNumCurrent,...
                        stm.Vis.TexSeq(1, stm.Vis.TexIdxCurrent),...
                        stm.Vis.TexSeq(2, stm.Vis.TexIdxCurrent) ),...
                    stm.Vis.TexSeq(1, stm.Vis.TexIdxCurrent)        ),...
                [], (stm.Vis.MonitorCenter([1 2 1 2]) + stm.Vis.MonitorTexPatchPos)');
%             Screen('DrawingFinished', stm.Vis.windowPtr);                           % hold
            vbl = Screen('Flip', stm.Vis.windowPtr, vbl + 0.5*stm.Vis.TrialIFI);    % flip
        end
    end
	pause(0.01);
end

Screen('FillRect', stm.Vis.windowPtr, bgcolor);
Screen('Flip', stm.Vis.windowPtr);

%save(['D:\XINTRINSIC\VisSeqData\' stm.Vis.SesTime '_VisSeqData.mat'], 'stm', '-v7.3');
%save([stm.DataDir, filesep, stm.Vis.SesTimeStr, '_', stm.Vis.Name, '_VisSeqData.mat'], 'stm', '-v7.3');
save([stm.DataDir, filesep, stm.Vis.SesTimeStr, '_Stimulus_', stm.Vis.Name, '_SequenceData.mat'], 'stm', '-v7.3');

%% Clean up
pause(0.5);
try sys.MsgBox.delete();
catch
end
Screen('Close');
sca;
if isfile(sys.TempStimFile)
    delete(sys.TempStimFile);
end
