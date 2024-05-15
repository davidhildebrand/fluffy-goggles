function XinStimEx_Vis_MT_Localizer_Freiwald
% MT localizer with multiple options

clearvars stm sys
global stm sys

ScanImagePath = 'E:\FreiwaldSync\MarmoScope\ScanImage\SI-Basic_2021.1.0_(2021-12-22)_2af5d7cfec';
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
    stm.DataDir = 'D:\XINTRINSIC\';
    stm.StimDir = 'E:\FreiwaldSync\XINTRINSIC\Stimuli\';
end

%% Switch multi-display mode
Screen('Close') 
if max(Screen('Screens')) ~= 2
    opts = struct(  'WindowStyle',  'modal',... 
                    'Interpreter',  'tex');
    errordlg(...
        [   '\fontsize{20} The monitors are not in all extended mode \newline ',...
            'Close thecurrent Matlab \newline',...
            'Extend all screens in windows \newline' ,...
            'And restart Matlab'],...
        '', opts)
    return
end
dos('C:\Windows\System32\DisplaySwitch.exe /extend');
sca;       
pause(0.5);

%% Specify Session Parameters
stm.Vis.SesTime = now;
stm.Vis.SesTimeStr = [datestr(stm.Vis.SesTime, 'yyyymmdd'), ...
    'd', datestr(stm.Vis.SesTime, 'HHMMSS'), 't'];
stm.Vis.Name = 'MTLocalizerTrial';
sys.screenNumber = max(Screen('Screens'));
for i = 1:max(Screen('Screens'))
    info = Screen('Resolution', i);
    if info.hz == 144
        sys.screenNumber = i;
        break
    end
end

% Session Timer 
% stm.TimerOption =       'simulated';
stm.TimerOption =       'NI-DAQ';

% Session Type Options
%stm.SesOption =         'Cali';
stm.SesOption =         'DLCL';  % Dot Localizer
% stm.SesOption =         'DCPS';  % Dot Center vs Periphery, Sinusoidal
% stm.SesOption =         'DRCWF';  % Dot Rotating Quarter, Clockwise, w/ Face;
% stm.SesOption =         'DRCCF';  % Dot Rotating Quarter, CounterClockwise';
% stm.SesOption =         'DRCW';  % Dot Rotating Quarter, Clockwise; SOC.
stm.DotAngleDivide =      10;

stm.SR =                    100e3;
% stm.SesCycleTime =          stm.CycleDur; %    10;     % in second
% sys.SesCycleNumTotal =      stm.CycleNum; %    40;     % rep # total
stm.SesCycleTime =              18;     % in second %SOC.
sys.SesCycleNumTotal =          22;     % rep # total
if strcmp(stm.SesOption, 'Cali')
	stm.SesCycleTime =              1.5;	% in second
    sys.SesCycleNumTotal =          80;     % rep # total
end
    
% Session Parameters
% stm.SesCycleTime =              1.25;     % in second
% sys.SesCycleNumTotal =          80;     % rep # total
% stm.SesCycleTime =              5;     % in second
% sys.SesCycleNumTotal =          80;     % rep # total

sys.SesCycleNumCurrent =        0;
stm.SesCycleTimeInitial =       tic;

% Display Device
% stm.MonitorName =       'LG 32GK850F-B';
% stm.MonitorDistance =   75;             % in cm
% stm.MonitorHeight =     0.02724*1440;	% in cm
% stm.MonitorWidth =      0.02724*2560;	% in cm
stm.MonitorName =       'ASUS_XG16AHPE';
stm.MonitorDistance =   30; % cm
stm.MonitorHeight =     193.5 / 10;	% cm
stm.MonitorWidth =      344 / 10; % cm

%% Prepare the Psychtoolbox window
%% Prepare the Psychtoolbox window
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% Define shades
white = WhiteIndex(sys.screenNumber);
black = BlackIndex(sys.screenNumber);   
gray =  GrayIndex(sys.screenNumber, 0.5);
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);

% switch stm.Vis.PicBackground
%     case 'white'
%         bgcolor = white;
%     case 'black'
%         bgcolor = black;
%     case 'gray'
%         bgcolor = gray;
% end

if strcmp(stm.SesOption, 'Cali')
    bgcolor = gray;
else
    bgcolor = black;
end
[stm.windowPtr, ~] = PsychImaging('OpenWindow', sys.screenNumber, bgcolor);

% Query: Get the size of the on screen window
[stm.MonitorPixelNumX, stm.MonitorPixelNumY] = Screen('WindowSize', stm.windowPtr);
% Query: the frame duration
stm.TrialIFI = Screen('GetFlipInterval', stm.windowPtr);
% Querry the max dot size allowed
[~, stm.DotDiameterInPixelMax, ~, ~] = Screen('DrawDots', stm.windowPtr);

%% Initialize parameters
switch stm.SesOption(1)
    case 'C'
%         cali_init(stm.windowPtr, @generate_half_sequence, 'dot', 0.5);  % 80 x 5s
%         cali_init(stm.windowPtr, @partition_small, 'face', 1);  % 20 x 20s
%         cali_init(stm.windowPtr, @partition_large, 'face', 1);  % 20 x 20s   
%         cali_init(stm.windowPtr, @partition_small, 'face', 0.5);  % 20 x 20s
%         cali_init(stm.windowPtr, @generate_half_sequence, 'dot', 0.5);   % 80 x 5s
%         cali_init(stm.windowPtr, @generate_center80_sequence,	'dot',  0.5, 'center'); 
        cali_init(stm.windowPtr, @generate_center80_sequence,	'face', 1, 'center');
%         cali_init(stm.windowPtr, @generate_center80_sequence,	'face', 0.5, 'center');
%         cali_init(stm.windowPtr, @generate_center78_sequence,	'face', 1); % 78x   1, 2s
%         cali_init(stm.windowPtr, @generate_full_sequence,       'face', 1); % 144x  1, 2s 
    case 'D'
        stm.DotDiameter =       0.4;        % in degree
        stm.DotMotionSpeedMax =	16;         % in degree / second
        stm.DotDensityAngle =   180;        % as how many dots are in the field

        stm.MonitorAngleX =         2*atan(stm.MonitorWidth/2/stm.MonitorDistance)/pi*180;  
        stm.MonitorAngleY =         2*atan(stm.MonitorHeight/2/stm.MonitorDistance)/pi*180;
        stm.MonitorPixelAngleX =    stm.MonitorAngleX/stm.MonitorPixelNumX;
        stm.MonitorPixelAngleY =    stm.MonitorAngleY/stm.MonitorPixelNumY;
        stm.MonitorPixelAngle =     mean([stm.MonitorPixelAngleX stm.MonitorPixelAngleY]);
        stm.MonitorCenter =         [stm.MonitorPixelNumX/2 stm.MonitorPixelNumY/2];
        stm.DotDiameterInPixel =    stm.DotDiameter/stm.MonitorPixelAngle;
        if stm.DotDiameterInPixel > stm.DotDiameterInPixelMax 
            errordlg('Dot diameter set too big!')
        end
        stm.DotCenterRadiusMax =    stm.MonitorAngleY/2;
        stm.DotVecLength =          stm.DotDensityAngle*1;
        stm.DotVecAngle =           rand(1, stm.DotDensityAngle)*360;
        stm.DotVecRadius =          rand(1, stm.DotDensityAngle)*stm.DotCenterRadiusMax;
        stm.DotVecAngleInit =       stm.DotVecAngle;
        stm.DotVecRadiusInit =      stm.DotVecRadius;      
        stm.DotVecPositionX =       cos(stm.DotVecAngle/180*pi).*...
                                    (stm.DotVecRadius/stm.MonitorPixelAngleX);
        stm.DotVecPositionY =       sin(stm.DotVecAngle/180*pi).*...
                                    (stm.DotVecRadius/stm.MonitorPixelAngleY);
        stm.DotMotionStepMax =      stm.DotMotionSpeedMax*stm.TrialIFI;
        stm.TrialFrameSeq =         1:round(stm.SesCycleTime/stm.TrialIFI);
        switch stm.SesOption(2:4)
            case 'LCL'
                stm.DotMotionPeakTime =     4;    
                slope =                     1/(stm.SesCycleTime/2 - stm.DotMotionPeakTime);
                shift1 =                    slope*stm.SesCycleTime/2;
                shift2 =                    slope*(stm.SesCycleTime/2 - stm.DotMotionPeakTime/2);
            case 'RCW'
                stm.DotAngleMinInitial =    -135;
                if length(stm.SesOption)>4
                    if stm.SesOption(5) == 'F'
                        for i = 1:9
                            stm.FaceLib{i} =    imread(['E:\FreiwaldSync\XINTRINSIC\Stimuli\Visual\profthecopyright_calibration\images\m' num2str(i) '.png']);
                            stm.FaceTexLib(i) = Screen('MakeTexture', stm.windowPtr, stm.FaceLib{i});
                            stm.FaceTic =       tic;
                            stm.FaceTocTime =   2;
                            stm.FaceIdxCurrent = randi([1 length(stm.FaceTexLib)]);
                        end
                    end
                end
            case 'RCC'
                stm.DotAngleMinInitial =    -135;
                if length(stm.SesOption)>4
                    if stm.SesOption(5) == 'F'
                        for i = 1:9
                            stm.FaceLib{i} =    imread(['E:\FreiwaldSync\XINTRINSIC\Stimuli\Visual\profthecopyright_calibration\images\m' num2str(i) '.png']);
                            stm.FaceTexLib(i) = Screen('MakeTexture', stm.windowPtr, stm.FaceLib{i});
                            stm.FaceTic =       tic;
                            stm.FaceTocTime =   2;
                            stm.FaceIdxCurrent = randi([1 length(stm.FaceTexLib)]);
                        end
                    end
                end
            case 'CPS'
                stm.DotAngleDivide =        stm.DotAngleDivide;
        end
        Screen('DrawDots', stm.windowPtr, [stm.DotVecPositionX; stm.DotVecPositionY],...
            stm.DotDiameterInPixel, white,	stm.MonitorCenter, 2);
        vbl = Screen('Flip', stm.windowPtr);
    otherwise
end

%% Setup Timer: NI-DAQ or PC timing
switch stm.TimerOption
    case 'NI-DAQ'
        sys.NIDAQ.D.InTimebaseRate =    stm.SR;
        sys.NIDAQ.D.CycleSmplHigh =     2;
        sys.NIDAQ.D.CycleSmplLow =      sys.NIDAQ.D.InTimebaseRate * stm.SesCycleTime - ...
                                        sys.NIDAQ.D.CycleSmplHigh;
        import dabs.ni.daqmx.*
        sys.NIDAQ.TaskCO = Task('Recording Session Cycle Switcher');
        sys.NIDAQ.TaskCO.createCOPulseChanTicks(...
            'Intrinsic_PCIe6323', 1, 'Cycle Counter', '100kHzTimebase', ...
            sys.NIDAQ.D.CycleSmplLow, sys.NIDAQ.D.CycleSmplHigh,...
            0, 'DAQmx_Val_Low');
        sys.NIDAQ.TaskCO.cfgImplicitTiming(...
            'DAQmx_Val_FiniteSamps',	sys.SesCycleNumTotal+1);
        sys.NIDAQ.TaskCO.cfgDigEdgeStartTrig(...
            'RTSI6',            'DAQmx_Val_Rising');
        sys.NIDAQ.TaskCO.registerSignalEvent(...
            @XinStimEx_Vis_MT_Localizer_Callback, 'DAQmx_Val_CounterOutputEvent');
        sys.NIDAQ.TaskCO.start();
        sys.MsgBox =                msgbox('Click to terminate after current cycle');
        stm.Running =               1;
    case 'simulated'
        sys.TimerH =                timer;
        sys.TimerH.TimerFcn =       @XinStimEx_Vis_MT_Localizer_Callback;
        sys.TimerH.Period =         stm.SesCycleTime;
        sys.TimerH.TasksToExecute = sys.SesCycleNumTotal+1;
        sys.TimerH.ExecutionMode =  'fixedRate';
        sys.MsgBox =                msgbox('Click to terminate the session after current cycle');
        stm.Running =               1;
        %pause; SOC.
        sys.TimerH.start;   
    otherwise
end

%% Play

while stm.Running
    % Session timing
    stm.SesOn = (sys.SesCycleNumCurrent > 0) && (sys.SesCycleNumCurrent <= sys.SesCycleNumTotal);
    stm.SesCycleTimeCurrent = toc(stm.SesCycleTimeInitial); 
    
    % Estimate current dot step parameters   
    if stm.SesOption(1) == 'D'
        switch stm.SesOption(2:4)
            case 'LCL'
                stm.DotMotionSpeedNorm =	stm.SesOn*(...
                        min(1, max(0, shift2-abs(stm.SesCycleTimeCurrent*slope-shift1) ) ));
                stm.DotRadiusMinCurrent =   0;
                stm.DotRadiusMaxCurrent =   stm.MonitorAngleY/2;
                stm.DotAngleMinCurrent =    0;
                stm.DotAngleMaxCurrent =    360;
                stm.DotVecMovingIdx =       1+0*stm.DotVecRadius;
            case 'RCW'
                stm.DotMotionSpeedNorm =    1;
                stm.DotRadiusMinCurrent =   0;
                stm.DotRadiusMaxCurrent =   stm.MonitorAngleY/2;
                stm.DotAngleMinCurrent =    stm.DotAngleMinInitial+stm.SesOn*360/stm.SesCycleTime* stm.SesCycleTimeCurrent;
                stm.DotAngleMaxCurrent =    stm.DotAngleMinCurrent+90;
                stm.DotVecMovingIdx =       mod(stm.DotVecAngle-stm.DotAngleMinCurrent, 360)<90;
                %disp(num2str(stm.DotAngleMinCurrent))
            case 'RCC'
                stm.DotMotionSpeedNorm =    1;
                stm.DotRadiusMinCurrent =   0;
                stm.DotRadiusMaxCurrent =   stm.MonitorAngleY/2;
                stm.DotAngleMinCurrent =    stm.DotAngleMinInitial-stm.SesOn*360/stm.SesCycleTime* stm.SesCycleTimeCurrent;
                stm.DotAngleMaxCurrent =    stm.DotAngleMinCurrent+90;
                stm.DotVecMovingIdx =       mod(stm.DotVecAngle-stm.DotAngleMinCurrent, 360)<90;
            case 'CPS'
                stm.DotMotionSpeedNorm =    1;
    %             stm.DotRadiusMinCurrent =   (cos(2*pi*stm.SesCycleTimeCurrent/stm.SesCycleTime)+1)/2 * stm.DotAngleDivide/2;
    %             stm.DotRadiusMaxCurrent =   (cos(2*pi*stm.SesCycleTimeCurrent/stm.SesCycleTime)+1)/2 * stm.MonitorAngleY/2 + ... 
    %                                         stm.DotAngleDivide/2;
                stm.DotRadiusMinCurrent =   0;
                stm.DotRadiusMaxCurrent =   stm.MonitorAngleY/2;
                stm.DotAngleMinCurrent =    0;
                stm.DotAngleMaxCurrent =    360;
                stm.DotVecMovingIdx =       stm.SesOn + 0*stm.DotVecRadius;
            otherwise
        end
        % Move the dots radially: move the ones on the DocVecMovingIdx, with DotMotionStepCurrent
        stm.DotVecIndexOut =        stm.DotVecRadius < stm.DotRadiusMaxCurrent;
        stm.DotMotionStepCurrent =  stm.DotMotionSpeedNorm*stm.DotMotionStepMax*stm.DotVecMovingIdx;
        stm.DotVecRadius =          stm.DotVecRadius + stm.DotMotionStepCurrent;
 
        % Replace the missing dots
        stm.DotVecIndexOut =        stm.DotVecIndexOut & (stm.DotVecRadius > stm.DotRadiusMaxCurrent);  
        stm.DotVecIndexIn =         ~stm.DotVecIndexOut;                      
        stm.DotVecRadius(stm.DotVecIndexOut) = ...
                                    stm.DotRadiusMinCurrent; 
        stm.DotVecAngle =           stm.DotVecIndexIn.*stm.DotVecAngle + ...
                                    stm.DotVecIndexOut.*...
                                    (  ( stm.DotAngleMinCurrent + ...
                                    (stm.DotAngleMaxCurrent-stm.DotAngleMinCurrent)*...
                                    rand(1,stm.DotVecLength) ).*stm.DotVecIndexOut  );    
        % Masking
        if strcmp(stm.SesOption(2:4), 'CPS')
            stm.DotRadiusMaskerMinCurrent =   (cos(2*pi*stm.SesCycleTimeCurrent/stm.SesCycleTime)+1)/2 * stm.DotAngleDivide/2;
            stm.DotRadiusMaskerMaxCurrent =   (cos(2*pi*stm.SesCycleTimeCurrent/stm.SesCycleTime)+1)/2 * stm.MonitorAngleY/2 + ... 
                                                    stm.DotAngleDivide/2;
            stm.DotMaskerIdx =      (stm.DotVecRadiusInit   >stm.DotRadiusMaskerMaxCurrent) | ...
                                    (stm.DotVecRadiusInit   <stm.DotRadiusMaskerMinCurrent);
            stm.DotMotionIdx =      (stm.DotVecRadius       <stm.DotRadiusMaskerMaxCurrent) & ...
                                    (stm.DotVecRadius       >stm.DotRadiusMaskerMinCurrent);   
            stm.DotVecAngleOut =	[stm.DotVecAngle(stm.DotMotionIdx)  stm.DotVecAngleInit(stm.DotMaskerIdx)];
            stm.DotVecRadiusOut =	[stm.DotVecRadius(stm.DotMotionIdx) stm.DotVecRadiusInit(stm.DotMaskerIdx)];            
        else
            stm.DotVecAngleOut =	stm.DotVecAngle;
            stm.DotVecRadiusOut =	stm.DotVecRadius;
        end

        % Translate the display coordinates
        stm.DotVecPositionX =       cos(stm.DotVecAngleOut/180*pi).*...
                                    (stm.DotVecRadiusOut/stm.MonitorPixelAngle);
        stm.DotVecPositionY =       sin(stm.DotVecAngleOut/180*pi).*...
                                    (stm.DotVecRadiusOut/stm.MonitorPixelAngle);

        %% Display 
        Screen('DrawDots', stm.windowPtr, [stm.DotVecPositionX; stm.DotVecPositionY],...
            stm.DotDiameterInPixel, white, stm.MonitorCenter, 2);      % dots
        if length(stm.SesOption) > 4
            if stm.SesOption(5) == 'F'
                if toc(stm.FaceTic) > stm.FaceTocTime
                    stm.FaceIdxCurrent = randi([1 length(stm.FaceTexLib)]);
                    stm.FaceTic = tic;
                end
                Screen('DrawTextures', stm.windowPtr, stm.FaceTexLib(stm.FaceIdxCurrent),...
                    [], (stm.MonitorCenter([1 2 1 2]) + [80 -80 -80 80])');
            end
        end                                                             % face (center attractor)
        Screen('DrawingFinished', stm.windowPtr);                       % hold
        vbl = Screen('Flip', stm.windowPtr, vbl + 0.5*stm.TrialIFI);    % flip
    else
        pause(0.1);
    end
end

%% Clear the display

Screen('FillRect', stm.windowPtr, bgcolor);
Screen('Flip', stm.windowPtr);

save([stm.DataDir, filesep, stm.Vis.SesTimeStr, '_Stimulus_', stm.Vis.Name, '_StimulusData.mat'], 'stm', '-v7.3');

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
disp('Finished.')