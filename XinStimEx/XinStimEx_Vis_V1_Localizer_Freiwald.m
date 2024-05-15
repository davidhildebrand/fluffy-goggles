function XinStimEx_Vis_V1_Localizer_Freiwald
% V1 localizer attempt

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
    opts = struct('WindowStyle', 'modal',... 
                  'Interpreter', 'tex');
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
stm.Vis.Name = 'V1LocalizerTrial';
sys.screenNumber = max(Screen('Screens'));
for i = 1:max(Screen('Screens'))
    info = Screen('Resolution', i);
    if info.hz == 144
        sys.screenNumber = i;
        break
    end
end

stm.SR = 100e3;  % DAQ sampling rate
stm.SesDurTotal = stm.CycleNum * stm.CycleDur;  % rep x sec
stm.Vis.Background = 'gray';
stm.Vis.GratingSF = 2;  % cpd
stm.Vis.GratingTF = 4;  % 2.85;  % Hz, cps
nd = 4;
stm.Vis.GratingDirs = linspace(0, 360-(360/nd), nd)';
% stm.Vis.GratingDirs = [0 90 180 270]
stm.Vis.GratingDur = 10;  % sec

stm.Vis.SesOptionContrast = 'LURD'; % 'URDL';
    % 'U';  Up
    % 'R';  Right
    % 'D';  Down               
    % 'L';  Left
stm.Vis.TrlDurTotal = stm.CycleDur;
stm.Vis.TrlDurPreStim = 2;
stm.Vis.TrlDurStim = stm.Vis.GratingDur;
stm.Vis.TrlDurPostStim = stm.Vis.TrlDurTotal - stm.Vis.TrlDurPreStim - stm.Vis.TrlDurStim;
stm.Vis.TrlNumTotal = length(stm.Vis.GratingDirs);
stm.Vis.TrlNames = cellstr(stm.Vis.SesOptionContrast');     
% stm.Vis.TrlIndexSoundNum =  1:stm.Vis.TrlNumTotal;
% stm.Vis.TrlIndexAddAttNum = ones(1, stm.Vis.TrlNumTotal);

%% Visual CO Parameters
% Session Timer 
stm.Vis.SesCycleDurTotal = stm.Vis.TrlDurTotal * stm.Vis.TrlNumTotal;  % sec
stm.Vis.SesCycleNumTotal = floor(stm.SesDurTotal / stm.Vis.SesCycleDurTotal);  % rep # total
if stm.Vis.SesCycleNumTotal == 0
    error('Provided session duration does not include enough time to cycle through stimuli.')
end
stm.Vis.SesDurTotal = stm.Vis.SesCycleNumTotal * stm.Vis.SesCycleDurTotal;  % sec
% stm.Vis.SesSoundDurTotal = stm.Vis.SesCycleDurTotal;

stm.Vis.SesTrlOrder = 'Randomized';
stm.Vis.SesTrlOrderMat = [];
for i = 1:stm.Vis.SesCycleNumTotal
    stm.Vis.SesTrlOrderMat = [stm.Vis.SesTrlOrderMat; randperm(stm.Vis.TrlNumTotal)];
end

stm.Vis.SesTrlOrderVec = reshape(stm.Vis.SesTrlOrderMat',1,[]);
disp(stm.Vis.SesTrlOrderVec)
% stm.Vis.SesTrlOrderSoundVec = stm.Vis.TrlIndexSoundNum(stm.Vis.SesTrlOrderVec);

stm.Vis.CtrlTrlNumTotal = stm.Vis.SesCycleNumTotal * stm.Vis.TrlNumTotal;  % sec

%% Prepare the Psychtoolbox window

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% Define shades
white = WhiteIndex(sys.screenNumber);
black = BlackIndex(sys.screenNumber);   
gray = GrayIndex(sys.screenNumber, 0.5);
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);

switch stm.Vis.Background
    case 'white'
        bgcolor = white;
    case 'black'
        bgcolor = black;
    case 'gray'
        bgcolor = gray;
end
[stm.Vis.windowPtr, ~] = PsychImaging('OpenWindow', sys.screenNumber, bgcolor);
AssertGLSL;


% Check frame duration and size
stm.Vis.TrialIFI = Screen('GetFlipInterval', stm.Vis.windowPtr);
if (abs(stm.Vis.TrialIFI - 1/144) / (1/144)) > 0.05 % 1/59)/(1/59) > 0.05 
    errordlg('screen is not at right refresh rate!');
    return;
end
vbl = Screen('Flip', stm.Vis.windowPtr);

%% Compare stimulus size in Wang and Freiwald labs and adjust accordingly

%Wang
% stm.Vis.MonitorName =       'LG 32GK850F-B';
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
% MonitorPixelAngle_Wang = stm.Vis.MonitorPixelAngle_Wang;

%Freiwald
stm.Vis.MonitorName =       'ASUS_XG16AHPE';
stm.Vis.MonitorDistance =   30;         % in cm
stm.Vis.MonitorHeight =     193.5 / 10; % in cm  % or 0.026875 cm/px
stm.Vis.MonitorWidth =      344 / 10;   % in cm  % or 0.026875 cm/px
[stm.Vis.MonitorPixelNumX, stm.Vis.MonitorPixelNumY] = ...
                                                Screen('WindowSize', stm.Vis.windowPtr);

stm.Vis.MonitorAngleX =         2 * atan((stm.Vis.MonitorWidth/2) / stm.Vis.MonitorDistance)/pi * 180;  
stm.Vis.MonitorAngleY =         2 * atan((stm.Vis.MonitorHeight/2) / stm.Vis.MonitorDistance)/pi * 180;
stm.Vis.MonitorPixelAngleX =    stm.Vis.MonitorAngleX / stm.Vis.MonitorPixelNumX;
stm.Vis.MonitorPixelAngleY =    stm.Vis.MonitorAngleY / stm.Vis.MonitorPixelNumY;

stm.Vis.MonitorPixelAngle_Freiwald = mean([stm.Vis.MonitorPixelAngleX stm.Vis.MonitorPixelAngleY]);

% MonitorPixelAngle_Freiwald = stm.Vis.MonitorPixelAngle_Freiwald;
% MonitorPixelAngle_Ratio = stm.Vis.MonitorPixelAngle_Wang / stm.Vis.MonitorPixelAngle_Freiwald;


%% Initialize parameters

stm.Vis.MonitorPixelAngle = stm.Vis.MonitorPixelAngle_Freiwald;
stm.Vis.MonitorCenter = [stm.Vis.MonitorPixelNumX/2 stm.Vis.MonitorPixelNumY/2];
% texture patch display size
stm.Vis.MonitorTexPatchPos = round([1280 -720 -1280 720] * ...
    stm.Vis.MonitorPixelAngle_Wang / stm.Vis.MonitorPixelAngle_Freiwald);

stm.SesOption = 'Gratings';
switch stm.SesOption
    case 'Gratings'
        gratingtex = CreateProceduralSquareWaveGrating(stm.Vis.windowPtr, ...
                                                       stm.Vis.MonitorPixelNumX, stm.Vis.MonitorPixelNumY, ...
                                                       [0.5 0.5 0.5 0.0]);
    otherwise
end

%% Setup Timer: NI-DAQ or PC timing

% stm.TimerOption = 'simulated';
stm.TimerOption = 'NI-DAQ';

%% NI-DAQ

switch stm.TimerOption
    case 'NI-DAQ'
        sys.NIDAQ.D.InTimebaseRate =    stm.SR;
        % sys.NIDAQ.D.CycleSmplHigh =     2;
        % sys.NIDAQ.D.CycleSmplLow =      sys.NIDAQ.D.InTimebaseRate * stm.SesCycleTime - ...
        %                                 sys.NIDAQ.D.CycleSmplHigh;
        sys.NIDAQ.D.TrialSmplHigh = 2;
        sys.NIDAQ.D.TrialSmplLow = sys.NIDAQ.D.InTimebaseRate * stm.Vis.TrlDurTotal - ...
                                   sys.NIDAQ.D.TrialSmplHigh;
        import dabs.ni.daqmx.*
        sys.NIDAQ.TaskCO = Task('Recording Session Cycle Switcher');
        sys.NIDAQ.TaskCO.createCOPulseChanTicks(...
            'Intrinsic_PCIe6323', 1, 'Cycle Counter', '100kHzTimebase', ...
            sys.NIDAQ.D.TrialSmplLow, sys.NIDAQ.D.TrialSmplHigh,...
            0, 'DAQmx_Val_Low');
        sys.NIDAQ.TaskCO.cfgImplicitTiming(...
            'DAQmx_Val_FiniteSamps', stm.Vis.CtrlTrlNumTotal+1);
        sys.NIDAQ.TaskCO.cfgDigEdgeStartTrig(...
            'RTSI6',            'DAQmx_Val_Rising');
        sys.NIDAQ.TaskCO.registerSignalEvent(...
            @XinStimEx_Vis_V1_Localizer_Callback, 'DAQmx_Val_CounterOutputEvent');
        sys.NIDAQ.TaskCO.start();
        sys.MsgBox =                msgbox('Click to terminate after current cycle');
        stm.Vis.Running =           1;
    case 'simulated'
        sys.TimerH =                timer;
        sys.TimerH.TimerFcn =       @XinStimEx_Vis_V1_Localizer_Callback;
        sys.TimerH.Period =         stm.Vis.SesCycleDurTotal;
        sys.TimerH.TasksToExecute = stm.Vis.SesCycleNumTotal + 1;
        sys.TimerH.ExecutionMode =  'fixedRate';
        sys.MsgBox =                msgbox('Click to terminate the session after current cycle');
        stm.Vis.Running =           1;
        sys.TimerH.start;   
    otherwise
end

%% Play

% Amplitude of the grating in units of absolute display intensity range: A
% setting of 0.5 means that the grating will extend over a range from -0.5
% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
% displayable range. As we select a background color and offset for the
% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
% will extend the sinewaves values from 0 = total black in the minima of
% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
% than 0.5 don't make sense, as parts of the grating would lie outside the
% displayable range for your computers displays:
amplitude = 0.5;

% freq = Frequency of sine grating in cycles per pixel.
degrees_per_cycle = 1 / stm.Vis.GratingSF;
monitor_width_pix = stm.Vis.MonitorPixelNumX;
monitor_width_mm = stm.Vis.MonitorWidth * 10;
mm_per_degree = pi * stm.Vis.MonitorDistance * 10 / 180;
mm_per_cycle = mm_per_degree * degrees_per_cycle;
monitor_pix_per_mm = monitor_width_pix / monitor_width_mm;
pix_per_cycle = monitor_pix_per_mm * mm_per_cycle;
cycle_per_pix = 1 / pix_per_cycle;
sfreq = cycle_per_pix

% Phase is the phase shift in degrees (0-360 etc.) applied to the grating:
phase = 0;
% Increment of phase shift per redraw
phaseincrement = (stm.Vis.GratingTF * 360) * stm.Vis.TrialIFI;

stm.Vis.CtrlTrlNumCurrent = 0;
stm.Vis.CtrlTrlDurCurrentTimer = tic;

while stm.Vis.Running
    % Session timing
    stm.Vis.SesOn = (stm.Vis.CtrlTrlNumCurrent > 0 && ...
                     stm.Vis.CtrlTrlNumCurrent <= stm.Vis.CtrlTrlNumTotal);
    stm.Vis.CtrlTrlDurCurrent = toc(stm.Vis.CtrlTrlDurCurrentTimer); 

    %% Display
    if stm.Vis.SesOn
        
        if stm.Vis.CtrlTrlDurCurrent <= stm.Vis.TrlDurPreStim
            Screen('FillRect', stm.Vis.windowPtr, bgcolor);
            Screen('Flip', stm.Vis.windowPtr);

        elseif stm.Vis.CtrlTrlDurCurrent >= stm.Vis.TrlDurPreStim && ...
               stm.Vis.CtrlTrlDurCurrent < (stm.Vis.TrlDurPreStim + stm.Vis.TrlDurStim)

           angle = stm.Vis.GratingDirs(stm.Vis.SesTrlOrderVec(stm.Vis.CtrlTrlNumCurrent));
           phase = phase + phaseincrement;
           
           % Draw the grating, centered on the screen, with given rotation 'angle',
           % sine grating 'phase' shift and amplitude, rotating via set
           % 'rotateMode'. Note that we pad the last argument with a 4th
           % component, which is 0. This is required, as this argument must be a
           % vector with a number of components that is an integral multiple of 4,
           % i.e. in our case it must have 4 components:
           %dest = (stm.Vis.MonitorCenter([1 2 1 2]) + stm.Vis.MonitorTexPatchPos)'; 
           dest = [];
           % Initial stimulus parameters for the grating:
           internalRotation = angle;
           if internalRotation
               rotateMode = kPsychUseTextureMatrixForRotation;
           else
               rotateMode = [];
           end
          
           Screen('DrawTexture', stm.Vis.windowPtr, gratingtex, ...
               [], dest, angle, [], [], [], [], rotateMode, [phase, sfreq, amplitude, 0]);
           % Screen('DrawTextures', stm.Vis.windowPtr, stm.Vis.TexIdxNow,...
           %    [], (stm.Vis.MonitorCenter([1 2 1 2]) + stm.Vis.MonitorTexPatchPos)');
           
           vbl = Screen('Flip', stm.Vis.windowPtr, vbl + 0.5*stm.Vis.TrialIFI);
           
        elseif stm.Vis.CtrlTrlDurCurrent >= (stm.Vis.TrlDurPreStim + stm.Vis.TrlDurStim) && ...
               stm.Vis.CtrlTrlDurCurrent < stm.Vis.TrlDurTotal
           Screen('FillRect', stm.Vis.windowPtr, bgcolor);
           Screen('Flip', stm.Vis.windowPtr);
        end
        %         cond = stm.Vis.CtrlTrlNumCurrent;
        %         if stm.Vis.GratingAngleCondCurrent ~= cond
        %
        %         end
    end
	pause(0.01);
end

%% Clear the display

Screen('FillRect', stm.Vis.windowPtr, bgcolor);
Screen('Flip', stm.Vis.windowPtr);

disp('Saving StimulusData file.')
save([stm.DataDir, filesep, stm.Vis.SesTimeStr, '_Stimulus_', stm.Vis.Name, '_StimulusData.mat'], 'stm', '-v7.3');
disp('Saved StimulusData file.')

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
