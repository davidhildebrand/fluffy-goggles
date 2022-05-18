
% function DriftDemo4([angle=0][, cyclespersecond=1][, freq=1/360][, gratingsize=360][, internalRotation=0])
% ___________________________________________________________________
%
% Display an animated grating, using the new Screen('DrawTexture') command.
% This demo demonstrates fast drawing of such a grating via use of procedural
% texture mapping. It only works on hardware with support for the GLSL
% shading language, vertex- and fragmentshaders. The demo ends if you press
% any key on the keyboard.
%
% The grating is not encoded into a texture, but instead a little algorithm - a
% procedural texture shader - is executed on the graphics processor (GPU)
% to compute the grating on-the-fly during drawing.
%
% This is very fast and efficient! All parameters of the grating can be
% changed dynamically. For a similar approach wrt. Gabors, check out
% ProceduralGaborDemo. For an extremely fast aproach for drawing many Gabor
% patches at once, check out ProceduralGarboriumDemo. That demo could be
% easily customized to draw many sine gratings by mixing code from that
% demo with setup code from this demo.
%
% Optional Parameters:
% 'angle' = Rotation angle of grating in degrees.
% 'internalRotation' = Shall the rectangular image patch be rotated
% (default), or the grating within the rectangular patch?
% gratingsize = Size of 2D grating patch in pixels.
% freq = Frequency of sine grating in cycles per pixel.
% cyclespersecond = Drift speed in cycles per second.
%

% History:
% 3/1/9  mk   Written.


% dos('C:\Windows\System32\DisplaySwitch.exe /extend');
% sca;                    % Clear the screen       
% pause(2);
% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;
clearvars;              % Clear the workspace
global stm sys

%stm.TimerOption = 'simulated';
stm.TimerOption = 'NI-DAQ';

stm.SesCycleTime =              20;% in seconds
sys.SesCycleNumTotal =          2;
TransitionDuration =       1;


%%% 2 - 15 deg/sec
%%% 0.5 - 2 cycle/deg
%%% 0, 45, 90, 135 deg orientations
%FREQUENCY

degrees_per_second = 10;
cycles_per_degree = 1;
angle = 0;

distance_monitor_mm = 300;
monitor_width_mm = 294;
monitor_width_pix = 1920;
res = [1920 1200];


mm_per_degree = pi * distance_monitor_mm / 180;
mm_per_cycle = mm_per_degree / cycles_per_degree;
monitor_pix_per_mm = monitor_width_pix / monitor_width_mm;
pix_per_cycle = monitor_pix_per_mm * mm_per_cycle;
cycle_per_pix = 1/pix_per_cycle;

%freq = 100/360; % Frequency of the grating in cycles per pixel: Here 1 cycle per pixel
freq = cycle_per_pix;% cycle_per_pix;

cyclespersecond = cycles_per_degree * degrees_per_second;

% Initial stimulus parameters for the grating patch:
internalRotation = 0;


if internalRotation
    rotateMode = kPsychUseTextureMatrixForRotation;
else
    rotateMode = [];
end  

% res is the total size of the patch in x- and y- direction, i.e., the
% width and height of the mathematical support:
%gratingsize = 2000;



% Choose screen with maximum id - the secondary display on a dual-display
% setup for display:
Screen('Screens')
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 128 = gray, i.e. 50% max intensity:
win = Screen('OpenWindow', screenid, 128);

% Make sure the GLSL shading language is supported:
AssertGLSL;

% Retrieve video redraw interval for later control of our animation timing:
ifi = Screen('GetFlipInterval', win);

% Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
phase = 0;

% Compute increment of phase shift per redraw 
phaseincrement = (cyclespersecond * 360) * ifi;

HalfPhaseDuration = 0.5 * ( stm.SesCycleTime /2 - TransitionDuration ) ;
 

% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]);

% Wait for release of all keys on keyboard, then sync us to retrace:
KbReleaseWait;
vbl = Screen('Flip', win);


% Amplitude of the grating in units of absolute display intensity range: A
% setting of 0.5 means that the grating will extend over a range from -0.5
% up to 0.5, i.e., it will cover a total range of 1.0 == 100% of the total
% displayable range. As we select a background color and offset for the
% grating of 0.5 (== 50% nominal intensity == a nice neutral gray), this
% will extend the sinewaves values from 0 = total black in the minima of
% the sine wave up to 1 = maximum white in the maxima. Amplitudes of more
% than 0.5 don't make sense, as parts of the grating would lie outside the
% displayable range for your computers displays:


sys.SesCycleNumCurrent =        0;
stm.SesCycleTimeInitial =       tic;


%% Setup Timer: NI-DAQ or PC timing
switch stm.TimerOption
    case 'NI-DAQ'
        sys.NIDAQ.D.InTimebaseRate =    100e3;
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
            @Grating_Intrinsic_Callback, 'DAQmx_Val_CounterOutputEvent');
        sys.NIDAQ.TaskCO.start();
        stm.Running =               1;
    case 'simulated'
        sys.TimerH =                timer;
        sys.TimerH.TimerFcn =       @Grating_Intrinsic_Callback;
        sys.TimerH.Period =         stm.SesCycleTime;
        sys.TimerH.TasksToExecute = sys.SesCycleNumTotal+1;
        sys.TimerH.ExecutionMode =  'fixedRate';
        sys.MsgBox =                msgbox('Click to terminate the session after current cycle');
        stm.Running =               1;
        %pause; SOC.
        sys.TimerH.start;   
    otherwise
end


while stm.Running
    % Session timing
    stm.SesOn =           (   sys.SesCycleNumCurrent>0 && ...
                                    sys.SesCycleNumCurrent<=sys.SesCycleNumTotal);
    stm.SesCycleTimeCurrent =   toc(stm.SesCycleTimeInitial);
    AbsTimeToHalfCycle = abs( stm.SesCycleTimeCurrent - stm.SesCycleTime/2 );
    
    %if stm.SesOn == 1
        %Amplitude modulation based on timing
        if AbsTimeToHalfCycle < HalfPhaseDuration %Stimulus ON
            amplitude = 0.5;
        elseif AbsTimeToHalfCycle > HalfPhaseDuration + TransitionDuration %Stimulus OFF
            amplitude = 0;
        else %Transition
            amplitude = 0.5 * (1 - (AbsTimeToHalfCycle - HalfPhaseDuration) / TransitionDuration);
        end
        %amplitude = 0.5;
        amplitude = amplitude * stm.SesOn;

        % Increment phase by 1 degree:
        phase = phase + phaseincrement;

        % Draw the grating, centered on the screen, with given rotation 'angle',
        % sine grating 'phase' shift and amplitude, rotating via set
        % 'rotateMode'. Note that we pad the last argument with a 4th
        % component, which is 0. This is required, as this argument must be a
        % vector with a number of components that is an integral multiple of 4,
        % i.e. in our case it must have 4 components:
        Screen('DrawTexture', win, gratingtex, [], [], angle, [], [], [], [], rotateMode, [phase, freq, amplitude, 0]);

        % Show it at next retrace:
        vbl = Screen('Flip', win, vbl + 0.5 * ifi);
    %end
   
end

% We're done. Close the window. This will also release all other ressources:
sca;