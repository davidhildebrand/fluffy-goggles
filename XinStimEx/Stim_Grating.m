
clearvars;              % Clear the workspace
global stm sys
% Session Timer 
stm.TimerOption =       'simulated'; %SOC.
% stm.TimerOption =       'NI-DAQ';


% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

% Tilt angle of the grating:
angle = 0;

% Frequency of the grating in cycles per pixel: Here 0.01 cycles per pixel:
freq = 1/360;
freq = freq;


degrees_per_cycle = 5;
distance_monitor_mm = 100;
mm_per_degree = pi * distance_monitor_mm / 180;
mm_per_cycle = mm_per_degree * degrees_per_cycle;



monitor_width_mm = 294;


cyclespersecond = 1;
gratingsize = 360;
% res is the total size of the patch in x- and y- direction, i.e., the
% width and height of the mathematical support:
res = [gratingsize gratingsize];


% Initial stimulus parameters for the grating patch:
internalRotation = 0;

if internalRotation
    rotateMode = kPsychUseTextureMatrixForRotation;
else
    rotateMode = [];
end





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

% Choose screen with maximum id - the secondary display on a dual-display
% setup for display:
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

% Compute increment of phase shift per redraw:
phaseincrement = (cyclespersecond * 360) * ifi;

% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.
gratingtex = CreateProceduralSineGrating(win, res(1), res(2), [0.5 0.5 0.5 0.0]);

% Wait for release of all keys on keyb oard, then sync us to retrace:
KbReleaseWait;
vbl = Screen('Flip', win);

stm.SesCycleTime =              20;% in seconds
sys.SesCycleNumTotal =          4;

stm.DotMotionPeakTime =     4;    
slope =                     1/(stm.SesCycleTime/2 - stm.DotMotionPeakTime);
shift1 =                    slope*stm.SesCycleTime/2;
shift2 =                    slope*(stm.SesCycleTime/2 - stm.DotMotionPeakTime/2);


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
            @XinStimEx_Vis_MT_Localizer_Callback, 'DAQmx_Val_CounterOutputEvent');
        sys.NIDAQ.TaskCO.start();
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


while stm.Running
    
    % Session timing
    stm.SesOn =           (   sys.SesCycleNumCurrent>0 && ...
                                    sys.SesCycleNumCurrent<=sys.SesCycleNumTotal);
    stm.SesCycleTimeCurrent =    toc(stm.SesCycleTimeInitial); 
    
    % Update some grating animation parameters:
    % DriftDemo4(45,2,36/360,1000)
    % Increment phase by 1 degree:
    amplitude =	0.5 * stm.SesOn*(...
                        min(1, max(0, shift2-abs(stm.SesCycleTimeCurrent*slope-shift1) ) )...
                                                            );
                                                        
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
end

% We're done. Close the window. This will also release all other ressources:
sca;
