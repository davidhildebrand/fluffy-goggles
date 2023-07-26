function XinStimEx_Vis_FacePatch_Trial_Callback(~,~)

global stm sys

%% Start a cycle
% This function is called everytime 
%   stm.Vis.TimerOption = 'simulated':  an event happens on the software timer
%   stm.Vis.TimerOption = 'NI-DAQ':     a signal event happens on the hardware CO

if toc(stm.Vis.CtrlTrlDurCurrentTimer) > 0.1    
    % This "if" is necessary for 'NI-DAQ' to require
    % at least 0.1s between real rising-edge triggering events to avoid 
    % pseudo-triggering by the falling-edge events
    stm.Vis.CtrlTrlNumCurrent = stm.Vis.CtrlTrlNumCurrent + 1;
    stm.Vis.CtrlTrlDurCurrentTimer = tic;
    tt = datestr(now, 'HH:MM:SS.FFF');
    try
        disp([sprintf('trial: #%d type ''%s'' ', stm.Vis.CtrlTrlNumCurrent, ...
            stm.Vis.TrlNames{stm.Vis.SesTrlOrderVec(stm.Vis.CtrlTrlNumCurrent)}), ...
            tt]);
    end
%     if strcmp(stm.Vis.SesOption, 'Cali')
%         cali_run(stm.Vis.CtrlTrlNumCurrent);
%     end
    drawnow;
    if isfield(sys, 'MsgBox')
        if ~ishandle(sys.MsgBox)
            stm.Vis.CtrlTrlNumCurrent = stm.Vis.CtrlTrlNumTotal + 1;
        end
    end
end

%% Stop session
if (stm.Vis.CtrlTrlNumCurrent == stm.Vis.CtrlTrlNumTotal + 1) && (stm.Vis.Running == 1)
    stm.Vis.Running = 0;    
    try
        sys.NIDAQ.TaskCO.abort();
        sys.NIDAQ.TaskCO.delete;
    catch
        disp('Encountered an error when trying to stop the internal stimulus trigger.')
        disp('Attempting to stop the software timer.');
        sys.TimerH.stop;
        sys.TimerH.delete;
        try sys.MsgBox.delete();    catch
        end
    end    
end