function XinRanAnalysis2_Sweep_ButtonDown(varargin)

if nargin==0                % call from GUI
    H =         gcbo;
else
    H =         varargin{1};
end

ButtonTag =     get(H,              'tag');
hFig =          gcf;
    
%% Update "Spec", "Tune", or "Temp", "Frme", "Pixl"
switch ButtonTag(1:4)
    case 'Spec';    SpecUpdt=1; TuneUpdt=1; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=0;
    case 'Tune';    SpecUpdt=0; TuneUpdt=1; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=0;
    case 'Frme';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=1;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=0;
    case 'Temp';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=1; Play=0; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=0;
    case 'Play';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=0; Play=1; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=0;
    case 'Pixl';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=1; AmpSUpdt=0; TitlUpdt=0;
    case 'AmpS';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=1; TitlUpdt=0;
    case 'Titl';    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=0; TitlUpdt=1;
    otherwise
                    SpecUpdt=0; TuneUpdt=0; TempAmpUpdt=0;  FrmeUpdt=0; Play=0; PixlUpdt=0; AmpSUpdt=0;
end
%% Amp(litude of the) S(pectrum)
if AmpSUpdt
    hAxesSpec =     getappdata(hFig,	'hAxesSpec');
    hDotSpec =      getappdata(hFig,	'hDotSpec');
    YLimRange =     get(hAxesSpec,      'YTick');
    CurYLim =       get(hAxesSpec,      'YLim');
    YLimRange =     YLimRange(3:end);
    CurYLim =       CurYLim(2);
    CurYLimIdx =    find(YLimRange == CurYLim, 1);
    if CurYLimIdx == length(YLimRange)
        CurYLimIdx = 1;
    else
        CurYLimIdx = CurYLimIdx +1;
    end        
    set(hAxesSpec,	'YLim',     [0 YLimRange(CurYLimIdx)]);
    set(hDotSpec,   'YData',	YLimRange(CurYLimIdx));
end

%% Spec(trum sample number selection)
if SpecUpdt
    hDotSpec =          getappdata(hFig,	'hDotSpec');
    hTextSpecRepNum =	getappdata(hFig,	'hTextSpecRepNum');
    hHistAll =          getappdata(hFig,	'hHistAll');
    hHistIn =           getappdata(hFig,	'hHistIn');
    dTrlDurPreStim =    getappdata(hFig,	'dTrlDurPreStim');
    dTrlDurStim =       getappdata(hFig,	'dTrlDurStim');
    dTrlDurTotal =      getappdata(hFig,	'dTrlDurTotal');
    dSesDurTotal =      getappdata(hFig,	'dSesDurTotal');
    dQCycNum =          getappdata(hFig,	'dQCycNum');
    dPtQt_FFTAmp =      getappdata(hFig,	'dPtQt_FFTAmp');
    dPtQt_FFTAgl =      getappdata(hFig,	'dPtQt_FFTAgl');
    dPtIndex_ROIin =	getappdata(hFig,	'dPtIndex_ROIin');
    cQRepNum =          getappdata(hFig,	'cQRepNum');
    SpecUpDown =        getappdata(H,       'UpDown');
    % update cQRepNum
    cQRepNum = cQRepNum + SpecUpDown;
    if cQRepNum<1;                      cQRepNum = 1;                               end
    if cQRepNum>size(dPtQt_FFTAmp,2)/2;	cQRepNum = floor(size(dPtQt_FFTAmp,2)/2);	end    
    % update AxesSpec
    setappdata(hFig,        'cQRepNum',	cQRepNum);
    set(hDotSpec,           'XData',    double((cQRepNum-1))/dSesDurTotal);
    set(hTextSpecRepNum,	'string',   sprintf('rep#:%d', cQRepNum-1));
    % update AxesTune
        dRawHue =	dPtQt_FFTAgl(:,cQRepNum)/(2*pi);
    if cQRepNum == dQCycNum
        dRawHue =	(dRawHue -   dTrlDurPreStim  /dTrlDurTotal) / ...
                            (   dTrlDurStim     /dTrlDurTotal);
    end
        dRawVal =    dPtQt_FFTAmp(:,cQRepNum);
    setappdata(hFig,     'dRawHue',   dRawHue);
    setappdata(hFig,     'dRawVal',   dRawVal);       
    % update AxesHist
    hHistAll.Data =	dPtQt_FFTAgl(:,cQRepNum);   
    hHistIn.Data =	dPtQt_FFTAgl(dPtIndex_ROIin,cQRepNum);     
    drawnow;
end

%% Tuning
if TuneUpdt
    hAxesTuneHueHid =	getappdata(hFig,	'hAxesTuneHueHid');
    hImageTune =		getappdata(hFig,    'hImageTune');
    dTrlDurStim =       getappdata(hFig,	'dTrlDurStim');
    dRawHue =	getappdata(hFig,	'dRawHue');
    dRawSat =	getappdata(hFig,	'dRawSat');
    dRawVal =	getappdata(hFig,	'dRawVal');

    HueMapOptions = getappdata(hImageTune,	'HueMapOptions');
    HueMapOrder =   getappdata(hImageTune,	'HueMapOrder');
    HueTempolate =	getappdata(hImageTune,	'HueTempolate');

    SatParaRange =  getappdata(hImageTune,	'SatParaRange');
    SatParaOrder =  getappdata(hImageTune,	'SatParaOrder');
    SatValSync =    getappdata(hImageTune,	'SatValSync');
    SatGroundOut =  getappdata(hImageTune,	'SatGroundOut');
    SatGroundTime =	getappdata(hImageTune,	'SatGroundTime');

    ValLimRange =   getappdata(hImageTune,	'ValLimRange');
    ValLimOrder =   getappdata(hImageTune,	'ValLimOrder');
    ValLim =        getappdata(hImageTune,	'ValLim');
    
    % What TUNE parameter has changed, if any
    switch ButtonTag(5:end)
        case 'ScalebarHue'  
            HueMapOptions = {...
                'HSLuv',        'HSLuv',        'HSVadjusted';
                'Circular',     'Linear',       'Circular'};
            HueMapOrder =   HueMapOrder + 1;
            if HueMapOrder> length(HueMapOptions)
                HueMapOrder = 1;
            end
            HueMap =        [   HueMapOptions{1, HueMapOrder}, '_',...
                                HueMapOptions{2, HueMapOrder}   ];
            HueTempolate =	HueRedist(  HueMapOptions{1, HueMapOrder},...
                                        HueMapOptions{2, HueMapOrder});
            HueColorMap =   hsv2rgb([HueTempolate ones(length(HueTempolate),2)]);
            setappdata(hImageTune,	'HueMapOptions',    HueMapOptions);
            setappdata(hImageTune,	'HueMapOrder',      HueMapOrder);
            setappdata(hImageTune,	'HueMap',           HueMap);
            setappdata(hImageTune,	'HueTempolate',     HueTempolate);
            setappdata(hImageTune,	'HueColorMap',      HueColorMap);
            colormap(hAxesTuneHueHid, HueColorMap);     caxis(hAxesTuneHueHid, [0 1]);   
            ylabel(H,	'\uparrow Onset                                                    Offset \uparrow',...
                        'FontSize',             8,...
                        'VerticalAlignment',    'Cap');    
        case 'ScalebarSat'
            SatParaOrder =  SatParaOrder + 1;
            if SatParaOrder> length(SatParaRange)
                SatParaOrder = 1;
            end
            SatValSync =    SatParaRange(1, SatParaOrder);
            SatGroundOut =  SatParaRange(2, SatParaOrder);
            SatGroundTime =	SatParaRange(3, SatParaOrder);
            setappdata(hImageTune,	'SatParaOrder',     SatParaOrder);
            setappdata(hImageTune,	'SatValSync',       SatValSync);
            setappdata(hImageTune,  'SatGroundOut',     SatGroundOut);
            setappdata(hImageTune,	'SatGroundTime',	SatGroundTime);
            switch sprintf('%d',[SatGroundOut SatValSync]) 
                case '00'; str = '                                                          everything \uparrow';
                case '01'; str = '\rightarrow sync w/ value \leftarrow';
                case '10'; str = ['\uparrow out of stim\pm', sprintf('%4.1fs range', abs(SatGroundTime)) '       wihtin the range \uparrow']; 
                case '11'; str = ['\uparrow out of stim\pm', sprintf('%4.1fs range', abs(SatGroundTime)) '  \rightarrow sync w/ value \leftarrow   ']; 
                otherwise
            end
            ylabel(H,	str,...
                        'FontSize',             8,...
                        'VerticalAlignment',    'Cap');
        case 'ScalebarVal'
            ValLimOrder =   ValLimOrder + 1;
            if ValLimOrder> length(ValLimRange)
                ValLimOrder = 1;
            end
            ValLim =        ValLimRange(ValLimOrder);
            setappdata(hImageTune,	'ValLimOrder',  ValLimOrder);
            setappdata(hImageTune,	'ValLim',       ValLim);
            ylabel(H,	sprintf('0                       (old plots: %5.3f %%)        %5.3f %%',...
                        ValLim*sqrt(2)*100, ValLim*100), ...
                        'FontSize',             8,...
                        'VerticalAlignment',    'Cap');
        otherwise
            % no TUNE parameter change, change should be on SPEC
    end
    % Calculate and Update TUNE image
    PtOne_Hue =             min(max(dRawHue, 0), 1);            % HUE
    InHueTempolate =        linspace(0, 1, length(HueTempolate))';
    PtOne_Hue =             interp1q(InHueTempolate, HueTempolate, PtOne_Hue);
    if SatValSync                                               % SATURATION
        PtOne_Sat =         min(dRawVal/ValLim,1);
    else
        PtOne_Sat =         dRawSat;
    end
    if SatGroundOut
        PtOne_Sat(dRawHue>(1+SatGroundTime/dTrlDurStim)) =	0;
        PtOne_Sat(dRawHue<(0-SatGroundTime/dTrlDurStim)) =	0;
    end
    PtOne_Val =             min(dRawVal/ValLim,1);              % VALUE
    PtThree_TuneMap =       hsv2rgb([	PtOne_Hue,...           % IMAGE
                                        PtOne_Sat,...
                                        PtOne_Val]);
	SizeImageTune =         size(get(hImageTune, 'CData'));
    PhPwThree_TuneMap =     uint8(255*reshape(PtThree_TuneMap,...
                                SizeImageTune(1),...
                                SizeImageTune(2),...
                                SizeImageTune(3) ) );
    set(hImageTune,    'CData',        PhPwThree_TuneMap);
end
%% Temp(oral) Amp(litude)      
if TempAmpUpdt  
    % callback is from "hFig2AxesFrmeScalebar"
    hAxesFrme =     getappdata(hFig,	'hAxesFrme');
    hAxesTemp =     getappdata(hFig,	'hAxesTemp');
    hAxesPixl =     getappdata(hFig,	'hAxesPixl');
    dAxesTempYLims= getappdata(hFig,	'dAxesTempYLims'); 
    Tmax =          get(hAxesTemp,      'YLim');
    Tmax =          Tmax(2);
    TmaxIdx =       find(dAxesTempYLims==Tmax);
    TmaxIdx =       TmaxIdx + 1;
    if TmaxIdx > length(dAxesTempYLims)
        TmaxIdx =   1;
    end
    Tmax =          dAxesTempYLims(TmaxIdx);
    TmaxLabels =    {	sprintf('-%3.1f%%',Tmax*100),...
                        '0%',...
                        sprintf(' %3.1f%%',Tmax*100) };
    caxis(hAxesFrme,                Tmax*[-1 1]);
    set(H,          'Limits',       Tmax*[-1 1]);
    set(hAxesTemp,	'YLim',         Tmax*[-1 1]);
    set(hAxesPixl,	'YLim',         Tmax*[-1 1]);
    set(H,          'Ticks',        Tmax*[-1 0 1]);
    set(hAxesTemp,	'YTick',        Tmax*[-1 0 1]);
    set(hAxesPixl,	'YTick',        Tmax*[-1 0 1]);
    set(H,          'TickLabels',   TmaxLabels);
    set(hAxesTemp,	'YTickLabels',  TmaxLabels);
    set(hAxesPixl,	'YTickLabels',  TmaxLabels);
end
%% Frame
if FrmeUpdt
    % callback is either from	"T.H.hFig2AxesTempTextRight" or
    %                           "T.H.hFig2AxesTempTextLeft"
    hImageFrme =        getappdata(hFig,	'hImageFrme');
    hTextTemp =         getappdata(hFig,	'hTextTemp');
    hTextFrmeXLabel =   getappdata(hFig,	'hTextFrmeXLabel');
    hDotTemp =          getappdata(hFig,	'hDotTemp');
	dFrmeTime =         getappdata(hFig,	'dFrmeTime');
    dFptPhPw_NormTrlAdj =   getappdata(hFig,	'dFptPhPw_NormTrlAdj');
    cFrmeNum =      getappdata(hFig,	'cFrmeNum');   
    cPlayingNow =	getappdata(hFig,	'cPlayingNow');   
    FrmeUpDown =    getappdata(H,       'UpDown');
    % update cFrmeNum
    TrialTotalFrmeNum =	size(dFptPhPw_NormTrlAdj,1);
    cFrmeNum =          cFrmeNum + FrmeUpDown;
    if cFrmeNum<1;                   cFrmeNum = 1;                  end
    if cFrmeNum>TrialTotalFrmeNum;   cFrmeNum = TrialTotalFrmeNum;	end   
    FrmeTime =          cFrmeNum * dFrmeTime;
    if cPlayingNow
        XLabelString =  sprintf( 'Amplitude @ %4.1f s, in adjusted polarity, playing',...
                            cFrmeNum*dFrmeTime );
    else
        XLabelString =  sprintf('Amplitude @ %4.1f s, in adjusted polarity, click to play',...
                            cFrmeNum*dFrmeTime );
    end
    % update AxesFrme GUI & data
    setappdata(hFig,        'cFrmeNum',	cFrmeNum);
    set(hImageFrme,         'CData',    squeeze(dFptPhPw_NormTrlAdj(cFrmeNum,:,:)));
    set(hTextFrmeXLabel,	'String',   XLabelString);
    set(hDotTemp,           'XData',    cFrmeNum);
    set(hTextTemp,          'string',   sprintf('Time: %4.1f s', FrmeTime));  
end
%% Playing
if Play
    % Video Playback
    cPlayingNow =	getappdata(hFig,	'cPlayingNow');
    cPlayingNow =	1 - cPlayingNow;
    setappdata(hFig,	'cPlayingNow',   cPlayingNow);
    if cPlayingNow
        hTextTempRight =	getappdata(hFig,	'hTextTempRight');
        dFrmeTime =         getappdata(hFig,	'dFrmeTime');
        dTrlDurTotal =      getappdata(hFig,    'dTrlDurTotal');
        dSoundWave =        getappdata(hFig,	'dSoundWave');
        hAP = audioplayer(dSoundWave, 100e3);        
        cFrmeNum =          0;  
        setappdata(hFig,	'hAudioPlayer',     hAP); 
        setappdata(hFig,	'cFrmeNum',         cFrmeNum);    
        pause(0.5);
        timebase = tic;
        hAP.play;
        while toc(timebase)<dTrlDurTotal && getappdata(hFig, 'cPlayingNow')
            if toc(timebase)>(cFrmeNum*dFrmeTime)
                XinRanAnalysis2_Sweep_ButtonDown(hTextTempRight);
                cFrmeNum = cFrmeNum + 1;
            end       
            pause(0.05);
        end
        hAP.pause;     
            cPlayingNow = 0;
            setappdata(hFig,	'cPlayingNow',   cPlayingNow);
                XinRanAnalysis2_Sweep_ButtonDown(hTextTempRight);
    else
        hAP = getappdata(hFig,	'hAudioPlayer');
        hAP.pause;        
    end
end
%% Pixel
if PixlUpdt
    PixelPosi =         get(gca,            'CurrentPoint');
    PixelPosi =         round(PixelPosi(1,1:2));    
    PhInd =             PixelPosi(2);    
    PwInd =             PixelPosi(1);
    hAxesPixl =         getappdata(hFig,	'hAxesPixl');
    hLinePixlAll =      getappdata(hFig,	'hLinePixlAll');
    hLinePixlMean =     getappdata(hFig,	'hLinePixlMean');
    hLineSpecPixl =     getappdata(hFig,	'hLineSpecPixl');
    hTextPixl =         getappdata(hFig,	'hTextPixl');
    dFptCtPhPw_NormTrl= getappdata(hFig,	'dFptCtPhPw_NormTrl');
	dPtQt_FFTAmp =      getappdata(hFig,	'dPtQt_FFTAmp');
    PtInd =             PhInd + (PwInd-1)*size(dFptCtPhPw_NormTrl,3);
    for i = 1:length(hLinePixlAll)
        hLinePixlAll(i).YData =	squeeze(dFptCtPhPw_NormTrl(:,i,PhInd,PwInd));
    end
    hLinePixlMean.YData =       mean(squeeze(dFptCtPhPw_NormTrl(:,:,PhInd,PwInd)),2);
    hAxesPixl.Title.String =    sprintf('Pixel(%dH,%dW): temporal traces across reps',...
        PhInd,PwInd);
    hTextPixl.String =	{...
        sprintf('Across-all:      %5.2f%%', 100*std(reshape(dFptCtPhPw_NormTrl(:,:,PhInd,PwInd), 1, []))),...
        sprintf('Across-rep:    %5.2f%%',	100*std(mean(dFptCtPhPw_NormTrl(:,:,PhInd,PwInd), 2), 0, 1)),...
        sprintf('Across-frame:%5.2f%%',     100*std(mean(dFptCtPhPw_NormTrl(:,:,PhInd,PwInd), 1), 0, 2)),...
        'on STD'};
    hLineSpecPixl.YData =       dPtQt_FFTAmp(PtInd, :);
    drawnow;
end
%% Title clipping
if TitlUpdt
    AxesH =             get(H,      'Parent');
    AxesPosiMatlab =    get(AxesH,  'Position');
    FigPosiMatlab =     get(gcf,    'Position');
    FigName =           get(gcf,    'Name');
    MoniPosiMatlab =    get(0,      'MonitorPosition');
    MoniNumMain =       find(MoniPosiMatlab(:,1) == 1);
    % -x -y are counted from the lowerleft corner of the WINDOWS MAIN DISPLAY
    %   MATLAB� sets the display information values for this property at startup.
    %   The values are static. If your system display settings change, 
    %   for example, if you plug in a new monitor, then the values do not update. 
    %   To refresh the values, restart MATLAB.
    switch ButtonTag(5:8)
        case 'Var1';    AxesPlus = [    -2  -56 54  17  ];  figclip = 1;
        case 'Var2';    AxesPlus = [    -2  -56 64  17  ];  figclip = 1;
        case 'Var3';    AxesPlus = [    -2  -56 64  17  ];  figclip = 1;
        case 'Temp';    AxesPlus = [    -58 -35 10  17  ];  figclip = 1;
        case 'Pixl';    AxesPlus = [    -58 -35 10  17  ];  figclip = 1;
        case 'Spec';    AxesPlus = [    -70 -35 10  17  ];  figclip = 1;
        case 'Tune';    AxesPlus = [    -5  -25 100 17  ];  figclip = 1;
        case 'TuIm';    AxesPlus = [    -3  -3  0   0   ];  figclip = 1;
        case 'Frme';    AxesPlus = [    -67 -25 0   17  ];  figclip = 1;
        case 'Hist';    AxesPlus = [    -55 -57 5   17  ];  figclip = 1;
        case 'Tex1';                                        figclip = 0;
        case 'Tex2';    AxesPlus = [    -1  -1  -1  -1  ];  figclip = 1;
                        AxesPosiMatlab = [ 0 0 FigPosiMatlab(3:4)];
        case 'Tex3';    SesName =  strtok(FigName, '"');
                        clipboard('copy',   SesName(1:end-7));...
                                                            figclip = 0;  
        case 'Tex4';    [~, SoundName] =    strtok(FigName, '"');
                        [~, SoundName] =	strtok(SoundName, '"');
                        clipboard('copy',   SoundName(2:end-1));       
                                                            figclip = 0;     
        otherwise
    end
    if figclip
        AxesAdd =       AxesPlus;
        AxesAdd(3:4) =  AxesAdd(3:4)-AxesAdd(1:2);
        FigCutMatlab =  AxesPosiMatlab + AxesAdd;
        MoniCutMatlab = FigCutMatlab + [FigPosiMatlab(1:2) 0 0];
        MoniCutWin =    MoniCutMatlab;
        MoniCutWin(2) = MoniPosiMatlab(MoniNumMain,4) - MoniCutWin(2) - MoniCutWin(4);
        % The following function requires NirCmd
        % http://www.nirsoft.net/utils/nircmd2.html#using
        % download NirCmd, unzip, and put the .exe files into Windows/System32
        dos([ 'C:\Windows\System32\Nircmd.exe savescreenshot *clipboard* ', ...
            sprintf('%d ', MoniCutWin) ]);
        % The coordinates for Nircmd is: -x, -y, width, height
        %   -x, -y are counted from the upperleft corner of the WINDOWS MAIN DISPLAY
        
    else
    end
end
