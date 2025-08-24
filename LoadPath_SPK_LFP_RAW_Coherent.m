function [LFP_SPK_PSD_RMS, LFP_SPK_Coherence, varargout] = LoadPath_SPK_LFP_RAW_Coherent(filepath, varargin) 

% varargin{1}: electrode (RE1,RE2,LE1,LE2);
% varargin{2}: Regions = 'Pre' or 'Dos' or 'Ven';
% varargin{3}: Colors;
% varargin{4}: axis(4 vectors: LFP_Signal,LFP_PSD,SPK_Signal,SPK_PSD);
%{
'CLFP_01','CLFP_02','CLFP_01_KHz','CLFP_02_KHz',...
    'CSPK_01','CSPK_02','CSPK_01_KHz','CSPK_02_KHz',...
    'CLFP_01_TimeBegin','CLFP_01_TimeEnd','CLFP_02_TimeBegin','CLFP_02_TimeEnd',...
    'CSPK_01_TimeBegin','CSPK_01_TimeEnd','CSPK_02_TimeBegin','CSPK_02_TimeEnd',...
    'CLFP_01_BitResolution','CLFP_01_Gain','CLFP_02_BitResolution','CLFP_02_Gain',...
    'CSPK_01_BitResolution','CSPK_01_Gain','CSPK_02_BitResolution','CSPK_02_Gain',...
%}
load(filepath, 'CRAW_01','CRAW_01_KHz','CRAW_02','CRAW_02_KHz',...
    'CRAW_01_TimeBegin','CRAW_01_TimeEnd','CRAW_02_TimeBegin','CRAW_02_TimeEnd',...
    'CRAW_01_BitResolution','CRAW_01_Gain','CRAW_02_BitResolution','CRAW_02_Gain');
% CRAW_01: raw signal acquired in A/D(analog/digital) value.
% CRAW_01_BitResolution: 1 A/D = CRAW_01_BitResolution uV (1 A/D value equivalent in how much uV).
% CRAW_01_Gain: how much the multiple the raw signal is amplified.
% CRAW_01_KHz: sampling rate in KHz.
% CRAW_01_TimeBegin: Time of the first sample starting to be recorded (Sec unit).
% CRAW_01_TimeEnd: Time of the last sample saved (Sec unit).

%% Transfer the A/D value(artificial unit) to truth value(uV unit)
%
CRAW_01 = double(CRAW_01);
CRAW_02 = double(CRAW_02);
RAW_01 = CRAW_01 .* (CRAW_01_BitResolution ./ CRAW_01_Gain);
RAW_02 = CRAW_02 .* (CRAW_02_BitResolution ./ CRAW_02_Gain);
% RAW_01 = CRAW_01 .* (CRAW_01_Gain./CRAW_01_BitResolution);
% RAW_02 = CRAW_02 .* (CRAW_02_Gain./CRAW_02_BitResolution);

%% Transfer time of recording and sampling rate
t_1 = linspace(0, CRAW_01_TimeEnd-CRAW_01_TimeBegin, length(RAW_01));
t_2 = linspace(0, CRAW_02_TimeEnd-CRAW_02_TimeBegin, length(RAW_02));
%}
%{
% note: remember to comment this part after getting the time
t_1 = CRAW_01_TimeEnd-CRAW_01_TimeBegin;
t_2 = CRAW_02_TimeEnd-CRAW_02_TimeBegin;
LFP_SPK_PSD_RMS = t_1;
LFP_SPK_Coherence = t_2;
if nargout==3
    varargout{1} = FH;
end
%}
%
Sample_Rate_1 = CRAW_01_KHz*1000; %(KHz to Hz)
Sample_Rate_2 = CRAW_02_KHz*1000;

%% this method just can help delete artificial noise in time domain. needn't to delete in frequency domain again
%RAW_01 = Delete_Artificial_Niose_Temporary(RAW_01, Sample_Rate_1);
%RAW_02 = Delete_Artificial_Niose_Temporary(RAW_02, Sample_Rate_2);

%% Set parameters for pwelch and filter
[Win_1,Win_2,Overlap_1,Overlap_2,f] = Set_Pwech_Parameters(Sample_Rate_1,Sample_Rate_2);

%% Do filter to get LFP and SPK, then do pwelch
[LFP_PSD_1,LFP_fre_1,LFP_RMS_1,LFP_Raw_1,~] = Do_Filter_Pwelch(RAW_01,Win_1,Overlap_1,f,Sample_Rate_1,'LFP');
[LFP_PSD_2,LFP_fre_2,LFP_RMS_2,LFP_Raw_2,~] = Do_Filter_Pwelch(RAW_02,Win_2,Overlap_2,f,Sample_Rate_2,'LFP');
[SPK_PSD_1,SPK_fre_1,SPK_RMS_1,SPK_Raw_1,SPK_Abs_Mean_1] = Do_Filter_Pwelch(RAW_01,Win_1,Overlap_1,f,Sample_Rate_1,'SPK');
[SPK_PSD_2,SPK_fre_2,SPK_RMS_2,SPK_Raw_2,SPK_Abs_Mean_2] = Do_Filter_Pwelch(RAW_02,Win_2,Overlap_2,f,Sample_Rate_2,'SPK');

%% Do coherence of LFP and SPK signal 
if ~all(isnan([LFP_Raw_1,LFP_Raw_2,SPK_Abs_Mean_1,SPK_Abs_Mean_2]))
    LFP_Signal_1 = LFP_Raw_1;         LFP_Signal_2 = LFP_Raw_2;
    SPK_Signal_1 = SPK_Abs_Mean_1;    SPK_Signal_2 = SPK_Abs_Mean_2;
    [LFP_SPK_Coherence_1, Fre] = Do_Signal_Coherence(LFP_Signal_1,SPK_Signal_1,Win_1,Overlap_1,f,Sample_Rate_1,Sample_Rate_1);
    [LFP_SPK_Coherence_2, Fre] = Do_Signal_Coherence(LFP_Signal_2,SPK_Signal_2,Win_2,Overlap_2,f,Sample_Rate_2,Sample_Rate_2);
else
    LFP_Signal_1 = NaN;         LFP_Signal_2 = NaN;
    SPK_Signal_1 = NaN;         SPK_Signal_2 = NaN;
    LFP_SPK_Coherence_1 = NaN(1, length(f));
    LFP_SPK_Coherence_2 = NaN(1, length(f));
    Fre = NaN(1, length(f));
end
%% Delete the artifical Hz (50 or n*50 Hz): this method just can help delete artificial noise in frequency domain
if 0
    LFP_PSD_1 = Delete_Artifact_Hz(LFP_PSD_1, LFP_fre_1);
    LFP_PSD_2 = Delete_Artifact_Hz(LFP_PSD_2, LFP_fre_2);
    SPK_PSD_1 = Delete_Artifact_Hz(SPK_PSD_1, SPK_fre_1);
    SPK_PSD_2 = Delete_Artifact_Hz(SPK_PSD_2, SPK_fre_2);
end
if 1
    LFP_PSD_1 = Linear_Interpolate_Each_Site(LFP_fre_1,LFP_PSD_1);
    LFP_PSD_2 = Linear_Interpolate_Each_Site(LFP_fre_2,LFP_PSD_2);
    SPK_PSD_1 = Linear_Interpolate_Each_Site(SPK_fre_1,SPK_PSD_1);
    SPK_PSD_2 = Linear_Interpolate_Each_Site(SPK_fre_2,SPK_PSD_2);
    LFP_SPK_Coherence_1 = Linear_Interpolate_Each_Site(Fre,LFP_SPK_Coherence_1);
    LFP_SPK_Coherence_2 = Linear_Interpolate_Each_Site(Fre,LFP_SPK_Coherence_2);
end

% test plot
if 0
    fh = figure; for ii = 1:6, ax{ii} = subplot(3,2,ii,'parent',fh); end
    plot(ax{1}, t_1,RAW_01);    plot(ax{2}, t_2,RAW_02);
    plot(ax{3}, t_1,LFP_Raw_1); plot(ax{4}, t_2,LFP_Raw_2);
    plot(ax{5}, t_1,SPK_Raw_1); plot(ax{6},t_2,SPK_Raw_2);
    close all;
end
%% plot signal and PSD when nargin>1
if nargin>1
    FH = figure('name','LFP and SPK Original signal','units','normalized','outerposition', [0.01, 0.04, 0.99, 0.95]);
    T = tiledlayout(5,6,'parent',FH,'Padding','compact');%,'Padding','compact'
    ax = nexttile(3,[1,2]);  Plot_Signal_PSD(ax,varargin,t_1(4*Sample_Rate_1:6*Sample_Rate_1),t_2(4*Sample_Rate_2:6*Sample_Rate_2),RAW_01(4*Sample_Rate_1:6*Sample_Rate_1),RAW_02(4*Sample_Rate_2:6*Sample_Rate_2),[],'RAW Signal (0.07~9000Hz)');
    ax = nexttile(7,[1,2]);  Plot_Signal_PSD(ax,varargin,t_1(4*Sample_Rate_1:6*Sample_Rate_1),t_2(4*Sample_Rate_2:6*Sample_Rate_2),LFP_Raw_1(4*Sample_Rate_1:6*Sample_Rate_1),LFP_Raw_2(4*Sample_Rate_2:6*Sample_Rate_2),[],'LFP Signal (2~300Hz)');
    ax = nexttile(9,[1,2]);  Plot_Signal_PSD(ax,varargin,[],[],LFP_SPK_Coherence_1,LFP_SPK_Coherence_2,Fre,'The Choherence of LFP and SPK signals (3~200HZ)');
    ax = nexttile(11,[1,2]); Plot_Signal_PSD(ax,varargin,t_1(4*Sample_Rate_1:6*Sample_Rate_1),t_2(4*Sample_Rate_2:6*Sample_Rate_2),SPK_Raw_1(4*Sample_Rate_1:6*Sample_Rate_1),SPK_Raw_2(4*Sample_Rate_2:6*Sample_Rate_2),[],'SPK Signal (300~6000Hz)');
    ax = nexttile(17,[1,2]); Plot_Signal_PSD(ax,varargin,t_1(4*Sample_Rate_1:6*Sample_Rate_1),t_2(4*Sample_Rate_2:6*Sample_Rate_2),SPK_Abs_Mean_1(4*Sample_Rate_1:6*Sample_Rate_1),SPK_Abs_Mean_2(4*Sample_Rate_2:6*Sample_Rate_2),[],'SPK Signal with Absolute Operator');
    ax = nexttile(19,[1,2]); Plot_Signal_PSD(ax,varargin,[],[],LFP_PSD_1,LFP_PSD_2,f,'LFP PSD (3~200HZ)');
    ax = nexttile(23,[1,2]); Plot_Signal_PSD(ax,varargin,[],[],SPK_PSD_1,SPK_PSD_2,f,'SPK PSD (3~200Hz)');
end

%% output arguments
LFP_SPK_PSD_RMS.LFP_PSD_1 = LFP_PSD_1; LFP_SPK_PSD_RMS.LFP_fre_1 = LFP_fre_1; LFP_SPK_PSD_RMS.LFP_RMS_1 = LFP_RMS_1;
LFP_SPK_PSD_RMS.LFP_PSD_2 = LFP_PSD_2; LFP_SPK_PSD_RMS.LFP_fre_2 = LFP_fre_2; LFP_SPK_PSD_RMS.LFP_RMS_2 = LFP_RMS_2;
LFP_SPK_PSD_RMS.SPK_PSD_1 = SPK_PSD_1; LFP_SPK_PSD_RMS.SPK_fre_1 = SPK_fre_1; LFP_SPK_PSD_RMS.SPK_RMS_1 = SPK_RMS_1;
LFP_SPK_PSD_RMS.SPK_PSD_2 = SPK_PSD_2; LFP_SPK_PSD_RMS.SPK_fre_2 = SPK_fre_2; LFP_SPK_PSD_RMS.SPK_RMS_2 = SPK_RMS_2;

LFP_SPK_Coherence.LFP_SPK_Coherence_1 = LFP_SPK_Coherence_1;
LFP_SPK_Coherence.LFP_SPK_Coherence_2 = LFP_SPK_Coherence_2;
LFP_SPK_Coherence.Fre = Fre;

if nargout==3
    varargout{1} = FH;
end


if nargout==4
    LS_Signal.LFP_Signal_1 = LFP_Signal_1;  LS_Signal.LFP_Signal_2 = LFP_Signal_2;
    LS_Signal.SPK_Signal_1 = SPK_Signal_1;  LS_Signal.SPK_Signal_2 = SPK_Signal_2;
    Pwelch_Parameters.Win_1 = Win_1;  Pwelch_Parameters.Overlap_1 = Overlap_1; Pwelch_Parameters.Sample_Rate_1 = Sample_Rate_1;
    Pwelch_Parameters.Win_2 = Win_2;  Pwelch_Parameters.Overlap_2 = Overlap_2; Pwelch_Parameters.Sample_Rate_2 = Sample_Rate_2;
    Pwelch_Parameters.f = f;
    varargout{1} = LS_Signal;
    varargout{2} = Pwelch_Parameters;
end
if nargout==5
    LS_Signal.LFP_Signal_1 = LFP_Signal_1;  LS_Signal.LFP_Signal_2 = LFP_Signal_2;
    LS_Signal.SPK_Signal_1 = SPK_Signal_1;  LS_Signal.SPK_Signal_2 = SPK_Signal_2;
    Pwelch_Parameters.Win_1 = Win_1;  Pwelch_Parameters.Overlap_1 = Overlap_1; Pwelch_Parameters.Sample_Rate_1 = Sample_Rate_1;
    Pwelch_Parameters.Win_2 = Win_2;  Pwelch_Parameters.Overlap_2 = Overlap_2; Pwelch_Parameters.Sample_Rate_2 = Sample_Rate_2;
    Pwelch_Parameters.f = f;
    t.t_1 = t_1;  t.t_2 = t_2;
    varargout{1} = LS_Signal;
    varargout{2} = Pwelch_Parameters;
    varargout{3} = t;
end

%}
end

%% Subfunction: set parameters foe pwelch and filter
function [Window_Length_1,Window_Length_2,Overlap_Length_1,Overlap_Length_2,f] = Set_Pwech_Parameters(Sample_Rate_1,Sample_Rate_2)

fre_res = 1/2;
t_window = 1/fre_res; % s
overlap_rate = 0.5;
n = floor((200-3+fre_res)/fre_res)+1; % (200.5 - 3)/(n-1) = LFP_fre_res; %3
f = linspace(3, 200+fre_res, n); % we choose the frequency from 3 to 200 Hz, because there is a big noise below 3 HZ. That's why we start from 3 Hz.
f = f(1:end-1); % we delete the last point, because we want to avoid the artifacte noise from the last point.
Window_Length_1 = t_window * Sample_Rate_1;
Window_Length_2 = t_window * Sample_Rate_2;
Overlap_Length_1 = Window_Length_1 * overlap_rate;
Overlap_Length_2 = Window_Length_2 * overlap_rate;

end

%% Subfunction: DO filter and pwelch
function [PSD, fre, RMS, Fitered_Sig, Abs_Mean_Sig] = Do_Filter_Pwelch(RAW, Win, Overlap, f, Sample_Rate,Tag)
Fre_Reso = f(2)-f(1);
if strcmp(Tag,'LFP')
    Low_pass = f(1);% - 2*Fre_Reso; %2
    High_pass = f(end)*1.5; %300
elseif strcmp(Tag,'SPK')
    Low_pass = 300;
    High_pass = 6000;
end
if length(RAW)>12
    %[b1,a1] = butter(2, [Low_pass,High_pass]/(Sample_Rate/2), 'bandpass'); %4
    [z,p,k] = butter(2, [Low_pass,High_pass]/(Sample_Rate/2), 'bandpass'); % in fact, it's 4-pole
    [sos,g] = zp2sos(z,p,k);
    %fvtool(b1,a1); % this built-in function can be used to plot the filter
    %RAW = filtfilt(b1, a1, RAW);
    RAW = filtfilt(sos, g, RAW);
    Fitered_Sig = RAW;
    % if just do the rectification for SPK, keep this if-condition
    % if do the rectification both for the LFP and SPK, remove this
    % if-condition
    %if strcmp(Tag,'SPK')% Absolute the SPK, and remove the mean of the absolute SPK
        RAW = abs(RAW);
        RAW = RAW - mean(RAW); % it doesn't affect the 1/f slope
        Abs_Mean_Sig = RAW;
    %else
    %    Abs_Mean_Sig = [];
    %end
    if length(RAW) > 1.5*Win
        [PSD, fre] = pwelch(RAW, Win, Overlap, f, Sample_Rate);
        RMS = rms(RAW);
        if 0 % to plot the LFP and rectified LFP
            if strcmp(Tag,'LFP')
                Rectified_LFP = abs(RAW);
                Rectified_LFP = Rectified_LFP - mean(Rectified_LFP);
                [Rectified_PSD, fre] = pwelch(Rectified_LFP, Win, Overlap, f, Sample_Rate);
                figure; lg(1) = plot(fre, log10(PSD), 'r', 'displayname', 'LFP'); hold on; 
                lg(2) = plot(fre, log10(Rectified_PSD), 'g', 'displayname', 'Rectified LFP'); legend(lg);
                set(gca, 'xscale', 'log'); set(gca, 'xlim', [fre(1), fre(end)]);
                close all;
            end
        end
    else
        PSD = NaN(1, length(f));
        fre = NaN(1, length(f));
        RMS = NaN;
    end
else
    PSD = NaN(1, length(f));
    fre = NaN(1, length(f));
    RMS = NaN;
    Fitered_Sig = NaN;
    Abs_Mean_Sig = NaN;
end
end

%% Subfunction: plot raw signal and PSD
function  Plot_Signal_PSD(ax,Ipt_C,t_1,t_2,Signal_PSD_1,Signal_PSD_2,Fre,Tag)
if ~isempty(Fre)
    if strcmp(Ipt_C{1}(end), '1')
        if strcmp(Tag(1:3),'LFP')||strcmp(Tag(1:3),'SPK')
            plot(ax,Fre,log10(Signal_PSD_1),'color','k','LineWidth',2);
        else
            plot(ax,Fre,Signal_PSD_1,'color','k','LineWidth',2);
        end
    elseif strcmp(Ipt_C{1}(end), '2')
        if strcmp(Tag(1:3),'LFP')||strcmp(Tag(1:3),'SPK')
            plot(ax,Fre,log10(Signal_PSD_2),'color','k','LineWidth',2);
        else
            plot(ax,Fre,Signal_PSD_2,'color','k','LineWidth',2);
        end
    end
    if strcmp(Tag(1:3),'LFP')||strcmp(Tag(1:3),'SPK')
        xlim(ax,[Fre(1),Fre(end)]);
    else
        xlim(ax,[Fre(1),70]);
        ylim(ax,[0,0.8]);
    end
    
    xlabel(ax, 'Frequency(Hz)');
    if strcmp(Tag(1:3),'LFP')||strcmp(Tag(1:3),'SPK')
        ylabel(ax, 'log10(Power)(uV^2/Hz)');
    end
    xticks(ax, [4, 10, 20, 40, 100, 200]);
else
    if strcmp(Ipt_C{1}(end), '1')
        plot(ax, t_1, Signal_PSD_1,'color','k','LineWidth',2);
        xlim(ax,[t_1(1),t_1(end)]);
    elseif strcmp(Ipt_C{1}(end), '2')
        plot(ax, t_2, Signal_PSD_2,'color','k','LineWidth',2);
        xlim(ax,[t_2(1),t_2(end)]);
    end
    xlabel(ax, 'Time(s)');
    ylabel(ax, 'Power(uV)');
end
title(ax, Tag);
%{
if strcmp(Tag, 'SPK Signal')
    ylim(ax,[-150,150])
end

if strcmp(Ipt_C{2}, 'Pre')
    ylabel(ax, Tag);
end
if strcmp(Tag, 'LFP PSD')
    ylim(ax,[-2.5,3])
    yticks(ax, [-2, 0, 2]);
end
if strcmp(Tag, 'SPK PSD')
    ylim(ax,[-3.2,-1.85])
    yticks(ax, [-3,-2.5,-2]);
end
%if strcmp(Tag(end-2:end), 'PSD')
%    set(ax,'YScale','log');
%end
%}
end
