function FigH = Fooof_Plot(Patient_Fooof, Iden_Tag) %Patient_Fooof: single fooof results of every depth

if ~exist('Iden_Tag', 'var') || isempty(Iden_Tag)
    Iden_Tag = ' ';
end

Fre = Patient_Fooof.Fre;

Raw_PSD = Patient_Fooof.PSD;
Aperiodic_FirstFit_PSD = Patient_Fooof.Aperiodic_First_BetterFitted_PSD;
Aperiodic_ReFit_PSD = Patient_Fooof.Aperiodic_PeakRemoved_ReFitPSD;
Gn_Single_Fit_PSD = Patient_Fooof.Gn_Single_Fit_PSD;
Fooof_PSD = Patient_Fooof.Aperiodic_Gn_Fit_PSD;
Diff_PSD = Raw_PSD - Fooof_PSD;
Raw_Peak = Raw_PSD - Aperiodic_ReFit_PSD;
Fooof_Peak = Fooof_PSD - Aperiodic_ReFit_PSD;

Aperiodic_ReFit_Coefficients = Patient_Fooof.Aperiodic_PeakRemoved_ReFitCoefficients; %(offset, knee, exp)
Gn_Fit_Coefficients = Patient_Fooof.Gn_Fit_Coefficients; %(h, c, s)

Correlation_Coefficient = Patient_Fooof.Correlation_Coefficient;
Eorror_Metric = Patient_Fooof.Initial_Parameters.error_metric;
Eorror = Patient_Fooof.Eorror;

Gn_Num = size(Gn_Single_Fit_PSD, 1);

FigH = figure('name', Iden_Tag, 'units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);

%% Raw_PSD + Ap_Fit_PSD + Ap_Refit_PSD + Every Single Gn_fit_PSD + Fooof_PSD (Linear-Log)
subplot(2, 2, 1, 'parent', FigH);

Raw_Obj = plot(Fre, Raw_PSD, 'color', [0 0.4470 0.7410], 'displayname', 'Raw PSD');
hold on;
ApFit_Obj = plot(Fre, Aperiodic_FirstFit_PSD, 'color', [0.9290 0.6940 0.1250], 'displayname', 'Aperiodic Fit PSD');
ApReFit_Obj = plot(Fre, Aperiodic_ReFit_PSD, 'color', [0.6350 0.0780 0.1840], 'displayname', 'Aperiodic Refit PSD');

for ii = Gn_Num: -1: 1
    Color_percent = ii/Gn_Num;
    SingleGn_ApReFit_PSD = Gn_Single_Fit_PSD(ii, :) + Aperiodic_ReFit_PSD';
    SiglGn_Obj(ii) = plot(Fre, SingleGn_ApReFit_PSD, 'color', [1-Color_percent, 1, 1-Color_percent], 'displayname', ['Gn No.', num2str(ii)], 'linestyle', ':', 'linewidth',3);
end

Fooof_Obj = plot(Fre, Fooof_PSD, 'color', [0.4940 0.1840 0.5560], 'displayname', 'Fooof PSD');

xlabel('Frequency-Hz');
ylabel('Log10(PSD)');
legend([Raw_Obj, ApFit_Obj, ApReFit_Obj, Fooof_Obj]);
title('Linear-Log PSD in 5 stages');

%% Diff_ PSD between raw and fooof PSD + Raw_Peak + Fooof_Peak + Every Single Gn_fit_PSD
subplot(2,2,2, 'parent', FigH);

plot(Fre([1, end]), [0, 0], 'color', [0.5, 0.5,0.5]);
hold on;
RawFooof_Obj = plot(Fre, Diff_PSD, 'color', [0 0.4470 0.7410], 'displayname', 'Raw-Fooof Residual PSD');
RawPeak_Obj = plot(Fre, Raw_Peak, 'color', [0.8500 0.3250 0.0980], 'displayname', 'Raw-ApRefit Peak PSD');
FooofPeak_Obj = plot(Fre, Fooof_Peak, 'color', [0.4940 0.1840 0.5560], 'displayname', 'Fooof-ApRefit Peak PSD');

for ii = Gn_Num: -1: 1
    Color_percent = ii/Gn_Num;
    SiglGn_Obj(ii) = plot(Fre, Gn_Single_Fit_PSD(ii, :), 'color', [1-Color_percent, 1, 1-Color_percent], 'displayname', ['Gn No.', num2str(ii)], 'linestyle', ':', 'linewidth', 2, 'marker', '*');
end

xlabel('Frequency-Hz');
ylabel('Log10(PSD)');
legend([RawFooof_Obj, RawPeak_Obj, FooofPeak_Obj]);
title('Different periodic PSD');

%% Raw_PSD + Ap_Fit_PSD + Ap_Refit_PSD + Every Single Gn_fit_PSD + Fooof_PSD (Log-Log)
subplot(2, 2, 3, 'parent', FigH);

Raw_Obj = plot(Fre, Raw_PSD, 'color', [0 0.4470 0.7410], 'displayname', 'Raw PSD');
hold on;
ApFit_Obj = plot(Fre, Aperiodic_FirstFit_PSD, 'color', [0.9290 0.6940 0.1250], 'displayname', 'Aperiodic Fit PSD');
ApReFit_Obj = plot(Fre, Aperiodic_ReFit_PSD, 'color', [0.6350 0.0780 0.1840], 'displayname', 'Aperiodic Refit PSD');

for ii = Gn_Num: -1: 1
    Color_percent = ii/Gn_Num;
    SingleGn_ApReFit_PSD = Gn_Single_Fit_PSD(ii, :) + Aperiodic_ReFit_PSD';
    SiglGn_Obj(ii) = plot(Fre, SingleGn_ApReFit_PSD, 'color', [1-Color_percent, 1, 1-Color_percent], 'displayname', ['Gn No.', num2str(ii)], 'linestyle', ':', 'linewidth', 3);
end

Fooof_Obj = plot(Fre, Fooof_PSD, 'color', [0.4940 0.1840 0.5560], 'displayname', 'Fooof PSD');

xlabel('Log(Frequency-Hz)');
ylabel('Log10(PSD)');
if ~exist('SiglGn_Obj', 'var')
    legend('There is not Gaussian', 'location', 'southwest');
else
    legend(SiglGn_Obj, 'location', 'southwest');
end
set(gca, 'xscale', 'log');
title('Log-Log PSD in 5 stages');

%% Ap_coefficients + R^2 + Error + Gn_NO. +Gn_Coefficiens
subplot(2, 2, 4, 'parent', FigH);

Aperiodic_ReFit_Coefficients = round(Aperiodic_ReFit_Coefficients.*1000)./1000;
Gn_Fit_Coefficients = round(Gn_Fit_Coefficients.*1000)./1000;

if length(Aperiodic_ReFit_Coefficients) == 2
    Ap_Text = [ 'Aperiodic Coefficients: offset: ', num2str(Aperiodic_ReFit_Coefficients(1)), 10, '                                       exp: ', num2str(Aperiodic_ReFit_Coefficients(2)) ];
elseif length(Aperiodic_ReFit_Coefficients) == 3
    Ap_Text = [ 'Aperiodic Coefficients: offset: ', num2str(Aperiodic_ReFit_Coefficients(1)), 10, '                                       knee: ', num2str(Aperiodic_ReFit_Coefficients(2)), 10, '                                       exp: ', num2str(Aperiodic_ReFit_Coefficients(3)) ];
end

Fooof_Text = ['Goodness of Fooof: R^2: ', num2str(Correlation_Coefficient), 10, '                                 ', Eorror_Metric, ': ', num2str(Eorror) ];

Gn_Text = 'Gaussian';
Gn_Hight_Text = 'Hight';
Gn_Mean_Text = 'Center-Fre';
Gn_Std_Text = 'Std';
for ii = 1: Gn_Num
    Gn_Text = [Gn_Text, 10, 'NO.', num2str(ii)];
    Gn_Hight_Text = [Gn_Hight_Text, 10, num2str(Gn_Fit_Coefficients(ii, 1))];
    Gn_Mean_Text = [Gn_Mean_Text, 10, num2str(Gn_Fit_Coefficients(ii, 2))];
    Gn_Std_Text = [Gn_Std_Text, 10, num2str(Gn_Fit_Coefficients(ii, 3))];
end

text(0.1, 0.75, Ap_Text);
text(0.1, 0.25, Fooof_Text);
text(0.55, 0.5, Gn_Text);
text(0.65, 0.5, Gn_Hight_Text);
text(0.75, 0.5, Gn_Mean_Text);
text(0.85, 0.5, Gn_Std_Text);


end

