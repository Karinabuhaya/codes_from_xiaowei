function Parameters = Fooof_Parameters(aperiodic_mode, peak_width_limits, peak_width_limits_Per,...
    max_n_peaks, min_peak_height, peak_threshold)

% aperiodic_mode: Whether to fit the aperiodic component with a knee  parameter ( 'fixed'[default], 'knee' ).
% aperiodic_percentile_threshold: Power percentile to select points to fit aperiodic component, unit: 2.5%.

% aperiodic_offset_knee_exponent = [offset, knee, exponent]: Seed values for aperiodic fitting, 
%     if offset is none, the first value of PSD is used as the offset guess, means offset_guess = PSD(1),
%     knee is decided by aperiodic_mode:  
%     if aperiodic_mode = 'fixed', knee_guess = [], means aperiodic_offset_knee_exponent = [offset, exponent];
%     if speriodic_mode = 'knee', knee_guess = 0, means aperiodic_offset_knee_exponent = [offset, knee, exponent], '0' is the initial value.    
%     if exponent is none, the abs(log-log slop) of the first and last point is used as the exponent guess, means exponent_ guess = abs( (PSD(end) - PSD(1)) / (log10(Fre(end)) - log10(Fre(1))) ). 

% aperiodic_fit_bounds: bound on offset, knee and exponent for aperiodic fit, aperiodic_fit_bounds is decided by aperiodic_mode:
%     if aperiodic_mode = 'fixed', knee_low_bound = [], knee_high_bound = [], means aperiodic_fit_bounds = [-Inf, Inf; -Inf, Inf];
%     if speriodic_mode = 'knee', knee_low_bound = -Inf, knee_high_bound = Inf, means aperiodic_fit_bounds = [-Inf, Inf; -Inf, Inf; -Inf, Inf]

% peak_width_limits: decide the limits of the peak bandwidth based on the frequency resolution, ... 
% means the minimum limit of the peak bandwidth and the maxmum limit of the peak bandwidth,  ... 
% the limits decide the accuracy of the Gaussian form.  unit: Hz, default: [0.5, 12].

% max_n_peaks: Maximum number of peaks to fit, default: Inf.
% min_peak_height: Absolute minimum power threshold for peaks, unit:  log10(PSD), default: 0.
% peak_threshold: Relative minimum power threshold for peaks, unit: std,  default: 2.

% bandwith_edge_threshold: Threshold for dropping peaks close to the edge ( <= 1.0std ); 
% means the distance between center frequency and two side edges should be larger than 1*std, otherwise, this peak should be dropped.

% gaussian_overlap_threshold: Threshold for dropping overlapping peaks ( <= 0.75std ).
% means the distance between the center frequency of two peaks should be larger than 0.75*std, otherwise, one peak should be dropped.

% gaussian_fit_bounds: Bound on centre-frequency for multi-gaussian fit.

% maxfunevals: Maximum number of evaluations of model allowed.
% error_metric: estimate the goodness of fit, the smaller of the value, the better of the fit, which is incontrast with R^2.

% the following values of the five parameters can be changed or be setted by ourseles
if ~exist('aperiodic_mode', 'var'),         aperiodic_mode = 'fixed';             end
if ~exist('peak_width_limits', 'var'),      peak_width_limits = [0.5, 12];        end
if ~exist('max_n_peaks', 'var'),              max_n_peaks = Inf;                        end
if ~exist('min_peak_height', 'var'),       min_peak_height = 0;                     end
if ~exist('peak_threshold', 'var'),          peak_threshold = 2;                        end

%% Store the values of those variables in the output struct (Parameters)

% used for the aperiodic fit
Parameters.aperiodic_mode = aperiodic_mode; % the percentile threshold value can be adjusted if needed, but in practice rarey needs to be
Parameters.aperiodic_percentile_threshold = 0.025; % can not be setted
Parameters.aperiodic_offset_knee_exponent = [NaN, 0, NaN]; % [offset, knee, exponent], can not be setted, get from the data or fit result
Parameters.aperiodic_fit_bounds = [-Inf, Inf; -Inf, Inf; -Inf, Inf];%[-Inf, Inf; -Inf, Inf; -Inf, Inf]; % [offset_low_bound, offset_high_bound; knee_low_bound, knee_high_bound; exponent_low_bound,   exponent_high_bound]


% used for the periodic fit (Gaussian)
Parameters.peak_width_limits = peak_width_limits;
Parameters.std_limits = peak_width_limits;
Parameters.std_limits_Per = peak_width_limits_Per;
Parameters.max_n_peaks = max_n_peaks;
Parameters.min_peak_height = min_peak_height;
Parameters.peak_threshold = peak_threshold;
Parameters.bandwith_edge_threshold = 1.0; % unit: std
Parameters.gaussian_overlap_threshold = 0.75; % unit: std
Parameters.gaussian_fit_bounds = 1.5; % unit: std [low_bound, high_bound] = [centre_Fre - 1.5std, centre_Fre + 1.5std]

% used for curve fit
Parameters.maxfunevals = 5000;
Parameters.error_metric = 'MAE'; % mae: Mean absolute error performance function, mae together with R^2 to estimate the goodness of fit

end