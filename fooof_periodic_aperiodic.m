function PSD_cpn = fooof_periodic_aperiodic( f, psd, plotflag, PSD_cpn, showdetail )
% FOOOF in        https://fooof-tools.github.io/fooof/index.html
% modified from   https://github.com/fooof-tools/fooof/
% tutorial        https://fooof-tools.github.io/fooof/auto_tutorials/index.html

if ~exist('plotflag','var')
    plotflag = true;
elseif isempty(plotflag)
    plotflag = false;
end
if ~exist('PSD_cpn','var') || isempty(PSD_cpn)
    % iniPSD parameters
    peak_width_limits = [1, 15.0];
    peak_width_limits_pro = [0.02, 0]; % used the max one with peak_width_limits % GJ
    max_n_peaks = Inf;
    min_peak_height = 0;
    peak_threshold = 2;
    aperiodic_mode = 'fixed'; % fixed knee
    verbose = true;
    PSD_cpn = fooof_ini( peak_width_limits, peak_width_limits_pro, max_n_peaks, ...
        min_peak_height, peak_threshold, aperiodic_mode, verbose);
else  %% for get input from    fooof_SPG_periodic_aperiodic
end
if ~exist('showdetail','var')
    showdetail = false;
end

% bug from  interplot  %% changed to linear insert in log-space
% psd(psd<=0) = 0.001;

% fooof.objs.fit -> add_data -> _prepare_data
if size(f,1)==1, f=f'; end
if size(psd,1)==1, psd=psd'; end
PSD_cpn.freqs = f;
PSD_cpn.power_spectrum = log10( psd );
PSD_cpn.freq_range = [min(f),max(f)];
PSD_cpn.freq_res = f(2)-f(1);

% fooof.objs.fit -> _check_width_limits
if PSD_cpn.verbose && 1.5 * PSD_cpn.freq_res>=PSD_cpn.peak_width_limits(1)
    msgbox(['low low-bounds limits comparing with freq resolution.',10, ...
        'recommend a lower bound of approximately 2x the frequency resolution.']);
end

if any(isnan(psd))
    psdNaN = nan(size(psd));
    PSD_cpn.aperiodic_params = [NaN,NaN];
    if strcmp(PSD_cpn.aperiodic_mode, 'knee'), PSD_cpn.aperiodic_params = [NaN,NaN,NaN]; end
    PSD_cpn.ap_fit_1st_posi = false(size(psd));
    PSD_cpn.ap_fit = psdNaN;
    PSD_cpn.ap_fit_1st = psdNaN;
    PSD_cpn.spectrum_flat = psdNaN;
    PSD_cpn.gaussian_params = zeros([0,3]);
    if showdetail
        PSD_cpn.gaussian_guess = zeros([0,3]);
        PSD_cpn.gaussian_guess_raw = zeros([0,3]);
    end
    PSD_cpn.peak_fit = psdNaN;
    PSD_cpn.spectrum_peak_rm = psdNaN;
    PSD_cpn.fooofed_spectrum = psdNaN;
    PSD_cpn.peak_params = zeros([0,3]);
    PSD_cpn.GNum = 0;
    PSD_cpn.r_squared = nan;
    PSD_cpn.error = nan;
    return;
end

%% Fit the aperiodic component
[PSD_cpn.aperiodic_params, PSD_cpn.ap_fit_1st_posi] = robust_ap_fit(PSD_cpn);
PSD_cpn.ap_fit = gen_aperiodic(PSD_cpn.freqs, PSD_cpn.aperiodic_params);
PSD_cpn.ap_fit_1st = PSD_cpn.ap_fit;

% Flatten the power spectrum using fit aperiodic fit
PSD_cpn.spectrum_flat = PSD_cpn.power_spectrum - PSD_cpn.ap_fit;

%% Find peaks, and fit them with gaussians
if showdetail
    [PSD_cpn.gaussian_params, PSD_cpn.gaussian_guess, PSD_cpn.gaussian_guess_raw] = fit_peaks( PSD_cpn );
else
    PSD_cpn.gaussian_params = fit_peaks( PSD_cpn );
end

% Calculate the peak fit
%   Note: if no peaks are found, this creates a flat (all zero) peak fit
PSD_cpn.peak_fit = gen_periodic( PSD_cpn.freqs, PSD_cpn.gaussian_params );

% Create peak-removed (but not flattened) power spectrum
PSD_cpn.spectrum_peak_rm = PSD_cpn.power_spectrum - PSD_cpn.peak_fit;

%% Run final aperiodic fit on peak-removed power spectrum
%   This overwrites previous aperiodic fit, and recomputes the flattened spectrum
PSD_cpn.aperiodic_params = simple_ap_fit(PSD_cpn, PSD_cpn.spectrum_peak_rm);
PSD_cpn.ap_fit = gen_aperiodic(PSD_cpn.freqs, PSD_cpn.aperiodic_params);
PSD_cpn.spectrum_flat = PSD_cpn.power_spectrum - PSD_cpn.ap_fit;

% Create full power_spectrum model fit
PSD_cpn.fooofed_spectrum = PSD_cpn.peak_fit + PSD_cpn.ap_fit;

% Convert gaussian definitions to peak parameters
PSD_cpn.peak_params = create_peak_params( PSD_cpn );
PSD_cpn.GNum = size(PSD_cpn.peak_params,1);

%% Calculate R^2 and error of the model fit
PSD_cpn = calc_r_squared( PSD_cpn );
PSD_cpn = calc_error( PSD_cpn );

%% draw
if plotflag
    fooof_plot_PSD_cpn( PSD_cpn );
end

end

function aperiodic_params = simple_ap_fit( PSD_cpn, power_spectrum )
%======  fooof.objs.fit -> _simple_ap_fit  ======
% Get the guess parameters and/or calculate from the data, as needed
%   Note that these are collected as lists, to concatenate with or without knee later

freqs = PSD_cpn.freqs;

if ~isnan(PSD_cpn.ap_guess(1))
    off_guess = PSD_cpn.ap_guess(1);
else
    off_guess = power_spectrum(1);
end

kne_guess = [];
if strcmp(PSD_cpn.aperiodic_mode,'knee'), kne_guess = PSD_cpn.ap_guess(2); end

if ~isnan(PSD_cpn.ap_guess(3))
    exp_guess = PSD_cpn.ap_guess(3);
else
    exp_guess = abs( (PSD_cpn.power_spectrum(end)-PSD_cpn.power_spectrum(1)) / ...
        (log10(PSD_cpn.freqs(end))-log10(PSD_cpn.freqs(1))) );
end

% Get bounds for aperiodic fitting, dropping knee bound if not set to fit knee
if strcmp(PSD_cpn.aperiodic_mode, 'knee')
    ap_bounds = PSD_cpn.ap_bounds;
else
    ap_bounds = PSD_cpn.ap_bounds(:,1:2:end);
end

% Collect together guess parameters
guess = [off_guess, kne_guess, exp_guess];

% Ignore warnings that are raised in curve_fit / fit
%   A runtime warning can occur while exploring parameters in curve fitting
%     This doesn't effect outcome - it won't settle on an answer that does this
%   It happens if / when b < 0 & |b| > x**2, as it leads to log of a negative number

%{
aperiodic_params, _ = curve_fit(get_ap_func(self.aperiodic_mode),
                                freqs, power_spectrum, p0=guess,
                                maxfev=self._maxfev, bounds=ap_bounds)
%}
aperiodic_params = fit( freqs, power_spectrum, get_ap_func(PSD_cpn.aperiodic_mode), ...
    'StartPoint', guess, 'MaxFunEvals', PSD_cpn.maxfev, ...
    'Lower', ap_bounds(1,:), 'Upper', ap_bounds(2,:) );
if length(guess) == 2
    aperiodic_params = [aperiodic_params.off, aperiodic_params.ep];
elseif length(guess) ==3
    aperiodic_params = [aperiodic_params.off, aperiodic_params.knee, aperiodic_params.ep];
end
end

function [aperiodic_params,perc_mask] = robust_ap_fit(PSD_cpn)
%======  fooof.objs.fit -> robust_ap_fit  ======
% freqs : 1d array
% Frequency values for the power spectrum, in linear scale.
% power_spectrum : 1d array
% Power values, in log10 scale.

% aperiodic_params : 1d array
% Parameter estimates for aperiodic fit.

freqs = PSD_cpn.freqs;
power_spectrum = PSD_cpn.power_spectrum;

% Do a quick, initial aperiodic fit
popt = simple_ap_fit(PSD_cpn,power_spectrum);
initial_fit = gen_aperiodic(freqs, popt);

% Flatten power_spectrum based on initial aperiodic fit
flatspec = power_spectrum - initial_fit;

% Flatten outliers, defined as any points that drop below 0
flatspec(flatspec < 0) = 0;

% Use percentile threshold, in terms of % of points, to extract and re-fit
perc_thresh = prctile(flatspec, PSD_cpn.ap_percentile_thresh); % np.percentile
perc_mask = flatspec <= perc_thresh;
freqs_ignore = freqs(perc_mask);
spectrum_ignore = power_spectrum(perc_mask);

% Get bounds for aperiodic fitting, dropping knee bound if not set to fit knee
if strcmp(PSD_cpn.aperiodic_mode, 'knee')
    ap_bounds = PSD_cpn.ap_bounds;
else
    ap_bounds = PSD_cpn.ap_bounds(:,1:2:end);
end

% Second aperiodic fit - using results of first fit as guess parameters
%  See note in simple_ap_fit about warnings
%{
aperiodic_params, _ = curve_fit(get_ap_func(self.aperiodic_mode),
                                freqs_ignore, spectrum_ignore, p0=popt,
                                maxfev=self._maxfev, bounds=ap_bounds)
%}
aperiodic_params = fit( freqs_ignore, spectrum_ignore, get_ap_func(PSD_cpn.aperiodic_mode), ...
    'StartPoint', popt, 'MaxFunEvals', PSD_cpn.maxfev, ...
    'Lower', ap_bounds(1,:), 'Upper', ap_bounds(2,:) );
if length(popt) == 2
    aperiodic_params = [aperiodic_params.off, aperiodic_params.ep];
elseif length(popt) ==3
    aperiodic_params = [aperiodic_params.off, aperiodic_params.knee, aperiodic_params.ep];
end
end

function [gaussian_params, guess, guess_raw] = fit_peaks(PSD_cpn)
%======  fooof.objs.fit -> _fit_peaks  ======
% flat_iter : 1d array
% Flattened power spectrum values.

% gaussian_params : 2d array
% Parameters that define the gaussian fit(s).
% Each row is a gaussian, as [mean, height, standard deviation].

flat_iter = PSD_cpn.spectrum_flat;

% Initialize matrix of guess parameters for gaussian fitting
guess = zeros([0, 3]);

% Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian
%   Stopping procedures: limit on % of peaks, or relative or absolute height thresholds
while size(guess,1) < PSD_cpn.max_n_peaks
    % Find candidate peak - the maximum point of the flattened spectrum
    [~,max_ind] = max(flat_iter);
    max_height = flat_iter(max_ind);
    
    % Stop searching for peaks once height drops below height threshold
    if max_height <= PSD_cpn.peak_threshold * std(flat_iter), break; end
    
    % Set the guess parameters for gaussian fitting, specifying the mean and height
    guess_freq = PSD_cpn.freqs(max_ind);
    guess_height = max_height;
    
    % Halt fitting process if candidate peak drops below minimum height
    if guess_height <= PSD_cpn.min_peak_height, break; end
    
    % Data-driven first guess at standard deviation
    %   Find half height index on each side of the center frequency
    half_height = 0.5 * max_height;
    le_ind = max_ind - find( flat_iter(max_ind-1:-1:1)<=half_height, 1, 'first' );
    ri_ind = find( flat_iter(max_ind+1:end)<=half_height, 1, 'first' ) + max_ind;
    
    % Guess bandwidth procedure: estimate the width of the peak
    % Get an estimated width from the shortest side of the peak
    %   We grab shortest to avoid estimating very large values from overlapping peaks
    % Grab the shortest side, ignoring a side if the half max was not found
    short_side = min(abs( [le_ind,ri_ind]-max_ind ));
    
    % Use the shortest side to estimate full-width, half max (converted to Hz)
    %   and use this to estimate that guess for gaussian standard deviation
    fwhm = short_side * 2 * PSD_cpn.freq_res;
    guess_std = fwhm / (2 * sqrt( 2*log(2) )); % compute_gauss_std
    
    % Check that guess value isn't outside preset limits - restrict if so
    %   Note: without this, curve_fitting fails if given guess > or < bounds
    crtstdlim = max( [PSD_cpn.gauss_std_limits; PSD_cpn.gauss_std_limits_pro.*guess_freq] );
    guess_std( guess_std<crtstdlim(1) ) = crtstdlim(1);
    guess_std( guess_std>crtstdlim(2) ) = crtstdlim(2);
    
    % Collect guess parameters and subtract this guess gaussian from the data
    guess = [guess; guess_freq, guess_height, guess_std]; %#ok<AGROW>
    peak_gauss = fooof_gaussian_function(guess_freq, guess_height, guess_std, PSD_cpn.freqs);
    flat_iter = flat_iter - peak_gauss;
end

% Check peaks based on edges, and on overlap, dropping any that violate requirements
guess_raw = guess;
guess = drop_peak_cf(PSD_cpn,guess);
guess = drop_peak_overlap(PSD_cpn,guess);

% If there are peak guesses, fit the peaks, and sort results
if size(guess,1) > 0
    gaussian_params = fit_peak_guess( PSD_cpn, guess );
    [~,fIX] = sort(gaussian_params(:,1));
    gaussian_params = gaussian_params( fIX,: );
else
    gaussian_params = zeros([0, 3]);
end
end

function gaussian_params = fit_peak_guess(PSD_cpn, guess)
%======  fooof.objs.fit -> _fit_peak_guess  ======
% guess,gaussian_params : 2d array, shape=[n_peaks, 3]
% as gaussian parameters.

% Set the bounds for CF, enforce positive height value, and set bandwidth limits
%   Note that 'guess' is in terms of gaussian std, so +/- BW is 2 * the guess_gauss_std
%   This set of list comprehensions is a way to end up with bounds in the form:
%     ((cf_low_peak1, height_low_peak1, bw_low_peak1, *repeated for n_peaks*),
%      (cf_high_peak1, height_high_peak1, bw_high_peak, *repeated for n_peaks*))
%     ^where each value sets the bound on the specified parameter

GNum = size(guess,1);

lo_bound = [guess(:,1) - 2 * PSD_cpn.cf_bound * guess(:,3), repmat([0,PSD_cpn.gauss_std_limits(:,1)],[GNum,1])];
hi_bound = [guess(:,1) + 2 * PSD_cpn.cf_bound * guess(:,3), repmat([Inf,PSD_cpn.gauss_std_limits(:,2)],[GNum,1])];

% Check that CF bounds are within frequency range
%   If they are  not, update them to be restricted to frequency range
lo_bound( lo_bound(:,1)<PSD_cpn.freq_range(1),1 ) = PSD_cpn.freq_range(1);
hi_bound( hi_bound(:,1)>PSD_cpn.freq_range(2),1 ) = PSD_cpn.freq_range(2);

% modify std by center frequency % add GJ   % likely not useful
stdlim = repmat(PSD_cpn.gauss_std_limits_pro,[GNum,1]).*repmat(guess(:,1),[1,2]);
lop = lo_bound(:,3)<stdlim(:,1);
lo_bound(lop,3) = stdlim(lop,1);
hip = hi_bound(:,3)<stdlim(:,2);
hi_bound(hip,3) = stdlim(hip,2);

% Unpacks the embedded lists into flat tuples
%   This is what the fit function requires as input
gaus_param_bounds = [reshape(lo_bound',1,[]);reshape(hi_bound',1,[])];

% Flatten guess, for use with curve fit
guess = reshape(guess',1,[]);

% Fit the peaks
%{
gaussian_params, _ = curve_fit(gaussian_function, self.freqs, self._spectrum_flat,
                                p0=guess, maxfev=self._maxfev, bounds=gaus_param_bounds)
%}
cG_func = [];
varStr = [];
varfitStr = [];
for ii = 1:GNum
    cnS = num2str(ii);
    varStr = [varStr,'c',cnS,',h',cnS,',s',cnS,',']; %#ok<AGROW>
    varfitStr = [varfitStr,'mGfit.c',cnS,',mGfit.h',cnS,',mGfit.s',cnS,',']; %#ok<AGROW>
end
eval(['cG_func = fittype( @(',varStr,'x) fooof_gaussian_function(',varStr,'x) );']);
mGfit = fit( PSD_cpn.freqs, PSD_cpn.spectrum_flat, cG_func, ...
    'StartPoint', guess, 'MaxFunEvals', PSD_cpn.maxfev, ...
    'Lower', gaus_param_bounds(1,:), 'Upper', gaus_param_bounds(2,:) );

% Re-organize params into 2d matrix
gaussian_params = mGfit; % not actural useful, for cancel matlab warning
eval(['gaussian_params = [',varfitStr(1:end-1),'];']);
gaussian_params = reshape(gaussian_params,3,[])';
end

function peak_params = create_peak_params( PSD_cpn )
%======  fooof.objs.fit -> _create_peak_params  ======
% gaus_params : 2d array
% Parameters that define the gaussian fit(s), as gaussian parameters.

% peak_params : 2d array
% Fitted parameter values for the peaks, with each row as [CF, PW, BW].

% The gaussian height is updated to reflect the height of the peak above
%       the aperiodic fit. This is returned instead of the gaussian height, as
%       the gaussian height is harder to interpret, due to peak overlaps.

gaus_params = PSD_cpn.gaussian_params;
GNum = size(gaus_params,1);

peak_params = nan([GNum, 3]);
for ii = 1:GNum
    % Gets the index of the power_spectrum at the frequency closest to the CF of the peak
    [~,ind] = min( abs(PSD_cpn.freqs - gaus_params(ii,1)) );
    
    % Collect peak parameter data
    peak_params(ii,1:3) = [ ...
        gaus_params(ii,1), PSD_cpn.fooofed_spectrum(ind)-PSD_cpn.ap_fit(ind), gaus_params(ii,3)*2];
end
end

function ap_vals = gen_aperiodic(freqs, ap_paras, aperiodic_mode)
%======  fooof.sim.gen -> gen_aperiodic  ======
% freqs : 1d array
%     Frequency vector to create aperiodic component for.
% aperiodic_params : list of float
%     Parameters that define the aperiodic component.
% aperiodic_mode : {'fixed', 'knee'}, optional
%     Which kind of aperiodic component to generate.
%     If not provided, is inferred from the parameters.

% ap_vals : 1d array
%     Aperiodic values, in log10 spacing.
if ~exist('aperiodic_mode','var')
    if length(ap_paras) == 2
        aperiodic_mode = 'fixed';
    elseif length(ap_paras) == 3
        aperiodic_mode = 'knee';
    end
end
ap_func = get_ap_func( aperiodic_mode );
if length(ap_paras) == 2
    ap_vals = ap_func( ap_paras(1), ap_paras(2), freqs );
elseif length(ap_paras) == 3
    ap_vals = ap_func( ap_paras(1), ap_paras(2), ap_paras(3), freqs );
end
end
function ap_func = get_ap_func( aperiodic_mode )
%======  fooof.core.funcs -> get_ap_func  ======
if strcmp( aperiodic_mode,'fixed')
    ap_func = fittype( @(off,ep,x) expo_nk_function(off,ep,x) );
elseif strcmp( aperiodic_mode,'knee')
    ap_func = fittype( @(off,knee,ep,x) expo_function(off,knee,ep,x) );
end
end

function pe_vals = gen_periodic(freqs, periodic_params, periodic_mode)
%======  fooof.sim.gen -> gen_periodic  ======
% freqs : 1d array
%     Frequency vector to create peak values for.
% periodic_params : list of float
%     Parameters to create the periodic component.
% periodic_mode : {'gaussian'}, optional
%     Which kind of periodic component to generate.

% peak_vals : 1d array
%     Peak values, in log10 spacing.
if ~exist('periodic_mode','var')
    periodic_mode = 'gaussian';
end
if strcmp( periodic_mode,'gaussian'), pe_func = @fooof_gaussian_function; end
pe_vals = pe_func(freqs, periodic_params);
end

function y = expo_nk_function(varargin)
x = varargin{end}; offset=varargin{1}; ep=varargin{2};
y = zeros(size(x));
y = y + offset-log10(x.^ep); % expo_nk_function;
end
function y = expo_function(varargin)
x = varargin{end}; offset=varargin{1}; knee=varargin{2}; ep=varargin{3};
y = zeros(size(x));
y = y + offset-log10(knee+x.^ep); % expo_function;
end

function guess = drop_peak_cf(PSD_cpn, guess)
%======  fooof.objs.fit -> _drop_peak_cf  ======
cf_params = guess(:,1);
bw_params = guess(:,3) .* PSD_cpn.bw_std_edge;

% Check if peaks within drop threshold from the edge of the frequency range
keepflag = abs(cf_params-PSD_cpn.freq_range(1))>bw_params & ...
    abs(cf_params-PSD_cpn.freq_range(end))>bw_params;

% Drop peaks that fail the center frequency edge criterion
guess = guess(keepflag,:);
end

function guess = drop_peak_overlap(PSD_cpn, guess)
%======  fooof.objs.fit -> _drop_peak_overlap  ======
% Checks whether to drop gaussians based on amount of overlap.
if size(guess,1)<=1, return; end

% Sort the peak guesses by increasing frequency
%   This is so adjacent peaks can be compared from right to left
[~,fIX] = sort(guess(:,1));
guess = guess(fIX,:);

% Calculate standard deviation bounds for checking amount of overlap
%   The bounds are the gaussian frequency +/- gaussian standard deviation
lo_bounds = guess(:,1) - guess(:,3) .* PSD_cpn.gauss_overlap_thresh;
hi_bounds = guess(:,1) + guess(:,3) .* PSD_cpn.gauss_overlap_thresh;

% Loop through peak bounds, comparing current bound to that of next peak
%   If the left peak's upper bound extends pass the right peaks lower bound,
%   then drop the Gaussian with the lower height
% comparing to all the right peak?? (haven't add) % GJ
%   (in case peak A,B,C. peak C cover A,B. comparing all to cancel A) % GJ
dropflag = false(size(lo_bounds));
for ii = 1:length(dropflag)-1
    if hi_bounds(ii)>=lo_bounds(ii+1) % overlap
        dropflag( ii + (guess(ii,2)>guess(ii+1,2)) ) = true;
    end
end

% Drop any peaks guesses that overlap too much, based on threshold
guess = guess(~dropflag,:);

end


function PSD_cpn = calc_r_squared(PSD_cpn)
%======  fooof.objs.fit -> _calc_r_squared  ======
% Calculate the r-squared goodness of fit of the model, compared to the original data.
r_val = corrcoef( PSD_cpn.power_spectrum, PSD_cpn.fooofed_spectrum );
PSD_cpn.r_squared = r_val(1,2) .^ 2;
end

function PSD_cpn = calc_error(PSD_cpn, metric)
%======  fooof.objs.fit -> _calc_error  ======
% Calculate the overall error of the model fit, compared to the original data.
% metric : {'MAE', 'MSE', 'RMSE'}, optional

% If metric is not specified, use the default approach
if ~exist('metric','var') || isempty(metric)
    metric = PSD_cpn.error_metric;
end

diff = PSD_cpn.power_spectrum-PSD_cpn.fooofed_spectrum;
switch metric
    case 'MAE'
        PSD_cpn.error = nanmean( abs(diff) );
    case 'MSE'
        PSD_cpn.error = nanmean( (diff.^2) );
    case 'RMSE'
        PSD_cpn.error = sqrt( nanmean(diff.^2) );
    otherwise
        error('unkonw method for calc_error');
end
end






