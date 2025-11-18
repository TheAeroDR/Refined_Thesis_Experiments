%filename constructor
position = 1;
particle = "poly";
airflow = "air";
magdir = "z";

filename = "pos_" + num2str(position) + "_" + particle + "_" + airflow + "_mag" + magdir + ".h5";

filename = "test_ni_output.h5";

a = h5totable(filename);

Fs = h5readatt(filename,'/','sampling_rate');

%%
%electrometer
test_b = movmean(a.BE1,Fs/10);
test_t = movmean(a.BE4,Fs/10);
figure(1)
plot(a.Time,test_t - test_b)
xlabel('Time [s]')
ylabel('$\Delta$Voltage [V]')
%%
%magnetic spectrogram
figure(2)
tiledlayout(3,1)
nexttile(1)
b = spectrogram_mag(a.IFG3, Fs, 'reduced', false, true);
plotspectrogram(b,100)
xticklabels([])
xlabel([])
nexttile(2)
b = spectrogram_mag(a.IFG2, Fs, 'reduced', false, true);
plotspectrogram(b,100)
xticklabels([])
xlabel([])
nexttile(3)
b = spectrogram_mag(a.IFG1, Fs, 'reduced', false, true);
plotspectrogram(b,100)
%%
%anemometer
TIn_zw = 1.2185;
Vert_zw = 1.1848;

TInW = movmean(a.ANTI_wind_speed,floor(Fs/2));
TInT = movmean(a.ANTI_temperature,floor(Fs/2));
VertW = movmean(a.ANV_wind_speed,floor(Fs/2));
VertT = movmean(a.ANV_temperature,floor(Fs/2));
T_Temp = convert_to_temperature(TInT);
T_Winds = convert_to_winsdpeed(TInW,TIn_zw,T_Temp,'ms');
V_Temp = convert_to_temperature(VertT);
V_Winds = convert_to_winsdpeed(VertW,Vert_zw,V_Temp,'ms');

figure(3)
tiledlayout(2,2)
nexttile(1)
plot(a.Time, T_Temp)
xticklabels([])
nexttile(2)
plot(a.Time, V_Temp)
xticklabels([])
nexttile(3)
plot(a.Time, T_Winds)
nexttile(4)
plot(a.Time,V_Winds)

%%
%photodiode
figure(4)
tiledlayout(3,1)
nexttile(1)
plot(a.Time,a.PD3)
xticklabels([])
nexttile(2)
plot(a.Time,a.PD2)
nexttile(3)
plot(a.Time,a.PD1)
xlabel('Time [s]')

%%
[PD1y,PD1x] = photodiode_demodulate(a.PD1,Fs,[100,200,300],5);
[PD2y,PD2x] = photodiode_demodulate(a.PD2,Fs,[100,200,300],5);
[PD3y,PD3x] = photodiode_demodulate(a.PD3,Fs,[100,200,300],5);

figure(5)
tiledlayout(3,1)
nexttile(1)
plot(PD3x,PD3y)
xticklabels([])
nexttile(2)
plot(PD2x,PD2y)
nexttile(3)
plot(PD1x,PD1y)
xlabel('Time [s]')
%%
function out = h5totable(filename)
    time = h5read(filename ,'/time');
    data = h5read(filename,'/data');
    column_names = h5read(filename,'/channel_names');
    column_names = ['Time'; column_names];

    out = table(time,data');

    out = splitvars(out);
    
    out.Properties.VariableNames = column_names;
end

function output = spectrogram_mag(data, Fs, notch_mode, apply_highpass, apply_notch)
    % Defaults for optional arguments
    if nargin < 2 || isempty(notch_mode), notch_mode = 'comb'; end
    if nargin < 3 || isempty(apply_highpass), apply_highpass = true; end
    if nargin < 4 || isempty(apply_notch), apply_notch = true; end

    data = data * 1000;
    data = data - mean(data);
  
    % High-pass filter
    if apply_highpass
        data = highpass_filter(data, Fs, 2);
    end

    % Notch filter
    if apply_notch
        switch lower(notch_mode)
            case 'butter'
                data = comb_notch_filter(data, Fs, 50, 5, 1);
            case 'comb'
                data = notch_50Hz(data, Fs);
            case 'reduced'
                Hd = notch_50_100(Fs, 1);
                data = Hd(data);
            otherwise
                warning('Unknown notch mode. Skipping notch filtering.');
        end
    end


    % Spectrogram settings
    M = Fs;
    L = Fs / 2;
    g = tukeywin(M,0.1);
    Ndft = 2^nextpow2(M);
    Ts = 1 / Fs;
    
    fft_single(data,Fs);

    [s, f, t] = spectrogram(data, g, L, Ndft, Fs, 'onesided', 'yaxis');
    s = s ./ M;
    s(2:end-1, :) = 2 * s(2:end-1, :);
    
    sc = s .* 35;

    output.s = sc;
    output.f = f;
    output.t = t;
end

function plotspectrogram(input,cutoff)
    [x,y] = meshgrid(input.t,input.f(input.f<cutoff));
    z = abs(input.s(input.f<cutoff,:));
    z(z>50)=50;
    h = surf(x,y,z);
    set(h,'linestyle','none');
    view(0,90)
    ylim([0,cutoff])
    ylabel('Frequency [Hz]')
    xlim([0,input.t(end)])
    xlabel('Time [s]')
    c = colorbar;
    c.Label.String = ('Magnetic Flux Density [nT]');
    c.Label.Interpreter = 'Latex';
end

function filtered_data = highpass_filter(data, Fs, cutoff_freq)
    [b, a] = butter(4, cutoff_freq / (Fs / 2), 'high');  % 4th order Butterworth
    filtered_data = filtfilt(b, a, data);               % zero-phase filtering
end

function filtered_data = low_pass(data, Fs, cutoff_freq)
    [b, a] = butter(4, cutoff_freq / (Fs / 2), 'low');  % 4th order Butterworth
    filtered_data = filtfilt(b, a, data);               % zero-phase filtering
end

function filtered_data = comb_notch_filter(data, Fs, base_freq, num_harmonics, bandwidth)
    filtered_data = data;
    for k = 1:num_harmonics
        notch_freq = k * base_freq;
        Wn = [notch_freq - bandwidth/2, notch_freq + bandwidth/2] / (Fs/2);
        [b, a] = butter(1, Wn, 'stop');
        filtered_data = filtfilt(b, a, filtered_data);
    end
end

function filtered_data = notch_50Hz(input, Fs)
    % Design a comb notch filter for 50 Hz and harmonics - this uses a
    % legacy function
    L = round(Fs / 50);           % spacing between notches
    BW = 1;                       % Bandwidth of each notch (Hz)
    GBW = -5;                     % Gain at the bandwidth edges (dB)
    Nsh = floor(Fs/(2*50)) - 1;   % Number of harmonics to suppress

    d = fdesign.comb('notch','L,BW,GBW,Nsh', L, BW, GBW, Nsh, Fs);
    h = design(d,'SystemObject',true);

    % Apply zero-phase filtering
    filtered_data = filtfilt(h.Numerator, h.Denominator, real(input));
end

function filtered_data = notch_50Hz_iirnotch(data, Fs)
    filtered_data = data;
    base_freq = 50;
    num_harmonics = floor(Fs/(2*base_freq)); % Nyquist limited
    Q = 35; % Quality factor - adjust to control bandwidth
    
    for k = 1:num_harmonics
        f0 = k * base_freq;
        wo = f0 / (Fs/2);  % Normalized frequency
        bw = wo / Q;       % Bandwidth
        [b,a] = iirnotch(wo, bw);
        filtered_data = filtfilt(b, a, filtered_data);
    end
end

function Hd = notch_50_100(Fs, BW)
    % Design notch filter at 50 Hz
    W0_50 = 50 / (Fs/2);
    Q_50 = 50 / BW;
    [b50, a50] = iirnotch(W0_50, W0_50/Q_50);

    % Design notch filter at 100 Hz
    W0_100 = 100 / (Fs/2);
    Q_100 = 100 / BW;
    [b100, a100] = iirnotch(W0_100, W0_100/Q_100);

    % Combine both filters in cascade
    Hd = dsp.FilterCascade( ...
        dsp.IIRFilter('Numerator', b50, 'Denominator', a50), ...
        dsp.IIRFilter('Numerator', b100, 'Denominator', a100));
end

function output = fft_single(data,Fs)

    L = length(data);
    f = Fs*(0:(L/2))/L;
    Y = fft(data);
    %matlab fft doesn't normalise by size
    %make double sided spectrum
    P2 = abs(Y/L);
    %make single sided spectrum
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    output.data = data;
    output.Y = Y;
    output.P1 = P1;
    output.f = f;
end

function out = convert_to_temperature(input)
 out = (input - 0.400)./(0.0195);
end

function out = convert_to_winsdpeed(input_w, cal_f, temp,unit)
 out = (((input_w - cal_f)./(3.038517 .* (temp.^(0.115157))))./(0.087288)).^(3.009364);
 if strcmp(unit,'mph')
    out(out<0) = 0;
    out = out;
 elseif strcmp(unit,'ms')
    out(out<0) = 0;
    out = out * 0.44704; %mph to m/s
 elseif strcmp(unit,'kph')
    out(out<0) = 0;
    out = out * 1.609344; %mph to kph
 else
     out = out;
 end
end

function [amp,t_out] = photodiode_demodulate(data,Fs,carriers,lpf)
    N = length(data);
    t_out = (0:N-1)' ./ Fs;

    numF = numel(carriers);

    data = data - mean(data);

    amp = zeros(N,numF);

    for i = 1:numF
        f = carriers(i);
        ref_sin = sin(2*pi*f*t_out);
        ref_cos = cos(2*pi*f*t_out);

        I = data .* ref_cos;
        Q = data .* ref_sin;

        %I = I - mean(I);
        %Q = Q - mean(Q);

        I_lp = low_pass(I,Fs,lpf);
        Q_lp = low_pass(Q,Fs,lpf);

        amp(:,i) = sqrt(I_lp.^2 + Q_lp.^2);
    end
end