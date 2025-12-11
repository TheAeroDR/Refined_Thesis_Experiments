%% photodiode

data = readtable("photodiode_test.csv");

figure
hold on
plot(data.time_s,data.PD1)
plot(data.time_s,data.PD2)
plot(data.time_s,data.PD3)
legend('PD1','PD2','PD3')

%%
LED_dat = readtable("LED_RLI.csv"); %wavelength vs RLI (%)
PD_dat = readtable("OPT101_PR.csv"); %wavelendth vs Photodiode Responsivity (A/W)
CIE_dat = readtable("CIE_sle_photopic.csv"); %https://doi.org/10.25039/CIE.DS.dktna2s3

%interpolate onto same range
dl = 1;
lambd = (200:dl:1100)';
LED_dat = interp1(LED_dat.Var1,LED_dat.Var2./sum(LED_dat.Var2),lambd,'linear',1e-12);
PD_dat = interp1(PD_dat.Var1,PD_dat.Var2,lambd,'linear',1e-12);
CIE_dat = interp1(CIE_dat.Var1,CIE_dat.Var2,lambd,'linear',1e-12);

LED_intensity = 38000; %mcd
viewing_angle = 15; %deg

LED_lm = LED_intensity * (2*pi*(1-cosd(viewing_angle/2)));

srf_wm = (LED_lm * LED_dat)./(683 * CIE_dat);

R_eff = sum(srf_wm .* PD_dat) * dl / (sum(srf_wm) * dl);

figure
plot(data.time_s,(data.PD1/(1e6*R_eff)))
hold on
plot(data.time_s,(data.PD2/(1e6*R_eff)))
plot(data.time_s,(data.PD3/(1e6*R_eff)))
legend('PD1','PD2','PD3')

%% electrometer
data = readtable("electrometer_test.csv",'NumHeaderLines',1);
firstRow = readcell('electrometer_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

R_fb = 1.01e9;
eps_0 = 8.8541878188e-12;
r_e = 6e-3;


BE1 = movmean(data.BE1,floor(Fs/10));
BE2 = movmean(data.BE2,floor(Fs/10));
BE3 = movmean(data.BE3,floor(Fs/10));
BE4 = movmean(data.BE4,floor(Fs/10));

BE1_ER = - BE1 / (4 * pi * r_e^2 * R_fb * eps_0);
BE2_ER = - BE2 / (4 * pi * r_e^2 * R_fb * eps_0);
BE3_ER = - BE3 / (4 * pi * r_e^2 * R_fb * eps_0);
BE4_ER = - BE4 / (4 * pi * r_e^2 * R_fb * eps_0);

figure
tiledlayout(4,1)
nexttile
yyaxis left
plot(data.time_s,data.BE4,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE4,'-')
yyaxis right
plot(data.time_s,BE4_ER)
%legend('BE4 Voltage','0.1 s Moving Average','$\frac{dE}{dt}$','Location','east')
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE3,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE3,'-')
yyaxis right
plot(data.time_s,BE3_ER)
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE2,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE2,'-')
yyaxis right
plot(data.time_s,BE2_ER)
xticklabels('')

nexttile
yyaxis left
plot(data.time_s,data.BE1,'Color',[0.30,0.75,0.93])
hold on
plot(data.time_s,BE1,'-')
ylabel('Voltage [V]')
yyaxis right
plot(data.time_s,BE1_ER)
xlabel('Time [s]')
ylabel('$\frac{dE}{dt}$ [V/m/s]')

a4 = decimate(BE4_ER,floor(Fs/10));
a3 = decimate(BE3_ER,floor(Fs/10));
a2 = decimate(BE2_ER,floor(Fs/10));
a1 = decimate(BE1_ER,floor(Fs/10));

t = decimate(data.time_s,floor(Fs/10));

E4 = 0.5 * (a4(1:end-1)+a4(2:end)) * (t(2)-t(1));
E4 = [NaN;E4];
E3 = 0.5 * (a3(1:end-1)+a3(2:end)) * (t(2)-t(1));
E3 = [NaN;E3];
E2 = 0.5 * (a2(1:end-1)+a2(2:end)) * (t(2)-t(1));
E2 = [NaN;E2];
E1 = 0.5 * (a1(1:end-1)+a1(2:end)) * (t(2)-t(1));
E1 = [NaN;E1];
figure
tiledlayout(4,1)
nexttile(1)
plot(t,E4)
nexttile(2)
plot(t,E3)
nexttile(3)
plot(t,E2)
nexttile(4)
plot(t,E1)
%% windspeed
%data = readtable("anemometer_zeros.csv",'NumHeaderLines',1);

%TIN_zw = mean(data.TInW);
%Vert_zw =  mean(data.VertW);

%data = readtable("anemometer_zeros2.csv",'NumHeaderLines',1);

%TIN_zw(2) = mean(data.TInW);
%Vert_zw(2) =  mean(data.VertW);

TIN_zw = 1.1848;
Vert_zw = 1.2185;

data = readtable("anemometer_test.csv",'NumHeaderLines',1);
firstRow = readcell('anemometer_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

TInW = movmean(data.TInW,floor(Fs/2));
TInT = movmean(data.TInT,floor(Fs/2));
VertW = movmean(data.VertW,floor(Fs/2));
VertT = movmean(data.VertT,floor(Fs/2));

TInW = decimate(TInW,floor(Fs/2));
TInT = decimate(TInT,floor(Fs/2));
VertW = decimate(VertW,floor(Fs/2));
VertT = decimate(VertT,floor(Fs/2));
time = decimate(data.time_s,floor(Fs/2));

temp = convert_to_temperature(TInT);
winds = convert_to_winsdpeed(TInW,TIN_zw,temp,'kph');

figure(1)
tiledlayout(2,1)
nexttile(1)
yyaxis left
plot(time,winds)
hold on
yyaxis right
plot(time,temp)

temp = convert_to_temperature(VertT);
winds = convert_to_winsdpeed(VertW,Vert_zw,temp,'kph');

nexttile(2)
yyaxis left
plot(time,winds)
ylabel('Windspeed [kph]')
hold on
yyaxis right
plot(time,temp)
ylabel('Temperature [deg C]')
xlabel('Time [s]')
%%
data = readtable("anemometer_test2.csv",'NumHeaderLines',1);
firstRow = readcell('anemometer_test2.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

TInW = movmean(data.TInW,floor(Fs/2));
TInT = movmean(data.TInT,floor(Fs/2));
VertW = movmean(data.VertW,floor(Fs/2));
VertT = movmean(data.VertT,floor(Fs/2));

TInW = decimate(TInW,floor(Fs/2));
TInT = decimate(TInT,floor(Fs/2));
VertW = decimate(VertW,floor(Fs/2));
VertT = decimate(VertT,floor(Fs/2));
time = decimate(data.time_s,floor(Fs/2));

temp = convert_to_temperature(TInT);
winds = convert_to_winsdpeed(TInW,TIN_zw,temp,'kph');

figure(2)
tiledlayout(2,1)
nexttile(1)
yyaxis left
plot(time,winds)
hold on
yyaxis right
plot(time,temp)

temp = convert_to_temperature(VertT);
winds = convert_to_winsdpeed(VertW,Vert_zw,temp,'kph');

nexttile(2)
yyaxis left
plot(time,winds)
ylabel('Windspeed [kph]')
hold on
yyaxis right
plot(time,temp)
ylabel('Temperature [deg C]')
xlabel('Time [s]')


%%
data = readtable("psd_from_image.csv");
histogram(data.size_um,200,'Normalization','pdf')
xlabel("Particle Size [$\mu$m]")
ylabel("Probability Density [-]")
set(gca,"XScale","log")
pd = fitdist(data{:,1}, 'Lognormal');

x_values = linspace(0, 5000, 1000);
y_values = pdf(pd, x_values);

hold on
plot(x_values, y_values, 'r-', 'LineWidth', 2);


%%
%magnetometers
data = readtable("mag1_test.csv");
firstRow = readcell('mag1_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

figure
tiledlayout(3,3)

nexttile(1)
a = spectrogram_mag(data.FGx, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
xticklabels([])

nexttile(4)
a = spectrogram_mag(data.FGy, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
xticklabels([])

nexttile(7)
a = spectrogram_mag(data.FGz, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])

data = readtable("mag2_test.csv");
firstRow = readcell('mag1_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

nexttile(2)
a = spectrogram_mag(data.FGx, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
xticklabels([])
yticklabels([])

nexttile(5)
a = spectrogram_mag(data.FGy, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
xticklabels([])
yticklabels([])

nexttile(8)
a = spectrogram_mag(data.FGz, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
ylabel([])
yticklabels([])

data = readtable("mag3_test.csv");
firstRow = readcell('mag1_test.csv', 'Range', '1:1');
Fs = sscanf(firstRow{1,1}, 'Sample rate: %f S/s');

nexttile(3)
a = spectrogram_mag(data.FGx, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
xticklabels([])
yticklabels([])

nexttile(6)
a = spectrogram_mag(data.FGy, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
xticklabels([])
yticklabels([])

nexttile(9)
a = spectrogram_mag(data.FGz, Fs, 'reduced', 0, 1);
plotspectrogram(a,50)
xlabel([])
ylabel([])
yticklabels([])
%%

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
    M = floor(Fs);
    L = M / 2;
    g = tukeywin(M,0.1);
    Ndft = 2^nextpow2(M);
    
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
    clim([0,8000])
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