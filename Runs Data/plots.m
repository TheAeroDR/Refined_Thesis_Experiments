%TERIE state
particle = "empty";
airflow = "noair";

figure(1) %electrometer
tiledlayout(4,3)
figure(2) %electrometer
tiledlayout(4,3)
figure(3) %photodiode
tiledlayout(3,3)
figure(4) %anemometer
tiledlayout(2,3)
figure(5) %magnetometer
tiledlayout(3,3)
figure(6) %magnetometer
tiledlayout(3,3)

sensors = {'BE1','BE2','BE3','BE4','PD1','PD2','PD3','ANTI_wind_speed',...
    'ANTI_temperature','ANTO_wind_speed','ANTO_temperature','ANV_wind_speed',...
    'ANV_temperature','IFG1_Z','IFG2_Z','IFG3_Z','IFG1_Y','IFG2_Y','IFG3_Y'};
%%
for i = 1:3
    filename = "pos_" + num2str(i) + "_" + particle + "_" + airflow + ".h5";

    data = h5totable(filename);

    Fs = h5readatt(filename,'/','sampling_rate');
    
    %electrometer
    R_fb = 1.01e9;
    eps_0 = 8.8541878188e-12;
    r_e = 6e-3;
    sensor_idx_offset = 0;
    for s = 1:4
        sensor_name = sensors{s+sensor_idx_offset};
        sensor_data = data.(sensor_name);
        
        sensor_mm = movmean(sensor_data,floor(Fs/10));
        
        sensor_ER = - sensor_mm / (4 * pi * r_e^2 * R_fb * eps_0);
        
        figure(1)
        nexttile(12 - 3*s + i)
        yyaxis left
        plot(data.Time,sensor_mm,'-')
        yyaxis right
        plot(data.Time,sensor_ER)
        if s ~= 1
            xticklabels('')
        else
            xlabel('Time [s]')
        end
        if i == 1
            yyaxis left
            ylabel('Voltage [V]')
        elseif i == 3
            yyaxis right
            ylabel('$\frac{dE}{dt}$ [V/m/s]')
        end

        if s == 4
            title(sprintf('Position %d',i))
        end

        sensor_E = decimate(sensor_ER,Fs/10);
        t = decimate(data.Time,Fs/10);
        sensor_E = 0.5 * (sensor_E(1:end-1)+sensor_E(2:end)) * (t(2)-t(1));
        sensor_E = [NaN;sensor_E];
    
        figure(2)
        nexttile(12 - 3*s + i)
        plot(t,sensor_E)

        if s ~= 1
            xticklabels('')
        else
            xlabel('Time [s]')
        end
        if i == 1
            ylabel('Electric Field [V/m]')
        end

        if s == 4
            title(sprintf('Position %d',i))
        end
    end
end
%%
for i = 1:3
    filename = "pos_" + num2str(i) + "_" + particle + "_" + airflow + ".h5";

    data = h5totable(filename);

    Fs = h5readatt(filename,'/','sampling_rate');
    %photodiodes
    sensor_idx_offset = 4;
    for s = 1:3
        sensor_name = sensors{s+sensor_idx_offset};
        sensor_data = data.(sensor_name);

        [PDy,PDx] = photodiode_demodulate(sensor_data,Fs,[100,200,300],5);
    
        figure(3)
        nexttile(9 - 3*s + i)
        plot(PDx,PDy)

        if s ~= 1
            xticklabels('')
        else
            xlabel('Time [s]')
        end

        if s == 3
            title(sprintf('Position %d',i))
        end
        
    end
end
%%
for i = 1:3
    filename = "pos_" + num2str(i) + "_" + particle + "_" + airflow + ".h5";

    data = h5totable(filename);

    Fs = h5readatt(filename,'/','sampling_rate');

    %anemometers
    zero_wind = [1.2185, NaN, 1.1848];
    sensor_idx_offset = 7;
    for s = 1:3
        if s ~= 2
            sensor_t_name = sensors{(2*s)+sensor_idx_offset};
            sensor_t_data = data.(sensor_t_name);

            sensor_w_name = sensors{(2*s-1)+sensor_idx_offset};
            sensor_w_data = data.(sensor_w_name);
    
            sensorW = movmean(sensor_w_data,floor(Fs/2));
            sensorT = movmean(sensor_t_data,floor(Fs/2));

            sensorT = convert_to_temperature(sensorT);
            sensorW = convert_to_winsdpeed(sensorW,zero_wind(s),sensorT,'ms');

            figure(4)
            sm = s;
            if s == 3
                sm = 2;
            end
            nexttile(6 - 3*sm + i)
            yyaxis left
            plot(data.Time,sensorW)
            
            yyaxis right
            plot(data.Time,sensorS)

            if s ~= 1
                xticklabels('')
            else
                xlabel('Time [s]')
            end
            if i == 1
                yyaxis left
                ylabel('Windspeed [m/s]')
            elseif i == 3
                yyaxis right
                ylabel('Temperature [K]')
            end
        
            if s == 3
                title(sprintf('Position %d',i))
            end
        end
    end
end
%%
for i = 1:3
    filename = "pos_" + num2str(i) + "_" + particle + "_" + airflow + ".h5";

    data = h5totable(filename);

    Fs = h5readatt(filename,'/','sampling_rate');
    %magnetometer z
    sensor_idx_offset = 13;
    for s = 1:3
        sensor_name = sensors{(s)+sensor_idx_offset};
        sensor_data = data.(sensor_name);
        sensor_b = spectrogram_mag(sensor_data, Fs, 'reduced', true, true);
        figure(5)
        nexttile(9 - 3*s + i)
        plotspectrogram(sensor_b,100)
        if s ~= 1
            xticklabels('')
        else
            xlabel('Time [s]')
        end
        
        if s == 3
            title(sprintf('Position %d',i))
        end
    end
end
%%
for i = 1:3
    filename = "pos_" + num2str(i) + "_" + particle + "_" + airflow + ".h5";

    data = h5totable(filename);

    Fs = h5readatt(filename,'/','sampling_rate');
    %magnetometer y
    sensor_idx_offset = 16;
    for s = 1:3
                sensor_name = sensors{(s)+sensor_idx_offset};
        sensor_data = data.(sensor_name);
        sensor_b = spectrogram_mag(sensor_data, Fs, 'reduced', true, true);
        figure(6)
        nexttile(9 - 3*s + i)
        plotspectrogram(sensor_b,100)
        if s ~= 1
            xticklabels('')
        else
            xlabel('Time [s]')
        end
        
        if s == 3
            title(sprintf('Position %d',i))
        end
    end
end
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
    0
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