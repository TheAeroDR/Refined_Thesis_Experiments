import h5py
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, decimate, iirnotch, spectrogram
from scipy.signal import freqz
from scipy.signal.windows import tukey
from scipy.ndimage import uniform_filter1d
from scipy.optimize import minimize
from scipy.fft import fft
import gc


# Matplotlib settings
plt.rcParams.update({
    # Grid and box
    'axes.grid': True,
    'axes.grid.axis': 'both',
    'axes.xmargin': 0,
    'axes.ymargin': 0,
    'axes.titlepad': 6.0,
    'axes.edgecolor': 'black',
    'axes.linewidth': 1.0,

    # Ticks
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.top': True,
    'ytick.right': True,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,

    # Fonts & LaTeX
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    'axes.labelsize': 30,
    'axes.titlesize': 30,
    'xtick.labelsize': 30,
    'ytick.labelsize': 30,
    'legend.fontsize': 30,
    'axes.formatter.limits': [-3, 3],


    # Line widths and figure rendering
    'lines.linewidth': 1.5,
    'lines.markersize': 6,
    'figure.autolayout': True,
    'figure.dpi': 600.0,
    'figure.figsize': [16.0, 9.0],
})
# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------

def h5totable(filename):
    with h5py.File(filename, 'r') as f:
        time = np.array(f['/time'])
        data = np.array(f['/data']).T
        channel_names = ["Time"] + [name.decode('utf-8') for name in f['/channel_names']]
    out = {"Time": time}
    for idx, name in enumerate(channel_names[1:]):
        out[name] = data[idx, :]
    return out

def highpass_filter(data, Fs, cutoff):
    b, a = butter(4, cutoff / (Fs / 2), 'highpass')
    return filtfilt(b, a, data)

def low_pass(data, Fs, cutoff):
    b, a = butter(4, cutoff / (Fs / 2), 'lowpass')
    return filtfilt(b, a, data)

def fft_single(data, Fs):
    data = data * 1000
    
    L = len(data)
    Y = fft(data)
    P2 = np.abs(Y/L)
    P1 = P2[: L//2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    P1 = P1 * 35
    f = Fs * np.arange(0, L/2 + 1) / L
    return {"data": data, "Y": Y, "P1": P1, "f": f}

def notch_50_100(Fs, BW):
    W0_50 = 50 / (Fs / 2)
    Q_50 = 50 / BW
    b50, a50 = iirnotch(W0_50, Q_50)

    W0_100 = 100 / (Fs / 2)
    Q_100 = 100 / BW
    b100, a100 = iirnotch(W0_100, Q_100)

    def filt(x):
        y = filtfilt(b50, a50, x)
        y = filtfilt(b100, a100, y)
        return y

    return filt

def spectrogram_mag(data, Fs, notch_mode='reduced', apply_highpass=True, apply_notch=True):
    data = data * 1000

    data = data - np.mean(data)

    if apply_highpass:
        data = highpass_filter(data, Fs, 2)

    if apply_notch:
        if notch_mode.lower() == 'reduced':
            Hd = notch_50_100(Fs, 1)
            data = Hd(data)

    M = int(Fs)
    L = int(Fs) // 2
    g = tukey(M, 0.1)
    Ndft = 2 ** int(np.ceil(np.log2(M)))

    f, t, S = spectrogram(data, window=g, noverlap=L, nfft=Ndft, fs=Fs, scaling='spectrum', mode='complex')

    sc = S * 35

    return {"s": sc, "f": f, "t": t}

def plotspectrogram(input_obj, cutoff, apply_lims, clim):
    mask = input_obj["f"] < cutoff
    f = input_obj["f"][mask]
    t = input_obj["t"]
    z = np.abs(input_obj["s"][mask, :])
    z[z > 50] = 50

    if apply_lims:
        plt.pcolormesh(t, f, z, shading='auto',cmap='viridis',vmin=clim[0],vmax=clim[1],zorder=3, rasterized=True)
    else:
        plt.pcolormesh(t, f, z, shading='auto',cmap='viridis',zorder=3, rasterized=True)
    plt.ylim([0, cutoff])
    #c = plt.colorbar()

def fakefft(input_obj, cutoff, apply_lims, clim):
    mask = input_obj["f"] < cutoff
    f = input_obj["f"][mask]
    t = input_obj["t"]
    z = np.abs(input_obj["s"][mask, :])

    z_mean = np.mean(z, axis=1)

    plt.plot(f, z_mean, color='#0072BD', zorder=3)

    if apply_lims:
        plt.ylim([clim[0], clim[1]])

def convert_to_temperature(x):
    return (x - 0.400) / 0.0195

def convert_to_winsdpeed(input_w, cal_f, temp, unit):
    out = ((input_w - cal_f) / (3.038517 * (temp ** 0.115157))) / 0.087288
    out = out ** 3.009364
    out[out < 0] = 0
    if unit == 'ms':
        out *= 0.44704
    elif unit == 'kph':
        out *= 1.609344
    return out

def photodiode_demodulate(data, Fs, carriers, lpf):
    N = len(data)
    t_out = np.arange(N) / Fs
    data = data - np.mean(data)

    amp = np.zeros((N, len(carriers)))

    for i, f in enumerate(carriers):
        ref_sin = np.sin(2 * np.pi * f * t_out)
        ref_cos = np.cos(2 * np.pi * f * t_out)

        I = data * ref_cos
        Q = data * ref_sin

        I_lp = low_pass(I, Fs, lpf)
        Q_lp = low_pass(Q, Fs, lpf)

        amp[:, i] = np.sqrt(I_lp**2 + Q_lp**2)

    amp = decimate(amp, int(Fs / lpf), axis=0)
    t_out = decimate(t_out, int(Fs / lpf))

    return amp, t_out

# Electrometer plot function ------------------------------------------------
def plot_electrometer(figure1=True, figure2=True, apply_lims=True, background_red = True):
    fig1,_ = plt.subplots(4,3,num=1)
    fig1.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)

    fig2,_ = plt.subplots(4,3,num=2)
    fig2.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)

    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}_2.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        R_fb = 1.01e9
        eps_0 = 8.8541878188e-12
        r_e = 6e-3


        for s in range(1, 5):
            sensor_name = sensors[s-1]
            sensor_data = data[sensor_name]

            if background_red:
                sensor_data = sensor_data - electrometer_background[s]

            N = int(np.floor(Fs/10))
            sensor_mm = uniform_filter1d(sensor_data, size=N, mode='nearest')
            sensor_mm = decimate(sensor_mm, N)
            t = decimate(data['Time'], N)

            sensor_ER = -sensor_mm / (4 * np.pi * r_e**2 * R_fb * eps_0)

            idx = 12 - 3*(s) + i

            if figure1:
                plt.figure(fig1.number)
                
                ax_left = plt.subplot(4, 3, idx)
                ax_left.plot(t, sensor_mm,color='#0072BD')
                ax_right = ax_left.twinx()
                ax_right.plot(t, sensor_ER,color='#D85319')

                if i == 1:
                    ax_left.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax_left.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                else:
                    ax_left.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax_left.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                ax_left.tick_params(axis='y', colors='#0072BD')
                ax_left.spines['left'].set_color('#0072BD')
                ax_right.tick_params(axis='y', colors='#D85319')
                ax_right.spines['top'].set_visible(False)
                ax_right.spines['bottom'].set_visible(False)
                ax_right.spines['left'].set_visible(False)
                ax_right.spines['right'].set_color('#D85319')
                if apply_lims:
                    limits_voltage = axis_limits[particle][airflow]['bell']['voltage'][s]
                    limits_dedt = axis_limits[particle][airflow]['bell']['dedt'][s]

                    ax_left.set_ylim(limits_voltage['ylim'])
                    ax_right.set_ylim(limits_dedt['ylim'])
                ax_left.set_xlim([0, 120])

                if i == 1:
                    if s == 3:
                        fig1.supylabel('Voltage [V]', fontsize=30, color='#0072BD') 
                    if s == 4:
                        ax_left.set_title('Position 1', fontsize=30)
                    ax_right.set_yticklabels([])
                elif i == 2:
                    ax_left.set_yticklabels([])
                    ax_right.set_yticklabels([])
                    if s == 1:
                        fig1.supxlabel('Time [s]', fontsize=30)
                    if s == 4:
                        ax_left.set_title('Position 2', fontsize=30)
                elif i == 3:
                    ax_left.set_yticklabels([])
                    if s == 3:
                        fig1.text(0.96, 0.5, r'$\frac{dE}{dt}$ [V/m/s]', va='center', rotation='vertical', fontsize=30, color='#D85319')
                    if s == 4:
                        ax_left.set_title('Position 3', fontsize=30)
                
                if s !=1:
                    ax_left.set_xticklabels([])
                
                if s == 1 and i == 1 and spike and particle == 'empty' and airflow == 'air':
                    plt.figure(num=10)
                    t_mod = t
                    I_mod = -sensor_mm / R_fb
                    mask = (t_mod >= 25) & (t_mod <= 35)
                    t_mod = t_mod[mask] 
                    I_mod = I_mod[  mask]
                    plt.plot(t_mod, I_mod)
                    plt.ylabel('Current [A]')
                    plt.xlabel('Time [s]')
                    #plt.xticks(fontsize=8)
                    #plt.yticks(fontsize=8)
                    #plt.xlim([-1,1])

                    
            sensor_E = 0.5*(sensor_ER[:-1] + sensor_ER[1:]) * (t[1] - t[0])
            sensor_E = np.insert(sensor_E, 0, np.nan)

            if figure2:
                plt.figure(fig2.number)
                ax = plt.subplot(4, 3, idx)
                ax.plot(t, sensor_E)

                if i == 1:
                    ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                else:
                    ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                if apply_lims:
                    limits_E = axis_limits[particle][airflow]['bell']['E'][s]

                    ax.set_ylim(limits_E['ylim'])
                ax.set_xlim([0, 120])

                if i == 1:
                    if s == 3:
                        fig2.supylabel('Electric Field [V/m]', fontsize=30) 
                    if s == 4:
                        ax.set_title('Position 1', fontsize=30)
                elif i == 2:
                    ax.set_yticklabels([])
                    if s == 1:
                        fig2.supxlabel('Time [s]', fontsize=30)
                    if s == 4:
                        ax.set_title('Position 2', fontsize=30)
                elif i == 3:
                    ax.set_yticklabels([])
                    if s == 4:
                        ax.set_title('Position 3', fontsize=30)

                if s !=1:
                    ax.set_xticklabels([])
        del data
        gc.collect()

# photodiode plot function ------------------------------------------------------
def plot_photodiode(figure3=True, apply_lims=True):
    fig3,_ = plt.subplots(3,3,num=3)
    fig3.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}_2.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            sensor_name = sensors[3 + s]
            sensor_data = data[sensor_name]
            PDy, PDx = photodiode_demodulate(sensor_data, Fs, [100,200,300], 5)

            if figure3:
                plt.figure(fig3.number)
                idx = 9 - 3*s+ i
                ax =plt.subplot(3, 3, idx)
                ax.plot(PDx, PDy)

                if i == 1:
                    ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                else:
                    ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                if apply_lims:
                    limits_pd = axis_limits[particle][airflow]['pd'][s]
                    ax.set_ylim(limits_pd['ylim'])
                ax.set_xlim([0, 120])
                if i == 1:
                    if s == 3:
                        fig3.supylabel('Voltage [V]', fontsize=30) 
                    if s == 3:
                        ax.set_title('Position 1', fontsize=30)
                    
                elif i == 2:
                    ax.set_yticklabels([])
                    if s == 1:
                        fig3.supxlabel('Time [s]', fontsize=30)
                    if s == 3:
                        ax.set_title('Position 2', fontsize=30)
                elif i == 3:
                    ax.set_yticklabels([])
                    if s == 3:
                        ax.set_title('Position 3', fontsize=30)

                if s !=1:
                    ax.set_xticklabels([])


        del data
        gc.collect()

# anemometer plot function ------------------------------------------------------
def plot_anemometer(figure4=True, apply_lims=True):
    fig4,_ = plt.subplots(2,3,num=4)
    fig4.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)
    zero_wind = [1.2185, np.nan, 1.1848]


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}_2.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            if s == 2:
                continue


            sensor_t_name = sensors[(2*s) + 6]
            sensor_w_name = sensors[(2*s-1) + 6]

            N = int(np.floor(Fs/2))

            sensor_T = uniform_filter1d(data[sensor_t_name], size=N, mode='nearest')
            sensor_W = uniform_filter1d(data[sensor_w_name], size=N, mode='nearest')

            sensorT = convert_to_temperature(sensor_T)
            sensorW = convert_to_winsdpeed(sensor_W, zero_wind[s-1], sensorT, 'ms')

            sensorT = decimate(sensorT, N)
            sensorW = decimate(sensorW, N)
            t = decimate(data['Time'], N)

            sm = s if s != 3 else 2

            if figure4:
                plt.figure(fig4.number)
                idx = 6 - 3*sm + i
                ax_left =plt.subplot(2, 3, idx)
                plt.plot(t, sensorW)

                if i == 1:
                    ax_left.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax_left.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                else:
                    ax_left.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax_left.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                ax_left.plot(t, sensorW,color='#0072BD')
                ax_right = ax_left.twinx()
                ax_right.plot(t, sensorT,color='#D85319')

                ax_left.tick_params(axis='y', colors='#0072BD')
                ax_left.spines['left'].set_color('#0072BD')
                ax_right.tick_params(axis='y', colors='#D85319')
                ax_right.spines['top'].set_visible(False)
                ax_right.spines['bottom'].set_visible(False)
                ax_right.spines['left'].set_visible(False)
                ax_right.spines['right'].set_color('#D85319')

                if apply_lims:
                    limits_w= axis_limits[particle][airflow]['anem']['wind'][s]
                    limits_t= axis_limits[particle][airflow]['anem']['temp'][s]

                    ax_left.set_ylim(limits_w['ylim'])
                    ax_right.set_ylim(limits_t['ylim'])
                ax_left.set_xlim([0, 120])

                if i == 1:
                    if s == 3:
                        fig4.supylabel('Windspeed [m/s]', fontsize=30, color='#0072BD') 
                        ax_left.set_title('Position 1', fontsize=30)
                    ax_right.set_yticklabels([])
                elif i == 2:
                    ax_left.set_yticklabels([])
                    ax_right.set_yticklabels([])
                    if s == 1:
                        fig4.supxlabel('Time [s]', fontsize=30)
                    if s == 3:
                        ax_left.set_title('Position 2', fontsize=30)
                elif i == 3:
                    ax_left.set_yticklabels([])
                    if s == 3:
                        fig4.text(0.96, 0.5, 'Temperature [deg C]', va='center', rotation='vertical', fontsize=30, color='#D85319')
                        ax_left.set_title('Position 3', fontsize=30)
                
                if s !=1:
                    ax_left.set_xticklabels([])
        
        del data
        gc.collect()

# magnetometer plot function ------------------------------------------------------
def plot_magnetometer(figure5=True, figure6=True, apply_lims=True, plot_freq=100, apply_highpass=True, apply_notch=True,fft_flag=False):
    # Z (Figure 5)
    fig5, _ = plt.subplots(3,3,num=5)
    fig5.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}_2.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            sensor_name = sensors[12 + s]
            sensor_b = spectrogram_mag(data[sensor_name], Fs, 'reduced', apply_highpass, apply_notch)
            idx = 9 - 3*s + i
            if figure5:
                plt.figure(fig5.number)
                ax = plt.subplot(3, 3, idx)
                if fft_flag:
                    fakefft(sensor_b, plot_freq, apply_lims, axis_limits[particle][airflow]['mag']['z'][s]['clim'])
                    
                    if i == 1:
                        ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                    else:
                        ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right') 

                    fig5.supylabel('Magnetic Flux Density [nT]', fontsize=30)
                    fig5.supxlabel('Frequency [Hz]', fontsize=30)

                else:
                    plotspectrogram(sensor_b, plot_freq, apply_lims, axis_limits[particle][airflow]['mag']['z'][s]['clim'])

                    if i == 1:
                        ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                    else:
                        ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                    ax.set_xlim([0, 120])
                    if i == 1:
                        if s == 3:
                            fig5.supylabel('Frequency [Hz]', fontsize=30) 
                        if s == 3:
                            ax.set_title('Position 1', fontsize=30)
                        
                    elif i == 2:
                        ax.set_yticklabels([])
                        if s == 1:
                            fig5.supxlabel('Time [s]', fontsize=30)
                        if s == 3:
                            ax.set_title('Position 2', fontsize=30)
                    elif i == 3:
                        ax.set_yticklabels([])
                        c = plt.colorbar(ax=ax)
                        if s == 3:
                            ax.set_title('Position 3', fontsize=30)
                        if s == 2:
                            c.set_label('Magnetic Flux Density [nT]', fontsize=30)

                    if s !=1:
                        ax.set_xticklabels([])

        del data
        gc.collect()

    # Y (Figure 6)
    fig6, _ = plt.subplots(3,3,num=6)
    fig6.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)

    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}_2.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            sensor_name = sensors[15 + s]
            sensor_b = spectrogram_mag(data[sensor_name], Fs, 'reduced', True, True)
            idx = 9 - 3*s + i
            if figure6:
                plt.figure(fig6.number)
                ax = plt.subplot(3, 3, idx)
                if fft_flag:
                    fakefft(sensor_b, plot_freq, apply_lims, axis_limits[particle][airflow]['mag']['y'][s]['clim'])

                    if i == 1:
                        ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                    else:
                        ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
    
                    fig6.supylabel('Magnetic Flux Density [nT]', fontsize=30)
                    fig6.supxlabel('Frequency [Hz]', fontsize=30)
                    
                else:
                    plotspectrogram(sensor_b, plot_freq, apply_lims, axis_limits[particle][airflow]['mag']['y'][s]['clim'])

                    if i == 1:
                        ax.text(-0.1, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')
                    else:
                        ax.text(-0.01, 1.05, r'$\textbf{%s}$' % labels[idx-1], transform=ax.transAxes, fontsize=22,  verticalalignment='bottom', horizontalalignment='right')

                    if i == 1:
                        if s == 3:
                            fig6.supylabel('Frequency [Hz]', fontsize=30) 
                        if s == 3:
                            ax.set_title('Position 1', fontsize=30)
                        
                    elif i == 2:
                        ax.set_yticklabels([])
                        if s == 1:
                            fig6.supxlabel('Time [s]', fontsize=30)
                        if s == 3:
                            ax.set_title('Position 2', fontsize=30)
                    elif i == 3:
                        ax.set_yticklabels([])
                        c = plt.colorbar(ax=ax)
                        if s == 3:
                            ax.set_title('Position 3', fontsize=30)
                        if s == 2:
                            c.set_label('Magnetic Flux Density [nT]', fontsize=30)

                    if s !=1:
                        ax.set_xticklabels([])

        del data
        gc.collect()

# axis limits dictionary ------------------------------------------------------
axis_limits = {
    'empty': {
        'noair': {
            'bell': {
                'voltage': {
                    1:{'ylim': [-0.02, 0.02]},
                    2:{'ylim': [-0.02, 0.02]},
                    3:{'ylim': [0.01, 0.05]},
                    4:{'ylim': [-0.02, 0.02]},
                },
                'dedt': {
                    1:{'ylim': [-5000, 5000]},
                    2:{'ylim': [-5000, 5000]},
                    3:{'ylim': [-12500, -2500]},
                    4:{'ylim': [-5000, 5000]},
                },
                'E': {
                    1:{'ylim': [140, 320]},
                    2:{'ylim': [-220, -80]},
                    3:{'ylim': [-460, -340]},
                    4:{'ylim': [80, 220]},
                },
            },
            'mag': {
                'z': {
                    1:{'clim': [0, 3]},
                    2:{'clim': [0, 3]},
                    3:{'clim': [0, 3]},
                },
                'y': {
                    1:{'clim': [0, 5]},
                    2:{'clim': [0, 5]},
                    3:{'clim': [0, 5]},
                },
            },
            'pd': {
                1:{'ylim': [0, 0.08]},
                2:{'ylim': [0, 0.08]},
                3:{'ylim': [0, 0.08]},
            },
            'anem': {
                'wind': {
                    1:{'ylim': [-1, 1]},
                    2:{'ylim': [-1, 1]},
                    3:{'ylim': [-1, 1]},
                },
                'temp': {
                    1:{'ylim': [23, 27]},
                    2:{'ylim': [23, 27]},
                    3:{'ylim': [23, 27]},
                },
            },
        },
        'air': {
            'bell': {
                'voltage': {
                    1:{'ylim': [-1.2, 1.2]},
                    2:{'ylim': [-0.5, 0.5]},
                    3:{'ylim': [-0.3, 0.3]},
                    4:{'ylim': [-0.3, 0.3]},
                },
                'dedt': {
                    1:{'ylim': [-100000, 100000]},
                    2:{'ylim': [-50000, 50000]},
                    3:{'ylim': [-50000, 50000]},
                    4:{'ylim': [-50000, 50000]},
                },
                'E': {
                    1:{'ylim': [-20000, 20000]},
                    2:{'ylim': [-7000, 7000]},
                    3:{'ylim': [-7000, 7000]},
                    4:{'ylim': [-2600, 2600]},
                },
            },
            'mag': {
                'z': {
                    1:{'clim': [0, 3]},
                    2:{'clim': [0, 3]},
                    3:{'clim': [0, 3]},
                },
                'zfft': {
                    1:{'ylim': [0, 100]},
                    2:{'ylim': [0, 100]},
                    3:{'ylim': [0, 100]},
                },
                'y': {
                    1:{'clim': [0, 5]},
                    2:{'clim': [0, 5]},
                    3:{'clim': [0, 5]},
                },
                'yfft': {
                    1:{'ylim': [0, 100]},
                    2:{'ylim': [0, 100]},
                    3:{'ylim': [0, 100]},
                },
            },
            'pd': {
                1:{'ylim': [0, 0.08]},
                2:{'ylim': [0, 0.08]},
                3:{'ylim': [0, 0.08]},
            },
            'anem': {
                'wind': {
                    1:{'ylim': [-1, 1]},
                    3:{'ylim': [-1, 1]},
                },
                'temp': {
                    1:{'ylim': [22, 24]},
                    3:{'ylim': [22, 24]},
                },
            },
        },
    },
    'poly': {
        'air': {
            'bell': {
                'voltage': {
                    1:{'ylim': [-2, 2]},
                    2:{'ylim': [-2, 2]},
                    3:{'ylim': [-2, 2]},
                    4:{'ylim': [-2, 2]},
                },
                'dedt': {
                    1:{'ylim': [-500000, 500000]},
                    2:{'ylim': [-500000, 500000]},
                    3:{'ylim': [-500000, 500000]},
                    4:{'ylim': [-500000, 500000]},
                },
                'E': {
                    1:{'ylim': [-50000, 50000]},
                    2:{'ylim': [-50000, 50000]},
                    3:{'ylim': [-50000, 50000]},
                    4:{'ylim': [-50000, 50000]},
                },
            },
            'mag': {
                'z': {
                    1:{'clim': [0, 10]},
                    2:{'clim': [0, 10]},
                    3:{'clim': [0, 10]},
                },
                'y': {
                    1:{'clim': [0, 10]},
                    2:{'clim': [0, 10]},
                    3:{'clim': [0, 10]},
                },
            },
            'pd': {
                1:{'ylim': [0, 0.08]},
                2:{'ylim': [0, 0.08]},
                3:{'ylim': [0, 0.08]},
            },
            'anem': {
                'wind': {
                    1:{'ylim': [-1, 1]},
                    2:{'ylim': [-1, 1]},
                    3:{'ylim': [-1, 1]},
                },
                'temp': {
                    1:{'ylim': [20, 30]},
                    2:{'ylim': [20, 30]},
                    3:{'ylim': [20, 30]},
                },
            },
        },
    },
    'sand': {
        'air': {
            'bell': {
                'voltage': {
                    1:{'ylim': [-2, 2]},
                    2:{'ylim': [-2, 2]},
                    3:{'ylim': [-2, 2]},
                    4:{'ylim': [-2, 2]},
                },
                'dedt': {
                    1:{'ylim': [-500000, 500000]},
                    2:{'ylim': [-500000, 500000]},
                    3:{'ylim': [-500000, 500000]},
                    4:{'ylim': [-500000, 500000]},
                },
                'E': {
                    1:{'ylim': [-50000, 50000]},
                    2:{'ylim': [-50000, 50000]},
                    3:{'ylim': [-50000, 50000]},
                    4:{'ylim': [-50000, 50000]},
                },
            },
            'mag': {
                'z': {
                    1:{'clim': [0, 10]},
                    2:{'clim': [0, 10]},
                    3:{'clim': [0, 10]},
                },
                'y': {
                    1:{'clim': [0, 10]},
                    2:{'clim': [0, 10]},
                    3:{'clim': [0, 10]},
                },
            },
            'pd': {
                1:{'ylim': [0, 0.08]},
                2:{'ylim': [0, 0.08]},
                3:{'ylim': [0, 0.08]},
            },
            'anem': {
                'wind': {
                    1:{'ylim': [-1, 1]},
                    2:{'ylim': [-1, 1]},
                    3:{'ylim': [-1, 1]},
                },
                'temp': {
                    1:{'ylim': [20, 30]},
                    2:{'ylim': [20, 30]},
                    3:{'ylim': [20, 30]},
                },
            },
        },
    },
}

electrometer_background = { 
    1: -0.0106,
    2: 0.0069,
    3: 0.0200,
    4: -0.0075,
}
# -----------------------------------------------------------------------------

particle = "empty"
airflow = "air"

folder = "Runs Data/"

sensors = ['BE1','BE2','BE3','BE4','PD1','PD2','PD3','ANTI_wind_speed',
    'ANTI_temperature','ANTO_wind_speed','ANTO_temperature','ANV_wind_speed',
    'ANV_temperature','IFG1_Z','IFG2_Z','IFG3_Z','IFG1_Y','IFG2_Y','IFG3_Y']

bells1 = True
bells2 = True
spike = False

mag_z = True
mag_y = True

pds = True
anems = True


apply_lims = True
background_red=True
apply_highpass=False
apply_notch=True
fft_flag=False

save_string = f"{particle}_{airflow}"

labels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)', '(q)', '(r)']

if bells1 or bells2:
    plot_electrometer(bells1, bells2, apply_lims, background_red)
    if bells1:
        plt.figure(1)
        plt.savefig(folder + save_string + "_bell1.eps")
        plt.savefig(folder + save_string + "_bell1.png")
        plt.clf()

    if bells2:
        plt.figure(2)
        plt.savefig(folder + save_string + "_bell2.eps")
        plt.savefig(folder + save_string + "_bell2.png")
        plt.clf()

    if spike:
        plt.figure(10)
        plt.savefig(folder + save_string + "_spike.eps")
        plt.savefig(folder + save_string + "_spike.png")
        plt.clf()

if pds:
    plot_photodiode(pds, apply_lims)
    if pds:
        plt.figure(3)
        plt.savefig(folder + save_string + "_pd.eps")
        plt.savefig(folder + save_string + "_pd.png")
        plt.clf()

if anems:
    plot_anemometer(anems, apply_lims)
    if anems:
        plt.figure(4)
        plt.savefig(folder + save_string + "_wind.eps")
        plt.savefig(folder + save_string + "_wind.png")
        plt.clf()

if mag_z or mag_y:
    plot_magnetometer(mag_z, mag_y, apply_lims, 50, apply_highpass, apply_notch, fft_flag)
    if mag_z:
        plt.figure(5)
        plt.savefig(folder + save_string + "_mag_z.eps")
        plt.savefig(folder + save_string + "_mag_z.png")
        plt.clf()

    if mag_y:
        plt.figure(6)
        plt.savefig(folder + save_string + "_mag_y.eps")
        plt.savefig(folder + save_string + "_mag_y.png")
        plt.clf()
