import h5py
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, decimate, iirnotch, spectrogram
from scipy.signal.windows import tukey
from scipy.ndimage import uniform_filter1d
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

    # Minor ticks
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,

    # Fonts & LaTeX
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    'axes.labelsize': 30,
    'axes.titlesize': 30,
    'xtick.labelsize': 30,
    'ytick.labelsize': 30,
    'legend.fontsize': 30,

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
    b, a = butter(4, cutoff / (Fs / 2), 'high')
    return filtfilt(b, a, data)

def low_pass(data, Fs, cutoff):
    b, a = butter(4, cutoff / (Fs / 2), 'low')
    return filtfilt(b, a, data)

def fft_single(data, Fs):
    L = len(data)
    Y = fft(data)
    P2 = np.abs(Y / L)
    P1 = P2[: L//2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    f = Fs * np.arange(0, L/2 + 1) / L
    return {"data": data, "Y": Y, "P1": P1, "f": f}

def notch_50_100(Fs, BW):
    W0_50 = 50 / (Fs / 2)
    Q_50 = 50 / BW
    b50, a50 = iirnotch(W0_50, W0_50 / Q_50)

    W0_100 = 100 / (Fs / 2)
    Q_100 = 100 / BW
    b100, a100 = iirnotch(W0_100, W0_100 / Q_100)

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

    fft_single(data, Fs)

    f, t, S = spectrogram(data, window=g, noverlap=L, nfft=Ndft, fs=Fs, scaling='spectrum', mode='complex')

    sc = S * 35

    return {"s": sc, "f": f, "t": t}

def plotspectrogram(input_obj, cutoff,clim):
    mask = input_obj["f"] < cutoff
    f = input_obj["f"][mask]
    t = input_obj["t"]
    z = np.abs(input_obj["s"][mask, :])
    z[z > 50] = 50

    plt.pcolormesh(t, f, z, shading='auto',cmap='viridis',vmin=clim[0],vmax=clim[1])
    plt.ylim([0, cutoff])
    #c = plt.colorbar()

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
def plot_electrometer(figure1=True, figure2=True):
    fig1,_ = plt.subplots(4,3,num=1)
    fig1.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)

    fig2,_ = plt.subplots(4,3,num=2)
    fig2.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)



    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        R_fb = 1.01e9
        eps_0 = 8.8541878188e-12
        r_e = 6e-3


        for s in range(1, 5):
            sensor_name = sensors[s-1]
            sensor_data = data[sensor_name]

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

                ax_left.tick_params(axis='y', colors='#0072BD')
                ax_left.spines['left'].set_color('#0072BD')
                ax_right.tick_params(axis='y', colors='#D85319')
                ax_right.spines['top'].set_visible(False)
                ax_right.spines['bottom'].set_visible(False)
                ax_right.spines['left'].set_visible(False)
                ax_right.spines['right'].set_color('#D85319')

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

            sensor_E = 0.5*(sensor_ER[:-1] + sensor_ER[1:]) * (t[1] - t[0])
            sensor_E = np.insert(sensor_E, 0, np.nan)

            if figure2:
                plt.figure(fig2.number)
                ax = plt.subplot(4, 3, idx)
                ax.plot(t, sensor_E)

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
def plot_photodiode(figure3=True):
    fig3,_ = plt.subplots(3,3,num=3)
    fig3.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}.h5"
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
def plot_anemometer(figure4=True):
    fig4,_ = plt.subplots(2,3,num=4)
    fig4.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)
    zero_wind = [1.2185, np.nan, 1.1848]


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}.h5"
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
def plot_magnetometer(figure5=True, figure6=True):
    # Z (Figure 5)
    fig5, _ = plt.subplots(3,3,num=5)
    fig5.subplots_adjust(left=0.12, right=0.88, top=0.95, bottom=0.1, hspace=0.4, wspace=0.3)


    for i in range(1, 4):
        filename = f"pos_{i}_{particle}_{airflow}.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            sensor_name = sensors[12 + s]
            sensor_b = spectrogram_mag(data[sensor_name], Fs, 'reduced', True, True)
            idx = 9 - 3*s + i
            if figure5:
                plt.figure(fig5.number)
                ax = plt.subplot(3, 3, idx)
                plotspectrogram(sensor_b, 100, axis_limits[particle][airflow]['mag']['z'][s]['clim'])

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
        filename = f"pos_{i}_{particle}_{airflow}.h5"
        data = h5totable(folder + filename)
        Fs = h5py.File(folder + filename, 'r').attrs['sampling_rate']


        for s in range(1, 4):
            sensor_name = sensors[15 + s]
            sensor_b = spectrogram_mag(data[sensor_name], Fs, 'reduced', True, True)
            idx = 9 - 3*s + i
            if figure6:
                plt.figure(fig6.number)
                ax = plt.subplot(3, 3, idx)
                plotspectrogram(sensor_b, 100, axis_limits[particle][airflow]['mag']['y'][s]['clim'])

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
                    1:{'clim': [0, 20]},
                    2:{'clim': [0, 20]},
                    3:{'clim': [0, 20]},
                },
                'y': {
                    1:{'clim': [0, 20]},
                    2:{'clim': [0, 20]},
                    3:{'clim': [0, 20]},
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
# -----------------------------------------------------------------------------

particle = "empty"
airflow = "noair"

folder = "Runs Data/"

sensors = ['BE1','BE2','BE3','BE4','PD1','PD2','PD3','ANTI_wind_speed',
    'ANTI_temperature','ANTO_wind_speed','ANTO_temperature','ANV_wind_speed',
    'ANV_temperature','IFG1_Z','IFG2_Z','IFG3_Z','IFG1_Y','IFG2_Y','IFG3_Y']

bells1 = False
bells2 = False
mag_z = False
mag_y = False
pds = False
anems = True

if bells1 or bells2:
    plot_electrometer(bells1, bells2)
if pds:
    plot_photodiode(pds)
if anems:
    plot_anemometer(anems)
if mag_z or mag_y:
    plot_magnetometer(mag_z, mag_y)

save_string = f"{particle}_{airflow}"

if bells1:
    plt.figure(1)
    plt.savefig(folder + save_string + "_bell1.eps")
    plt.savefig(folder + save_string + "_bell1.png")

if bells2:
    plt.figure(2)
    plt.savefig(folder + save_string + "_bell2.eps")
    plt.savefig(folder + save_string + "_bell2.png")

if pds:
   plt.figure(3)
   plt.savefig(folder + save_string + "_pd.eps")
   plt.savefig(folder + save_string + "_pd.png")

if anems:
   plt.figure(4)
   plt.savefig(folder + save_string + "_wind.eps")
   plt.savefig(folder + save_string + "_wind.png")

if mag_z:
   plt.figure(5)
   plt.savefig(folder + save_string + "_mag_z.eps")
   plt.savefig(folder + save_string + "_mag_z.png")

if mag_y:
   plt.figure(6)
   plt.savefig(folder + save_string + "_mag_y.eps")
   plt.savefig(folder + save_string + "_mag_y.png")