import nidaqmx
import nidaqmx.stream_readers
import numpy as np  
from datetime import datetime
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Qt5Agg')
import h5py


sampling_rate = 34000 # S/s/channel
duration = 120 # s
#filename = 'pos_1_poly_air_magz.h5'

filename = 'test_ni_output.h5'


with nidaqmx.Task() as task:
    ## bells
    #BE1
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai0", min_val=-10.0, max_val=10.0)
    #BE2
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai1", min_val=-10.0, max_val=10.0)
    #BE3
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai2", min_val=-10.0, max_val=10.0)
    #BE4
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai3", min_val=-10.0, max_val=10.0)

    ## photodiodes
    #PD1
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai4", min_val=-10.0, max_val=10.0)
    #PD2
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai5", min_val=-10.0, max_val=10.0)
    #PD3
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai6", min_val=-10.0, max_val=10.0)

    ## anemometers
    #ANTI - wind speed
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai10", min_val=-10.0, max_val=10.0)
    #ANTI - temperature
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai11", min_val=-10.0, max_val=10.0)
    #ANTO - wind speed
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai12", min_val=-10.0, max_val=10.0)
    #ANTO - temperature
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai13", min_val=-10.0, max_val=10.0)
    #ANV - wind speed
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai14", min_val=-10.0, max_val=10.0)
    #ANV - temperature
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai15", min_val=-10.0, max_val=10.0)

    ## magnetometers
    #IFG1_Z
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD2/ai0", min_val=-10.0, max_val=10.0)
    #IFG2_Z
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD2/ai1", min_val=-10.0, max_val=10.0)
    #IFG3_Z
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD2/ai2", min_val=-10.0, max_val=10.0)
    #IFG1_Y
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai7", min_val=-10.0, max_val=10.0)
    #IFG2_Y
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai8", min_val=-10.0, max_val=10.0)
    #IFG3_Y
    task.ai_channels.add_ai_voltage_chan("cDAQ2MOD1/ai9", min_val=-10.0, max_val=10.0)
    
    task.timing.cfg_samp_clk_timing(rate=sampling_rate, sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)
    

    sampling_rate = task.timing.samp_clk_rate
    total_samples = int(sampling_rate * duration)
    num_channel = len(task.ai_channels)

    chunk_samps = int(sampling_rate)

    task.in_stream.input_buf_size = max(512, chunk_samps * 2)


    n_chunks = int(np.ceil(total_samples / chunk_samps))

    memmap_filename = "acq_data.dat"
    data_mm = np.memmap(memmap_filename, dtype='float64', mode='w+', shape=(num_channel, total_samples))

    reader = nidaqmx.stream_readers.AnalogMultiChannelReader(task.in_stream)
    chunk_buf = np.zeros((num_channel, chunk_samps), dtype='float64')

    start = datetime.now()
    print("Starting continuous acquisition at", start)
    task.start()
    try:
        for i in range(n_chunks):
            to_read = min(chunk_samps, total_samples - i * chunk_samps)
            # read into the front of chunk_buf (reader accepts larger buffer but we'll only use to_read)
            reader.read_many_sample(chunk_buf[:, :to_read],
                                    number_of_samples_per_channel=to_read,
                                    timeout=10.0)
            data_mm[:, i * chunk_samps : i * chunk_samps + to_read] = chunk_buf[:, :to_read]
    finally:
        # ensure task is stopped and memmap flushed
        task.stop()
        data_mm.flush()
        end = datetime.now()
        print("Acquisition finished at", end, "duration:", end - start)


data_array = np.array(data_mm.T)  # shape (samples, channels)

# Create time axis
time_axis = np.arange(total_samples) / sampling_rate

channel_names = ["BE1", "BE2", "BE3", "BE4", "PD1", "PD2", "PD3", 
                "ANTI_wind_speed", "ANTI_temperature", "ANTO_wind_speed", 
                "ANTO_temperature", "ANV_wind_speed", "ANV_temperature", 
                "IFG1_Z", "IFG2_Z", "IFG3_Z","IFG1_Y", "IFG2_Y", "IFG3_Y"]

with h5py.File(filename, 'w') as f:
    # Create datasets
    f.create_dataset('data', data=data_array, compression='gzip')
    f.create_dataset('time', data=time_axis, compression='gzip')
    
    # Store metadata as attributes
    f.attrs['sampling_rate'] = sampling_rate
    f.attrs['duration'] = duration
    f.attrs['total_samples'] = total_samples
    f.attrs['num_channels'] = num_channel
    f.attrs['start_time'] = start.isoformat()
    f.attrs['end_time'] = end.isoformat()
    f.attrs['acquisition_duration'] = str(end - start)
    
    # Store channel names
    f.create_dataset('channel_names', data=[name.encode() for name in channel_names])
    
    # Create individual channel datasets for easier access
    channels_group = f.create_group('channels')
    for i, name in enumerate(channel_names):
        channels_group.create_dataset(name, data=data_array[:, i], compression='gzip')

print(f"Data saved to HDF5: {filename}")

# Plot the data
plt.figure(figsize=(10, 5))
plt.plot(time_axis, data_array)
plt.xlabel("Time [s]")
plt.ylabel("Voltage [V]")
plt.grid(True)
plt.show()
