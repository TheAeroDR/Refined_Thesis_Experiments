from datetime import datetime
import sys
import numpy as np
import pandas as pd
from labjack import ljm
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

# configuration (inputs)
scanRate = 40000 # S/s
scanDuration = 30  # s
csv_path = 'mag3_test.csv'
#selected_channels = ['AIN0', 'AIN1', 'AIN2']
#selected_channels = ['AIN0', 'AIN1', 'AIN2', 'AIN3']
#selected_channels = ['AIN0', 'AIN1', 'AIN2', 'AIN3', 'AIN4', 'AIN5']
selected_channels = ['AIN0', 'AIN1', 'AIN2']
#channel_names = ['PD1', 'PD2', 'PD3']
#channel_names = ['BE1', 'BE2', 'BE3', 'BE4']
#channel_names = ['TInW', 'TInT', 'TOutW', 'TOutT','VertW', 'VertT']
channel_names = ['FGx', 'FGy', 'FGz']

def labmjack_stream(scanRate, scanDuration):
    handle = ljm.openS("T8", "ANY", "480010419")  #T8 device, any connection, my serial number

    info = ljm.getHandleInfo(handle)
    print("Opened a LabJack with Device type: %i, Connection type: %i,\n"
        "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
        (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))

    deviceType = info[0]

    # Stream Configuration
    aScanListNames = ["AIN0", "AIN1", "AIN2", "AIN3", "AIN4", "AIN5", "AIN6", "AIN7"]  # Scan list names to stream
    numAddresses = len(aScanListNames)
    aScanList = ljm.namesToAddresses(numAddresses, aScanListNames)[0]
    scansPerRead = int(scanRate / 2)

    data_list = []

    try:
        # Controls when stream scanning will start
        # 0 = Start when stream is enabled
        # 2000 = Start when DIO_EF0 detects an edge
        # 2001 = Start when DIO_EF1 detects an edge

        ljm.eWriteName(handle, "STREAM_TRIGGER_INDEX", 0)

        # Enabling internally-clocked stream.

        ljm.eWriteName(handle, "STREAM_CLOCK_SOURCE", 0)

        # AIN_RANGE for T8:
        # 0.0=Default → ±11V
        # 11.0 → ±11.0
        # 9.6 → ±9.6
        # 4.8 → ±4.8
        # 2.4 → ±2.4
        # 1.2 → ±1.2
        # 0.6 → ±0.6
        # 0.3 → ±0.3
        # 0.15 → ±0.15
        # 0.075 → ±0.075
        # 0.036 → ±0.036
        # 0.018 → ±0.018

        # STREAM_RESOLUTION_INDEX:
        # 0-16. A value of 0 will instruct the T8 to use the best resolution for the rate specified.

        aNames = ["AIN_ALL_RANGE", "STREAM_RESOLUTION_INDEX"]
        aValues = [11.0, 0]

        # Write the analog inputs ranges and stream resolution configuration.
        numFrames = len(aNames)
        ljm.eWriteNames(handle, numFrames, aNames, aValues)

        # Configure and start stream
        scanRate = ljm.eStreamStart(handle, scansPerRead, numAddresses, aScanList, scanRate)
        print("\nStream started with a scan rate of %0.0f Hz." % scanRate)

        MAX_REQUESTS = np.ceil((scanDuration * scanRate) / scansPerRead)

        print("\nPerforming %i stream reads." % MAX_REQUESTS)
        start = datetime.now()
        totScans = 0
        totSkip = 0  # Total skipped samples

        i = 1
        while i <= MAX_REQUESTS:
            ret = ljm.eStreamRead(handle)

            aData = ret[0]
            scans = len(aData) / numAddresses
            totScans += scans

            # Count the skipped samples which are indicated by -9999 values. Missed
            # samples occur after a device's stream buffer overflows and are
            # reported after auto-recover mode ends.
            curSkip = aData.count(-9999.0)
            totSkip += curSkip


            arr = np.array(aData, dtype=float).reshape(-1, numAddresses)

            data_list.append(arr)

            print("  Scans Skipped = %0.0f, Scan Backlogs: Device = %i, LJM = "
                "%i" % (curSkip/numAddresses, ret[1], ret[2]))
            i += 1

        end = datetime.now()

        print("\nTotal scans = %i" % (totScans))
        tt = (end - start).seconds + float((end - start).microseconds) / 1000000
        print("Time taken = %f seconds" % (tt))
        print("LJM Scan Rate = %f scans/second" % (scanRate))
        print("Timed Scan Rate = %f scans/second" % (totScans / tt))
        print("Timed Sample Rate = %f samples/second" % (totScans * numAddresses / tt))
        print("Skipped scans = %0.0f" % (totSkip / numAddresses))
    except ljm.LJMError:
        ljme = sys.exc_info()[1]
        print(ljme)
    except Exception:
        e = sys.exc_info()[1]
        print(e)

    try:
        print("\nStop Stream")
        ljm.eStreamStop(handle)
    except ljm.LJMError:
        ljme = sys.exc_info()[1]
        print(ljme)
    except Exception:
        e = sys.exc_info()[1]
        print(e)

    # Close handle
    ljm.close(handle)
    return data_list, aScanListNames, scanRate

def save_data_to_csv(data_list, aScanListNames, selected_channels, scanRate, csv_path, channel_names):
    if len(data_list) > 0:
        all_data = np.vstack(data_list)

        df = pd.DataFrame(all_data, columns=aScanListNames)

        df = df[selected_channels]

        try:
            rate = float(scanRate)
            df.insert(0, 'time_s', np.arange(len(df)) / rate)
        except Exception:
            pass

        extra_lines = ["Sample rate: {} S/s".format(scanRate)]
        with open(csv_path, "w", encoding="utf-8") as f:
            for line in extra_lines:
                f.write(line + "\n")
            # Append DataFrame
            df.to_csv(f, index=False, header=[df.columns.tolist()[0]] + channel_names)
            print(f"Saved streamed data to {csv_path} (rows={len(df)}, cols={len(df.columns)})")
    else:
        print("No data collected; DataFrame not created.")
    return df

def plot_data(df, selected_channels):
    if len(df) > 0:
        for i, channel in enumerate(selected_channels):
            fig = plt.subplot(len(selected_channels),1,i+1)
            plt.plot(df['time_s'], df[channel], label=channel)
        plt.show()

if __name__ == "__main__":
    data_list, aScanListNames, scanRate = labmjack_stream(scanRate, scanDuration)
    df = save_data_to_csv(data_list, aScanListNames, selected_channels, scanRate, csv_path, channel_names)
    plot_data(df, selected_channels)