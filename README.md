# Refined_Thesis_Experiment


[![DOI](https://zenodo.org/badge/1099001151.svg)](https://doi.org/10.5281/zenodo.17879197)


The repository contains the code and dataproducts for the final chapter of my thesis, the experimental TERIE apparatus.

It is split into folders for separate aspects of the thesis.

**./3_freq_LEDS** contains the code to make an arduino output 3 square waves on the hardware timer (used to drive photodiode LEDs)

**./Prelminary Sensor Testing** contains the python script to aquire data from labjack T8, as well as the microscope, raw data products and the matlab plotting code for the section outlining the experimental apparatus/

**./Runs Data** contains the data and plotting code (python) for the results of the experiment presented in the thesis chapter.

**./Sand_Camsizer** contains the camsizer data and plotting code for the variety of materials presented in the thesis chapter.

**./Tank_Sim** contains the results of running the field generation simulation model on the tank. Several overrides of parameterisations were used to generate these.

**./helmholtx_openscad** contains the .scad and .stl of the helmholtz coils for 3d printing. The .scad is parameteric so can make helmholtz of any size.

**text_ni.py** is the script used to interface to and collect data from the NI cDAQ chassis (without having to use labview). This saves .h5 files.

**Before committing changes to large files:**

```bash
python3 split.py
```

**After pulling or checking out changes:**

```bash
python3 recombine.py
```