Quick Start for HyosPy
=====


The use of HyosPy requires the user to create a new Python script for running simulations, in which all the parameters 
are defined. This run script allows the user to customize the simulation by specifying the parameters. This tutorial 
starts with creating the run script, [hyospy.py](https://github.com/UT-CWE/Hyospy/blob/Hyospy-develop/Hyospy_ensemble/hyospy.py). 

First import the `upper_wrapper` module: 

```python
    from upper_wrapper import upper_wrapper
```

Then initialize the `upper_wrapper` object:

```python
    UW = upper_wrapper()
```

Once the object is initialized, we can set the parameters for a customized Hyospy simulation. The start time and end time
for the hydrodynamic model SUNTANS is set as:

```python
    starttime = '2017-05-29-00'
    endtime = '2017-05-31-00'
```

Then we need to set the time information for the oil spill simulation. Note the start time for the oil spill model needs to 
be later than the start time of the Hydrodynamic model. The period of the oil spill simulations is specified by the user. If 
the running period exceeds the data length, the program will raise an error. The start time and period of the oil spill are

```python
    starttime2 = '2017-05-29-20'  # starttime for oil spill model
    period = 20  # unit: hours
```

The parameter for choosing the hydrodynamic model is called `hydro_model`. There are three options: 'SUNTANS', 'ROMS' and 'BLENDED'.
Each one refers to a hydrodynamic model. Note the 'BLENDED' option will run SUNTANS simulation first and dynamically link SUNTANS
and ROMS model. However, if the output from the hydrodynamic model is already available, running it again is not necessary. Then,
the user can choose whether he needs to run the hydrodynamic model first by modify the parameter `hydro_run`. These two parameters
are set as:

```python
    ## Hydrodynamic model
    UW.hydro_model = 'SUNTANS'  # hydrodynamic model options: 'SUNTANS', 'ROMS', 'BLENDED'
    UW.hydro_run = False  # choose whether or not to run a hydrodynamic model first
```

The output of the ROMS model is downloaded from a [Texas A&M server](http://barataria.tamu.edu:8080/thredds/catalog.html). When
the server gets unstable sometime, the downloading process will terminate. Besides, the download speed depends on the available 
network environment. As the data size is large, downloading ROMS data may cost several hours under a poor Internet speed. HyosPy 
provides two options for sourcing the ROMS data. The parameter `ROMS_datasource` can either be `offline` or `online`. 

```python
    ## ROMS datasource
    UW.ROMS_datasource = 'offline' # ROMS data source options: 'online', 'offline'
```

- By choosing the `online` option, ROMS data will be downloaded during a simulation. If you run an ensemble of simulations, the
data is downloaded multiple times. 

- The `offline` option allows the user to download ROMS data before any simulation starts. This ensures the simulation will not
terminate as a result of failure of accessing the online data. By choosing `offline`, the HyosPy system will assume the output 
of the hydrodynamic model is ready to use and assign `hydro_run = False`. An example of downloading the ROMS output is in 
[ROMS_offline.py](https://github.com/UT-CWE/Hyospy/blob/Hyospy-develop/Hyospy_ensemble/lib/ROMS_offline.py). Or you can write 
a new Python script something like:

```python
    import lib
    from ROMS_offline import offline_download

    starttime='2017-05-29-00'
    endtime='2017-06-02-00'
    roms = offline_download(starttime, endtime)
    outfile = 'DATA/ROMS_offline_data.nc'
    roms.Writefile(outfile)
    roms.Go()
```

There are two drifter models (GNOME and TracPy) integrated in HyosPy. A user can run both or just one of them by specifying the
parameters, `runGNOME` and `runTracPy`. Note GNOME is the primary oil spill model in our system, which means it has more funtions and is available for all three hydrodynamic model options.
The parameters for the GNOME simulation can be modified in the function `make_model` in 
[oilspill_wrapper.py](https://github.com/UT-CWE/Hyospy/blob/Hyospy-develop/Hyospy_ensemble/oilspill_wrapper.py). TracPy is compatiable with curvilinear grid and is only available for `ROMS` and `BLENDED` models.  

```python
    ## Oil spill model
    UW.runGNOME = True
    UW.runTracPy = False
```

The grid for the blended model has a very high resolution with about 1000000 cells. GNOME runs slow on this grid. A user can
choose to subset the grid for GNOME simulation by choosing `gnome_subset = True` and defining the subset domain. An example is given as

```python
    ## GNOME settings
    UW.gnome_subset = True       # For blended current product, "subset=True" makes GNOME run faster
    UW.gnome_bbox = [28.17,-95.53,30.0,-93.9]
```

The probability map is a new feature in this version of HyosPy and is only available for GNOME. The map shows the probability of the oil spill trajectory of one or multiple simulations. This function can be turned on by specifying `probability_map = True`. The parameter `google_earth` determines whether or not the probability map will be imposed on Google Earth. If `google_earth=True`, a kml file will be created and a user can open the Google Earth software and visualize the result. Otherwise, when `google_earth=False`, the map is plotted by a python library, basemap. Note Google Earth requires high performance of the video card. Turning off the `google_earth` usually yields faster visualizations. `mpl` is the parameter that adjust the speed of generating the probability map. Increasing this number can accelerate the process but will lower the accurancy of the map. 

```python
    ## Probability Map
    UW.probability_map = True
    UW.google_earth = False
    UW.mpl = 8 # probability calculation time interval
```

Next step is to define the ensemble parameters, `number` and `interval`. The parameter `number` determines how many simulations you need and `interval` specifies the time interval between each simulation. `interval` has the unit of seconds and a default value of 10800. This is because the wind data is updated every 3 hours. 

```python
    ## ensemble parameters
    UW.number = 1
    UW.interval = 500 # units: seconds
```

The user is encouraged to modify the boundary and initial conditions for running the hydrodynamic model, SUNTANS. When `hydro_model` is set as `'SUNTANS'` or `'BLENDED'`, SUNTANS is a primary model in HyosPy. There are multiple choices for SUNTANS open boundary condition `OBC_opt`, while the most recommended one is `'ROMSFILE'`, that is, using salt and temperature from the outer ROMS model and surface elevation from the real-time observational data. In this situation,  `IC_opt` should be `'ROMS'`. This allows SUNTANS to use the existing coarse-resolution ROMS simulation as its initial condition. 
However, if the ROMS output is not available and SUNTANS is the only hydrodynamic model, the user can always set `OBC_opt = 'file'`. This makes SUNTANS use the surface elevation from the observation only and the boundary tracers are 0. The initial condition is inherently constant, i.e. `IC_opt = 'constant'`. Below is an example. 

```python
    ## SUNTANS IC and BC
    UW.OBC_opt = 'ROMSFILE'  # Type 3 boundary condition option: 'constant',
    #'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'
    UW.IC_opt = 'ROMS' # initial condition options: 'constant', 'ROMS'
```

Finally, as all the parameters are set, a user can call the `upper_wrapper` function to run the system. The initial spill location is assigned to the variable `init_latlon`. Note that if you choose SUNTANS as the only hydrodynamic model, the initial location has to be inside or close to the Galveston Bay, that is, in the SUNTANS domain, otherwise the oil spill model will raise an error. 

```python
UW(starttime, endtime, starttime2, period, init_latlon=[28.353786, -95.315109])  #ROMS domain
```










