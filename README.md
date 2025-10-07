# indpy

## Introduction
This is a python package to calculate inductances and mutual inductances for superconductors. It uses the 3D-MLSI software from M. Khapaev. This software combines London and Maxwell equations.


## Installation
- Clone this folder on your computer.
- Copy the 3D-MLSI folder somewhere on your computer (it is in roger/QCMX/Softwares/3D_MLSI_VK).
- Put the 3D-MLSI folder path in the first line of indpy.py. It may cause issues if this path as blank spaces in it.
- Open a python console in the folder containing indpy.py.
- Dependencies: you need the following python packages
  - numpy
  - os
  - subprocess
  - gdspy
- Run ```import indpy```.
- Run ```indpy.get_inductance_from_gds('test')```.
- If you obtain an inductance matrix like this, the installation should be fine: ```{(1, 1): 39.0, (1, 2): 6.415, (2, 2): 24.03}```.


## Usage
- Put the gds file (file_name.gds) of the design you want to simulate in the data folder and run ```indpy.get_inductance_from_gds('file_name')```.
- You obtain a dictionnary whose keys are the inductance paths and the values are the inductance: for instance (1, 1) means the self inductance of path 1 and (1, 2) is the mutual inductance between 1 and 2.
- You can check who is who by running ```indpy.check_mlscs('file_name')``` and click on Do/Results from disk.
  
### Formatting the gds file
- Try to put as few elements as possible.
- Use one layer per conductor. You can specify the vertical position and height for each conductor when you run ```indpy.get_inductance_from_gds```. For that, use the ```conds``` parameters (see below).
- Merge all elements of each layer. The program will crash if two metal polygons (on the same conductor) touch.
- Use layer 99 to define terminals. They define the current paths. To define a terminal, draw a path in klayout. It has to match one side of a polygon. For the program to work, you need to have an even number of terminals: for each pair, the first terminal is the source of current and the second is the drain.


### Parameters for get_inductance_from_gds
- Type ```indpy.get_inductance_from_gds?``` to display help.
- ```name``` is the filename (without extension).
- ```check_conv``` can be set to ```True``` if you want to have a look at the gds -> mlscs conversion. (Close the window to continue)
- ```cell_name``` is the name of the cell to extract from the GDS file. Default is 'TOP'.
- ```lmbd``` is the London penetration depth (in µm). Default is 0.1 µm.
- ```conds``` is used to specify the vertical position and height of each conductor. The format is [[a0, a1], [b0, b1], [c0, c1], ...], where each pair (in µm) corresponds to the bottom and the top of the conductor. If you don't put anything, they will all be set to [0, 0.1]
- ```ah``` is the mesh size parameter for the MLSCS file. Default is 1.0
- ```ahb``` is similar to ```ah```. Default is 0.25. No idea what it does.
- ```verbose``` can be set to ```False``` if you want less text to be printed.


### Possibility to move conductor
The ```move_conductor_and_simulate``` function can be used to move a conductor directly in the MLSCS file and then run the simulation. Very useful if you want to see the effect of the distance between two conductors for instance.
The following code is an example of how do do that:

```
distances = np.linspace(0, 10, 11) # in µm
Ms = np.zeros_like(distances)
for i, d in enumerate(distances):
    indu = indpy.move_conductor_and_simulate('filename', 1, dx=d, dy=0, verbose=False, recalc=True, del_files=True)
    Ms[i] = indu[(1, 2)]
    time.sleep(0.1) # wait a bit to avoid file access conflicts
```
