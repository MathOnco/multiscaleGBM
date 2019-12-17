# multiscaleGBM

This is the code used to produce the results in [From cells to tissue: How cell scale heterogeneity impacts glioblastoma growth and treatment](https://www.biorxiv.org/content/early/2019/05/26/650150.1.full.pdf) and the [interactive website](http://jillagal.github.io/multiscaleGBM/) with representative movies from this code. In this work, we ask how heterogeneity impacts glioblastoma (GBM) tumor growth and response to treatment, using a hybrid agent-based model fit to multiscale data from a rat experimental model system. 

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
You will need Java to run this program. It is known to work with:

```
Java(TM) SE Runtime Environment (build 1.7.0_79-b15)
Command line or other IDE software (we use IntelliJ IDEA)
```

### Install & Run

First download a ZIP file of the repository contents to your local machine. Double-click to unZIP. 

We will show here two ways to get the program up and running.

#### 1) Command line

Go to the directory folder from the command line, and run the shell script*:

```
cd local-directory/multiscaleGBM-master/
./exe.sh
```

*You may need to change the file permissions before running the script in order to make it executable:

```
chmod +x exe.sh
```

#### 2) IntelliJ IDEA

From IntelliJ, open the program

```
Click Open and select the local multiscaleGBM-master folder.
```

Ensure that the BashSupport plugin is installed. 

```
Click on the IntelliJ IDEA tab and select Preferences.
Go to the Plugins tab and click on Plugins from the menu on the left.
Click the Browse repositories button toward the bottom of the panel.
Type BashSupport into the search bar. 
Click INSTALL, then click Restart IntelliJ IDEA.
```

Set up the bash script to run the program.

```
Click the Run tab, and click Edit Configurations. 
Click the + sign, and select Bash.
Fill out the following run configurations.
  Name: any name
  Script: find and select the local exe.sh file
  Interpreter path: you may need edit this to find your Java interpreter, we used /usr/bin/env
  Interpreter options: bash
  Working directory: browse to and select local-directory/multiscaleGBM-master/source
Click OK.
```

It should now be set up, and will run when you click the green arrow.

### Parameter file
The provided params.txt file (within the input folder) contains the 5 parameters sets used for each example in the corresponding paper. This is called from the exe.sh script, and the order follows the same as in Table 1 and Table S2. The parameters have been scaled to integers in the params.txt file with respect to their values in Table S2 (scaling can be found in the setPars() function of the World.class file). 

### Other input files
The input folder contains several files needed to execute the code. The noMatter.txt and whiteMatter.txt files contain arrays of 0's and 1's that represent whether there is any brain tissue for the former or whether there is white matter specfically for the latter. These are both used to define the tissue structure. The persGo.dat and persStop.dat files contain the persistence times for moving and stopping, respectively, that were found from the data.  

## Basic Code Structure
This section describes the general code structure to get you oriented, and each class file is expanded upon below. The general idea for the flow of operations is given in the diagram below.

![Alt text](flow.png?raw=true "Title")

This code uses a bash script (exe.sh) to compile and execute the main java program (Main.java). The Main class initializes the program and calss functions from the World class, where the tissue field is setup and defined (using functions from Field.class), cells are stored (as an arrayList containing instances of the Cell.class), the frameUpdate is defined, and the interactions between cells and field occur. The Cell class defines attributes and functions for a cell. The Field class defines attributes and functions for the gray/white matter and the PDGF field on a hexagonal grid. The Pars class defines the parameters used in the simulation. The Data class contains variables and calculations for the data output, and Track class defines variables and functions related to storing tracks of single cells. The Treatment class contains the functions needed for the given treatments, and the Functions class contains several generic functions used in the code.

### exe.sh
This is the bash script that deletes any old class files, compiles the java files, and passes the defined arguments and those imported from the params.txt file to run the java code. The tumor type (corresponding to 1-5 from the params.txt file) and the treatment type (0-3) can be changed. Also, the output for movies, tracks, and phenotypes can be toggled on and off. These can be modified in the exe1.sh file, and the numbers used in the paper are given there in the comments.

### Main.java
This is the main source file that sets the variables passed from the bash script, and initializes and runs the simulation. The output directories are created here, and the timing, parameters, field, and cells are set. The frame update is called with conditions for ending the simulation. The functions to write the graphics are also called from this class.

### Pars.java
This file contains most of the parameters for the simulation, which are labeled.

### World.java
This file contains most of the functions for updating the simulation. 
* **frameUpdate()**: This function is called in the Main class and proceeds as follows:
  1. At the top of the frame, at certain time points, data is written to files. 
  2. Treatment is applied at the specified time point.
  3. The simulation is killed when reaching the kill time.
  4. The attributes of the cell density and the indexes of active and quiescent cells are recorded.
  5. At a certain time point, cell tracks are setup and data from the previous step are finalized. 
  6. The cell loop is called, which only includes activated cells.
  7. The diffusion of PDGF is calculated over the time step.
  8. Any cells marked for death are removed.
  9. Metrics are updated and recorded.
* **cellLoop()**: This function loops through just the activated cells, keeping track of midX and midY, which helps determine the center position of the tumor to get the radius.
  1. If the cell is not already set to die, check to see if cell is done with cell cycle. If so, the cell will divide. Otherwise, the cell will update migration.
  2. Each cell updates the field by consuming and/or secreting PDGF.
  3. Find the (x,y) center of mass for the tumor by dividing the sum of the positions by the number of cells.
* **resetDensity()**: This resets the mesh and metrics and loops through the cells to calculate the density with respect to the mesh. Cells can be added to the dead list, the activated list, and set to be quiescent here. This function also calculates the maximum Ki67 activity mesh point.
* **division()**: This function calls findNewPosition() to create a new cell at an angle a diameter away from the parental cell. It also assigns new trait values and resets other variables using the divNewParams() function found in the Cell class.
* **killCells()**: This function first sorts the death list in ascending order and then calls shiftTrackIndexes() to shift the indexes so that the cell tracks list maintains the correct indexes once they are removed. Then the dead cells are removed.
* **initializeMets()**: This function defines which cells will be tracked. Only cells that are labeled, not quiescent, not on the death list, and within a density greater than that defined by the parameter rimPercent in Pars.class. 
* **setupTracksCollection()**: If time to collect tracks data, tracksOn is flagged true, and the frame number is recorded. If the individual data is to be collected, setTracks() is called, and Pars.collect is flagged true.
* **metricsCollection()**: This records tumor size data, and individual cell migration data. If reached the end of data collection time period, it is turned off and the collection time is updated to the next time point.
* **collectData()**: This calculates trait distributions and infected/recruited cell ratios.
* **finalizeData()**: This finalizes the data in terms of the data collected in the experiment and writes to files.
* **findMaxKi67()**: This finds the mesh area with the maximum activity of Ki67 and prints it to a file.
* **setPars()**: This function takes in the values from the params.txt file and defines the parameters with consistent units. It also imports the tissue field values and the persistence times.
* **setTracks()**: This takes all cells in the trackable list, shuffles the indexes, and takes up to Pars.maxTracks number of cells to track.
* **setField()**: The gray/white matter, the mesh density, and the carrying capacity (based on the gray/white matter) is defined here.
* **setCells()**: The cells are set up here, as far as position and starting values for traits.

### Field.java
The hexagonal field keeps track of the white/gray matter and PDGF concentration. The indexes are labeled according to the picture below (left). The grayed out regions are not included in any calculations. Neighborhoods for each lattice point for diffusion are defined as shown (right).

![Hexagonal field setup.](hex.png?raw=true "Hexagonal field setup.")

* **setField()**: This function defines the field according to the imported rat atlas attrributes (converting from the rectangular to hex mesh) and sets the initial PDGF concentration, which is concentrated at the initial tumor site.
* **update()**: updates field due to secretion of infected cells and consumption of all activated cells within the nearest hex lattice point of the cell's center.

### Cell.java
This file contains specific functions and attributes for the Cell class.

* **init()**: Defines initial atrributes of cells.
* **divNewParams()**: This function resets attributes for a newly divided cell.
* **envRespDiv()**: This function returns the proliferation rate of a cell dependent on its max potential value and the local environment.
* **envRespSp()**: This function returns the migration rate of a cell dependent on its max potential value and the local environment.
* **envRespWalk()**: This defines the persistence of a cell depending on whether it is stopped or moving and in gray or white matter.
* **findEnvironment()**: This defines the cell's local PDGF concentration and whether it is in white/gray matter.
* **reset()**: This resets a cell's velocity state (moving or stopped) and angle. It also resets its speed, proliferation rate, and persistence according to the local environment.
* **moveIntoTissue()**: If a cell happens to be placed outside of the tissue, this function will be called to look in the local neighborhood to find a position within the tissue.
* **move()**: This function moves the cell according to the defined angle and speed. If the cell moves outside the tissue, moveIntoTissue() is called. If the cell is set to move into an area that is above the carrying capacity and greater than the previous density, then it will remain in the previous place and change its direction of movement. The persistence time is updated.
* **divUpdate()**: If the cell is not quiescent, this updates the cell cycle counter for division. Quiescent cells remain the same.
* **getKi67()**: This defines cells as Ki67+ if it has less than 10 hours left in the cell cycle.

### Treatment.java
This file contains the treatment functions.
* **get()**: Gets each treatment by index.
* **antiPro()**: Anti-proliferative treatment, which targets cells below IMT threshold, not quiescent, and part of the mass.
* **cellSlow()**: Anti-migratory treatment, which slows all cells to 10% of the original value.

### Track.java
This class defines a track for a single cell. Positions, distance traveled (for migration collection), and times spent moving and stopped are collected. 

### Data.java
This file contains variables and functions for calculating averages and standard deviations for traits output into the data folder.
* **setMesh()**: This initializes the mesh to find whether it is mostly white or gary matter.
* **findGW()**: This takes average values from setMesh() and defines the point as gray or white, then sets the max carrying capacity accordingly.
* **populateMesh()**: This populates the mesh grid coordinates with each cells index and records the density. It also sums up the Ki67 activity (to later find the max). Finally, these mesh points will be added to the quiescent list, which just contains mesh points of labeled cells or cells that have at least divided once (cell attribute mass=true). This list is used later to check whether the density has exceeded the carrying capacity. 
* **findAngArea()**: This calculates the area of each radial mesh point from the tumor center.
* **findDistStats()**: This finds the proliferation rate and migration speeds within radial distances of the center of the tumor mass.
* **findIR()**: This calculates the infected to recruited ratio of cells within a section of the tumor core.
* **getDistTrav()**: This calculates the distance travelled for a cell over the tracking period.
* **getM()**: This finalizes the migration and proliferation data to be output.

### Functions.java
Several generic functions are stored here for writing to files and sampling from distributions. The hexagonal diffusion function is also found here, as well as some drawing functions. 

### multiscaleGBM.iml
This is a module file that saves settings from IntelliJ IDEA.

## Contributions & Feedback
Please contact me if you have any suggestions for improvement.

## Authors
* Code, Analysis, Methodology, & Visualization - Jill A. Gallaher 
* Investigation - Jill A. Gallaher, Susan C. Massey, Andrea Hawkins-Daarud, 
* Conceptualization, Writing, & Editing - Jill A. Gallaher, Susan C. Massey, Andrea Hawkins-Daarud, Russell C. Rockne, Orlando Gil, Joe Juliano, Peter Canoll, Kristin R. Swanson, & Alexander R. A. Anderson.
