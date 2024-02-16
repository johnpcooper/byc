# BYC (Budding Yeast Chemostat)

<p align="center">
<img src="https://www.dropbox.com/s/fsv2bvrxtlpkzte/byc_schematic.png?raw=true" width="300">
</p>

`byc` is a python library for processing, annotating, and analyzing time lapse microscopy of budding and fission yeast trapped in the BYC device illustrated above. Installation and usage are detailed below.

## Install  [`byc`](https://github.com/johnpcooper/byc) and [`imagejpc`](https://github.com/johnpcooper/imagejpc) libraries to set up your environment

### 1. Install and configure `byc`

1. Clone `byc` from repo into your projects directory:

    ```sh
    cd c:/Users/usrname/Projects
    git clone https://github.com/johnpcooper/byc.git
    ```

1. Create a python virtual environment, install required packages, and install byc

    ```sh
    # Install a virtual env manager
    pip install virtualenv
    # Create a directory for virual environments
    cd c/Users/usrname/Projects
    mkdir envs
    cd envs
    # Create the the virtual environment for byc
    python -m venv .byc
    # Activate the .byc virtual environment
    .byc/Scripts/activate # windows
    source .byc/bin/activate # macOS/linux
    # Install python packages required to run byc in your
    # .byc environment
    cd c/Users/usrname/Projects/byc
    # This may take a few minutes
    pip install -r requirements_minimal.txt
    # Add currently active .byc environment to the list
    # of environments accessible in ipykernels ( like jupyter
    # notebooks)
    python -m ipykernel install --user --name=.byc
    # Install the byc library in your currently active 
    # .byc environment using pip
    cd c/Users/usrname/Projects/byc
    pip install .
    # Or install the byc using virtual environment wrapper
    pip install virtualenvwrapper-win # windows
    pip install virtualenvwrapper # macOS/linux
    # Before running virtualenvwrapper CLI tools in your unix system,
    # you'll need to add some paths to your .bash_profile/.bashrc
    # See the docs here: https://virtualenvwrapper.readthedocs.io/en/latest/install.html
    # Add the byc source directory to the list of source paths
    # virtualenvwrapper will look for packages to import
    add2virtualenv .
    ```

### 2. Install and configure `imagejpc`

1. Download [fiji](https://imagej.net/software/fiji/downloads). If you're on windows, extract files and move `fiji.app` folder to your appdata directory (e.g. `c/Users/usrname/AppData/Local`). If you're on macOS, leave the downloaded Fiji shortcut wherever its convenient and access plugins folder etc. byc '^ctrl' + 'left click' on the Fiji shortcut and select 'Show package contents' 

1. Clone `imagejpc` into your projects directory:

    ```sh
    cd c/Users/usrname/Projects/
    git clone https://github.com/johnpcooper/imagejpc
    ``` 

1. Edit hardcoded script location references in `macros/addCell.ijm`. Change the values that `script` and `python` are set to to reflect your local `byc` and `.byc` environment paths
1. Copy all plugin files (`.py` and `.ijm`) from `imagejpc/plugins` and `imagejpc/macros` into the plugins folder of Fiji at `C:/Users/usrname/AppData/Local/Fiji.app/plugins`
1. Copy `imagejpc/macros/addCell.ijm` to your Fiji installation macros folder at `C:/Users/usrname/AppData/Local/Fiji.app/macros` (so that `IJ.runMacroFile('addCell.ijm', arg)` will be able to find it when run from a plugin)
1. Run Fiji with the executable at `C:/Users/usrname/AppData/Local/Fiji.app/ImageJ-win64.exe`. Fiji will automatically find the files you added to its plugins folder above and install them

## Process and analyze BYC data

### 1. Register timelapse output from &micro;manager to correct for stage drift between timepoints

1. If you don't have a dataset to analyze, download an example dataset from this [box link](https://utexas.box.com/s/tzfxumwhate722d4n3x8k273p65xuh9v)

1. Run the `byc` alignment script on raw &micro;manager output

    ```sh
    # activate your .byc environment
    cd C:/Users/usrname/Projects
    envs/.byc/Scripts/activate
    # Run the alignment script. Once run, a GUI window will ask you
    # to choose the directory holding data. In the example dataset,
    # this is 20230126_byc_1
    python byc/bin/align_byc_expt.py # this will take a few minutes per xy position
    ```

1. The above aligned channel stacks are saved in a folder called `output`. You'll now want to copy these stacks into your data directory as detailed below before starting to annotate the data


### 2. Annotate raw data with cell location and bud/cell fission events

#### Create experiment and compartment directories and add channel stacks to compartment directories

1. Add a byc data directory to your `byc` installation location. This is where `byc` modules look for data during annotation and following analysis steps. If you installed `byc` using `add2virtualenv` above, the data directory that comes with the repo is already in the "installation location" and no new data directory needs to be created

    ```sh
    mkdir C:/Users/usrname/Projects/envs/.byc/Lib/site-packages/data
    # if you used add2virtualenv, the data directory will already be at 
    # C:/Users/usrname/Projects/byc/data
    ```

1. Create an experiment directory in the byc data directory created above

    ```sh
    mkdir C:/Users/usrname/Projects/envs/.byc/Lib/site-packages/data/data/20230126_byc
    # or if you used add2virtualenv
    mkdir C:/Users/usrname/Projects/byc/data/20230126_byc
    # 20230126_byc is the experiment name in the example dataset
    ```

1. Create a compartment directory in the experiment directory created above. This will hold all channel stacks corresponding to xy positions collected in the experiment that are within a certain flow compartment of the microfluidic device. The name of the compartment directory needs to identify the strain being imaged and potentially other conditions like small molecule concentration etc.

    ```sh
    # This compartment directory name typically includes the strain number
    # and any other information necessary for describing the strain and 
    # conditions of cells in this compartment
    mkdir C:/Users/usrname/Projects/envs/.byc/Lib/site-packages/data/20230126_byc/20230126_byc_JPC228_UBL-YFP-Su9_BY4741
    # Or if you used add2virtualenv
    mkdir C:/Users/usrname/Projects/byc/data/20230126_byc/20230126_byc_JPC228_UBL-YFP-Su9_BY4741
    ```

1. Copy the channel stacks from `output` to the appropriate compartment directory (`20230126_byc_JPC228_UBL-YFP-Su9_BY4741` in the example)
1. You are now ready start annotating data

### 3. Annotate cell location (crop ROIs) and budding events (bud ROIs)

1. Run Fiji (`C:/Users/usrname/AppData/Local/Fiji.app/ImageJ-win64.exe`)
1. Open the xy stack you want to analyze (`20230126_byc\20230126_byc_JPC228_UBL-YFP-Su9_BY4741\20230126_byc_xy09_bf_stack.tif` inside the byc data directory you created above)
1. For each cell you want to analyze:
    1. Annotate crop ROIs, which are essentially key frames, to track the cell of interest. 
        * Draw a rectangular box with its center somewhere inside the cell at the first frame you want to segment. Press 't' to add this box selection to `RoiManager`. 
        * Scroll through the stack until you need to adjust the box selection to keep its center within the cell of interest. When needed, move the box so its center is within the cell of interest and press 't' to add the box to `RoiManager`
        * Add another ROI like above at the last frame you want to segment
        * Make sure that your ROIs in `RoiManager` are sorted by position (select all ROIs, right click, click 'sort')
        * Press 'L' to focus the Fiji search bar and type 'save cell roi set` and press enter
        * When prompted, enter ROI set type as "crop". Enter other annotation information as relevant
    1. Annotate bud ROIs, which are rectangular selections marking the frame before a bud first becomes visible (for cerevisiae data) or the frame at which the vertical fission line first becomes visible (for pombe data). These frame of interest are referred to as the "bud frame":
        * For each bud frame (as described above) draw a rectangular box around the area where the budding or fission event occurs. If analyzing budding yeast data, this frame will be used to annotate the shape of the daughter cell that came before the daughter appearing at the bud frame. If the previous daughter was round, press "6" to add the bud frame ROI to `RoiManager`. If the previous daughter was elongated, press "7" to add the bud frame ROI to `RoiManager`. If you're annotating fission yeast or daughter isn't relevant, simply press 't' to add an unlablled bud ROI to `RoiManager`
        * Once you have created an ROI for each bud frame as described above, add one more frame annotating the end of our observation of the current cell. If the cell escapes or dies, create an ROI at the frame in which the cell was last seen alive/in its catch tube. If the cells is still alive in the device when the experiment ends, create an ROI at the last frame of the stack
        * Sort the ROIs in `RoiManager`, press 'L' to focus the search bar, and enter 'save cell roi set'
        * When prompted, enter 'bud' as ROI set type. For 'end_event_type', enter 'death' if the cell dies during data collection, enter 'escape' if the cell escaped before dying or the end of the experiment, or enter 'sen' if the cell was still alive at the end of the experiment

### 4. Segment cells and measure fluorescence using annotation data created above

This step is executed in the jupyter notebook at [byc/notebooks/BYC_pipeline.ipynb](https://github.com/johnpcooper/byc/blob/master/notebooks/BYC_pipeline.ipynb). This  notebook details the above installation steps and executes cell segmentation, curve fitting to fluorescence measurements, and writing of aggregated single cell data (fluorescence and bud timing) to .csv files.