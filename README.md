# Budding yeast chemostat (byc)

![img](https://www.dropbox.com/s/ii73daki5wtva2f/byc_schematic.png)

## Installation

### Windows

Use powershell (pwsh.exe, not cmd.exe)

```sh
# Change to C:\ directory and clone in repository. If this doesn't work,
# you made need to run powershell in administrator mode
cd C:\
git clone https://github.com/johnpcooper/byc.git
# Make sure virtualenv is installed
pip install virtualenv
pip --version virtualenv
# Move to the byc project folder and create the byc environment
# using virtualenv
cd .\byc\
python -m venv .byc
# Activate .byc environment and the install byc requirements to 
# the byc environment
.\.byc\Scripts\activate
pip install -r .\requirements.txt
# Add .byc to the list of envs that jupyter notebook can access
python -m ipykernel install --user --name=.byc
# Install the byc package to .byc environment
pip install .
```

### Linux

This is specifically tested in the Windows subsystem for linux, so there may be quirks. To enter linux in Windows 10, run `bash` in a powershell or cmd prompt.

```sh
# Go to your home directory and clone the byc repository there
cd ~
git clone https://github.com/johnpcooper/byc.git
# update apt and install pip3 and venv
sudo apt update
sudo apt install python3-pip
sudo apt install virtualenv
# Create and activate byc environment
python3 -m venv .byc_env
source .byc_env/bin/activate
# Install build dependencies to the .byc_env, this could take a while
pip3 install -r requirements.txt
# Add .byc_env to the list of environments accessible in ipykernels (jupyter notebook and ipython)
python -m ipykernel install --user --name=.byc_env
pip3 install .
# Deactivate environment if desired
deactivate
```

This package uses tkinter. To use tkinter with the WSL , you'll need to install [Xming for windows](https://sourceforge.net/projects/xming/) and run Xlaunch (search in start menu) before running byc functionality that uses tkinter. You'll also need to make sure you have a valid DISPLAY environment variable. Once Xlaunch is running with display number set to 0:

```bash
export DISPLAY=:0
```

Then you should be able to use tkinter.

## Test run_alignment on example data:

Example  `alignment_input` directory can be downloaded [here](https://utexas.box.com/s/wzkp7ijc9v3ksvhf9rj5sr564aoa1o3k). When you run the script below, you'll be prompted first to select an input directory (choose one holding just the data downloaded above) and an output directory.

### Windows (powershell)

```sh
cd C:\
# Activate .byc env if not already active
.\byc\.byc\Scripts\activate
# Run the alignment script and follow prompts
python .\byc\bin\run_alignment.py
```

### Linux

```sh
cd ~
cd byc
source .byc_env/bin/activate
python bin/run_alignment
```

Python will then ask how you want to name channels in the experiment (should be the same for linux and windows):

```sh
# Python will prompt you to name channels in the terminal.
# Enter names as shown below
Detected 3 channels
Enter names below:
Name for channel 0: bf
Name for channel 1: dsred
Name for channel 2: yfp
```

Once run_alignment is finished running, the aligned stacks for each channel will be saved in the output folder selected in the file dialogs above.
