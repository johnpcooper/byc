# Budding yeast chemostat (byc)

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

This is specifically the Windows subsystem for linux, so there may be quirks. To enter linux in Windows 10, run `bash` in a powershell or cmd prompt.

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

## Test run_alignment on example data:

```sh
cd C:\
# Activate .byc env if not already active
.\byc\.byc\Scripts\activate
# Run the alignment script and follow prompts
python .\byc\bin\run_alignment.py
```

A file dialogue will come up and prompt you to select the input directory. Select `...Finkelstein-Matouschek\byc_data\example_byc_expts\20200221_byc\tifs\alignment_input ` for input folder and then `...Finkelstein-Matouschek\byc_data\example_byc_expts\20200221_byc\tifs\alignment_output` for output folder.

Full paths of input and output directories will depend on where they are on your computer. These are just the paths to the folder within my Box Sync\Finkelstein-Matouschek backup

Python will then ask how you want to name channels in the experiment:

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
