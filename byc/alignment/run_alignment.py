import alignment as a
import tifffile as tf
import tkinter as tk
import tkinter.filedialog as dia
import os

def set_paths():
    
    """ Return input_path and output_path, two directories chosen by the user. 
        Input_path should contain the sorted xy directories for each fov. """
    
    # Choose the directory holding all the fields of view that you'll align
    root = tk.Tk()
    input_path = dia.askdirectory(parent=root,
                                 title='Choose the directory holding the experiment you want to align')
    root.destroy()

    # ask the user where they would like to save the output stacks
    root = tk.Tk()
    output_path = dia.askdirectory(parent=root,
                            title='Choose the directory where you want to save aligned images')
    root.destroy()
    
    return input_path, output_path

expt_path, save_path = set_paths()

# Generate a list of directories, one for each fov directory in the expt_path
def set_fov_dirs():

    fov_dirs = os.listdir(expt_path)
    fov_paths = []

    for directory in fov_dirs:
        fov_paths.append(expt_path + '/' + directory)
    return fov_dirs, fov_paths

def run():

    fov_dirs, fov_paths = set_fov_dirs()
    len_channels, channel_names = a.get_channel_names(fov_paths[0])

    for i in range(0, len(fov_dirs)):

        fov_name = fov_dirs[i]
        fov_path = fov_paths[i]
        print("Aligning %s" % fov_name)

        translated_images_dict = a.align_images(fov_path, channel_names)
        a.save_stacks(translated_images_dict, save_path, fov_name)

if __name__ == "__main__":
  run()