import tkinter as tk
import tkinter.filedialog as dia

def set_fp(prompt):

    """ Return the path to a file of interest. Call with the prompt you would like to 
        display in the dialog box."""

    # create the dialog box and set the fn
    root = tk.Tk()
    fp = dia.askopenfilename(parent=root, title=prompt)
    root.destroy() # very important to destroy the root object, otherwise python 
    # just keeps running in a dialog box

    return fp # return the path to the file you just selected