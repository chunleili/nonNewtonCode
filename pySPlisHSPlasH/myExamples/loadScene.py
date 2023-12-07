import pysplishsplash as sph
import tkinter as tk
from tkinter import filedialog

base = sph.Exec.SimulatorBase()

tk.Tk().withdraw() # Dont show main window
custom_scene = filedialog.askopenfilename()

base.init(sceneFile=custom_scene)

gui = sph.GUI.Simulator_GUI_imgui(base)
base.setGui(gui)
base.run()