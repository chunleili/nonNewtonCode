import pysplishsplash as sph
import numpy as np

counter = 0
function_list = ['command', 'command2']

def init(base):
    global counter
    print("init test")
    counter = 1
    
def step():
    global counter
    sim = sph.Simulation.getCurrent()
    fluid = sim.getFluidModel(0)
    tm = sph.TimeManager.getCurrent()
    print(fluid.getPosition(0))
    print(tm.getTime())
    print(counter)
    counter += 1
    print("---")
    
def reset():
    print("reset test")
    
def command():
    print("tst cmd")
        
def command2():
    print("tst cmd2")