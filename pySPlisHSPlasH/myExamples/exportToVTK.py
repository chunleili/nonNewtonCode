import pysplishsplash as sph
from pysplishsplash.Extras import Scenes
import os
from pathlib import Path

def main():
    base = sph.Exec.SimulatorBase()
    # output_dir = os.path.abspath("bin/output/pydata/")+ "_" + os.path.basename(__file__)
    output_dir = os.path.abspath("bin/output/py/") +"_"+ Path(__file__).stem
    base.init(useGui=False, outputDir=output_dir, sceneFile=Scenes.DoubleDamBreak)
    base.setValueFloat(base.STOP_AT, 20.0) # Important to have the dot to denote a float
    base.activateExporter("VTK Exporter", True)
    # Uncomment the next line to set the output FPS value (must be float)
    # base.setValueFloat(base.DATA_EXPORT_FPS, 10000.) 
    base.run()

if __name__ == "__main__":
    main()