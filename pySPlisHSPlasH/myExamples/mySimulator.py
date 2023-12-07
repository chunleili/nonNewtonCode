import pysplishsplash
import pysplishsplash.Utilities.SceneLoaderStructs as Scene

def main():
    base = pysplishsplash.Exec.SimulatorBase()
    args = base.init()
    gui = pysplishsplash.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    scene = base.getScene()
    add_block = Scene.FluidBlock('Fluid', Scene.Box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 0, [0.0, 0.0, 0.0])
    scene.fluidBlocks[1] = add_block # In Place construction not supported yet
    base.run()

if __name__ == "__main__":
    main()