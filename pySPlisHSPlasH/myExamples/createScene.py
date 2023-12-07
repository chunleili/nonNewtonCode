import pysplishsplash
import pysplishsplash.Utilities.SceneLoaderStructs as Scene

def main():
    base = pysplishsplash.Exec.SimulatorBase()
    args = base.init()
    gui = pysplishsplash.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    scene = pysplishsplash.Exec.SceneConfiguration.getCurrent().getScene()
    scene.fluidBlocks.append(Scene.FluidBlock(id='Fluid', box=Scene.Box([0.0, 0.0, 0.0], [0.5, 1.0, 1.0]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    base.run()

if __name__ == "__main__":
    main()