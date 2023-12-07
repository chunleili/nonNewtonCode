#include "Simulation.h"

#include "SPlisHSPlasH/Drag/DragForce_Macklin2014.h"
#include "SPlisHSPlasH/Drag/DragForce_Gissler2017.h"

#include "SPlisHSPlasH/Viscosity/Viscosity_Standard.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_XSPH.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Bender2017.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Peer2015.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Peer2016.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Takahashi2015.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Weiler2018.h"

#include "SPlisHSPlasH/Vorticity/VorticityConfinement.h"
#include "SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h"

#include "Elasticity/Elasticity_Becker2009.h"
#include "Elasticity/Elasticity_Peer2018.h"
#include "Elasticity/Elasticity_Kugelstadt2021.h"

#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#ifdef USE_THIRD_PARTY_METHODS
#include "SurfaceTension/SurfaceTension_ZorillaRitter2020.h"
#endif

//MYADD
#include "My/Viscosity_Casson/Viscosity_Casson.h"
#include "My/Diffusion/Coagulation.h"
#include "My/Diffusion/TemperatureDiffusion.h"
#include "My/Diffusion/Melting.h"
#include "My/ShapeMatching/ShapeMatching.h"
#include "My/Plastic/Plastic.h"

#include "My/NonNewton/NonNewton.h"
#include "My/NonNewton/NonNewton_XSPH.h"
#include "My/NonNewton/NonNewton_standard.h"
#include "My/NonNewton/NonNewton_Peer2015.h"
#include "My/NonNewton/NonNewton_Peer2016.h"
#include "My/NonNewton/NonNewton_Bender2017.h"
#include "My/NonNewton/NonNewton_Weiler2018.h"
#include "My/NonNewton/NonNewton_Takahashi2015.h"

using namespace SPH;

void Simulation::registerNonpressureForces()
{
	addDragMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addDragMethod("Macklin et al. 2014", DragForce_Macklin2014::creator);
	addDragMethod("Gissler et al. 2017", DragForce_Gissler2017::creator);

	addElasticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addElasticityMethod("Becker et al. 2009", Elasticity_Becker2009::creator);
	addElasticityMethod("Peer et al. 2018", Elasticity_Peer2018::creator);
	addElasticityMethod("Kugelstadt et al. 2021", Elasticity_Kugelstadt2021::creator);

	addSurfaceTensionMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addSurfaceTensionMethod("Becker & Teschner 2007", SurfaceTension_Becker2007::creator);
	addSurfaceTensionMethod("Akinci et al. 2013", SurfaceTension_Akinci2013::creator);
	addSurfaceTensionMethod("He et al. 2014", SurfaceTension_He2014::creator);

#ifdef USE_THIRD_PARTY_METHODS
	addSurfaceTensionMethod("Zorilla, Ritter, et al. 2020", SurfaceTension_ZorillaRitter2020::creator);
#endif

	addViscosityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addViscosityMethod("Standard", Viscosity_Standard::creator);
	addViscosityMethod("XSPH", Viscosity_XSPH::creator);
	addViscosityMethod("Bender and Koschier 2017", Viscosity_Bender2017::creator);
	addViscosityMethod("Peer et al. 2015", Viscosity_Peer2015::creator);
	addViscosityMethod("Peer et al. 2016", Viscosity_Peer2016::creator);
	addViscosityMethod("Takahashi et al. 2015 (improved)", Viscosity_Takahashi2015::creator);
	addViscosityMethod("Weiler et al. 2018", Viscosity_Weiler2018::creator);

	addVorticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addVorticityMethod("Micropolar model", MicropolarModel_Bender2017::creator);
	addVorticityMethod("Vorticity confinement", VorticityConfinement::creator);

	//MYADD
	addViscosityMethod("Viscosity_Casson", Viscosity_Casson::creator);
	addSurfaceTensionMethod("Coagulation method", Coagulation::creator);
	addSurfaceTensionMethod("TemperatureDiffusion", TemperatureDiffusion::creator);
	addVorticityMethod("ShapeMatching", ShapeMatching::creator);
	addElasticityMethod("Plastic", Plastic::creator);

	addSurfaceTensionMethod("NonNewton", NonNewton::creator);
	addViscosityMethod("NonNewton_XSPH", NonNewton_XSPH::creator);
	addViscosityMethod("NonNewton_Bender2017", NonNewton_Bender2017::creator);
	addViscosityMethod("NonNewton_Peer2015", NonNewton_Peer2015::creator);
	addViscosityMethod("NonNewton_Peer2016", NonNewton_Peer2016::creator);
	addViscosityMethod("NonNewton_Weiler2018", NonNewton_Weiler2018::creator);
	addViscosityMethod("NonNewton_Takahashi2015", NonNewton_Takahashi2015::creator);

	addSurfaceTensionMethod("Melting", Melting::creator);
}
