#include "SimulatorBase.h"
#include "Exporter/ExporterBase.h"
#include "Exporter/ParticleExporter_Partio.h"
#include "Exporter/ParticleExporter_VTK.h"
#include "Exporter/RigidBodyExporter_BIN.h"
#include "Exporter/RigidBodyExporter_OBJ.h"
#include "Exporter/RigidBodyExporter_VTK.h"
#include "Exporter/RigidBodyExporter_PLY.h"
#include "Exporter/ParticleExporter_xyz.h"
#include "Exporter/ParticleExporter_MyPartio.h"

using namespace SPH;

void SimulatorBase::createExporters()
{
	addParticleExporter("enablePartioExport", "Partio Exporter", "Enable/disable partio export.", new ParticleExporter_Partio(this));
	addParticleExporter("enableVTKExport", "VTK Exporter", "Enable/disable VTK export.", new ParticleExporter_VTK(this));
	addParticleExporter("enableXyzExport", "xyz Exporter", "Enable/disable xyz export.", new ParticleExporter_xyz(this));
	addParticleExporter("enableMyPartioExport", "Houdini Exporter(MyPartio)", "Enable/disable my partio export.", new ParticleExporter_MyPartio(this));


	addRigidBodyExporter("enableRigidBodyExport", "Rigid Body Exporter", "Enable/disable rigid body BIN export.", new RigidBodyExporter_BIN(this));
	addRigidBodyExporter("enableRigidBodyOBJExport", "Rigid Body OBJ Exporter", "Enable/disable rigid body OBJ export.", new RigidBodyExporter_OBJ(this));
	addRigidBodyExporter("enableRigidBodyVTKExport", "Rigid Body VTK Exporter", "Enable/disable rigid body VTK export.", new RigidBodyExporter_VTK(this));
	addRigidBodyExporter("enableRigidBodyPLYExport", "Rigid Body PLY Exporter", "Enable/disable rigid body PLY export.", new RigidBodyExporter_PLY(this));
}