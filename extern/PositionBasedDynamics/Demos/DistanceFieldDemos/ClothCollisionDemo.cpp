#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "GL/glut.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include <iostream>
#include "Demos/Visualization/Visualization.h"
#include "Simulation/DistanceFieldCollisionDetection.h"
#include "Utils/OBJLoader.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Demos/Common/DemoBase.h"
#include "Demos/Common/TweakBarParameters.h"
#include "Simulation/Simulation.h"


// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
	#define new DEBUG_NEW 
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

void timeStep ();
void buildModel ();
void createMesh();
void render ();
void reset();
void initParameters();
void TW_CALL setBendingMethod(const void *value, void *clientData);
void TW_CALL getBendingMethod(void *value, void *clientData);
void TW_CALL setSimulationMethod(const void *value, void *clientData);
void TW_CALL getSimulationMethod(void *value, void *clientData);


const int nRows = 50;
const int nCols = 50;
const Real width = 10.0;
const Real height = 10.0;
short simulationMethod = 2;
short bendingMethod = 2;
bool doPause = true;
DemoBase *base;
DistanceFieldCollisionDetection cd;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS

	base = new DemoBase();
	base->init(argc, argv, "Cloth demo");

	SimulationModel *model = new SimulationModel();
	model->init();
	Simulation::getCurrent()->setModel(model);

	buildModel();

	initParameters();

	Simulation::getCurrent()->setSimulationMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject(Simulation::getCurrent()->getTimeStep()); });

	// OpenGL
	MiniGL::setClientIdleFunc (50, timeStep);		
	MiniGL::setKeyFunc(0, 'r', reset);
	MiniGL::setClientSceneFunc(render);			
	MiniGL::setViewport (40.0f, 0.1f, 500.0f, Vector3r (5.0, 10.0, 30.0), Vector3r (5.0, 0.0, 0.0));

	TwType enumType2 = TwDefineEnum("SimulationMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "SimulationMethod", enumType2, setSimulationMethod, getSimulationMethod, &simulationMethod, " label='Simulation method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics}' group=Simulation");
	TwType enumType3 = TwDefineEnum("BendingMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "BendingMethod", enumType3, setBendingMethod, getBendingMethod, &bendingMethod, " label='Bending method' enum='0 {None}, 1 {Dihedral angle}, 2 {Isometric bending}' group=Bending");

	glutMainLoop ();	

	base->cleanup ();

	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	delete Simulation::getCurrent();
	delete base;
	delete model;

	return 0;
}

void initParameters()
{
	TwRemoveAllVars(MiniGL::getTweakBar());
	TweakBarParameters::cleanup();

	MiniGL::initTweakBarParameters();

	TweakBarParameters::createParameterGUI();
	TweakBarParameters::createParameterObjectGUI(base);
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent());
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getModel());
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getTimeStep());
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();
	base->getSelectedParticles().clear();

	Simulation::getCurrent()->getModel()->cleanup();

	buildModel();
}

void timeStep ()
{
	const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(DemoBase::PAUSE, true);

	if (base->getValue<bool>(DemoBase::PAUSE))
		return;

	// Simulation code
	SimulationModel *model = Simulation::getCurrent()->getModel();
	const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step(*model);
		STOP_TIMING_AVG;
	}

	for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
		model->getTriangleModels()[i]->updateMeshNormals(model->getParticles());
}

void loadObj(const std::string &filename, VertexData &vd, IndexedFaceMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<OBJLoader::Vec2f> texCoords;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	const unsigned int nTexCoords = (unsigned int)texCoords.size();
	mesh.initMesh(nPoints, nFaces * 2, nFaces);
	vd.reserve(nPoints);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		vd.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nTexCoords; i++)
	{
		mesh.addUV(texCoords[i][0], texCoords[i][1]);
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		int texIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
			if (nTexCoords > 0)
			{
				texIndices[j] = faces[i].texIndices[j] - 1;
				mesh.addUVIndex(texIndices[j]);
			}
		}

		mesh.addFace(&posIndices[0]);
	}
	mesh.buildNeighbors();

	mesh.updateNormals(vd, 0);
	mesh.updateVertexNormals(vd);

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}

void buildModel ()
{
	TimeManager::getCurrent ()->setTimeStepSize (static_cast<Real>(0.005));

	createMesh();

	// create static rigid body
	string fileName = FileSystem::normalizePath(base->getDataPath() + "/models/cube.obj");
	IndexedFaceMesh mesh;
	VertexData vd;
	loadObj(fileName, vd, mesh, Vector3r::Ones());

	string fileNameTorus = FileSystem::normalizePath(base->getDataPath() + "/models/torus.obj");
	IndexedFaceMesh meshTorus;
	VertexData vdTorus;
	loadObj(fileNameTorus, vdTorus, meshTorus, Vector3r::Ones());

	SimulationModel *model = Simulation::getCurrent()->getModel();
	SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
	rb.resize(2);

	// floor
	rb[0] = new RigidBody();
	rb[0]->initBody(1.0,
		Vector3r(0.0, -5.5, 0.0),
		Quaternionr(1.0, 0.0, 0.0, 0.0),
		vd, mesh,
		Vector3r(100.0, 1.0, 100.0));
	rb[0]->setMass(0.0);

	// torus
	rb[1] = new RigidBody();
	rb[1]->initBody(1.0,
		Vector3r(5.0, -1.5, 5.0),
		Quaternionr(1.0, 0.0, 0.0, 0.0),
		vdTorus, meshTorus,
		Vector3r(2.0, 2.0, 2.0));
	rb[1]->setMass(0.0);
	rb[1]->setFrictionCoeff(static_cast<Real>(0.1));

	Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);
	cd.setTolerance(static_cast<Real>(0.05));

	const std::vector<Vector3r> *vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert1 = static_cast<unsigned int>(vertices1->size());
	cd.addCollisionBox(0, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices1)[0], nVert1, Vector3r(100.0, 1.0, 100.0));

	const std::vector<Vector3r> *vertices2 = rb[1]->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert2 = static_cast<unsigned int>(vertices2->size());
	cd.addCollisionTorus(1, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices2)[0], nVert2, Vector2r(2.0, 1.0));
	
	SimulationModel::TriangleModelVector &tm = model->getTriangleModels();
	ParticleData &pd = model->getParticles();
	for (unsigned int i = 0; i < tm.size(); i++)
	{
		const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
		unsigned int offset = tm[i]->getIndexOffset();
		tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
		cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
	}
}


void render ()
{
	base->render();
}


/** Create a particle model mesh 
*/
void createMesh()
{
	TriangleModel::ParticleMesh::UVs uvs;
	uvs.resize(nRows*nCols);

	const Real dy = width / (Real)(nCols-1);
	const Real dx = height / (Real)(nRows-1);

	Vector3r points[nRows*nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			const Real y = (Real)dy*j;
			const Real x = (Real)dx*i;
			points[i*nCols + j] = Vector3r(x, 1.0, y);

			uvs[i*nCols + j][0] = x/width;
			uvs[i*nCols + j][1] = y/height;
		}
	}
	const int nIndices = 6 * (nRows - 1)*(nCols - 1);

	TriangleModel::ParticleMesh::UVIndices uvIndices;
	uvIndices.resize(nIndices);

	unsigned int indices[nIndices];
	int index = 0;
	for (int i = 0; i < nRows - 1; i++)
	{
		for (int j = 0; j < nCols - 1; j++)
		{
			int helper = 0;
			if (i % 2 == j % 2)
				helper = 1;

			indices[index] = i*nCols + j;
			indices[index + 1] = i*nCols + j + 1;
			indices[index + 2] = (i + 1)*nCols + j + helper;

			uvIndices[index] = i*nCols + j;
			uvIndices[index + 1] = i*nCols + j + 1;
			uvIndices[index + 2] = (i + 1)*nCols + j + helper;
			index += 3;

			indices[index] = (i + 1)*nCols + j + 1;
			indices[index + 1] = (i + 1)*nCols + j;
			indices[index + 2] = i*nCols + j + 1 - helper;

			uvIndices[index] = (i + 1)*nCols + j + 1;
			uvIndices[index + 1] = (i + 1)*nCols + j;
			uvIndices[index + 2] = i*nCols + j + 1 - helper;
			index += 3;
		}
	}

	SimulationModel *model = Simulation::getCurrent()->getModel();
	model->addTriangleModel(nRows*nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);
	
	ParticleData &pd = model->getParticles();
	for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
	{
		pd.setMass(i, 1.0);
	}

	// init constraints
	for (unsigned int cm = 0; cm < model->getTriangleModels().size(); cm++)
	{
		if (simulationMethod == 1)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			const unsigned int nEdges = model->getTriangleModels()[cm]->getParticleMesh().numEdges();
			const IndexedFaceMesh::Edge *edges = model->getTriangleModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;

				model->addDistanceConstraint(v1, v2);
			}
		}
		else if (simulationMethod == 2)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh &mesh = model->getTriangleModels()[cm]->getParticleMesh();
			const unsigned int *tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				model->addFEMTriangleConstraint(v1, v2, v3);
			}
		}
		else if (simulationMethod == 3)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh &mesh = model->getTriangleModels()[cm]->getParticleMesh();
			const unsigned int *tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				model->addStrainTriangleConstraint(v1, v2, v3);
			}
		}
		if (bendingMethod != 0)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh &mesh = model->getTriangleModels()[cm]->getParticleMesh();
			unsigned int nEdges = mesh.numEdges();
			const TriangleModel::ParticleMesh::Edge *edges = mesh.getEdges().data();
			const unsigned int *tris = mesh.getFaces().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const int tri1 = edges[i].m_face[0];
				const int tri2 = edges[i].m_face[1];
				if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff))
				{
					// Find the triangle points which do not lie on the axis
					const int axisPoint1 = edges[i].m_vert[0];
					const int axisPoint2 = edges[i].m_vert[1];
					int point1 = -1;
					int point2 = -1;
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2))
						{
							point1 = tris[3 * tri1 + j];
							break;
						}
					}
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2))
						{
							point2 = tris[3 * tri2 + j];
							break;
						}
					}
					if ((point1 != -1) && (point2 != -1))
					{
						const unsigned int vertex1 = point1 + offset;
						const unsigned int vertex2 = point2 + offset;
						const unsigned int vertex3 = edges[i].m_vert[0] + offset;
						const unsigned int vertex4 = edges[i].m_vert[1] + offset;
						if (bendingMethod == 1)
							model->addDihedralConstraint(vertex1, vertex2, vertex3, vertex4);
						else if (bendingMethod == 2)
							model->addIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
					}
				}
			}
		}
	}

	LOG_INFO << "Number of triangles: " << nIndices / 3;
	LOG_INFO << "Number of vertices: " << nRows*nCols;

}

void TW_CALL setBendingMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	*((short*)clientData) = val;
	reset();
}

void TW_CALL getBendingMethod(void *value, void *clientData)
{
	*(short *)(value) = *((short*)clientData);
}

void TW_CALL setSimulationMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	*((short*)clientData) = val;
	reset();
}

void TW_CALL getSimulationMethod(void *value, void *clientData)
{
	*(short *)(value) = *((short*)clientData);
}
