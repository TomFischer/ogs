/**
 * \file
 * \author
 * \date   2013-10-08
 * \brief  Implementation of simple ground water flow process
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstdlib>

// AssemblerLib
#include "AssemblerLib/SerialDenseSetup.h"
#include "AssemblerLib/VectorMatrixAssembler.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ThirdParty/tclap
#include "ThirdParty/tclap/CmdLine.h"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"
#include "Configure.h"

// FileIO
#include "XmlIO/Boost/BoostXmlCndInterface.h"
#include "Legacy/OGSIOVer4.h"
#include "readMeshFromFile.h"
#include "XmlIO/Qt/XmlGmlInterface.h"

// GeoLib
#include "GEOObjects.h"
#include "GeoObject.h"

// MathLib
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

// MeshGeoToolsLib
#include "MeshNodeSearcher.h"
#include "MeshNodesToPoints.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"

// OGS
#include "BoundaryCondition.h"
#include "ProjectData.h"

class LocalGWAssembler
{
public:
	LocalGWAssembler() :
			_m(4, 4)
	{
		_m(0, 0) = 4.0;
		_m(0, 1) = -1.0;
		_m(0, 2) = -2.0;
		_m(0, 3) = -1.0;
		_m(1, 1) = 4.0;
		_m(1, 2) = -1.0;
		_m(1, 3) = -2.0;
		_m(2, 2) = 4.0;
		_m(2, 3) = -1.0;
		_m(3, 3) = 4.0;

		// copy upper triangle to lower to make symmetric
		for (std::size_t i = 0; i < 4; i++)
			for (std::size_t j = 0; j < i; j++)
				_m(i, j) = _m(j, i);

		//_m *= 1.e-11/6.0;
		for (std::size_t i = 0; i < 4; i++)
			for (std::size_t j = 0; j < 4; j++)
				_m(i, j) *= 1.0;
	}

	void operator()(const MeshLib::Element & /*e*/, MathLib::DenseMatrix<double> &localA,
			MathLib::DenseVector<double> & /*rhs*/)
	{
		for (std::size_t i = 0; i < 4; i++)
			for (std::size_t j = 0; j < 4; j++)
				localA(i, j) = _m(i, j);
	}

private:
	MathLib::DenseMatrix<double> _m;
};


int main(int argc, char *argv[])
{
	// logog
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format(new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog_cout(new logog::Cout);
	logog_cout->SetFormatter(*custom_format);

	// tclap
	TCLAP::CmdLine cmd(
			"Simple ground water flow test, reading mesh (only 2d quad elements), geometry and bc and simulate ground water flow",
			' ', "0.1");

	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh", "file name of the mesh", true, "", "string");
	cmd.add(mesh_arg);

	TCLAP::ValueArg<std::string> geometry_arg("g", "geometry", "file name of the geometry", true,
			"", "string");
	cmd.add(geometry_arg);

	TCLAP::ValueArg<std::string> bc_arg("", "boundary_condition",
			"file name of the boundary condition", true, "", "string");
	cmd.add(bc_arg);

	cmd.parse(argc, argv);

	ProjectData project_data;

	// *** read geometry
//	const std::string schema_file(std::string(SOURCEPATH).append("/FileIO/OpenGeoSysGLI.xsd"));
//	FileIO::XmlGmlInterface geo_io(&project_data, schema_file);
//	geo_io.readFile(geometry_arg.getValue());
	std::string unique_name;
	std::vector < std::string > errors;
	if (!FileIO::readGLIFileV4(geometry_arg.getValue(), project_data.getGEOObjects(), unique_name,
			errors)) {
		ERR("Could not read geometry file \"%s\".", geometry_arg.getValue().c_str());
		abort();
	}

	// *** read mesh
	project_data.addMesh(FileIO::readMeshFromFile(mesh_arg.getValue()));

	// *** read boundary conditions
	FileIO::BoostXmlCndInterface xml_io(project_data);
	xml_io.readFile(bc_arg.getValue());

	std::vector<FEMCondition*> bcs(
			project_data.getConditions(FiniteElement::GROUNDWATER_FLOW, unique_name,
					FEMCondition::BOUNDARY_CONDITION));

	std::vector < std::size_t > bc_mesh_node_ids;
	const std::string mesh_name(BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue()));
	MeshGeoTools::MeshNodeSearcher searcher(*project_data.getMesh(mesh_name));
	for (auto it(bcs.cbegin()); it != bcs.cend(); it++) {
		// fetch geometry obj from condition obj
		GeoLib::GeoObject const* geom_obj((*it)->getGeoObj());
		if (dynamic_cast<GeoLib::Point const*>(geom_obj) != nullptr) {
			GeoLib::Point const& pnt(*dynamic_cast<GeoLib::Point const*>(geom_obj));
			bc_mesh_node_ids.push_back(searcher.getMeshNodeIDForPoint(pnt));
		} else {
			if (dynamic_cast<GeoLib::Polyline const*>(geom_obj) != nullptr) {
				GeoLib::Polyline const& ply(*dynamic_cast<GeoLib::Polyline const*>(geom_obj));
				std::vector<std::size_t> const& ids(searcher.getMeshNodeIDsAlongPolyline(ply));
				bc_mesh_node_ids.insert(bc_mesh_node_ids.end(), ids.cbegin(), ids.cend());
			}
		}
	}

	//--------------------------------------------------------------------------
	// Prepare mesh items where data are assigned
	//--------------------------------------------------------------------------
	MeshLib::Mesh const& mesh(*project_data.getMesh(mesh_name));
	const MeshLib::MeshSubset mesh_items_all_nodes(mesh, mesh.getNodes());

	//-------------------------------------------------------------------------
	// Allocate a (global) coefficient matrix, RHS and solution vectors
	//-------------------------------------------------------------------------
	// define a mesh item composition in a vector
	std::vector<MeshLib::MeshSubsets*> vec_comp_dis;
	vec_comp_dis.push_back(new MeshLib::MeshSubsets(&mesh_items_all_nodes));
	AssemblerLib::MeshComponentMap vec1_composition(vec_comp_dis,
			AssemblerLib::ComponentOrder::BY_COMPONENT);

	//--------------------------------------------------------------------------
	// Choose implementation type
	//--------------------------------------------------------------------------
	typedef AssemblerLib::SerialDenseSetup GlobalSetup;
	const GlobalSetup global_setup;

	// allocate a vector and matrix
	typedef GlobalSetup::VectorType GlobalVector;
	typedef GlobalSetup::MatrixType GlobalMatrix;
	std::unique_ptr < GlobalMatrix > A(global_setup.createMatrix(vec1_composition));
	A->setZero();
	std::unique_ptr < GlobalVector > rhs(global_setup.createVector(vec1_composition));
	std::unique_ptr < GlobalVector > x(global_setup.createVector(vec1_composition));

	//--------------------------------------------------------------------------
	// Construct a linear system
	//--------------------------------------------------------------------------
	// create a mapping table from element nodes to entries in the linear system
	std::vector<MeshLib::Element*> const& all_eles = mesh.getElements();
	std::vector < std::vector<std::size_t> > map_ele_nodes2vec_entries;
	map_ele_nodes2vec_entries.reserve(all_eles.size());
	for (auto e = all_eles.cbegin(); e != all_eles.cend(); ++e) {
		std::size_t const nnodes = (*e)->getNNodes();
		std::size_t const mesh_id = mesh.getID();
		std::vector<MeshLib::Location> vec_items;
		vec_items.reserve(nnodes);
		for (std::size_t j = 0; j < nnodes; j++)
			vec_items.emplace_back(mesh_id, MeshLib::MeshItemType::Node, (*e)->getNode(j)->getID());

		map_ele_nodes2vec_entries.push_back(
				vec1_composition.getGlobalIndices<AssemblerLib::ComponentOrder::BY_COMPONENT>(
						vec_items));
	}

	// Local and global assemblers.
	LocalGWAssembler local_gw_assembler;

	typedef AssemblerLib::VectorMatrixAssembler<GlobalMatrix, GlobalVector, MeshLib::Element,
			LocalGWAssembler, MathLib::DenseMatrix<double>, MathLib::DenseVector<double> > GlobalAssembler;

	GlobalAssembler global_assembler(*A.get(), *rhs.get(), local_gw_assembler,
			AssemblerLib::LocalToGlobalIndexMap(map_ele_nodes2vec_entries));

	// Call global assembler for each mesh element.
	global_setup.execute(global_assembler, mesh.getElements());

	//--------------------------------------------------------------------------
	// solve x=A^-1 rhs
	//--------------------------------------------------------------------------
	std::cout << "A=\n";
	for (std::size_t i = 0; i < 30; i++) {
		for (std::size_t j = 0; j < 30; j++)
			std::cout << (*A)(i, j) << " ";
		std::cout << std::endl;
	}

	MathLib::GaussAlgorithm<GlobalMatrix, GlobalVector> ls(*A);
	ls.solve(*rhs, *x);

	std::cout << "x=\n";
	for (std::size_t i = 0; i < 30; i++) {
		std::cout << (*x)[i] << " ";
	}
	std::cout << std::endl;
//	std::string geo_name("bc_as_geo_test");
//	MeshGeoTools::MeshNodesToPoints(*(project_data.getMesh(mesh_name)),
//			b_mesh_node_ids,
//			*(project_data.getGEOObjects()),
//			geo_name);
//
//	FileIO::XmlGmlInterface xml_out(*(project_data.getGEOObjects()));
//	xml_out.setNameForExport(geo_name);
//	xml_out.writeToFile(geo_name+".gml");

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
