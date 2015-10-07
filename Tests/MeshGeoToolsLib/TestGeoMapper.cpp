/**
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "../MathLib/AutoCheckTools.h"

#include "FileIO/VtkIO/VtuInterface.h"

namespace ac = autocheck;

struct MeshGeoToolsLibGeoMapper : public ::testing::Test
{
	ac::randomTupleGenerator<double, 3> tuple_gen;
	ac::cons_generator<GeoLib::Point, ac::randomTupleGenerator<double, 3>>
	    points_gen{tuple_gen};

	std::unique_ptr<MeshLib::Mesh> _surface_mesh{
		MeshLib::MeshGenerator::createSurfaceMesh(
			"Test", MathLib::Point3d{{0.0, 0.0, 0.0}},
			MathLib::Point3d{{1.0, 1.0, 0.0}}, {110,60},
				[](double x, double y) { return std::cos(x+y); })};

	ac::gtest_reporter gtest_reporter;
};

TEST_F(MeshGeoToolsLibGeoMapper, PointsOnSurfaceMesh)
{
	auto mapPointsOnMeshSurface = [this](std::vector<GeoLib::Point> & pnts) -> bool
	{
		GeoLib::GEOObjects geo_obj;
		std::string geo_name("TestGeoMapperPoints");
		std::vector<GeoLib::Point*> * points(new std::vector<GeoLib::Point*>);
		for (auto & p : pnts) {
			points->push_back(new GeoLib::Point(p));
		}
		geo_obj.addPointVec(points, geo_name);
		MeshGeoToolsLib::GeoMapper geo_mapper(geo_obj, geo_name);

		geo_mapper.mapOnMesh(_surface_mesh.get());

		double const eps(0.01);
		for (auto pnt : *points) {
			GeoLib::Point const& p(*pnt);
			if (0.0 <= p[0] && p[0] <= 1.0 && 0.0 <= p[1] && p[1] <= 1.0) {
				if (std::abs(std::cos(p[0]+p[1]) - p[2]) >= eps) {
					INFO("std::cos(%f + %f) = %f, %f",
						p[0], p[1], cos(p[0]+p[1]), p[2]);
					return false;
				}
			}
		}
		return true;
	};

	ac::check<std::vector<GeoLib::Point>>(
		mapPointsOnMeshSurface,
		10,
		ac::make_arbitrary(ac::fix(10000,list_of(points_gen))),
		gtest_reporter);
}

