/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <sstream>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/IO/readGeometryFromFile.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/AABB.h"

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Program computes the axis aligned bounding box of given points", ' ', "0.1");
    TCLAP::ValueArg<std::string> geo_fname_arg("g", "geometry-file",
                                               "non-empty file name", true, "",
                                               "file name string");
    cmd.add(geo_fname_arg);
    cmd.parse(argc, argv);

    GeoLib::GEOObjects geo_objs;
    GeoLib::IO::readGeometryFromFile(geo_fname_arg.getValue(), geo_objs);

    std::vector<std::string> geo_names;
    geo_objs.getGeometryNames(geo_names);
    std::vector<GeoLib::Point*> const* pnts(
        geo_objs.getPointVec(geo_names.back()));

    GeoLib::AABB aabb(pnts->begin(), pnts->end());
    MathLib::Point3d const& min(aabb.getMinPoint());
    MathLib::Point3d const& max(aabb.getMaxPoint());
    std::stringstream out;
    out.precision(std::numeric_limits<double>::digits10);
    out << "[" << min << ", " << max << "]";
    INFO("AABB: %s", out.str().c_str());

    return EXIT_SUCCESS;
}
