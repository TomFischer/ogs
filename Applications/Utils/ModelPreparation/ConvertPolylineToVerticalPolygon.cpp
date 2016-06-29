/**
 * \brief  Program to convert a polyline to a vertical polygon.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <chrono>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileFinder.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/IO/readGeometryFromFile.h"
#include "GeoLib/IO/writeGeometryToFile.h"
#include "GeoLib/GEOObjects.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

#include "ExtractPolygon.h"

void createSurfaceFromVerticalPolygon(GeoLib::GEOObjects & geo_objs,
    GeoLib::Polygon const* polygon, std::string & sfc_name)
{
    // create copy of polygon points
    std::unique_ptr<std::vector<GeoLib::Point*>> sfc_pnts(
        new std::vector<GeoLib::Point*>);
    const std::size_t n_polygon_pnts(polygon->getNumberOfPoints()-1);
    for (std::size_t k(0); k<n_polygon_pnts; k++) {
        sfc_pnts->push_back(new GeoLib::Point(*(polygon->getPoint(k))));
    }

    // add points to GEOObject object
    geo_objs.addPointVec(std::move(sfc_pnts), sfc_name);
    std::vector<GeoLib::Point*> const& pnts(*geo_objs.getPointVec(sfc_name));

    if (pnts.size() != n_polygon_pnts) {
        geo_objs.removePointVec(sfc_name);
        return;
    }

    // create the surface
    GeoLib::Surface *sfc(new GeoLib::Surface(pnts));
    // deploying the special structure of the polygon to create the surface
    sfc->addTriangle(0, 1, 2);
    sfc->addTriangle(0, 2, n_polygon_pnts-1);
    for (std::size_t k(2); k<n_polygon_pnts/2; k++) {
        // triangulate quadrilateral (n_polygon_pnts-k+1, k, k+1, n_polygon_pnts-k)
        sfc->addTriangle(n_polygon_pnts-k-1, k, k+1);
        sfc->addTriangle(n_polygon_pnts-k-1, k+1, n_polygon_pnts-k);
    }
    std::unique_ptr<std::vector<GeoLib::Surface*>> sfcs(new
        std::vector<GeoLib::Surface*>);
    sfcs->push_back(sfc);
    INFO("Added surface \"%s\" with %d triangles.", sfc_name.c_str(),
        sfc->getNumberOfTriangles());

    std::map<std::string, std::size_t>* sfc_name_map(new std::map<std::string, std::size_t>);
    sfc_name_map->insert(std::pair<std::string, std::size_t>(sfc_name, 0));

    geo_objs.addSurfaceVec(std::move(sfcs), sfc_name, sfc_name_map);
}

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog;

    TCLAP::CmdLine cmd("Program calculates a polygon (that is vertical located) "
        "out of a polyline embedded in a mesh. For this reason the polyline "
        "is projected to the bottom surface and top surface of the mesh. "
        "The resulting two polylines are linked to a closed polyline, i.e., "
        "a polygon.", ' ', "0.1");

    TCLAP::ValueArg<std::string> mesh_arg("m", "mesh",
            "*layered* mesh", true, "", "filename for layered mesh");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> geo_arg("", "geometry-file", "",
            true, "", "name of gli or gml file");
    cmd.add(geo_arg);

    TCLAP::ValueArg<std::size_t> ply_id_upper_arg("", "upper",
        "upper boundary of polyline range", true, 0, "");
    cmd.add(ply_id_upper_arg);

    TCLAP::ValueArg<std::size_t> ply_id_lower_arg("", "lower",
        "lower boundary of polyline range", true, 0, "");
    cmd.add(ply_id_lower_arg);

    TCLAP::ValueArg < std::string > sfc_out_arg(
        "s", "surface-file", "name of file the surface will be written to",
        false, "", "string");
    cmd.add(sfc_out_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    // *** read geometry
    GeoLib::GEOObjects geo_objs;
    GeoLib::IO::readGeometryFromFile(geo_arg.getValue(), geo_objs);

    std::vector<std::string> names;
    geo_objs.getGeometryNames(names);
    std::string const& unique_name(names.back());

    // *** get Polygon
    const std::vector<GeoLib::Polyline*>* plys(geo_objs.getPolylineVec(unique_name));
    INFO("Fetched %d polylines.", plys->size());
    if (!plys) {
        ERR("Could not get vector of polylines.");
        return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();

    MeshLib::ExtractPolygon extract_mesh_nodes(*(mesh.get()));

    // *** generate a polygon from polyline
    std::size_t const ply_id_upper(std::min(plys->size(), ply_id_upper_arg.getValue()+1));
    std::size_t const ply_id_lower(ply_id_lower_arg.getValue());

    std::vector<std::string> sfc_names;
    std::string const prefix_name(BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue()));
    for (std::size_t k(ply_id_lower); k < ply_id_upper; k++) {
        bool closed((*plys)[k]->isClosed());
        if (!closed) {
            std::string ply_name;
            geo_objs.getPolylineVecObj(unique_name)->getNameOfElement((*plys)[k], ply_name);
            if (ply_name.empty())
                ply_name = "ply-"+std::to_string(k);
            INFO("Converting polyline %d (%s) to polygon (closed polyline).", k, ply_name.c_str());
            GeoLib::Polygon* polygon(nullptr);
            extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo_objs,
                unique_name, polygon, 50);

            std::string sfc_name("sfc-" + std::to_string(k));
            if (polygon) {
                INFO("Creating surface \"%s\" from polygon.", sfc_name.c_str());
                createSurfaceFromVerticalPolygon(geo_objs, polygon, sfc_name);

                delete polygon;

                if (geo_objs.getPointVec(sfc_name)) {
                    sfc_names.push_back(sfc_name);
                }
            } else {
                INFO("Surface %s could not be created.", sfc_name.c_str());
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto t = std::chrono::duration<double>(end-start);

    INFO("generated %d surfaces in %f seconds", sfc_names.size(), t);
    for (auto name : sfc_names)
        INFO("\t\"%s\"", name.c_str());

    if (! sfc_names.empty()) {
        std::string sfc_project_name(prefix_name + "-Surfaces");
        if (sfc_names.size() > 1)
            geo_objs.mergeGeometries(sfc_names, sfc_project_name);
        else
            sfc_project_name = sfc_names[0];

        GeoLib::IO::writeGeometryToFile(sfc_project_name, geo_objs, sfc_project_name);
    }

    return EXIT_SUCCESS;
}
