/*
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>

#include "ExtractPolygon.h"

#include "GeoLib/Point.h"

#include "MathLib/MathTools.h"

#include "MeshLib/Node.h"

namespace MeshLib
{
ExtractPolygon::ExtractPolygon(MeshLib::Mesh const& mesh) :
    _mesh(mesh)
{
}

std::vector<GeoLib::Point>
ExtractPolygon::getOrthogonalProjectedMeshNodesAlongPolyline(
    GeoLib::Polyline const& polyline, double eps) const
{
    // get all nodes of mesh
    std::vector<MeshLib::Node*> const& nodes(_mesh.getNodes());

    std::size_t number_of_ply_pnts(polyline.getNumberOfPoints());
    if (polyline.isClosed())
        number_of_ply_pnts--;

    std::vector<GeoLib::Point> nodes_as_points;
    // This algorithm returns only mesh nodes that are close (eps) to the end
    // points of the line segments. It does not check mesh nodes along the line segment ;-(
    for (std::size_t k(0); k < number_of_ply_pnts; k++) {
        GeoLib::Point const& ply_pnt(*(polyline.getPoint(k)));
        for (auto pnode : nodes) {
            if (MathLib::sqrDist2d(*pnode, ply_pnt) < eps)
                nodes_as_points.push_back(GeoLib::Point(*pnode, pnode->getID()));
        }
    }

    std::sort(nodes_as_points.begin(), nodes_as_points.end());
    return nodes_as_points;
}

std::vector<GeoLib::Point*> ExtractPolygon::extractTopNodes(
    std::vector<GeoLib::Point> const& points) const
{
    double fp_tol(std::numeric_limits<double>::epsilon());
    std::vector<GeoLib::Point*> top_points;
    std::size_t upper_bound(points.size() - 1);
    for (std::size_t k(0); k < upper_bound; ++k)
    {
        if (MathLib::sqrDist2d(points[k], points[k+1]) > fp_tol)
            top_points.push_back(new GeoLib::Point(points[k]));
    }
    top_points.push_back(new GeoLib::Point(points[upper_bound]));
    return top_points;
}

std::vector<GeoLib::Point*> ExtractPolygon::extractBottomNodes(
    std::vector<GeoLib::Point> const& points) const
{
    double fp_tol(std::numeric_limits<double>::epsilon());
    std::vector<GeoLib::Point*> bottom_points;

    bottom_points.push_back(new GeoLib::Point(points[0]));
    std::size_t const upper_bound(points.size()-1);
    for (std::size_t k(0); k < upper_bound; ++k) {
        if (MathLib::sqrDist2d(points[k], points[k+1]) > fp_tol)
            bottom_points.push_back(new GeoLib::Point(points[k+1]));
    }
    return bottom_points;
}

void ExtractPolygon::getPolygonFromPolyline(const GeoLib::Polyline& polyline,
                                            GeoLib::GEOObjects& geo_obj,
                                            std::string const& name,
                                            GeoLib::Polygon*& polygon,
                                            double eps) const
{
    std::vector<GeoLib::Point> points =
        getOrthogonalProjectedMeshNodesAlongPolyline(polyline, eps);
    std::vector<GeoLib::Point*> top_polygon_pnts = extractTopNodes(points);
    std::vector<GeoLib::Point*> bottom_polygon_pnts = extractBottomNodes(points);

    if (top_polygon_pnts.size() != bottom_polygon_pnts.size() ||
        top_polygon_pnts.size() <= 1)
    {
        WARN(
            "Number of top polygon points (%d) does not match number of bottom "
            "polygon points (%d).",
            top_polygon_pnts.size(), bottom_polygon_pnts.size());
        return;
    }

    // append new points to the end of the points vector
    std::vector<std::size_t> top_ids;
    GeoLib::PointVec & pnt_vec(*(geo_obj.getPointVecObj(name)));
    for (auto p : top_polygon_pnts) {
        top_ids.push_back(pnt_vec.push_back(p));
    }
    std::vector<std::size_t> bottom_ids;
    for (auto p : bottom_polygon_pnts) {
        bottom_ids.push_back(pnt_vec.push_back(p));
    }

    // create (an empty) polygon
    GeoLib::Polyline ply(*(geo_obj.getPointVec(name)));

    std::vector<GeoLib::Point*> const* orig_pnts(pnt_vec.getVector());

    // *** add ids of points to polygon
    // for top polyline sort points along polyline
    std::size_t s (top_ids.size());
    for (std::size_t j(0); j < polyline.getNumberOfPoints(); j++) {
        for (std::size_t k(0); k < s; k++) {
            GeoLib::Point const& pnt(*(*orig_pnts)[top_ids[k]]);
            if (MathLib::sqrDist2d(*polyline.getPoint(j),pnt) < eps) {
                ply.addPoint(top_ids[k]);
                k = s;
            }
        }
    }

    // for bottom polyline sort points along polyline in reverse order
    s = bottom_ids.size();
    for (int j(polyline.getNumberOfPoints() - 1); j > -1; j--) {
        for (std::size_t k(0); k < s; k++) {
            GeoLib::Point const& pnt(*(*orig_pnts)[bottom_ids[k]]);
            if (MathLib::sqrDist2d(*polyline.getPoint(j), pnt) < eps) {
                ply.addPoint(bottom_ids[k]);
                k = s;
            }
        }
    }

    // close polygon
    polygon = new GeoLib::Polygon(ply);
    polygon->addPoint(polygon->getPointID(0));
    polygon->initialise();
}

} // end namespace MeshLib
