/**
 * \brief  Implements algorithms for extracting information from a mesh.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EXTRACTPOLYGON_H_
#define EXTRACTPOLYGON_H_

#include <iostream>

#include "MeshLib/Mesh.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"

namespace MeshLib
{
/**
 * This class implements an algorithm to extract mesh node ids from a given *extruded* mesh.
 */
class ExtractPolygon
{
public:
    /**
     * constructor - takes an *extruded* mesh
     * @param mesh an instance of class MeshLib::Mesh
     */
    ExtractPolygon(MeshLib::Mesh const& mesh);

    /**
     * Computes the polygon to a given polyline that is consisting of the projection
     * of this polyline to the bottom and the top surface of the mesh and the links between
     * these two polylines.
     * @param polyline the ("defining") polyline
     * @param geo_obj geometric objects manager
     * @param name the name of the group of geometric objects
     * @param polygon pointer to the resulting polygon
     *      warning: the pointer to an already existing polygon will be destroyed
     * @param eps tolerance for Euclidean distance
     */
    void getPolygonFromPolyline(const GeoLib::Polyline& polyline,
                                GeoLib::GEOObjects& geo_obj,
                                std::string const& name,
                                GeoLib::Polygon*& polygon, double eps) const;

private:
    /**
     * Extracts the nodes along a polyline belonging to the bottom surface
     * @param bottom_points the bottom mesh nodes as points
     */
    std::vector<GeoLib::Point*> extractBottomNodes(
        std::vector<GeoLib::Point> const& points) const;

    /**
     * Extracts the nodes along a polyline belonging to the top surface
     * @param points the orthogonal projected nodes along the polyline
     * @return vector of pointers to newly created GeoLib::Point object that
     * belongs to the top surface
     */
    std::vector<GeoLib::Point*> extractTopNodes(
        std::vector<GeoLib::Point> const& points) const;

    /**
     * This method searchs all mesh nodes with the same x and y coordinates
     * like the polyline points. It returns the mesh nodes as points and the ids
     * within objects of type Point. The mesh nodes / points are
     * sorted lexicographically.
     * @param polyline the "defining" polyline
     * @param nodes_as_points vector of GeoLib::Point objects
     * @param eps tolerance for Euclidean distance
     */
    std::vector<GeoLib::Point> getOrthogonalProjectedMeshNodesAlongPolyline(
        GeoLib::Polyline const& polyline, double eps) const;

    MeshLib::Mesh const&  _mesh;
};

}

#endif /* EXTRACTPOLYGON_H_ */
