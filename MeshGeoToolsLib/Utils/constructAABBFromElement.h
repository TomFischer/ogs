/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef CONSTRUCTAABBFROMELEMENT_H_
#define CONSTRUCTAABBFROMELEMENT_H_

#include "MeshLib/Elements/Element.h"
#include "GeoLib/AABB.h"

namespace MeshGeoToolsLib
{

GeoLib::AABB constructAABBFromElement(MeshLib::Element const& e)
{
    return GeoLib::AABB(e.getNodes(), e.getNodes() + e.getNumberOfBaseNodes());
}

}  // end namespace MeshGeoToolsLib
#endif
