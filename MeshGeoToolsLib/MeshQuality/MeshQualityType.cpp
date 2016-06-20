/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshQualityType.h"

namespace MeshGeoToolsLib
{

const std::string MeshQualityType2String(const MeshQualityType t)
{
    if (t == MeshQualityType::ELEMENTSIZE)
        return "ElementSize";
    if (t == MeshQualityType::EDGERATIO)
        return "EdgeRatio";
    if (t == MeshQualityType::EQUIANGLESKEW)
        return "EquiAngleSkew";
    if (t == MeshQualityType::RADIUSEDGERATIO)
        return "RadiusEdgeRatio";
    if (t == MeshQualityType::SIZEDIFFERENCE)
        return "SizeDifference";
    return "none";
}

}  // end namespace MeshGeoToolsLib
