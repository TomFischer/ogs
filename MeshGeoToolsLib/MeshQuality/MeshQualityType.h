/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHQUALITYTYPE_H
#define MESHQUALITYTYPE_H

#include <string>

namespace MeshGeoToolsLib
{

/**
 * \brief Describes a mesh quality metric.
 */
enum class MeshQualityType
{
    INVALID = 0,
    ELEMENTSIZE,
    SIZEDIFFERENCE,
    EDGERATIO,
    EQUIANGLESKEW,
    RADIUSEDGERATIO
};

const std::string MeshQualityType2String(const MeshQualityType t);

}  // end namespace MeshGeoToolsLib

#endif
