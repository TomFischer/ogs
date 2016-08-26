
/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <string>

namespace MeshLib
{

class Mesh;

void fixNoDataValuesByInterpolation(MeshLib::Mesh& mesh, double no_data_value,
                                    std::string const& original_property_name,
                                    std::string & interpolated_property_name);
}  // end namespace MeshLib
