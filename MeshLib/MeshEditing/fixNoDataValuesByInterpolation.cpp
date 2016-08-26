/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "fixNoDataValuesByInterpolation.h"

#include <numeric>
#include <vector>

#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
void fixNoDataValuesByInterpolation(MeshLib::Mesh& mesh, double no_data_value,
                                    std::string const& original_property_name,
                                    std::string & interpolated_property_name)
{
    std::vector<MeshLib::Element*> const& elements(mesh.getElements());

    MeshLib::PropertyVector<double> const* orig_prop(
        mesh.getProperties().getPropertyVector<double>(original_property_name));
    if (!orig_prop)
    {
        OGS_FATAL("PropertyVector \"%s\" doesn't exist in the mesh.",
                  original_property_name.c_str());
    }

    MeshLib::PropertyVector<double>* int_prop(
        mesh.getProperties().getPropertyVector<double>(
            interpolated_property_name));
    if (int_prop)
    {
        OGS_FATAL(
            "PropertyVector \"%s\" for interpolation of no data values already "
            "exist in the mesh.",
            interpolated_property_name.c_str());
    }

    int_prop = mesh.getProperties().createNewPropertyVector<double>(
        interpolated_property_name, orig_prop->getMeshItemType(),
        orig_prop->getNumberOfComponents());
    int_prop->resize(orig_prop->size());

    for (std::size_t k(0); k < orig_prop->size(); ++k)
    {
        if ((*orig_prop)[k] == no_data_value)
        {
            (*int_prop)[k] = (*orig_prop)[k];
            continue;
        }

        std::size_t const n_neighbors(elements[k]->getNumberOfNeighbors());
        double prop_value(0.0);
        std::size_t cnt(0);  // count neighbors with valid property values
        // simple average
        for (std::size_t j(0); j < n_neighbors; ++j)
        {
            MeshLib::Element const* const j_th_neighbor(
                elements[k]->getNeighbor(j));
            if (j_th_neighbor != nullptr)
            {
                double neighbor_val((*orig_prop)[j_th_neighbor->getID()]);
                if (std::abs(neighbor_val - no_data_value) >
                    std::numeric_limits<double>::epsilon())
                {
                    prop_value += neighbor_val;
                    cnt++;
                }
            }
        }
        prop_value /= cnt;
        (*int_prop)[k] = prop_value;
    }
}
}  // end namespace MeshLib

