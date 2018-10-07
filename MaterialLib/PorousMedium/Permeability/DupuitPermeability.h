/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>

#include "BaseLib/Error.h"
#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace PorousMedium
{
/// The class implements a basic permeability model that employs a parameter
/// (for instance a constant parameter or mesh cell dependend parameter) to fill
/// the intrinsic permeability tensor.
class DupuitPermeability final : public Permeability
{
public:
    explicit DupuitPermeability(
        ProcessLib::Parameter<double> const& permeability_parameter,
        int const dimension)
        : Permeability(permeability_parameter, dimension)
    {
    }

    ~DupuitPermeability() = default;

    /**
     *  Get property value.
     *  @param t point in time
     *  @param pos spatial position
     *  @param variable    A variable with any double type value.
     *  @param temperature Temperature with any double type value.
     */
    Eigen::MatrixXd const& getValue(
        const double t,
        ProcessLib::SpatialPosition const& pos,
        const double variable,
        const double temperature) const override
    {
        _intrinsic_permeability_tensor =
            variable * Permeability::getValue(t, pos, variable, temperature);

        return _intrinsic_permeability_tensor;
    }

private:
//    ProcessLib::Parameter<double> const& _thickness;
};

}  // end of namespace
}  // end of namespace
