/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTransportPorosityFromMassBalance.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"
#include "TransportPorosityFromMassBalance.h"

namespace MaterialPropertyLib
{
std::unique_ptr<TransportPorosityFromMassBalance>
createTransportPorosityFromMassBalance(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "TransportPorosityFromMassBalance");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create TransportPorosityFromMassBalance medium property {:s}.",
         property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__TransportPorosityFromMassBalance__initial_porosity}
        config.getConfigParameter<std::string>("initial_porosity");
    auto const& initial_porosity = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);

    //! \ogs_file_param{properties__property__TransportPorosityFromMassBalance__minimal_porosity}
    auto const& phi_min = config.getConfigParameter<double>("minimal_porosity");

    //! \ogs_file_param{properties__property__TransportPorosityFromMassBalance__maximal_porosity}
    auto const& phi_max = config.getConfigParameter<double>("maximal_porosity");

    return std::make_unique<TransportPorosityFromMassBalance>(
        std::move(property_name), initial_porosity, phi_min, phi_max);
}
}  // namespace MaterialPropertyLib
