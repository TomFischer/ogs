/**
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <memory>
#include <cmath>
#include <numeric>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

double computeInitialTemperatureValueDomain(MeshLib::Node const& node)
{
    const double z = node[2];
    // initial condition for the total domain
    return 20. - (150.0 / 5500.0) * z;
}

double computeInitialTemperatureValueFault(MeshLib::Node const& node)
{
    const double y = node[1];
    const double z = node[2];
    // initial condition for the fault domain
    return 20. - (150.0 / 5500.0) * z +
           sin(3.1415 * z / 5500.) * cos(3.1415 * y / 5500.);
}

double computeInitialPressureValue(MeshLib::Node const& node)
{
    const double z = node[2];
    return -9810 * z;
}

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logo_setup;

    TCLAP::CmdLine cmd(
        "Add a property to a mesh",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> out_mesh_arg(
        "o",
        "out-mesh",
        "the mesh is stored to a file of this name",
        true,
        "",
        "filename for mesh output");
    cmd.add(out_mesh_arg);

    TCLAP::ValueArg<std::string> mesh_arg("i",
                                          "in-mesh",
                                          "the mesh is read from this file",
                                          true,
                                          "",
                                          "file name");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    MeshLib::PropertyVector<double> *pv_t;
    if (!mesh->getProperties().existsPropertyVector<double>(
            "initial_temperature"))
    {
        // create new PropertyVector for initial temperature
        pv_t = mesh->getProperties().createNewPropertyVector<double>(
            "initial_temperature", MeshLib::MeshItemType::Node, 1);
        pv_t->resize(mesh->getNumberOfNodes());
    }
    else
    {
        pv_t = mesh->getProperties().getPropertyVector<double>(
            "initial_temperature");
    }

    MeshLib::PropertyVector<double> *pv_p;
    if (!mesh->getProperties().existsPropertyVector<double>(
            "initial_pressure"))
    {
        // create new PropertyVector for initial temperature
        pv_p = mesh->getProperties().createNewPropertyVector<double>(
            "initial_pressure", MeshLib::MeshItemType::Node, 1);
        pv_p->resize(mesh->getNumberOfNodes());
    }
    else
    {
        pv_p = mesh->getProperties().getPropertyVector<double>(
            "initial_pressure");
    }

    auto* material_ids =
        mesh->getProperties().getPropertyVector<int>("MaterialIDs");

    for (std::size_t k(0); k < mesh->getNumberOfElements(); ++k)
    {
        auto const* element(mesh->getElement(k));
        std::size_t n_nodes(element->getNumberOfNodes());
        for (std::size_t i(0); i < n_nodes; ++i)
        {
            (*pv_t)[element->getNode(i)->getID()] =
                computeInitialTemperatureValueDomain(*element->getNode(i));
            (*pv_p)[element->getNode(i)->getID()] =
                computeInitialPressureValue(*element->getNode(i));
        }
    }
    for (std::size_t k(0); k < mesh->getNumberOfElements(); ++k)
    {
        if ((*material_ids)[k] == 1)
        {
            auto const* element(mesh->getElement(k));
            std::size_t n_nodes(element->getNumberOfNodes());
            for (std::size_t i(0); i < n_nodes; ++i)
            {
                (*pv_t)[element->getNode(i)->getID()] =
                    computeInitialTemperatureValueFault(*element->getNode(i));
            }
        }
    }


    MeshLib::IO::writeMeshToFile(*mesh, out_mesh_arg.getValue());

    return EXIT_SUCCESS;
}
