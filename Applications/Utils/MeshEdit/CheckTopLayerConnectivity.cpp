/*
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"
#include "BaseLib/Histogram.h"

#include "MeshLib/IO/readMeshFromFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"

#include "MathLib/Vector3.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Checks the connectivity of the top layer material id.",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> mesh_arg("i", "input-mesh-file",
        "the name of the file containing the mesh", true,
        "", "file name");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
    auto subsfc_mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!subsfc_mesh) {
        ERR("Error reading mesh \"%s\".", mesh_arg.getValue().c_str());
        return EXIT_FAILURE;
    }
    INFO("done.");

    // creating statistic
    boost::optional<MeshLib::PropertyVector<int> const&> materials(
        subsfc_mesh->getProperties().getPropertyVector<int>("MaterialIDs")
    );
    if (!materials) {
        ERR("No MaterialIDs property in the mesh.");
        return EXIT_FAILURE;
    }

    int const max_layer_id(
        *std::max_element(materials->cbegin(), materials->cend()));
    MathLib::Vector3 const up{0, 0, 1};
    std::vector<std::size_t> neighbor_count;
    for (std::size_t k(0); k < materials->size(); ++k)
    {
        if ((*materials)[k] != max_layer_id)
            continue;

        auto const* elem(subsfc_mesh->getElement(k));
        neighbor_count.push_back(0);
        for (unsigned j(0); j < elem->getNumberOfFaces(); ++j)
        {
            std::unique_ptr<MeshLib::Element const> face{elem->getFace(j)};
            MathLib::Vector3 const normal{
                MeshLib::FaceRule::getSurfaceNormal(face.get())};
            if (MathLib::scalarProduct(normal, up) < 0.9 * normal.getLength())
            {
                auto const n = elem->getNeighbor(j);
                if (n != nullptr && (*materials)[n->getID()] == max_layer_id)
                    neighbor_count.back()++;
            }
        }
    }
    std::cout << "neighbor statistic: \n";
    std::size_t const max_value(
        *std::max_element(neighbor_count.cbegin(), neighbor_count.cend()));
        std::cout << "max_value: " << max_value << "\n";
    std::size_t const min_value(
        *std::min_element(neighbor_count.cbegin(), neighbor_count.cend()));
    BaseLib::Histogram<std::size_t> h(neighbor_count, max_value-min_value+1);
    h.prettyPrint(std::cout);
    INFO("size of neighbor count: %d.", neighbor_count.size());
    return EXIT_SUCCESS;
}
