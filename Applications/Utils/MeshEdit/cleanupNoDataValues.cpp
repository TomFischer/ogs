/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_arg(
        "m", "mesh", "input mesh file", true, "", "string");
    cmd.add( mesh_arg );
    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh", "output mesh file", false, "", "string");
    cmd.add(mesh_out_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    if (!mesh)
    {
        ERR("Could not read mesh from file \"%s\".",
            mesh_arg.getValue().c_str());
        return EXIT_FAILURE;
    }

    MeshLib::IO::writeMeshToFile(*mesh, out_fname);

    return EXIT_SUCCESS;
}
