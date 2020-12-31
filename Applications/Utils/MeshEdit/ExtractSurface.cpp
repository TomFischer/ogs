/**
 * \file
 * \brief Extracts the surface from the given mesh.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSurfaceExtraction.h"

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Extracts a 2D surface from a 3D input mesh by specifying a normal "
        "vector and an angle. (The surface normal (0, 0, 0) will extract the "
        "complete outer boundary of the 3D mesh.)\n"
        "An extensive documentation can be found at "
        "https://docs.opengeosys.org/docs/tools/meshing-submeshes/extract-surface\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg use_ascii_arg("", "ascii-output",
                                   "If the switch is set use ascii instead of "
                                   "binary format for data in the vtu output.",
                                   false);
    cmd.add(use_ascii_arg);
    TCLAP::ValueArg<std::string> face_prop_name(
        "f", "face-property-name",
        "the name of the data array the surface face id of the subsurface/bulk "
        "element will be stored to",
        false, "bulk_face_ids", "string");
    cmd.add(face_prop_name);
    TCLAP::ValueArg<std::string> element_prop_name(
        "e", "element-property-name",
        "the name of the data array the subsurface/bulk element id will be "
        "stored to",
        false, "bulk_element_ids", "string");
    cmd.add(element_prop_name);
    TCLAP::ValueArg<std::string> node_prop_name(
        "n", "node-property-name",
        "the name of the data array the subsurface/bulk node id will be stored "
        "to",
        false, "bulk_node_ids", "string");
    cmd.add(node_prop_name);
    TCLAP::ValueArg<double> angle_arg(
        "a", "angle",
        "tolerated angle (in degrees) between given normal and element normal",
        false, 90, "floating point value");
    cmd.add(angle_arg);
    TCLAP::ValueArg<double> z("z", "z-component", "z component of the normal",
                              false, -1.0, "floating point value");
    cmd.add(z);
    TCLAP::ValueArg<double> y("y", "y-component", "y component of the normal",
                              false, 0, "floating point value");
    cmd.add(y);
    TCLAP::ValueArg<double> x("x", "x-component", "x component of the normal",
                              false, 0, "floating point value");
    cmd.add(x);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the surface mesh should be written to", false, "",
        "file name of output mesh");
    cmd.add(mesh_out);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));

    if (!mesh)
    {
        ERR("Error reading mesh file.");
        return EXIT_FAILURE;
    }

    if (mesh->getDimension() != 3)
    {
        ERR("Surfaces can currently only be extracted from 3D meshes.");
        return EXIT_FAILURE;
    }

    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    // extract surface
    Eigen::Vector3d const dir({x.getValue(), y.getValue(), z.getValue()});
    double const angle(angle_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::MeshSurfaceExtraction::getMeshSurface(
            *mesh, dir, angle, node_prop_name.getValue(),
            element_prop_name.getValue(), face_prop_name.getValue()));

    std::string out_fname(mesh_out.getValue());
    if (out_fname.empty())
    {
        out_fname = BaseLib::dropFileExtension(mesh_in.getValue()) + "_sfc.vtu";
    }

    auto const data_mode =
        use_ascii_arg.getValue() ? vtkXMLWriter::Ascii : vtkXMLWriter::Binary;
    MeshLib::IO::writeVtu(*surface_mesh, out_fname, data_mode);

    return EXIT_SUCCESS;
}
