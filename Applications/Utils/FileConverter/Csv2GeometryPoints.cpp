/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "FileIO/CsvInterface.h"

#include "GeoLib/IO/writeGeometryToFile.h"
#include "GeoLib/GEOObjects.h"

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Program reads a CSV file and transforms the input in a geometry", ' ',
        "0.1");
    TCLAP::ValueArg<std::string> csv_fname_arg("i", "csv-input",
                                               "non-empty file name", true, "",
                                               "file name string");
    cmd.add(csv_fname_arg);
    TCLAP::ValueArg<unsigned> x_col_arg(
        "x", "x-col", "number of the column that contains x coordinate", false,
        0u, "integer number used as column to read the x coordinate");
    cmd.add(x_col_arg);
    TCLAP::ValueArg<unsigned> y_col_arg(
        "y", "y-col", "number of the column that contains y coordinate", false,
        1u, "integer number used as column to read the y coordinate");
    cmd.add(y_col_arg);
    TCLAP::ValueArg<unsigned> z_col_arg(
        "z", "z-col", "number of the column that contains z coordinate", false,
        2u, "integer number used as column to read the z coordinate");
    cmd.add(z_col_arg);
    TCLAP::ValueArg<std::string> geo_fname_arg("o", "output",
                                               "non-empty file name", true, "",
                                               "file name string");
    cmd.add(geo_fname_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<std::vector<GeoLib::Point*>> pnts(
        new std::vector<GeoLib::Point*>);
    FileIO::CsvInterface::readPoints(csv_fname_arg.getValue(), ' ', *pnts,
                                     x_col_arg.getValue(), y_col_arg.getValue(),
                                     z_col_arg.getValue());
    GeoLib::GEOObjects geo_objs;
    std::string geo_name(
        "CsvPoints-" + BaseLib::dropFileExtension(
                           BaseLib::extractBaseName(csv_fname_arg.getValue())));
    geo_objs.addPointVec(std::move(pnts), geo_name);
    GeoLib::IO::writeGeometryToFile(geo_name, geo_objs, geo_fname_arg.getValue());
}
