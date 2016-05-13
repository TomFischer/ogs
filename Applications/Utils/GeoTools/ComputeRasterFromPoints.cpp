/**
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>
#include <vector>
#include <string>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/IO/readGeometryFromFile.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Raster.h"

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Program computes a raster out of a given point set", ' ', "0.1");
    TCLAP::ValueArg<std::string> geo_fname_arg("g", "geometry-file",
                                               "non-empty file name", true, "",
                                               "file name string");
    cmd.add(geo_fname_arg);
    TCLAP::ValueArg<unsigned> n_rows_arg(
        "r", "nRows", "number of rows", true, 10u,
        "integer number used as number of rows");
    cmd.add(n_rows_arg);
    TCLAP::ValueArg<unsigned> n_cols_arg(
        "c", "nCols", "number of columns", true, 10u,
        "integer number used as number of columns");
    cmd.add(n_cols_arg);
    TCLAP::ValueArg<double> x_min_arg(
        "", "xmin", "minimal value of x ", false, 0.0,
        "floating point number for the minimal value for x");
    cmd.add(x_min_arg);
    TCLAP::ValueArg<double> y_min_arg(
        "", "ymin", "minimal value of y", false, 0.0,
        "floating point number for the minimal value for y");
    cmd.add(y_min_arg);
    TCLAP::ValueArg<double> edge_size_of_cell_arg(
        "e", "edge-size", "size of the edge of a cell", false, 1.0,
        "floating point number for the size of an edge of a cell");
    cmd.add(edge_size_of_cell_arg);
    TCLAP::ValueArg<std::size_t> dilate_arg(
        "d", "dilate-pixels",
        "number of neighbor pixels that have no-data values and that have a "
        "neighbor neighbor pixel with a valid value",
        false, 0, "positiv integer number (of neighbor pixels to dilate)");
    cmd.add(dilate_arg);
    cmd.parse(argc, argv);

    GeoLib::GEOObjects geo_objs;
    GeoLib::IO::readGeometryFromFile(geo_fname_arg.getValue(), geo_objs);

    std::vector<std::string> geo_names;
    geo_objs.getGeometryNames(geo_names);
    std::vector<GeoLib::Point*> const* pnts(
        geo_objs.getPointVec(geo_names.back()));

    // init raster data header
    GeoLib::RasterHeader rh{
        n_cols_arg.getValue(), n_rows_arg.getValue(),
        MathLib::Point3d{{{x_min_arg.getValue(), y_min_arg.getValue(), 0.0}}},
        edge_size_of_cell_arg.getValue(), -9999.0};

    std::vector<double> data(rh.n_cols*rh.n_rows, rh.no_data);

    for (std::size_t k(0); k<pnts->size(); k++) {
        GeoLib::Point const& pnt(*(*pnts)[k]);
        // compute pos in grid
        const double x((pnt[0] - rh.origin[0])/rh.cell_size - 0.5);
        const double y((pnt[1] - rh.origin[1])/rh.cell_size - 0.5);
        if (x<0 || y<0 || x>rh.n_cols || y>rh.n_rows)
            continue;
        const std::size_t xpix(static_cast<std::size_t>(x));
        const std::size_t ypix(static_cast<std::size_t>(y));
        data[xpix + ypix*rh.n_cols] = pnt[2];
    }

    const double eps (std::numeric_limits<double>::epsilon());
    for (std::size_t k(0); k<dilate_arg.getValue(); k++) {
        // dilate one layer
        for (std::size_t i(1); i<rh.n_cols; i++) {
            for (std::size_t j(1); j<rh.n_rows; j++) {
                if (data[j*rh.n_cols+i] != rh.no_data) {
                    if (std::abs(data[j*rh.n_cols+i-1] - rh.no_data) < eps)
                        data[j*rh.n_cols+i-1] = data[j*rh.n_cols+i];
                    if (std::abs(data[(j-1)*rh.n_cols+i] - rh.no_data) < eps)
                        data[(j-1)*rh.n_cols+i] = data[j*rh.n_cols+i];
                }
            }
        }
        for (std::size_t i(0); i<rh.n_cols-1; i++) {
            for (std::size_t j(0); j<rh.n_rows-1; j++) {
                if (data[(rh.n_rows - j - 1) * rh.n_cols +
                         (rh.n_cols - i - 1)] != rh.no_data)
                {
                    if (std::abs(data[(rh.n_rows - j) * rh.n_cols +
                                      (rh.n_cols - i - 1)] -
                                 rh.no_data) < eps)
                        data[(rh.n_rows - j) * rh.n_cols +
                             (rh.n_cols - i - 1)] =
                            data[(rh.n_rows - j - 1) * rh.n_cols + rh.n_cols -
                                 i - 1];
                }
            }
        }
    }

    GeoLib::Raster raster(rh, data.begin(), data.end());

    GeoLib::IO::AsciiRasterInterface::writeRasterAsASC(
        raster, geo_fname_arg.getValue() + ".asc");

    return EXIT_SUCCESS;
}
