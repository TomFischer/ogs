/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef TESTS_MESHLIB_AUTOCHECKTOOLS_H_
#define TESTS_MESHLIB_AUTOCHECKTOOLS_H_

#include "autocheck/autocheck.hpp"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "Tests/MathLib/AutoCheckTools.h"

namespace autocheck
{

template <typename T, typename Gen = ac::generator<T>>
struct RandomPairGenerator
{
    Gen source;
    using result_type = std::tuple<T, T>;

    result_type operator()(std::size_t)
    {
        auto p = std::make_tuple(source(), source());
        return p;
    }
};

template <typename T, std::size_t N, typename Gen = generator<T>>
struct RandomSurfaceMeshGenerator
{
    Gen source;

    using result_type = MeshLib::Mesh;

    result_type operator()(std::size_t size = 0)
    {
        result_type rv;
        return rv;
    }
};

}  // namespace autocheck
#endif  // TESTS_MESHLIB_AUTOCHECKTOOLS_H_
