/**
 * \file   AngleSkewMetric.h
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Definition of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ANGLESKEWMETRIC_H_
#define ANGLESKEWMETRIC_H_

#include "ElementQualityMetric.h"

namespace MeshGeoToolsLib
{

/**
 * Calculates the quality of mesh elements based on the EquiAngleSkew measure
 */
class AngleSkewMetric : public ElementQualityMetric
{
public:
    AngleSkewMetric(MeshLib::Mesh const& mesh);
    virtual ~AngleSkewMetric();

    virtual void calculateQuality ();

private:
    double checkTriangle(MeshLib::Element const& elem) const;
    double checkQuad(MeshLib::Element const& elem) const;
    double checkTetrahedron(MeshLib::Element const& elem) const;
    double checkHexahedron(MeshLib::Element const& elem) const;
    double checkPrism (MeshLib::Element const& elem) const;
    void getMinMaxAngleFromQuad(double const* const n0,
                                double const* const n1, double const* const n2,
                                double const* const n3, double &min_angle,
                                double &max_angle) const;
    void getMinMaxAngleFromTriangle(double const* const n0,
                                    double const* const n1, double const* const n2,
                                    double &min_angle, double &max_angle) const;
};
}

#endif /* ANGLESKEWMETRIC_H_ */
