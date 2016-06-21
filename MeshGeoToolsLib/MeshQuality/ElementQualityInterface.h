/**
 * \file   ElementQualityInterface.h
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the ElementQualityInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTQUALITYINTERFACE_H_
#define ELEMENTQUALITYINTERFACE_H_

#include <vector>

#include "BaseLib/Histogram.h"

#include "MeshLib/Mesh.h"
#include "MeshGeoToolsLib/MeshQuality/ElementQualityMetric.h"
#include "MeshGeoToolsLib/MeshQuality/EdgeRatioMetric.h"
#include "MeshGeoToolsLib/MeshQuality/ElementSizeMetric.h"
#include "MeshGeoToolsLib/MeshQuality/SizeDifferenceMetric.h"
#include "MeshGeoToolsLib/MeshQuality/AngleSkewMetric.h"
#include "MeshGeoToolsLib/MeshQuality/RadiusEdgeRatioMetric.h"

namespace MeshGeoToolsLib
{

/**
 * Interface class for handling mesh element quality metrics
 */
class ElementQualityInterface
{
public:
    /// Constructor
    ElementQualityInterface(MeshLib::Mesh const& mesh, MeshQualityType t)
    : _type(t), _mesh(mesh), _quality_tester(nullptr)
    {
        calculateElementQuality(_mesh, _type);
    }

    /// Destructor
    ~ElementQualityInterface()
    {
        delete _quality_tester;
    }

    /// Returns the vector containing a quality measure for each element.
    std::vector<double> const getQualityVector() const
    {
        if (_quality_tester)
            return _quality_tester->getElementQuality();

        std::vector<double> empty_quality_vec(0);
        return empty_quality_vec;
    }

    /// Returns a histogram of the quality vector seperated into the given number of bins.
    /// If no number of bins is specified, one will be calculated based on the Sturges criterium.
    BaseLib::Histogram<double> getHistogram(std::size_t n_bins = 0) const
    {
        if (_quality_tester)
            return _quality_tester->getHistogram(static_cast<std::size_t>(n_bins));

        std::vector<double> empty_quality_vec(0);
        return empty_quality_vec;
    }

    /// Writes a histogram of the quality vector to a specified file.
    int writeHistogram(std::string const& file_name, std::size_t n_bins = 0) const
    {
        if (_quality_tester == nullptr)
            return 1;

        BaseLib::Histogram<double> const histogram (_quality_tester->getHistogram(n_bins));
        histogram.write(file_name, _mesh.getName(), MeshQualityType2String(_type));
        return 0;
    }

private:
    /// Calculates the quality of each mesh element based on the specified metric
    void calculateElementQuality(MeshLib::Mesh const& mesh, MeshQualityType t)
    {
        if (t == MeshQualityType::EDGERATIO)
            _quality_tester = new MeshGeoToolsLib::EdgeRatioMetric(mesh);
        else if (t == MeshQualityType::ELEMENTSIZE)
            _quality_tester = new MeshGeoToolsLib::ElementSizeMetric(mesh);
        else if (t == MeshQualityType::SIZEDIFFERENCE)
            _quality_tester = new MeshGeoToolsLib::SizeDifferenceMetric(mesh);
        else if (t == MeshQualityType::EQUIANGLESKEW)
            _quality_tester = new MeshGeoToolsLib::AngleSkewMetric(mesh);
        else if (t == MeshQualityType::RADIUSEDGERATIO)
            _quality_tester = new MeshGeoToolsLib::RadiusEdgeRatioMetric(mesh);
        else
        {
            ERR("ElementQualityInterface::calculateElementQuality(): Unknown MeshQualityType.");
            return;
        }
        _quality_tester->calculateQuality();
    }

    MeshQualityType const _type;
    MeshLib::Mesh const& _mesh;
    MeshGeoToolsLib::ElementQualityMetric* _quality_tester;
};

}

#endif /* ELEMENTQUALITYINTERFACE_H_ */
