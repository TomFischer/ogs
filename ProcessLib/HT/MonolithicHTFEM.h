/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on October 11, 2017, 2:33 PM
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <boost/math/special_functions/pow.hpp>

#include "HTProcessData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HTFEM.h"

namespace ProcessLib
{
namespace HT
{
const unsigned NUM_NODAL_DOF = 2;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class MonolithicHTFEM
    : public HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    MonolithicHTFEM(MeshLib::Element const& element,
                    std::size_t const local_matrix_size,
                    bool is_axially_symmetric,
                    unsigned const integration_order,
                    HTProcessData const& process_data)
        : HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
              element, local_matrix_size, is_axially_symmetric,
              integration_order, process_data, NUM_NODAL_DOF)
    {
    }

    void assemble(double const t, double const /*dt*/,
                  std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        auto KTT = local_K.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);
        auto MTT = local_M.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);
        auto Kpp = local_K.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Mpp = local_M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto MpT = local_M.template block<pressure_size, temperature_size>(
            pressure_index, temperature_index);
        auto Bp = local_b.template segment<pressure_size>(pressure_index);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        auto p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);

        auto const& process_data = this->_process_data;
        auto const& medium =
            *process_data.media_map->getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& b = process_data.specific_body_force;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            auto const porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t);

            auto const intrinsic_permeability =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    medium
                        .property(
                            MaterialPropertyLib::PropertyType::permeability)
                        .value(vars, pos, t));

            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(vars, pos, t);

            // Use the fluid density model to compute the density
            auto const fluid_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t);
            auto const drho_dp =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const drho_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);

            // Use the viscosity model to compute the viscosity
            auto const viscosity =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t);
            GlobalDimMatrixType K_over_mu = intrinsic_permeability / viscosity;

            GlobalDimVectorType const velocity =
                process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu * (dNdx * p_nodal_values -
                                                        fluid_density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);

            // matrix assembly
            GlobalDimMatrixType const thermal_conductivity_dispersivity =
                this->getThermalConductivityDispersivity(
                    vars, porosity, fluid_density, specific_heat_capacity_fluid,
                    velocity, I, pos, t);
            KTT.noalias() +=
                (dNdx.transpose() * thermal_conductivity_dispersivity * dNdx +
                 N.transpose() * velocity.transpose() * dNdx * fluid_density *
                     specific_heat_capacity_fluid) *
                w;
            Kpp.noalias() +=
                (w * fluid_density) * dNdx.transpose() * K_over_mu * dNdx;
            MTT.noalias() += w *
                             this->getHeatEnergyCoefficient(
                                 vars, porosity, fluid_density,
                                 specific_heat_capacity_fluid, pos, t) *
                             N.transpose() * N;
            Mpp.noalias() += (w * porosity * drho_dp) * N.transpose() * N;
            MpT.noalias() += (w * drho_dT) * N.transpose() * porosity * N;
            if (process_data.has_gravity)
            {
                Bp += w * boost::math::pow<2>(fluid_density) *
                      dNdx.transpose() * K_over_mu * b;
            }
        }
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        assert(local_x.size() == pressure_size + temperature_size);

        auto p = Eigen::Map<typename ShapeMatricesType::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);
        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                     temperature_size);
        auto p_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            pressure_size> const>(local_xdot.data() + pressure_index,
                                  pressure_size);
        auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);
        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesType::template MatrixType<
                pressure_size + temperature_size,
                pressure_size + temperature_size>>(
            local_Jac_data, pressure_size + temperature_size,
            pressure_size + temperature_size);
        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesType::template VectorType<pressure_size +
                                                            temperature_size>>(
            local_rhs_data, pressure_size + temperature_size);

        auto const& process_data = this->_process_data;
        auto const& b = process_data.specific_body_force;
        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        auto const& medium =
            *process_data.media_map->getMedium(this->_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& solid_phase = medium.phase("Solid");

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));
        GlobalDimMatrixType const& Zero(
            GlobalDimMatrixType::Zero(GlobalDim, GlobalDim));

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);

            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // \todo the first argument has to be changed for non constant
            // porosity model
            auto const porosity =
                solid_phase.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t);

            auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t));
            auto const dK_dp = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::permeability)
                    .dValue(vars, MaterialPropertyLib::Variable::phase_pressure,
                            pos, t));
            auto const mu =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t);
            auto const dmu_dp =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const dmu_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);

            GlobalDimMatrixType const K_over_mu = K / mu;
            GlobalDimMatrixType const K_over_mu2 = K / boost::math::pow<2>(mu);

            auto const dspecific_heat_capacity_fluid_dp =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const dspecific_heat_capacity_fluid_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(vars, pos, t);

            auto const specific_heat_capacity_solid =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(vars, pos, t);
            auto const dspecific_heat_capacity_solid_dp =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const dspecific_heat_capacity_solid_dT =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);

            // fluid density and derivations
            auto const rho =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t);
            auto const drho_dp = liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure);
            auto const drho_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);

            auto const d2rho_dp2 =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template d2Value<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        MaterialPropertyLib::Variable::phase_pressure, pos, t);
            auto const d2rho_dp_dT =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template d2Value<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        MaterialPropertyLib::Variable::temperature, pos, t);
            auto const d2rho_dT2 =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template d2Value<double>(
                        vars, MaterialPropertyLib::Variable::temperature,
                        MaterialPropertyLib::Variable::temperature, pos, t);

            auto const rho_solid =
                solid_phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t);
            auto const drho_solid_dp =
                solid_phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const drho_solid_dT =
                solid_phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            auto const C_p =
                porosity * rho * specific_heat_capacity_fluid +
                rho_solid * (1 - porosity) * specific_heat_capacity_solid;
            auto const dC_p_dp =
                porosity * (drho_dp * specific_heat_capacity_fluid +
                            rho * dspecific_heat_capacity_fluid_dp) +
                (1 - porosity) * (drho_solid_dp * specific_heat_capacity_solid +
                                  rho_solid * dspecific_heat_capacity_solid_dp);
            auto const dC_p_dT =
                porosity * (drho_dT * specific_heat_capacity_fluid +
                            rho * dspecific_heat_capacity_fluid_dT) +
                (1 - porosity) * (drho_solid_dT * specific_heat_capacity_solid +
                                  rho_solid * dspecific_heat_capacity_solid_dT);
            // INFO("C_p: %e, dC_p_dp: %e, dC_p_dT: %e", C_p, dC_p_dp, dC_p_dT);

            // assemble first part of Jac_pp
            auto C_pp_jac =
                local_Jac.template block<pressure_size, pressure_size>(
                    pressure_index, pressure_index);
            C_pp_jac.noalias() +=
                N.transpose() * porosity * d2rho_dp2 * N * p_dot * N * w;
            C_pp_jac.noalias() +=
                1 / dt * N.transpose() * porosity * drho_dp * N * w;
            C_pp_jac.noalias() +=
                N.transpose() * porosity * d2rho_dp_dT * N * T_dot * N * w;

            auto C_pT_jac =
                local_Jac.template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index);
            // equation (4.33)
            C_pT_jac.noalias() +=
                N.transpose() * porosity * d2rho_dp_dT * N * p_dot * N * w;
            C_pT_jac.noalias() +=
                N.transpose() * porosity * d2rho_dT2 * N * T_dot * N * w;
            C_pT_jac.noalias() +=
                1 / dt * N.transpose() * porosity * drho_dT * N * w;

            auto C_TT_jac =
                local_Jac.template block<temperature_size, temperature_size>(
                    temperature_index, temperature_index);
            // equation (4.35)
            C_TT_jac.noalias() += N.transpose() * dC_p_dT * N * T_dot * N * w;
            C_TT_jac.noalias() += 1 / dt * N.transpose() * C_p * N * w;

            GlobalDimVectorType const darcy_velocity =
                process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu * (dNdx * p - rho * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p);

            // pressure equation - pressure derivations
            auto K_pp_jac =
                local_Jac.template block<pressure_size, pressure_size>(
                    pressure_index, pressure_index);
            // equation (4.31)
            K_pp_jac.noalias() +=
                dNdx.transpose() * drho_dp * K_over_mu * dNdx * p * N * w;
            K_pp_jac.noalias() += dNdx.transpose() * rho /
                                  boost::math::pow<2>(mu) *
                                  (dK_dp * mu - K * dmu_dp) * dNdx * p * N * w;
            K_pp_jac.noalias() += dNdx.transpose() * rho * K_over_mu * dNdx * w;
            if (process_data.has_gravity)
            {
                K_pp_jac.noalias() += dNdx.transpose() * 2 * rho * drho_dp *
                                      K_over_mu * b * N * w;
                K_pp_jac.noalias() -= dNdx.transpose() *
                                      boost::math::pow<2>(rho) * K_over_mu2 *
                                      dmu_dp * b * N * w;
            }

            auto const dlambda_fluid_dp =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure);
            auto const dlambda_solid_dp =
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const dlambda_cond_dp =
                porosity * dlambda_fluid_dp + (1 - porosity) * dlambda_solid_dp;

            auto const alpha_L =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_longitudinal_dispersivity)
                    .template value<double>(vars, pos, t);
            auto const alpha_T =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_transversal_dispersivity)
                    .template value<double>(vars, pos, t);
            auto const dalpha_L_dp =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_longitudinal_dispersivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            auto const dalpha_T_dp =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_transversal_dispersivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t);
            double const darcy_norm = darcy_velocity.norm();
            auto const dlambda_disp_dp =
                darcy_norm == 0
                    ? Zero
                    : (drho_dp * specific_heat_capacity_fluid +
                       rho * dspecific_heat_capacity_fluid_dp) *
                              (alpha_T * darcy_velocity.norm() * I +
                               ((alpha_L - alpha_T) / darcy_velocity.norm()) *
                                   darcy_velocity *
                                   darcy_velocity.transpose()) +
                          rho * specific_heat_capacity_fluid *
                              (dalpha_T_dp * darcy_velocity.norm() * I +
                               ((dalpha_L_dp - dalpha_T_dp) /
                                darcy_velocity.norm()) *
                                   darcy_velocity * darcy_velocity.transpose());
            auto const dlambda_dp = dlambda_cond_dp * I + dlambda_disp_dp;

            auto K_pT_jac =
                local_Jac.template block<pressure_size, temperature_size>(
                    pressure_index, temperature_index);
            // equation(4.37)
            if (process_data.has_gravity)
            {
                K_pT_jac.noalias() += dNdx.transpose() * drho_dT * K_over_mu *
                                      (dNdx * p + 2 * rho * b) * N * w;
                K_pT_jac.noalias() -= dNdx.transpose() *
                                      boost::math::pow<2>(rho) / mu *
                                      K_over_mu * dmu_dT * b * N * w;
                K_pT_jac.noalias() -= dNdx.transpose() * rho / mu * K_over_mu *
                                      dmu_dT * dNdx * p * N * w;
            }
            else
            {
                K_pT_jac.noalias() += dNdx.transpose() * drho_dT * K_over_mu *
                                      dNdx * p * N * w;
                K_pT_jac.noalias() -= dNdx.transpose() * rho / mu * K_over_mu *
                                      dmu_dT * dNdx * p * N * w;
            }

            auto const dlambda_fluid_dT =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            auto const dlambda_solid_dT =
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            auto const dlambda_cond_dT =
                porosity * dlambda_fluid_dT + (1 - porosity) * dlambda_solid_dT;

            auto const dalpha_L_dT =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_longitudinal_dispersivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            auto const dalpha_T_dT =
                medium
                    .property(MaterialPropertyLib::PropertyType::
                                  thermal_transversal_dispersivity)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::temperature, pos,
                        t);
            GlobalDimMatrixType const dlambda_disp_dT =
                darcy_norm == 0.0
                    ? Zero
                    : (drho_dT * specific_heat_capacity_fluid +
                       rho * dspecific_heat_capacity_fluid_dT) *
                              (alpha_T * darcy_velocity.norm() * I +
                               ((alpha_L - alpha_T) / darcy_velocity.norm()) *
                                   darcy_velocity *
                                   darcy_velocity.transpose()) +
                          rho * specific_heat_capacity_fluid *
                              (dalpha_T_dT * darcy_velocity.norm() * I +
                               ((dalpha_L_dT - dalpha_T_dT) /
                                darcy_velocity.norm()) *
                                   darcy_velocity * darcy_velocity.transpose());
            GlobalDimMatrixType const dlambda_dT =
                dlambda_cond_dT * I + dlambda_disp_dT;

            auto const lambda_fluid =
                liquid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(vars, pos, t);
            auto const lambda_solid =
                solid_phase
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .template value<double>(vars, pos, t);
            auto const lambda_cond =
                porosity * lambda_fluid + (1 - porosity) * lambda_solid;
            auto const lambda_dist =
                darcy_norm == 0
                    ? Zero
                    : rho * specific_heat_capacity_fluid *
                          (alpha_T * darcy_norm * I +
                           (alpha_L - alpha_T) / darcy_norm * darcy_velocity *
                               darcy_velocity.transpose());
            auto const lambda = lambda_cond * I + lambda_dist;

            // temperature equation - temperature derivations
            auto K_TT_jac =
                local_Jac.template block<temperature_size, temperature_size>(
                    temperature_index, temperature_index);
            K_TT_jac.noalias() +=
                dNdx.transpose() * dlambda_dT * dNdx * T * N * w;
            K_TT_jac.noalias() += dNdx.transpose() * lambda * dNdx * w;
            K_TT_jac.noalias() += drho_dT * specific_heat_capacity_fluid *
                                  N.transpose() * darcy_velocity.transpose() *
                                  dNdx * T * N * w;
            K_TT_jac.noalias() += rho * dspecific_heat_capacity_fluid_dT *
                                  N.transpose() * darcy_velocity.transpose() *
                                  dNdx * T * N * w;
            K_TT_jac.noalias() +=
                rho * specific_heat_capacity_fluid * N.transpose() *
                (K_over_mu2 * dNdx * p).transpose() * dmu_dT * dNdx * T * N * w;
            if (process_data.has_gravity)
            {
                K_TT_jac.noalias() += rho * specific_heat_capacity_fluid *
                                      N.transpose() *
                                      (K_over_mu2 * rho * b).transpose() *
                                      dmu_dT * dNdx * T * N * w;
                K_TT_jac.noalias() -=
                    rho * specific_heat_capacity_fluid * N.transpose() *
                    (K_over_mu * drho_dT * b).transpose() * dNdx * T * N * w;
            }
            K_TT_jac.noalias() += rho * specific_heat_capacity_fluid *
                                  N.transpose() * darcy_velocity.transpose() *
                                  dNdx * w;

            auto C_Tp_jac =
                local_Jac.template block<temperature_size, pressure_size>(
                    temperature_index, pressure_index);
            C_Tp_jac.noalias() += N.transpose() * dC_p_dp * N * T_dot * N * w;

            auto K_Tp_jac =
                local_Jac.template block<temperature_size, pressure_size>(
                    temperature_index, pressure_index);
            K_Tp_jac.noalias() +=
                dNdx.transpose() * dlambda_dp * dNdx * T * N * w;
            K_Tp_jac.noalias() += N.transpose() * drho_dp *
                                  specific_heat_capacity_fluid *
                                  darcy_velocity.transpose() * dNdx * T * N * w;
            K_Tp_jac.noalias() += N.transpose() * rho *
                                  dspecific_heat_capacity_fluid_dp *
                                  darcy_velocity.transpose() * dNdx * T * N * w;
            K_Tp_jac.noalias() +=
                N.transpose() * rho * specific_heat_capacity_fluid *
                (K_over_mu2 * dmu_dp * dNdx * p).transpose() * dNdx * T * N * w;
            K_Tp_jac.noalias() +=
                N.transpose() * rho * specific_heat_capacity_fluid *
                (K_over_mu * dNdx * p).transpose() * dNdx * T * N * w;
            if (process_data.has_gravity)
            {
                K_Tp_jac.noalias() +=
                    N.transpose() * rho * specific_heat_capacity_fluid *
                    (K_over_mu2 * dmu_dp * rho * b).transpose() *
                    dNdx * T * N * w;
                K_Tp_jac.noalias() +=
                    N.transpose() * rho * specific_heat_capacity_fluid *
                    (K_over_mu * drho_dp * b).transpose() * dNdx * T * N * w;
            }

            auto psi_p =
                local_rhs.template segment<pressure_size>(pressure_index);
            psi_p.noalias() -=
                dNdx.transpose() * K_over_mu * rho * dNdx * p * w;
            if (process_data.has_gravity)
            {
                psi_p.noalias() +=
                    dNdx.transpose() * K_over_mu * rho * rho * b * w;
            }
            psi_p.noalias() -=
                N.transpose() * porosity * drho_dp * N * p_dot * w;
            psi_p.noalias() -=
                N.transpose() * porosity * drho_dT * N * T_dot * w;

            auto psi_T =
                local_rhs.template segment<temperature_size>(temperature_index);
            psi_T.noalias() -= N.transpose() * C_p * N * T_dot * w;
            psi_T.noalias() -= dNdx.transpose() * lambda * dNdx * T * w;
            psi_T.noalias() -= (rho * specific_heat_capacity_fluid) *
                               (N.transpose() * darcy_velocity.transpose()) *
                               dNdx * T * w;
        }
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const indices =
            NumLib::getIndices(this->_element.getID(), dof_table);
        assert(!indices.empty());
        auto const& local_x = current_solution.get(indices);

        return this->getIntPtDarcyVelocityLocal(t, local_x, cache);
    }

private:
    using HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::pressure_index;
    using HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::pressure_size;
    using HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::temperature_index;
    using HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>::temperature_size;
};

}  // namespace HT
}  // namespace ProcessLib
