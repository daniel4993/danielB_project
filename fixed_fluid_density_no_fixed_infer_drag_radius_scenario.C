//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// DIP - Droplet Inverse Problem
//
// Copyright (C) 2017 Paul T. Bauman, Matthew Ringuette, David Salac
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

// This class
#include <dip/fixed_fluid_density_no_fixed_infer_drag_radius_scenario.h>

// DIP
#include <dip/droplet_data.h>
#include <dip/droplet_likelihood.h>

// QUESO
#include <queso/Defines.h>
#include <queso/getpot.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/BoxSubset.h>
#include <queso/UniformVectorRV.h>
#include <queso/1DQuadrature.h>
#include <queso/TensorProductQuadrature.h>

// C++
#include <sstream>

namespace DIP
{

  template<typename VectorType, typename MatrixType>
  FixedFluidDensityNoFixedInferDragRadiusScenario<VectorType,MatrixType>::
  FixedFluidDensityNoFixedInferDragRadiusScenario( const QUESO::GetPot & input,
                                                               const QUESO::BaseEnvironment & queso_env,
                                                               const DropletData<VectorType,MatrixType> & data )
    : FixedFluidDensityNoFixedScenarioBase<VectorType,MatrixType>(input,data),
    //_n_params(1+this->_data.n_data())
      _n_params(3)
    //_n_marg_params(1)
  {
    // Setup parameter space, domain, and RV
    this->_param_space.reset( new QUESO::VectorSpace<VectorType,MatrixType>(queso_env,
                                                                            "param_space_",
                                                                            this->_n_params,
                                                                            NULL) );

    std::unique_ptr<VectorType> param_mins( this->_param_space->newVector() );
    std::unique_ptr<VectorType> param_maxs( this->_param_space->newVector() );

    (*param_mins)[0] = this->_rho_d_min;
    (*param_mins)[1] = this->_cd_min;
    (*param_mins)[2] = this->_r_d_min;


    (*param_maxs)[0] = this->_rho_d_max;
    (*param_maxs)[1] = this->_cd_max;
    (*param_maxs)[2] = this->_r_d_max;

    this->_param_domain.reset( new QUESO::BoxSubset<VectorType,MatrixType>("param_domain_", *(this->_param_space),
                                                                           *param_mins,*param_maxs) );

    this->_prior_rv.reset( new QUESO::UniformVectorRV<VectorType,MatrixType>("prior_rv_",
                                                                             *(this->_param_domain)) );

    this->_likelihood = DropletLikelihood<VectorType,MatrixType>::build_no_marg(*this);
  }

  template<typename VectorType, typename MatrixType>
  void
  FixedFluidDensityNoFixedInferDragRadiusScenario<VectorType,MatrixType>::
  eval_model( const VectorType & domain_vector,
              VectorType & model_output ) const
  {
    queso_assert_equal_to(domain_vector.sizeGlobal(), this->_n_params);
   // queso_assert_equal_to(marg_vector.sizeGlobal(), this->_n_marg_params);
    queso_assert_equal_to(model_output.sizeGlobal(), this->_data.n_data());

    const double rho_d = this->_rho_d_nom*domain_vector[0];
    //const double _cd = this->_cd*domain_vector[0];


   for( unsigned int i = 0; i < this->_data.n_data(); i++ )
      {
        const double t = this->_data.time(i);

        const double r_d = this->_data.droplet_radius(i);

        model_output[i] = this->_model.position(this->_rho_f,
                                                rho_d,
                                                r_d,
                                                this->_cd,
                                                t);
      }
  }

} // end namespace DIP

// Instantiate
template class DIP::FixedFluidDensityNoFixedInferDragRadiusScenario<QUESO::GslVector,QUESO::GslMatrix>;
