/*
*  adapt_lif.cpp
*
*  This file is part of NEST.
*
*  Copyright (C) 2004 The NEST Initiative
*
*  NEST is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 2 of the License, or
*  (at your option) any later version.
*
*  NEST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
*
*  2019-04-28 19:32:29.743894
*/

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "adapt_lif.h"


/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<adapt_lif> adapt_lif::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <> void RecordablesMap<adapt_lif>::create(){
  // use standard names whereever you can for consistency!  

  insert_("refr_spikes_buffer", &adapt_lif::get_refr_spikes_buffer);

  insert_("V_abs", &adapt_lif::get_V_abs);

  insert_("w", &adapt_lif::get_w);
  }
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * Note: the implementation is empty. The initialization is of variables
 * is a part of the adapt_lif's constructor.
 * ---------------------------------------------------------------- */
adapt_lif::Parameters_::Parameters_(){}

adapt_lif::State_::State_(){}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

adapt_lif::Buffers_::Buffers_(adapt_lif &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

adapt_lif::Buffers_::Buffers_(const Buffers_ &, adapt_lif &n):
  logger_(n), __s( 0 ), __c( 0 ), __e( 0 ){
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */
adapt_lif::adapt_lif():Archiving_Node(), P_(), S_(), B_(*this)
{
  recordablesMap_.create();
  // use a default `good` enough value for the absolute error.
  // it cab be adjusted via `SetStatus`
  P_.__gsl_error_tol = 1e-3;
  
  P_.tau_m = 10*1.0; // as ms
  
  P_.C_m = 250*1.0; // as pF
  
  P_.t_ref = 2*1.0; // as ms
  
  P_.tau_syn = 2*1.0; // as ms
  
  P_.E_L = 0*1.0; // as mV
  
  P_.I_e = 0*1.0; // as pA
  
  P_.V_reset = 10*1.0; // as mV
  
  P_.Theta = 20*1.0; // as mV
  
  P_.tau_w = 100*1.0; // as ms
  
  P_.a = 1*1.0; // as nS
  
  P_.b = 2*1.0; // as pA
  
  P_.with_refr_input = false; // as boolean
  
  S_.refr_spikes_buffer = 0*1.0; // as mV
  
  S_.r = 0; // as integer
  
  S_.ode_state[State_::V_abs] = 0*1.0; // as mV
  
  S_.ode_state[State_::w] = 0*1.0; // as pA
}

adapt_lif::adapt_lif(const adapt_lif& __n):
  Archiving_Node(), P_(__n.P_), S_(__n.S_), B_(__n.B_, *this){
  P_.tau_m = __n.P_.tau_m;
  P_.C_m = __n.P_.C_m;
  P_.t_ref = __n.P_.t_ref;
  P_.tau_syn = __n.P_.tau_syn;
  P_.E_L = __n.P_.E_L;
  P_.I_e = __n.P_.I_e;
  P_.V_reset = __n.P_.V_reset;
  P_.Theta = __n.P_.Theta;
  P_.tau_w = __n.P_.tau_w;
  P_.a = __n.P_.a;
  P_.b = __n.P_.b;
  P_.with_refr_input = __n.P_.with_refr_input;
  
  S_.refr_spikes_buffer = __n.S_.refr_spikes_buffer;
  S_.r = __n.S_.r;
  
  S_.ode_state[State_::V_abs] = __n.S_.ode_state[State_::V_abs];
  S_.ode_state[State_::w] = __n.S_.ode_state[State_::w];
  
  V_.h = __n.V_.h;
  V_.RefractoryCounts = __n.V_.RefractoryCounts;
  V_.__h = __n.V_.__h;
  
}

adapt_lif::~adapt_lif(){ 
  // GSL structs may not have been allocated, so we need to protect destruction
  if (B_.__s)
    gsl_odeiv_step_free( B_.__s );
  if (B_.__c)
    gsl_odeiv_control_free( B_.__c );
  if (B_.__e)
    gsl_odeiv_evolve_free( B_.__e );
}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void adapt_lif::init_state_(const Node& proto){
  const adapt_lif& pr = downcast<adapt_lif>(proto);
  S_ = pr.S_;
}



extern "C" inline int adapt_lif_dynamics(double, const double ode_state[], double f[], void* pnode){
  typedef adapt_lif::State_ State_;
  // get access to node so we can almost work as in a member function
  assert( pnode );
  const adapt_lif& node = *( reinterpret_cast< adapt_lif* >( pnode ) );

  // ode_state[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.ode_state[].
  
  f[State_::V_abs] = -(ode_state[State_::V_abs]) / node.get_tau_m() - ( ode_state[State_::w]) / node.get_C_m();
  f[State_::w] = -(ode_state[State_::w]) / node.get_tau_w();
  return GSL_SUCCESS;
}



void adapt_lif::init_buffers_(){
  get_spikes().clear(); //includes resize
  get_currents().clear(); //includes resize
  
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
  
  if ( B_.__s == 0 ){
    B_.__s = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, 2 );
  } else {
    gsl_odeiv_step_reset( B_.__s );
  }

  if ( B_.__c == 0 ){
    B_.__c = gsl_odeiv_control_y_new( P_.__gsl_error_tol, 0.0 );
  } else {
    gsl_odeiv_control_init( B_.__c, P_.__gsl_error_tol, 0.0, 1.0, 0.0 );
  }

  if ( B_.__e == 0 ){
    B_.__e = gsl_odeiv_evolve_alloc( 2 );
  } else {
    gsl_odeiv_evolve_reset( B_.__e );
  }

  B_.__sys.function = adapt_lif_dynamics;
  B_.__sys.jacobian = NULL;
  B_.__sys.dimension = 2;
  B_.__sys.params = reinterpret_cast< void* >( this );
  B_.__step = nest::Time::get_resolution().get_ms();
  B_.__integration_step = nest::Time::get_resolution().get_ms();
}

void adapt_lif::calibrate(){
  B_.logger_.init();
  
  
  V_.h =nest::Time::get_resolution().get_ms();
  
  
  V_.RefractoryCounts =nest::Time(nest::Time::ms((double) P_.t_ref)).get_steps();
  
  
  V_.__h =nest::Time::get_resolution().get_ms();
}

/* ----------------------------------------------------------------
* Update and spike handling functions
* ---------------------------------------------------------------- */

/*
 *
 */
void adapt_lif::update(nest::Time const & origin,const long from, const long to){
  double __t = 0;

  for ( long lag = from ; lag < to ; ++lag ) {
    B_.spikes_grid_sum_ = get_spikes().get_value(lag);
    B_.currents_grid_sum_ = get_currents().get_value(lag);
      
    
    __t = 0;
    while ( __t < B_.__step )
    {
      const int status = gsl_odeiv_evolve_apply(B_.__e,
                                                B_.__c,
                                                B_.__s,
                                                &B_.__sys,              // system of ODE
                                                &__t,                   // from t
                                                B_.__step,              // to t <= step
                                                &B_.__integration_step, // integration step size
                                                S_.ode_state);          // neuronal state

      if ( status != GSL_SUCCESS ) {
        throw nest::GSLSolverFailure( get_name(), status );
      }
    }
    



    if (S_.r>0) {
      

      S_.r = S_.r-1;
      

      S_.ode_state[State_::V_abs] = P_.V_reset;
    }else
    {
	  S_.ode_state[State_::V_abs] += B_.spikes_grid_sum_;

    if(S_.ode_state[State_::V_abs]>=P_.Theta) {
      

      S_.ode_state[State_::w] += P_.b;
      

      S_.r = V_.RefractoryCounts;
      

      S_.ode_state[State_::V_abs] = P_.V_reset;
      
      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
    } /* if end */
    }

    // voltage logging
    B_.logger_.record_data(origin.get_steps()+lag);
  }

}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void adapt_lif::handle(nest::DataLoggingRequest& e){
  B_.logger_.handle(e);
}


void adapt_lif::handle(nest::SpikeEvent &e){
  assert(e.get_delay_steps() > 0);
  
  const double weight = e.get_weight();
  const double multiplicity = e.get_multiplicity();
  
  if ( weight >= 0.0 ){ // excitatory
    get_spikes().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                       weight * multiplicity );
  }
  if ( weight < 0.0 ){ // inhibitory
    get_spikes().
        add_value(e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
                      
                       weight * multiplicity );
  }
}

void adapt_lif::handle(nest::CurrentEvent& e){
  assert(e.get_delay_steps() > 0);

  const double current=e.get_current();
  const double weight=e.get_weight();

  // add weighted current; HEP 2002-10-04
  get_currents().add_value(
               e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin()),
               weight * current );
  
}
