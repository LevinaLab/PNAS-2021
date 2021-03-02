
/*
*  adapt_lif.h
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
*  2019-04-28 19:32:29.708206
*/
#ifndef ADAPT_LIF
#define ADAPT_LIF

#include "config.h"


#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// forwards the declaration of the function
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" inline int adapt_lif_dynamics( double, const double y[], double f[], void* pnode );


// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


// Includes from sli:
#include "dictdatum.h"

/* BeginDocumentation
  Name: adapt_lif.

  Description:  
    
  Name: iaf_psc_delta - Adaptive leaky integrate-and-fire neuron with 
  delta-shaped PSCs.

  Description:

  based on iaf_psc_delta and aeif_cond_alpha
  iaf_psc_delta is an implementation of a leaky integrate-and-fire model
  where the potential jumps on each spike arrival.

  The threshold crossing is followed by an absolute refractory period
  during which the membrane potential is clamped to the resting potential.

  Spikes arriving while the neuron is refractory, are discarded by
  default. If the property "refractory_input" is set to true, such
  spikes are added to the membrane potential at the end of the
  refractory period, dampened according to the interval between
  arrival and end of refractoriness.

  The linear subthresold dynamics is integrated by the Exact
  Integration scheme [1]. The neuron dynamics is solved on the time
  grid given by the computation step size. Incoming as well as emitted
  spikes are forced to that grid.

  An additional state variable and the corresponding differential
  equation represents a piecewise constant external current.

  The general framework for the consistent formulation of systems with
  neuron like dynamics interacting by point events is described in
  [1].  A flow chart can be found in [2].

  Critical tests for the formulation of the neuron model are the
  comparisons of simulation results for different computation step
  sizes. sli/testsuite/nest contains a number of such tests.

  The iaf_psc_delta is the standard model used to check the consistency
  of the nest simulation kernel because it is at the same time complex
  enough to exhibit non-trivial dynamics and simple enough compute
  relevant measures analytically.

  Remarks:

  The present implementation uses individual variables for the
  components of the state vector and the non-zero matrix elements of
  the propagator.  Because the propagator is a lower triangular matrix
  no full matrix multiplication needs to be carried out and the
  computation can be done "in place" i.e. no temporary state vector
  object is required.

  The template support of recent C++ compilers enables a more succinct
  formulation without loss of runtime performance already at minimal
  optimization levels. A future version of iaf_psc_delta will probably
  address the problem of efficient usage of appropriate vector and
  matrix objects.

  References:
  [1] Rotter S & Diesmann M (1999) Exact digital simulation of time-invariant
  linear systems with applications to neuronal modeling. Biologial Cybernetics
  81:381-402.
  [2] Diesmann M, Gewaltig M-O, Rotter S, & Aertsen A (2001) State space
  analysis of synchronous spiking in cortical neural networks.
  Neurocomputing 38-40:565-571.

  Sends: SpikeEvent

  Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

  Author:  September 1999, Diesmann, Gewaltig
  SeeAlso: iaf_psc_alpha, iaf_psc_exp, iaf, iaf_psc_delta_canon



  Parameters:
  The following parameters can be set in the status dictionary.
  tau_m [ms]  Membrane time constant.
  C_m [pF]  Capacity of the membrane
  t_ref [ms]  Duration of refractory period.
  tau_syn [ms]  Time constant of synaptic current.
  E_L [mV]  Resting membrane potential.
  I_e [pA]  Constant input current.
  V_reset [mV]  Reset potential of the membrane.
  Theta [mV]  Spike threshold.
  b [pA] V_min mV = -inf * 1 mV            Absolute lower value for the membrane potential
  with_refr_input [boolean] V_min mV = -inf * 1 mV            Absolute lower value for the membrane potential
 If true, do not discard input during  refractory period. Default: false.
  

  Dynamic state variables:
  r [integer]  counts number of tick during the refractory period
  

  Initial values:
  V_abs [mV] function V_m mV = V_abs + E_L  Membrane potential.
  w [pA] function V_m mV = V_abs + E_L  Membrane potential.
  

  References: Empty

  Sends: nest::SpikeEvent

  Receives: Spike, Current, DataLoggingRequest
*/
class adapt_lif : public nest::Archiving_Node{
public:
  /**
  * The constructor is only used to create the model prototype in the model manager.
  */
  adapt_lif();

  /**
  * The copy constructor is used to create model copies and instances of the model.
  * @node The copy constructor needs to initialize the parameters and the state.
  *       Initialization of buffers and interal variables is deferred to
  *       @c init_buffers_() and @c calibrate().
  */
  adapt_lif(const adapt_lif &);

  /**
  * Releases resources.
  */
  ~adapt_lif();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handles_test_event;
  using nest::Node::handle;

  /**
  * Used to validate that we can send nest::SpikeEvent to desired target:port.
  */
  nest::port send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool);

  /**
  * @defgroup mynest_handle Functions handling incoming events.
  * We tell nest that we can handle incoming events of various types by
  * defining @c handle() and @c connect_sender() for the given event.
  * @{
  */
  void handle(nest::SpikeEvent &);        //! accept spikes
  void handle(nest::CurrentEvent &);      //! accept input current
  void handle(nest::DataLoggingRequest &);//! allow recording with multimeter

  nest::port handles_test_event(nest::SpikeEvent&, nest::port);
  nest::port handles_test_event(nest::CurrentEvent&, nest::port);
  nest::port handles_test_event(nest::DataLoggingRequest&, nest::port);
  /** @} */

  // SLI communication functions:
  void get_status(DictionaryDatum &) const;
  void set_status(const DictionaryDatum &);

private:
  //! Reset parameters and state of neuron.

  //! Reset state of neuron.
  void init_state_(const Node& proto);

  //! Reset internal buffers of neuron.
  void init_buffers_();

  //! Initialize auxiliary quantities, leave parameters and state untouched.
  void calibrate();

  //! Take neuron through given time interval
  void update(nest::Time const &, const long, const long);

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap<adapt_lif>;
  friend class nest::UniversalDataLogger<adapt_lif>;

  /**
  * Free parameters of the neuron.
  *
  *
  *
  * These are the parameters that can be set by the user through @c SetStatus.
  * They are initialized from the model prototype when the node is created.
  * Parameters do not change during calls to @c update() and are not reset by
  * @c ResetNetwork.
  *
  * @note Parameters_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If Parameters_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If Parameters_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct Parameters_{
        
        

    //!  Membrane time constant.
    double tau_m;

    //!  Capacity of the membrane
    double C_m;

    //!  Duration of refractory period.
    double t_ref;

    //!  Time constant of synaptic current.
    double tau_syn;

    //!  Resting membrane potential.
    double E_L;

    //!  Constant input current.
    double I_e;

    //!  Reset potential of the membrane.
    double V_reset;

    //!  Spike threshold.
    double Theta;

    double tau_w;

    double a;

    //! V_min mV = -inf * 1 mV            Absolute lower value for the membrane potential
    double b;

    //! V_min mV = -inf * 1 mV            Absolute lower value for the membrane potential
    //!  If true, do not discard input during  refractory period. Default: false.
    bool with_refr_input;

    double __gsl_error_tol;
    /** Initialize parameters to their default values. */
    Parameters_();
  };

  /**
  * Dynamic state of the neuron.
  *
  *
  *
  * These are the state variables that are advanced in time by calls to
  * @c update(). In many models, some or all of them can be set by the user
  * through @c SetStatus. The state variables are initialized from the model
  * prototype when the node is created. State variables are reset by @c ResetNetwork.
  *
  * @note State_ need neither copy constructor nor @c operator=(), since
  *       all its members are copied properly by the default copy constructor
  *       and assignment operator. Important:
  *       - If State_ contained @c Time members, you need to define the
  *         assignment operator to recalibrate all members of type @c Time . You
  *         may also want to define the assignment operator.
  *       - If State_ contained members that cannot copy themselves, such
  *         as C-style arrays, you need to define the copy constructor and
  *         assignment operator to copy those members.
  */
  struct State_{
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems{
    // function V_m mV = V_abs + E_L  Membrane potential.
      V_abs,
      // function V_m mV = V_abs + E_L  Membrane potential.
      w,
      STATE_VEC_SIZE
    };
    //! state vector, must be C-array for GSL solver
    double ode_state[STATE_VEC_SIZE];    

    double refr_spikes_buffer;

    //!  counts number of tick during the refractory period
    long r;    
        State_();
  };

  /**
  * Internal variables of the neuron.
  *
  *
  *
  * These variables must be initialized by @c calibrate, which is called before
  * the first call to @c update() upon each call to @c Simulate.
  * @node Variables_ needs neither constructor, copy constructor or assignment operator,
  *       since it is initialized by @c calibrate(). If Variables_ has members that
  *       cannot destroy themselves, Variables_ will need a destructor.
  */
  struct Variables_ {    

    double h;
        

    //!  refractory time in steps
    long RefractoryCounts;
        

    double __h;
    
  };

  /**
    * Buffers of the neuron.
    * Ususally buffers for incoming spikes and data logged for analog recorders.
    * Buffers must be initialized by @c init_buffers_(), which is called before
    * @c calibrate() on the first call to @c Simulate after the start of NEST,
    * ResetKernel or ResetNetwork.
    * @node Buffers_ needs neither constructor, copy constructor or assignment operator,
    *       since it is initialized by @c init_nodes_(). If Buffers_ has members that
    *       cannot destroy themselves, Buffers_ will need a destructor.
    */
  struct Buffers_ {
    Buffers_(adapt_lif &);
    Buffers_(const Buffers_ &, adapt_lif &);

    /** Logger for all analog data */
    nest::UniversalDataLogger<adapt_lif> logger_;
    
    inline nest::RingBuffer& get_spikes() {return spikes;}
    //!< Buffer incoming pAs through delay, as sum
    nest::RingBuffer spikes;
    double spikes_grid_sum_;
    
    //!< Buffer incoming pAs through delay, as sum
    nest::RingBuffer currents;
    inline nest::RingBuffer& get_currents() {return currents;}
    double currents_grid_sum_;
    /** GSL ODE stuff */
    gsl_odeiv_step* __s;    //!< stepping function
    gsl_odeiv_control* __c; //!< adaptive stepsize control function
    gsl_odeiv_evolve* __e;  //!< evolution function
    gsl_odeiv_system __sys; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double __step;             //!< step size in ms
    double __integration_step; //!< current integration time step, updated by GSL
    };
  inline double get_refr_spikes_buffer() const {
    return S_.refr_spikes_buffer;
  }
  inline void set_refr_spikes_buffer(const double __v) {
    S_.refr_spikes_buffer = __v;
  }

  inline long get_r() const {
    return S_.r;
  }
  inline void set_r(const long __v) {
    S_.r = __v;
  }

  inline double get_V_abs() const {
    return S_.ode_state[State_::V_abs];
  }
  inline void set_V_abs(const double __v) {
    S_.ode_state[State_::V_abs] = __v;
  }

  inline double get_w() const {
    return S_.ode_state[State_::w];
  }
  inline void set_w(const double __v) {
    S_.ode_state[State_::w] = __v;
  }

  inline double get_tau_m() const {
    return P_.tau_m;
  }
  inline void set_tau_m(const double __v) {
    P_.tau_m = __v;
  }

  inline double get_C_m() const {
    return P_.C_m;
  }
  inline void set_C_m(const double __v) {
    P_.C_m = __v;
  }

  inline double get_t_ref() const {
    return P_.t_ref;
  }
  inline void set_t_ref(const double __v) {
    P_.t_ref = __v;
  }

  inline double get_tau_syn() const {
    return P_.tau_syn;
  }
  inline void set_tau_syn(const double __v) {
    P_.tau_syn = __v;
  }

  inline double get_E_L() const {
    return P_.E_L;
  }
  inline void set_E_L(const double __v) {
    P_.E_L = __v;
  }

  inline double get_I_e() const {
    return P_.I_e;
  }
  inline void set_I_e(const double __v) {
    P_.I_e = __v;
  }

  inline double get_V_reset() const {
    return P_.V_reset;
  }
  inline void set_V_reset(const double __v) {
    P_.V_reset = __v;
  }

  inline double get_Theta() const {
    return P_.Theta;
  }
  inline void set_Theta(const double __v) {
    P_.Theta = __v;
  }

  inline double get_tau_w() const {
    return P_.tau_w;
  }
  inline void set_tau_w(const double __v) {
    P_.tau_w = __v;
  }

  inline double get_a() const {
    return P_.a;
  }
  inline void set_a(const double __v) {
    P_.a = __v;
  }

  inline double get_b() const {
    return P_.b;
  }
  inline void set_b(const double __v) {
    P_.b = __v;
  }

  inline bool get_with_refr_input() const {
    return P_.with_refr_input;
  }
  inline void set_with_refr_input(const bool __v) {
    P_.with_refr_input = __v;
  }

  inline double get_h() const {
    return V_.h;
  }
  inline void set_h(const double __v) {
    V_.h = __v;
  }

  inline long get_RefractoryCounts() const {
    return V_.RefractoryCounts;
  }
  inline void set_RefractoryCounts(const long __v) {
    V_.RefractoryCounts = __v;
  }

  inline double get___h() const {
    return V_.__h;
  }
  inline void set___h(const double __v) {
    V_.__h = __v;
  }


  
  inline nest::RingBuffer& get_spikes() {return B_.get_spikes();};
  
  inline nest::RingBuffer& get_currents() {return B_.get_currents();};
  

  // Generate function header
  
  /**
  * @defgroup pif_members Member variables of neuron model.
  * Each model neuron should have precisely the following four data members,
  * which are one instance each of the parameters, state, buffers and variables
  * structures. Experience indicates that the state and variables member should
  * be next to each other to achieve good efficiency (caching).
  * @note Devices require one additional data member, an instance of the @c Device
  *       child class they belong to.
  * @{
  */
  Parameters_ P_;  //!< Free parameters.
  State_      S_;  //!< Dynamic state.
  Variables_  V_;  //!< Internal Variables
  Buffers_    B_;  //!< Buffers.

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap<adapt_lif> recordablesMap_;

  friend int adapt_lif_dynamics( double, const double y[], double f[], void* pnode );
  
/** @} */
}; /* neuron adapt_lif */

inline nest::port adapt_lif::send_test_event(
    nest::Node& target, nest::rport receptor_type, nest::synindex, bool){
  // You should usually not change the code in this function.
  // It confirms that the target of connection @c c accepts @c nest::SpikeEvent on
  // the given @c receptor_type.
  nest::SpikeEvent e;
  e.set_sender(*this);
  return target.handles_test_event(e, receptor_type);
}

inline nest::port adapt_lif::handles_test_event(nest::SpikeEvent&, nest::port receptor_type){
  
    // You should usually not change the code in this function.
    // It confirms to the connection management system that we are able
    // to handle @c SpikeEvent on port 0. You need to extend the function
    // if you want to differentiate between input ports.
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
}



inline nest::port adapt_lif::handles_test_event(
    nest::CurrentEvent&, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c CurrentEvent on port 0. You need to extend the function
  // if you want to differentiate between input ports.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());
  return 0;
}

inline nest::port adapt_lif::handles_test_event(
    nest::DataLoggingRequest& dlr, nest::port receptor_type){
  // You should usually not change the code in this function.
  // It confirms to the connection management system that we are able
  // to handle @c DataLoggingRequest on port 0.
  // The function also tells the built-in UniversalDataLogger that this node
  // is recorded from and that it thus needs to collect data during simulation.
  if (receptor_type != 0)
  throw nest::UnknownReceptorType(receptor_type, get_name());

  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

// TODO call get_status on used or internal components
inline void adapt_lif::get_status(DictionaryDatum &__d) const{  
  def<double>(__d, "tau_m", get_tau_m());
      
  def<double>(__d, "C_m", get_C_m());
      
  def<double>(__d, "t_ref", get_t_ref());
      
  def<double>(__d, "tau_syn", get_tau_syn());
      
  def<double>(__d, "E_L", get_E_L());
      
  def<double>(__d, "I_e", get_I_e());
      
  def<double>(__d, "V_reset", get_V_reset());
      
  def<double>(__d, "Theta", get_Theta());
      
  def<double>(__d, "tau_w", get_tau_w());
      
  def<double>(__d, "a", get_a());
      
  def<double>(__d, "b", get_b());
      
  def<bool>(__d, "with_refr_input", get_with_refr_input());
      
  def<double>(__d, "refr_spikes_buffer", get_refr_spikes_buffer());
      
  def<long>(__d, "r", get_r());
      
  def<double>(__d, "V_abs", get_V_abs());
      
  def<double>(__d, "w", get_w());
    

  (*__d)[nest::names::recordables] = recordablesMap_.get_list();
  
  def< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  

}

inline void adapt_lif::set_status(const DictionaryDatum &__d){

  double tmp_tau_m = get_tau_m();
  updateValue<double>(__d, "tau_m", tmp_tau_m);


  double tmp_C_m = get_C_m();
  updateValue<double>(__d, "C_m", tmp_C_m);


  double tmp_t_ref = get_t_ref();
  updateValue<double>(__d, "t_ref", tmp_t_ref);


  double tmp_tau_syn = get_tau_syn();
  updateValue<double>(__d, "tau_syn", tmp_tau_syn);


  double tmp_E_L = get_E_L();
  updateValue<double>(__d, "E_L", tmp_E_L);


  double tmp_I_e = get_I_e();
  updateValue<double>(__d, "I_e", tmp_I_e);


  double tmp_V_reset = get_V_reset();
  updateValue<double>(__d, "V_reset", tmp_V_reset);


  double tmp_Theta = get_Theta();
  updateValue<double>(__d, "Theta", tmp_Theta);


  double tmp_tau_w = get_tau_w();
  updateValue<double>(__d, "tau_w", tmp_tau_w);


  double tmp_a = get_a();
  updateValue<double>(__d, "a", tmp_a);


  double tmp_b = get_b();
  updateValue<double>(__d, "b", tmp_b);


  bool tmp_with_refr_input = get_with_refr_input();
  updateValue<bool>(__d, "with_refr_input", tmp_with_refr_input);


  double tmp_refr_spikes_buffer = get_refr_spikes_buffer();
  updateValue<double>(__d, "refr_spikes_buffer", tmp_refr_spikes_buffer);

  

  long tmp_r = get_r();
  updateValue<long>(__d, "r", tmp_r);

  

  double tmp_V_abs = get_V_abs();
  updateValue<double>(__d, "V_abs", tmp_V_abs);

  

  double tmp_w = get_w();
  updateValue<double>(__d, "w", tmp_w);

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status(__d);

  // if we get here, temporaries contain consistent set of properties


  set_tau_m(tmp_tau_m);



  set_C_m(tmp_C_m);



  set_t_ref(tmp_t_ref);



  set_tau_syn(tmp_tau_syn);



  set_E_L(tmp_E_L);



  set_I_e(tmp_I_e);



  set_V_reset(tmp_V_reset);



  set_Theta(tmp_Theta);



  set_tau_w(tmp_tau_w);



  set_a(tmp_a);



  set_b(tmp_b);



  set_with_refr_input(tmp_with_refr_input);



  set_refr_spikes_buffer(tmp_refr_spikes_buffer);



  set_r(tmp_r);



  set_V_abs(tmp_V_abs);



  set_w(tmp_w);


  
  updateValue< double >(__d, nest::names::gsl_error_tol, P_.__gsl_error_tol);
  if ( P_.__gsl_error_tol <= 0. ){
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
  
};

#endif /* #ifndef ADAPT_LIF */
#endif /* HAVE GSL */
