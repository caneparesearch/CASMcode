#ifndef CASM_ChargeNeutralGrandCanonical_HH
#define CASM_ChargeNeutralGrandCanonical_HH
#include <vector>
#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"
#include "casm/monte_carlo/DoFMod.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalEvent.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"

namespace CASM{

/// \brief Data structure for storing information regarding a proposed charge neutral grand canonical Monte Carlo event
/// Jerry: this is used as a data framework to store all informations related to charge neutral GCMC
/// Proposing a ChargeNeutralGrandCanonicalEvent will propose 3 GrandCanonicalEvent:
/// The 1st event is a random swap, i.e. swap either Na/Va or Si/P; 
/// The 2nd event will swap the other site to make charge neutral, i.e. if in 1st event Na/Va was swapped, swap Si/P.
/// Then check if the charge is neutral, if so, propose the 3rd event.
/// The 3rd event is a combination of 1st and 2nd events. The energy will be calculated based on this event for GCMC.

class ChargeNeutralGrandCanonicalEvent {
	public:
		//Construct the event
		ChargeNeutralGrandCanonicalEvent(){
		};
	
		//Calculates the change in energy due to application of the event
		calculate_deltaE(EventType &event);

  	private:

    	/// \brief Change in (extensive) formation energy due to this event
    	double m_dEf;

    	/// \brief Change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    	double m_dEpot;

    	/// \brief Change in number of each species in supercell due to this event.
    	///        The order is determined by primclex.get_param_comp().get_components()
    	Eigen::VectorXl m_dN;

    	/// \brief The ConfigDoF modification performed by this event
    	OccMod m_occ_mod;
};

  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
	ChargeNeutralGrandCanonicalEvent::ChargeNeutralGrandCanonicalEvent(){
	// propose the 1st event
		GrandCanonicalEvent event_1 = GrandCanonical::propose();
	// if the site is Na/Va
		if (event_1.occupational_change().sublat < 10){
			// then keep on trying proposing new event to meet the charge balance condition
			do {
				GrandCanonicalEvent event_2 = GrandCanonical::propose();
			}
			while(event_2.occupational_change().sublat < 10);
			// do the 3rd event
			
			}
		else{
			do{
				GrandCanonicalEvent event_2 = GrandCanonical::propose();
			}
			while(event_2.occupational_change().sublat < 10);
			// do the 3rd event
			
		}

	}

	bool check_charge(EventType &event1,EventType &event2){
		if(event_1.occupational_change().sublat)
	}

	double calculate_deltaE(EventType &event){
		return m_deltaE;
	}

	  /// \brief Set the change in total (formation) energy associated with this event
	  inline void GrandCanonicalEvent::set_dEf(double dEf) {
	    m_dEf = dEf;
	  }

	  /// \brief Return change in total (formation) energy associated with this event
	  inline double GrandCanonicalEvent::dEf() const {
	    return m_dEf;
	  }

	  /// \brief Access change in number of all species (extensive). Order as in CompositionConverter::components().
	  inline Eigen::VectorXl &GrandCanonicalEvent::dN() {
	    return m_dN;
	  }

	  /// \brief const Access change in number of all species (extensive). Order as in CompositionConverter::components().
	  inline const Eigen::VectorXl &GrandCanonicalEvent::dN() const {
	    return m_dN;
	  }

	  /// \brief Set the change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
	  inline void GrandCanonicalEvent::set_dN(size_type species_type_index, long int dNi) {
	    m_dN(species_type_index) = dNi;
	  }

	  /// \brief Return change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
	  inline long int GrandCanonicalEvent::dN(size_type species_type_index) const {
	    return m_dN(species_type_index);
	  }

	  /// \brief Set the change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
	  inline void GrandCanonicalEvent::set_dEpot(double dEpot) {
	    m_dEpot = dEpot;
	  }

	  /// \brief Return change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
	  inline double GrandCanonicalEvent::dEpot() const {
	    return m_dEpot;
	  }


	  }
}


#endif
