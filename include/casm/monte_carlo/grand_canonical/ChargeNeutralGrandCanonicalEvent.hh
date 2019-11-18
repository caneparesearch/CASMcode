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
		double m_deltaE;

};

  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
ChargeNeutralGrandCanonicalEvent::ChargeNeutralGrandCanonicalEvent(){
// propose the 1st event
	GrandCanonicalEvent event_1 = GrandCanonicalEvent();
	if (event_1.occupational_change().sublat < 10){
		GrandCanonicalEvent event_2 = GrandCanonicalEvent();
		if (event_2.occupational_change().sublat > 10) {
			
		}
	}

}

double calculate_deltaE(EventType &event){
	&event.calculateE()
	return m_deltaE;
}

}


#endif
