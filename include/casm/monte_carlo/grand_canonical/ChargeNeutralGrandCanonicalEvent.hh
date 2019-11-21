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
/// Zeyu: this is used as a data framework to store all informations related to charge neutral GCMC
/// Proposing a ChargeNeutralGrandCanonicalEvent will propose 3 GrandCanonicalEvent:
/// A event will be proposed: pick two sites, one Na/Va site and one Si/P site, and apply the same to_value() value in OccMod
/// in this case the charge is always balanced

class ChargeNeutralGrandCanonicalEvent {
	public:
		//Construct the event
		ChargeNeutralGrandCanonicalEvent(){
		};
		
    	/// \brief Set the change in (extensive) formation energy associated with this event
    	void set_dEf(double dE);

    	/// \brief Return change in (extensive) formation energy associated with this event
    	double dEf() const;


    	/// \brief Access change in number of species per supercell. Order as in CompositionConverter::components().
    	Eigen::VectorXl &dN();

    	/// \brief const Access change in number of species per supercell. Order as in CompositionConverter::components().
    	const Eigen::VectorXl &dN() const;

    	/// \brief Set the change in number of species in supercell. Order as in CompositionConverter::components().
    	void set_dN(size_type species_type_index, long int dn);

    	/// \brief Return change in number of species in supercell. Order as in CompositionConverter::components().
    	long int dN(size_type species_type_index) const;


    	/// \brief Set change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    	void set_dEpot(double dpot_nrg);

    	/// \brief Return change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    	double dEpot() const;

    	/// \brief Access the changes in (extensive) correlations associated with this event
    	Eigen::VectorXd &dCorr();

    	/// \brief const Access the changes in (extensive) correlations associated with this event
    	const Eigen::VectorXd &dCorr() const;

    	/// \brief  Access the occupational modification for this event
    	std::pair<OccMod,OccMod> &occupational_change();

    	/// \brief const Access the occupational modification for this event
    	const std::pair<OccMod,OccMod> &occupational_change();

  	private:
    	/// \brief Change in (extensive) correlations due to this event
    	Eigen::VectorXd m_dCorr;
    	/// \brief Change in (extensive) formation energy due to this event
    	double m_dEf;

    	/// \brief Change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    	double m_dEpot;

    	/// \brief Change in number of each species in supercell due to this event.
    	///        The order is determined by primclex.get_param_comp().get_components()
    	Eigen::VectorXl m_dN;

    	/// \brief The ConfigDoF modification performed by this event , Pairs
    	std::pair <OccMod,OccMod> m_occ_mod;

		/// dEpot for two swaps
		double m_dEpot_swapped_twice

};

  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
	ChargeNeutralGrandCanonicalEvent::ChargeNeutralGrandCanonicalEvent(){
	}

	  /// \brief Return change in total (formation) energy associated with this event
	  inline double ChargeNeutralGrandCanonicalEvent::dEf() const {
	    return m_dEf;
	  }

	  /// \brief Access change in number of all species (extensive). Order as in CompositionConverter::components().
	  inline Eigen::VectorXl &ChargeNeutralGrandCanonicalEvent::dN() {
	    return m_dN;
	  }

	  /// \brief const Access change in number of all species (extensive). Order as in CompositionConverter::components().
	  inline const Eigen::VectorXl &ChargeNeutralGrandCanonicalEvent::dN() const {
	    return m_dN;
	  }

	  /// \brief Set the change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
	  inline void ChargeNeutralGrandCanonicalEvent::set_dN(size_type species_type_index, long int dNi) {
	    m_dN(species_type_index) = dNi;
	  }

	  /// \brief Return change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
	  inline long int ChargeNeutralGrandCanonicalEvent::dN(size_type species_type_index) const {
	    return m_dN(species_type_index);
	  }

	  /// \brief Set the change in total (formation) energy associated with this event
	  inline void ChargeNeutralGrandCanonicalEvent::set_dEf(double dEf) {
	    m_dEf = dEf;
	  }

	  /// \brief Set the change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
	  inline void ChargeNeutralGrandCanonicalEvent::set_dEpot(double dEpot) {
	    m_dEpot = dEpot;
	  }
	  inline void ChargeNeutralGrandCanonicalEvent::set_dEpot_swapped_twice(double dEpot_swapped_twice) {
	    m_dEpot_swapped_twice = dEpot_swapped_twice;
	  }
	  inline void ChargeNeutralGrandCanonicalEvent::dEpot_swapped_twice() {
	    return m_dEpot_swapped_twice;
	  }
	  /// \brief Return change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
	  inline double ChargeNeutralGrandCanonicalEvent::dEpot() const {
	    return m_dEpot;
	
	  inline std::pair<OccMod,OccMod> &ChargeNeutralGrandCanonicalEvent::occupational_change(){
		  return m_occ_mod;
	  }

	  inline const std::pair<OccMod,OccMod> &ChargeNeutralGrandCanonicalEvent::occupational_change(){
		  return m_occ_mod;
	  }
	  }

	  }
}


#endif
