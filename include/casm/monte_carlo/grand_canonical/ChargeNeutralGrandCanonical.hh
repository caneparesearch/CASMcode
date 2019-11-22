#ifndef CASM_ChargeNeutralGrandCanonical_HH
#define CASM_ChargeNeutralGrandCanonical_HH

#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings.hh"
#include "casm/monte_carlo/grand_canonical/ChargeNeutralGrandCanonicalEvent.hh"
#include "casm/clex/Clex.hh"
#include "casm/monte_carlo/SiteExchanger.hh"
namespace CASM {
  ///
  /// Derives from base MonteCarlo class, to be used for simulations at constant
  /// temperature and chemical potential and accounts for charge neutral swaps.
  ///
  /// As with all the other derived Monte Carlo classes, member functions must
  /// follow a specific naming convention to be used with templated routines currently
  /// defined in MonteDriver.hh:
  ///      -conditions (should be able to steal from GrandCanonical)
  ///      -set_conditions (should be able to steal from GrandCanonical)
  ///      -propose (we have to write this for our specialized event)
  ///      -check  (we have to write this for our specialized event)
  ///      -accept(we have to write this for our specialized event)
  ///      -reject(we have to write this for our specialized event)
  ///      -write_results
  ///
  //


class ChargeNeutralGrandCanonical : public MonteCarlo{
	public:
    static const Monte::ENSEMBLE ensemble;
    typedef ChargeNeutralGrandCanonicalEvent EventType;
    typedef GrandCanonicalConditions CondType;
    typedef GrandCanonicalSettings SettingsType;

	ChargeNeutralGrandCanonical(PrimClex &primclex, const SettingsType &settings, Log &_log);


    /// \brief Return current conditions
    const CondType &conditions() const;

    ///Keeps track of what sites can change to what
    const SiteExchanger m_site_swaps;

	/// \brief Set conditions and clear previously collected data
    void set_conditions(const CondType &new_conditions);

	/// This function needs to do all the math for energy and correlation deltas and store
	/// the results inside the containers hosted by event.
	void _update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<int,int> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const;

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    const EventType &propose();
    
	/// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
    bool check(const EventType &event);
    
	/// \brief Accept proposed event. Change configuration accordingly and update energies etc.
    void accept(const EventType &event);

    /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
    void reject(const EventType &event);
    
	/// \brief Write results to files
    void write_results(Index cond_index) const;

    /// \brief Formation energy, normalized per primitive cell
    const double &formation_energy() const {
      return *m_formation_energy;
    }

    /// \brief Potential energy, normalized per primitive cell
    const double &potential_energy() const {
      return *m_potential_energy;
    }

    /// \brief Correlations, normalized per primitive cell
    const Eigen::VectorXd &corr() const {
      return *m_corr;
    }

    /// \brief Number of atoms of each type, normalized per primitive cell
    const Eigen::VectorXd &comp_n() const {
      return *m_comp_n;
    }

    /// \brief Get potential energy
    double potential_energy(const Configuration &config) const;
    
    private:

    /// \brief Formation energy, normalized per primitive cell
    double &_formation_energy() {
      return *m_formation_energy;
    }

    /// \brief Potential energy, normalized per primitive cell
    double &_potential_energy() {
      return *m_potential_energy;
    }

    /// \brief Correlations, normalized per primitive cell
    Eigen::VectorXd &_corr() {
      return *m_corr;
    }

    /// \brief Number of atoms of each type, normalized per primitive cell
    Eigen::VectorXd &_comp_n() {
      return *m_comp_n;
    }

    Clexulator &_clexulator() const {
      return m_formation_energy_clex.clexulator();
    }

    
   /// Conditions (T, mu). Initially determined by m_settings, but can be changed halfway through the run
    GrandCanonicalConditions m_condition;

    /// Holds Clexulator and ECI references
    Clex m_formation_energy_clex;

    /// If true, calculate all correlations; if false, calculate correlations with non-zero eci
    bool m_all_correlations;

    /// Event to propose, check, accept/reject:
    EventType m_event;

    /// \brief If the supercell is large enough, calculate delta correlations directly
    bool m_use_deltas;


    // ---- Pointers to properties for faster access

    /// \brief Formation energy, normalized per primitive cell
    double *m_formation_energy;

    /// \brief Potential energy, normalized per primitive cell
    double *m_potential_energy;

    /// \brief Correlations, normalized per primitive cell
    Eigen::VectorXd *m_corr;

    /// \brief Number of atoms of each type, normalized per primitive cell
    Eigen::VectorXd *m_comp_n;
    
};




}
#endif
