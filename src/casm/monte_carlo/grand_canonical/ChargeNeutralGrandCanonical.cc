#include "casm/monte_carlo/grand_canonical/ChargeNeutralGrandCanonical.hh"

namespace CASM {
    // <- Jerry: same as GrandCanonical.cc
    ChargeNeutralGrandCanonical::ChargeNeutralGrandCanonical(PrimClex &primclex, const SettingsType &settings, Log &_log):
    MonteCarlo(primclex, settings, log),
    m_site_swaps(supercell()),
    m_formation_energy_clex(primclex, settings.formation_energy(primclex)),
    m_all_correlations(settings.all_correlations()),
    m_event(primclex.composition_axes().components().size(), _clexulator().corr_size()) {
        const auto &desc = m_formation_energy_clex.desc();

        // set the SuperNeighborList...
        set_nlist();

        // If the simulation is big enough, use delta cluster functions;
        // else, calculate all cluster functions
        m_use_deltas = !nlist().overlaps();        

        _log().construct("Charge Neutral Grand Canonical Monte Carlo");
        _log() << "project: " << this->primclex().get_path() << "\n";
        _log() << "formation_energy cluster expansion: " << desc.name << "\n";
        _log() << std::setw(16) << "property: " << desc.property << "\n";
        _log() << std::setw(16) << "calctype: " << desc.calctype << "\n";
        _log() << std::setw(16) << "ref: " << desc.ref << "\n";
        _log() << std::setw(16) << "bset: " << desc.bset << "\n";
        _log() << std::setw(16) << "eci: " << desc.eci << "\n";
        _log() << "supercell: \n" << supercell().get_transf_mat() << "\n";
        _log() << "use_deltas: " << std::boolalpha << m_use_deltas << "\n";
        _log() << "\nSampling: \n";
        _log() << std::setw(24) << "quantity" << std::setw(24) << "requested_precision" << "\n";
        for(auto it = samplers().begin(); it != samplers().end(); ++it) {
          _log() << std::setw(24) << it->first;
          if(it->second->must_converge()) {
            _log() << std::setw(24) << it->second->requested_precision() << std::endl;
          }
          else {
            _log() << std::setw(24) << "none" << std::endl;
          }
        }
        _log() << "\nautomatic convergence mode?: " << std::boolalpha << must_converge() << std::endl;
        _log() << std::endl;
    }

    /// \brief Return current conditions 
    /// <- Jerry: same as GrandCanonical.cc
    const ChargeNeutralGrandCanonical::CondType &conditions() const{
        return m_condition;
    }
    
	/// \brief Set conditions and clear previously collected data 
    /// <- Jerry: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::set_conditions(const CondType &new_conditions){
        _log().set("Conditions");
        _log() << new_conditions << std::endl << std::endl;

        m_condition = new_conditions;

        clear_samples();
        _update_properties();

        return;        
    }

	/// This function needs to do all the math for energy and correlation deltas and store
	/// the results inside the containers hosted by event. 
    /// <- Jerry: This is different from GrandCanonical.cc, not finished
	void ChargeNeutralGrandCanonical::_update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<int,int> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const{
        // ---- set OccMod --------------

        event.occupational_change().set(mutating_site, sublat, new_occupant);

        // ---- set dspecies --------------

        for(int i = 0; i < event.dN().size(); ++i) {
          event.set_dN(i, 0);
        }
        Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
        Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];
        event.set_dN(curr_species, -1);
        event.set_dN(new_species, 1);


        // ---- set dcorr --------------

        _set_dCorr(event, mutating_site, sublat, current_occupant, new_occupant, m_use_deltas, m_all_correlations);

        // ---- set dformation_energy --------------

        event.set_dEf(_eci() * event.dCorr().data());


        // ---- set dpotential_energy --------------

        event.set_dEpot(event.dEf() - m_condition.exchange_chem_pot(new_species, curr_species));
    }

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    /// <- Jerry: This is different from GrandCanonical.cc, under construction......
    const ChargeNeutralGrandCanonical::EventType &propose(){

        GrandCanonicalEvent event_1 = GrandCanonical::propose();

        if (event_1.occupational_change().sublat < 10){
            
        // Randomly pick a site that's allowed more than one occupant
        Index random_variable_site = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

        // Determine what that site's linear index is and what the sublattice index is
        Index mutating_site = m_site_swaps.variable_sites()[random_variable_site];
        Index sublat = m_site_swaps.sublat()[random_variable_site];

        // Determine the current occupant of the mutating site
        int current_occupant = configdof().occ(mutating_site);

        // Randomly pick a new occupant for the mutating site
        const std::vector<int> &possible_mutation = m_site_swaps.possible_swap()[sublat][current_occupant];
        int new_occupant = possible_mutation[_mtrand().randInt(possible_mutation.size() - 1)];

        if(debug()) {
          const auto &site_occ = primclex().get_prim().basis[sublat].site_occupant();
          _log().custom("Propose event");

          _log()  << "  Mutating site (linear index): " << mutating_site << "\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site) << "\n"
                  << "  Current occupant: " << current_occupant << " (" << site_occ[current_occupant].name << ")\n"
                  << "  Proposed occupant: " << new_occupant << " (" << site_occ[new_occupant].name << ")\n\n"
                  << "  beta: " << m_condition.beta() << "\n"
                  << "  T: " << m_condition.temperature() << std::endl;
        }

        // Update delta properties in m_event
        _update_deltas(m_event, mutating_site, sublat, current_occupant, new_occupant);

        if(debug()) {

          auto origin = primclex().composition_axes().origin();
          auto exchange_chem_pot = m_condition.exchange_chem_pot();
          auto param_chem_pot = m_condition.param_chem_pot();
          auto Mpinv = primclex().composition_axes().dparam_dmol();
          // auto V = supercell().volume();
          Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
          Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];

          _log() << "  components: " << jsonParser(primclex().composition_axes().components()) << "\n"
                 << "  d(N): " << m_event.dN().transpose() << "\n"
                 << "    dx_dn: \n" << Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn: \n" << param_chem_pot.transpose()*Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn * dN: " << param_chem_pot.transpose()*Mpinv *m_event.dN().cast<double>() << "\n"
                 << "  d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species, curr_species) << "\n"
                 << "  d(Ef): " << m_event.dEf() << "\n"
                 << "  d(Epot): " << m_event.dEf() - exchange_chem_pot(new_species, curr_species) << "\n"
                 << std::endl;


        }

        return m_event;

    }
    
	/// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
    bool ChargeNeutralGrandCanonical::check(const EventType &event){

    }
    
	/// \brief Accept proposed event. Change configuration accordingly and update energies etc.
    /// <- Jerry: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::accept(const EventType &event){
        if(debug()) {
          _log().custom("Accept Event");
          _log() << std::endl;
        }

        // First apply changes to configuration (just a single occupant change)
        _configdof().occ(event.occupational_change().site_index()) = event.occupational_change().to_value();

        // Next update all properties that changed from the event
        _formation_energy() += event.dEf() / supercell().volume();
        _potential_energy() += event.dEpot() / supercell().volume();
        _corr() += event.dCorr() / supercell().volume();
        _comp_n() += event.dN().cast<double>() / supercell().volume();

        return;
    }

    /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
    /// <- Jerry: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::reject(const EventType &event){
        if(debug()) {
          _log().custom("Reject Event");
          _log() << std::endl;
        }
        return;
    }
    
	/// \brief Write results to files
    /// <- Jerry: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::write_results(Index cond_index) const{
        CASM::write_results(settings(), *this, _log());
        write_conditions_json(settings(), *this, cond_index, _log());
        write_observations(settings(), *this, cond_index, _log());
        write_trajectory(settings(), *this, cond_index, _log());
        //write_pos_trajectory(settings(), *this, cond_index);
    }
}