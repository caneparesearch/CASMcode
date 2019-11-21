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
    /// do this site by site and then calculate total dEpot and store in ChargeNeutralGrandCanonicalEvent
    /// and use it to for check()
	void ChargeNeutralGrandCanonical::_update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<int,int> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const{
        // Site Na/Va
        // ---- set OccMod --------------
        event.occupational_change().set(mutating_sites.first, sublats.first, new_occs.first);

        // ---- set dspecies --------------

        for(int i = 0; i < event.dN().size(); ++i) {
          event.set_dN(i, 0);
        }
        Index curr_species = m_site_swaps.sublat_to_mol()[sublats.first][curr_occs.first];
        Index new_species = m_site_swaps.sublat_to_mol()[sublats.first][new_occs.first];
        event.set_dN(curr_species, -1);
        event.set_dN(new_species, 1);

        // ---- set dcorr --------------

        _set_dCorr(event, mutating_sites.first, sublats.first, current_occupants.first, new_occupants.first, m_use_delta, m_all_correlation);

        // ---- set dformation_energy --------------

        event.set_dEf(_eci() * event.dCorr().data());


        // ---- set dpotential_energy --------------
        dEpot_NaVa = event.dEf() - m_condition.exchange_chem_pot(new_species, curr_species);
        event.set_dEpot(dEpot_NaVa);

        // Site Si/P
        // ---- set OccMod --------------

        event.occupational_change().set(mutating_sites.second, sublats.second, new_occs.second);

        // ---- set dspecies --------------

        for(int i = 0; i < event.dN().size(); ++i) {
          event.set_dN(i, 0);
        }
        Index curr_species = m_site_swaps.sublat_to_mol()[sublats.second][curr_occs.second];
        Index new_species = m_site_swaps.sublat_to_mol()[sublat.second][new_occs.second];
        event.set_dN(curr_species, -1);
        event.set_dN(new_species, 1);

        // ---- set dcorr --------------

        _set_dCorr(event, mutating_sites.second, sublats.second, current_occupants.second, new_occupants.second, m_use_delta, m_all_correlation);
       
        // ---- set dformation_energy --------------

        event.set_dEf(_eci() * event.dCorr().data());

        // ---- set dpotential_energy --------------
        dEpot_SiP = event.dEf() - m_condition.exchange_chem_pot(new_species, curr_species);
        event.set_dEpot(dEpot_SiP);
        
        // Calculate dEpot after two swaps
        event.set_dEpot_swapped_twice(dEpot_NaVa+dEpot_SiP)
        
    }

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    /// <- Jerry: This is different from GrandCanonical.cc, under construction......
    const ChargeNeutralGrandCanonical::EventType &propose(){

        // Jerry: 2 mutations at the same time; pick one Na/Va and one Si/P with the same occupancy and flip them together
        do{
          // Randomly pick a site that's allowed more than one occupant
          Index random_variable_site_NaVa = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);
          Index random_variable_site_SiP = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

        // Determine what that site's linear index is and what the sublattice index is
          Index mutating_site_NaVa = m_site_swaps.variable_sites()[random_variable_site_NaVa];
          Index mutating_site_SiP = m_site_swaps.variable_sites()[random_variable_site_SiP];

          Index sublat_NaVa = m_site_swaps.sublat()[random_variable_site_NaVa];
          Index sublat_SiP = m_site_swaps.sublat()[random_variable_site_SiP];
      
          // Determine the current occupant of the mutating site
          int current_occupant_NaVa = configdof().occ(mutating_site_NaVa);
          int current_occupant_SiP = configdof().occ(mutating_site_SiP);
        }
        while (sublat_NaVa < 10 && sublat_SiP > 10 && current_occupant_NaVa == current_occupant_SiP);

        // Randomly pick a new occupant for the mutating site
        const std::vector<int> &possible_mutation_NaVa = m_site_swaps.possible_swap()[sublat_NaVa][current_occupant_NaVa];
        int new_occupant_NaVa = possible_mutation_NaVa[_mtrand().randInt(possible_mutation_NaVa.size() - 1)];
        const std::vector<int> &possible_mutation_SiP = m_site_swaps.possible_swap()[sublat_SiP][current_occupant_SiP];
        int new_occupant_SiP = possible_mutation_SiP[_mtrand().randInt(possible_mutation_SiP.size() - 1)];

        if(debug()) {
          const auto &site_occ = primclex().get_prim().basis[sublat].site_occupant();
          _log().custom("Propose charge neutral grand canonical event");

          _log()  << "  Mutating site Na/Va (linear index): " << mutating_site_NaVa << "\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_NaVa) << "\n"
                  << "  Current occupant: " << current_occupant_NaVa << " (" << site_occ[current_occupant_NaVa].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_NaVa << " (" << site_occ[new_occupant_NaVa].name << ")\n\n"

                  << "  Mutating site Si/P (linear index): " << mutating_site_SiP << "\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_SiP) << "\n"
                  << "  Current occupant: " << current_occupant_SiP << " (" << site_occ[current_occupant_SiP].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_SiP << " (" << site_occ[new_occupant_SiP].name << ")\n\n"

                  << "  beta: " << m_condition.beta() << "\n"
                  << "  T: " << m_condition.temperature() << std::endl;
        }

        // creating pairs
        std::pair<Index,Index> mutating_sites (mutating_site_NaVa,mutating_site_SiP);
        std::pair<Index,Index> sublats (sublat_NaVa,sublat_SiP);
        std::pair<int,int> current_occupants (current_occupant_NaVa,current_occupant_NaVa);
        std::pair<int,int> new_occupants (new_occupant_NaVa,new_occupant_SiP)

        // Update delta properties in m_event
        // Pairs are passing into _update_deltas()
        _update_deltas(m_event, mutating_sites, sublats, current_occupants, new_occupants);

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
      if(event.dEpot_swapped_twice() < 0.0) {

        if(debug()) {
          _log().custom("Check event");
          _log() << "Probability to accept: 1.0\n" << std::endl;
        }
        return true;
      }

      double rand = _mtrand().rand53();
      double prob = exp(-event.dEpot_swapped_twice() * m_condition.beta()); // is this beta correct????

      if(debug()) {
        _log().custom("Check event");
        _log() << "Probability to accept: " << prob << "\n"
               << "Random number: " << rand << "\n" << std::endl;
      }

      return rand < prob;
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