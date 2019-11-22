#include "casm/monte_carlo/grand_canonical/ChargeNeutralGrandCanonical.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {
    const Monte::ENSEMBLE ChargeNeutralGrandCanonical::ensemble = Monte::ENSEMBLE::GrandCanonical; //Zeyu : not sure if this is correct??

    // <- Zeyu: same as GrandCanonical.cc
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
    /// <- Zeyu: same as GrandCanonical.cc
    const ChargeNeutralGrandCanonical::CondType &conditions() const{
        return m_condition;
    }
    
	/// \brief Set conditions and clear previously collected data 
    /// <- Zeyu: same as GrandCanonical.cc
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
    /// <- Zeyu: This is different from GrandCanonical.cc, not finished
    /// do this site by site and then calculate total dEpot and store in ChargeNeutralGrandCanonicalEvent
    /// and use it to for check()
	void ChargeNeutralGrandCanonical::_update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<int,int> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const{
        // Site 1
        // ---- set OccMod --------------
        event.occupational_change().first.set(mutating_sites.first, sublats.first, new_occs.first);

        // ---- set dspecies --------------

        for(int i = 0; i < event.dN().first.size(); ++i) {
          event.set_dN1(i, 0);
        }
        Index curr_species = m_site_swaps.sublat_to_mol()[sublats.first][curr_occs.first];
        Index new_species = m_site_swaps.sublat_to_mol()[sublats.first][new_occs.first];
        event.set_dN1(curr_species, -1);
        event.set_dN1(new_species, 1);

        // ---- set dcorr --------------

        _set_dCorr(event, mutating_sites.first, sublats.first, current_occupants.first, new_occupants.first, m_use_delta, m_all_correlation, 1); // Zeyu: Shall we rewrite _set_dCorr?

        // ---- set dformation_energy --------------

        event.set_dEf1(_eci() * event.dCorr().data());


        // ---- set dpotential_energy --------------
        double dEpot_1 = event.dEf().first - m_condition.exchange_chem_pot(new_species, curr_species);
        event.set_dEpot1(dEpot_1);

        // // Site 1 modification finished, update configuration ....
        _configdof().occ(event.occupational_change().first.site_index()) = event.occupational_change().first.to_value();

        // Site 2
        // ---- set OccMod --------------

        event.occupational_change().second.set(mutating_sites.second, sublats.second, new_occs.second);

        // ---- set dspecies --------------

        for(int i = 0; i < event.dN().first.size(); ++i) {
          event.set_dN2(i, 0);
        }
        Index curr_species = m_site_swaps.sublat_to_mol()[sublats.second][curr_occs.second];
        Index new_species = m_site_swaps.sublat_to_mol()[sublat.second][new_occs.second];
        event.set_dN2(curr_species, -1);
        event.set_dN2(new_species, 1);

        // ---- set dcorr --------------

        _set_dCorr(event, mutating_sites.second, sublats.second, current_occupants.second, new_occupants.second, m_use_delta, m_all_correlation, 2);
       
        // ---- set dformation_energy --------------

        event.set_dEf2(_eci() * event.dCorr().second.data());

        // ---- set dpotential_energy --------------
        double dEpot_2 = event.dEf().second - m_condition.exchange_chem_pot(new_species, curr_species);
        event.set_dEpot2(dEpot_2);
        
        // Calculate dEpot after two swaps
        event.set_dEpot_swapped_twice(dEpot_1+dEpot_2)

        // Zeyu: after get dEpot_swapped_twice, change configuration back to origin....
        _configdof().occ(event.occupational_change().first.site_index()) = event.occupational_change().first.from_value();

        
    }
  /// \brief Calculate delta correlations for an event

  void ChargeNeutralGrandCanonical::_set_dCorr(ChargeNeutralGrandCanonicalEvent &event,
                                  Index mutating_site,
                                  int sublat,
                                  int current_occupant,
                                  int new_occupant,
                                  bool use_deltas,
                                  bool all_correlations, int this_site) const {

    // uses _clexulator(), nlist(), _configdof()

    // Point the Clexulator to the right neighborhood and right ConfigDoF
    _clexulator().set_config_occ(_configdof().occupation().begin());
    _clexulator().set_nlist(nlist().sites(nlist().unitcell_index(mutating_site)).data());

    if(use_deltas) {
      // Calculate the change in correlations due to this event
      if(all_correlations) {
        if (this_site == 1) {
        _clexulator().calc_delta_point_corr(sublat,
                                            current_occupant,
                                            new_occupant,
                                            event.dCorr().first.data());
                                            }
        if (this_site == 2) {
        _clexulator().calc_delta_point_corr(sublat,
                                            current_occupant,
                                            new_occupant,
                                            event.dCorr().second.data());
                                            }
      }
      else {
        auto begin = _eci().index().data();
        auto end = begin + _eci().index().size();
        if (this_site == 1){
        _clexulator().calc_restricted_delta_point_corr(sublat,
                                                       current_occupant,
                                                       new_occupant,
                                                       event.dCorr().first.data(),
                                                       begin,
                                                       end);
        }
        if (this_site == 2){
        _clexulator().calc_restricted_delta_point_corr(sublat,
                                                       current_occupant,
                                                       new_occupant,
                                                       event.dCorr().second.data(),
                                                       begin,
                                                       end);          
        }
      }
    }
    else {
      if (this_site == 1){
      Eigen::VectorXd before { Eigen::VectorXd::Zero(event.dCorr().first.size()) };
      Eigen::VectorXd after { Eigen::VectorXd::Zero(event.dCorr().first.size()) };
      }
      if (this_site == 2){
      Eigen::VectorXd before { Eigen::VectorXd::Zero(event.dCorr().second.size()) };
      Eigen::VectorXd after { Eigen::VectorXd::Zero(event.dCorr().second.size()) };        
      }

      // Calculate the change in points correlations due to this event
      if(all_correlations) {

        // Calculate before
        _clexulator().calc_point_corr(sublat, before.data());

        // Apply change
        _configdof().occ(mutating_site) = new_occupant;

        // Calculate after
        _clexulator().calc_point_corr(sublat, after.data());
      }
      else {
        auto begin = _eci().index().data();
        auto end = begin + _eci().index().size();

        // Calculate before
        _clexulator().calc_restricted_point_corr(sublat, before.data(), begin, end);

        // Apply change
        _configdof().occ(mutating_site) = new_occupant;

        // Calculate after
        _clexulator().calc_restricted_point_corr(sublat, after.data(), begin, end);

      }

      // Calculate the change in correlations due to this event
      if (this_site == 1){
      event.dCorr().first = after - before;
      }
      if (this_site == 2){
      event.dCorr().second = after - before;  
      }

      // Unapply changes
      _configdof().occ(mutating_site) = current_occupant;
    }

    if(debug()) {
      if (this_site == 1){
      _print_correlations(event.dCorr().first, "delta correlations", "dCorr", all_correlations);
      }
      if (this_site == 2){
      _print_correlations(event.dCorr().second, "delta correlations", "dCorr", all_correlations);
      }
    }
  }

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    /// <- Zeyu: This is different from GrandCanonical.cc, under construction......
    const ChargeNeutralGrandCanonical::EventType &propose(){

        // Zeyu: 2 mutations at the same time; pick one Na/Va and one Si/P with the same occupancy and flip them together
        do{
          // Randomly pick a site that's allowed more than one occupant
          Index random_variable_site_1 = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);
          Index random_variable_site_2 = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

        // Determine what that site's linear index is and what the sublattice index is
          Index mutating_site_1 = m_site_swaps.variable_sites()[random_variable_site_1];
          Index mutating_site_2 = m_site_swaps.variable_sites()[random_variable_site_2];

          Index sublat_1 = m_site_swaps.sublat()[random_variable_site_1];
          Index sublat_2 = m_site_swaps.sublat()[random_variable_site_2];
      
          // Determine the current occupant of the mutating site
          int current_occupant_1 = configdof().occ(mutating_site_1);
          int current_occupant_2 = configdof().occ(mutating_site_2);
        }
        while (((sublat_1 <= 8 && sublat_2 > 8) || (sublat_1 > 8 && sublat_2 <= 8)) && (current_occupant_1 == current_occupant_2));

        // Randomly pick a new occupant for the mutating site
        const std::vector<int> &possible_mutation_1 = m_site_swaps.possible_swap()[sublat_1][current_occupant_1];
        int new_occupant_1 = possible_mutation_1[_mtrand().randInt(possible_mutation_1.size() - 1)];
        const std::vector<int> &possible_mutation_2 = m_site_swaps.possible_swap()[sublat_2][current_occupant_2];
        int new_occupant_2 = possible_mutation_2[_mtrand().randInt(possible_mutation_2.size() - 1)];

        if(debug()) {
          const auto &site_occ = primclex().get_prim().basis[sublat].site_occupant();
          _log().custom("Propose charge neutral grand canonical event");

          _log()  << "  Mutating site 1 (linear index): " << mutating_site_1 << "\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_1) << "\n"
                  << "  Current occupant: " << current_occupant_1 << " (" << site_occ[current_occupant_1].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_1 << " (" << site_occ[new_occupant_1].name << ")\n\n"

                  << "  Mutating site 2 (linear index): " << mutating_site_2 << "\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_2) << "\n"
                  << "  Current occupant: " << current_occupant_2 << " (" << site_occ[current_occupant_2].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_2 << " (" << site_occ[new_occupant_2].name << ")\n\n"

                  << "  beta: " << m_condition.beta() << "\n"
                  << "  T: " << m_condition.temperature() << std::endl;
        }

        // Zeyu: creating pairs
        std::pair<Index,Index> mutating_sites (mutating_site_1,mutating_site_2);
        std::pair<Index,Index> sublats (sublat_1,sublat_2);
        std::pair<int,int> current_occupants (current_occupant_1,current_occupant_1);
        std::pair<int,int> new_occupants (new_occupant_1,new_occupant_2)

        // Update delta properties in m_event
        // Zeyu: Pairs are passing into _update_deltas()
        _update_deltas(m_event, mutating_sites, sublats, current_occupants, new_occupants);

        if(debug()) {

          auto origin = primclex().composition_axes().origin();
          auto exchange_chem_pot = m_condition.exchange_chem_pot();
          auto param_chem_pot = m_condition.param_chem_pot();
          auto Mpinv = primclex().composition_axes().dparam_dmol();
          // auto V = supercell().volume();
          Index curr_species_1 = m_site_swaps.sublat_to_mol()[sublat_1][current_occupant_1];
          Index new_species_1 = m_site_swaps.sublat_to_mol()[sublat_1][new_occupant_1];
          Index curr_species_2 = m_site_swaps.sublat_to_mol()[sublat_2][current_occupant_2];
          Index new_species_2 = m_site_swaps.sublat_to_mol()[sublat_2][new_occupant_2];

          _log() << "  components: " << jsonParser(primclex().composition_axes().components()) << "\n"
                 << "  d(N): " << m_event.dN().transpose() << "\n"
                 << "    dx_dn: \n" << Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn: \n" << param_chem_pot.transpose()*Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn * dN: " << param_chem_pot.transpose()*Mpinv *m_event.dN().cast<double>() << "\n"
                 << "Swap step 1: d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species_1, curr_species_1) << "\n"
                 << "Swap step 2: d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species_2, curr_species_2) << "\n"

                 << "  d(Ef): " << m_event.dEf() << "\n"
                 << "  d(Epot): " << m_event.dEf() - exchange_chem_pot(new_species_2, curr_species_2) << "\n"
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
      double prob = exp(-event.dEpot_swapped_twice() * m_condition.beta()); 

      if(debug()) {
        _log().custom("Check event");
        _log() << "Probability to accept: " << prob << "\n"
               << "Random number: " << rand << "\n" << std::endl;
      }

      return rand < prob;
    }
    
	/// \brief Accept proposed event. Change configuration accordingly and update energies etc.
    /// <- Zeyu: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::accept(const EventType &event){
        if(debug()) {
          _log().custom("Accept Event");
          _log() << std::endl;
        }

        // First apply changes to configuration (just a single occupant change)
        _configdof().occ(event.occupational_change().first.site_index()) = event.occupational_change().first.to_value();
        _configdof().occ(event.occupational_change().second.site_index()) = event.occupational_change().second.to_value();

        // Next update all properties that changed from the event // Zeyu: update twice, the volume does not change throughout the simulation
        _formation_energy() += event.dEf().first / supercell().volume();
        _formation_energy() += event.dEf().second / supercell().volume();
        _potential_energy() += event.dEpot().first / supercell().volume();
        _potential_energy() += event.dEpot().second / supercell().volume();
        _corr() += event.dCorr().first / supercell().volume();
        _corr() += event.dCorr().second / supercell().volume();
        _comp_n() += event.dN().first.cast<double>() / supercell().volume();
        _comp_n() += event.dN().second.cast<double>() / supercell().volume();

        return;
    }

    /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
    /// <- Zeyu: same as GrandCanonical.cc
    void ChargeNeutralGrandCanonical::reject(const EventType &event){
        if(debug()) {
          _log().custom("Reject Event");
          _log() << std::endl;
        }
        return;
    }
    
	/// \brief Write results to files
    void ChargeNeutralGrandCanonical::write_results(Index cond_index) const{
        CASM::write_results(settings(), *this, _log());
        write_conditions_json(settings(), *this, cond_index, _log());
        write_observations(settings(), *this, cond_index, _log());
        write_trajectory(settings(), *this, cond_index, _log());
        //write_pos_trajectory(settings(), *this, cond_index);
    }
}