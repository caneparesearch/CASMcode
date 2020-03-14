#include "casm/monte_carlo/grand_canonical/ChargeNeutralGrandCanonical.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/grand_canonical/ChargeNeutralGrandCanonicalIO.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"


namespace CASM {
    const Monte::ENSEMBLE ChargeNeutralGrandCanonical::ensemble = Monte::ENSEMBLE::ChargeNeutralGrandCanonical; //Zeyu : not sure if this is correct??

    ChargeNeutralGrandCanonical::ChargeNeutralGrandCanonical(PrimClex &primclex, const SettingsType &settings, Log &log):
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

  /// \brief Return number of steps per pass. Equals number of sites with variable occupation.
  Index ChargeNeutralGrandCanonical::steps_per_pass() const {
    return m_site_swaps.variable_sites().size();
  }

    /// \brief Return current conditions 
    /// <- Zeyu: same as GrandCanonical.cc
    const ChargeNeutralGrandCanonical::CondType &ChargeNeutralGrandCanonical::conditions() const{
        return m_condition;
    }
  
  
  /// \brief Set conditions and clear previously collected data
  void ChargeNeutralGrandCanonical::set_conditions(const GrandCanonicalConditions &new_conditions) {
    _log().set("Conditions");
    _log() << new_conditions << std::endl << std::endl;

    m_condition = new_conditions;

    clear_samples();
    _update_properties();

    return;
  }

  /// \brief Set configdof and clear previously collected data
  void ChargeNeutralGrandCanonical::set_configdof(const ConfigDoF &configdof, const std::string &msg) {
    _log().set("DoF");
    if(!msg.empty()) {
      _log() << msg << "\n";
    }
    _log() << std::endl;

    reset(configdof);
    _update_properties();
  }


  /// \brief Set configdof and conditions and clear previously collected data
  ///
  /// \returns Specified ConfigDoF and configname (or configdof path)
  ///
  std::pair<ConfigDoF, std::string> ChargeNeutralGrandCanonical::set_state(
    const GrandCanonicalConditions &new_conditions,
    const GrandCanonicalSettings &settings) {

    _log().set("Conditions");
    _log() << new_conditions << std::endl;

    m_condition = new_conditions;

    ConfigDoF configdof;
    std::string configname;

    if(settings.is_motif_configname()) {

      configname = settings.motif_configname();

      if(configname == "default") {
        configdof = _default_motif();
      }
      else if(configname == "auto") {
        std::tie(configdof, configname) = _auto_motif(new_conditions);
      }
      else if(configname == "restricted_auto") {
        std::tie(configdof, configname) = _restricted_auto_motif(new_conditions);
      }
      else {
        configdof = _configname_motif(configname);
      }

    }
    else if(settings.is_motif_configdof()) {
      _log().set("DoF");
      _log() << "motif configdof: " << settings.motif_configdof_path() << "\n";
      _log() << "using configdof: " << settings.motif_configdof_path() << "\n" << std::endl;
      configdof = settings.motif_configdof();
      configname = settings.motif_configdof_path().string();
    }
    else {
      throw std::runtime_error("Error: Must specify motif \"configname\" or \"configdof\"");
    }

    reset(configdof);
    _update_properties();

    return std::make_pair(configdof, configname);
  }

  /// \brief Set configdof and conditions and clear previously collected data
  void ChargeNeutralGrandCanonical::set_state(const GrandCanonicalConditions &new_conditions,
                                 const ConfigDoF &configdof,
                                 const std::string &msg) {
    _log().set("Conditions");
    _log() << new_conditions << std::endl << std::endl;

    m_condition = new_conditions;

    _log().set("DoF");
    if(!msg.empty()) {
      _log() << msg << "\n";
    }
    _log() << std::endl;

    reset(configdof);
    _update_properties();

    return;
  }

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    /// <- Zeyu: This is different from GrandCanonical.cc, under construction......
    const ChargeNeutralGrandCanonical::EventType &ChargeNeutralGrandCanonical::propose(){
        Index random_variable_site_1,random_variable_site_2;
        Index mutating_site_1,mutating_site_2;
        Index sublat_1,sublat_2;
        int current_occupant_1,current_occupant_2;
	      int n_Na = 8;
        
        // Zeyu: 2 mutations at the same time; pick one Na/Va and one Si/P with the same occupancy and flip them together
        do{
          // Randomly pick a site that's allowed more than one occupant
          random_variable_site_1 = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);
          random_variable_site_2 = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

        // Determine what that site's linear index is and what the sublattice index is
          mutating_site_1 = m_site_swaps.variable_sites()[random_variable_site_1];
          mutating_site_2 = m_site_swaps.variable_sites()[random_variable_site_2];

          sublat_1 = m_site_swaps.sublat()[random_variable_site_1];
          sublat_2 = m_site_swaps.sublat()[random_variable_site_2];
      
          // Determine the current occupant of the mutating site
          current_occupant_1 = configdof().occ(mutating_site_1);
          current_occupant_2 = configdof().occ(mutating_site_2);
        }
        while (!(((sublat_1 <= n_Na && sublat_2 > n_Na) || (sublat_1 > n_Na && sublat_2 <= n_Na)) && (current_occupant_1 == current_occupant_2)));

        // Randomly pick a new occupant for the mutating site
        const std::vector<int> &possible_mutation_1 = m_site_swaps.possible_swap()[sublat_1][current_occupant_1];
        int new_occupant_1 = possible_mutation_1[_mtrand().randInt(possible_mutation_1.size() - 1)];
        const std::vector<int> &possible_mutation_2 = m_site_swaps.possible_swap()[sublat_2][current_occupant_2];
        int new_occupant_2 = possible_mutation_2[_mtrand().randInt(possible_mutation_2.size() - 1)];

        if(debug()) {
          const auto &site_occ_1 = primclex().get_prim().basis[sublat_1].site_occupant();
          const auto &site_occ_2 = primclex().get_prim().basis[sublat_2].site_occupant();
          _log().custom("Propose charge neutral grand canonical event");

          _log()  << "  Mutating site 1 (linear index): " << mutating_site_1 << "\n"
                  << "  Sublattice: "<< sublat_1<<"\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_1) << "\n"
                  << "  Current occupant: " << current_occupant_1 << " (" << site_occ_1[current_occupant_1].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_1 << " (" << site_occ_1[new_occupant_1].name << ")\n\n"

                  << "  Mutating site 2 (linear index): " << mutating_site_2 << "\n"
                  << "  Sublattice: "<< sublat_2<<"\n"
                  << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site_2) << "\n"
                  << "  Current occupant: " << current_occupant_2 << " (" << site_occ_2[current_occupant_2].name << ")\n"
                  << "  Proposed occupant: " << new_occupant_2 << " (" << site_occ_2[new_occupant_2].name << ")\n\n"

                  << "  beta: " << m_condition.beta() << "\n"
                  << "  T: " << m_condition.temperature() << std::endl;
        }

        // Zeyu: creating pairs
        std::pair<Index,Index> mutating_sites (mutating_site_1,mutating_site_2);
        std::pair<Index,Index> sublats (sublat_1,sublat_2);
        std::pair<int,int> current_occupants (current_occupant_1,current_occupant_1);
        std::pair<int,int> new_occupants (new_occupant_1,new_occupant_2);

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
                 << "  d(N1): " << m_event.dN().first.transpose() << "\n"
                 << "  d(N2): " << m_event.dN().second.transpose() << "\n"
                 << "    dx_dn: \n" << Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn: \n" << param_chem_pot.transpose()*Mpinv << "\n"
                 << "    param_chem_pot.transpose() * dx_dn * dN: " << param_chem_pot.transpose()*Mpinv *m_event.dN().first.cast<double>() << "\n"
                 << "Swap step 1: d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species_1, curr_species_1) << "\n"


                 << "Swap step 2: d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species_2, curr_species_2) << "\n"

                 << "  d(Ef2): " << m_event.dEf().second << "\n"
                 << "  d(Epot2): " << m_event.dEf().second - exchange_chem_pot(new_species_2, curr_species_2) << "\n"
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

  /// \brief Calculate the single spin flip low temperature expansion of the grand canonical potential
  ///
  /// Returns low temperature expansion estimate of the grand canonical free energy.
  /// Works with the current ConfigDoF as groundstate.
  ///
  /// Quick derivation:
  /// Z: partition function
  /// boltz(x): exp(-x/kBT)
  /// \Omega: (E-SUM(chem_pot*comp_n))*N
  /// N: number of unit cells in supercell
  ///
  /// The partition function is
  /// Z=SUM(boltz(\Omega_s))    summing over all microstates s
  ///
  /// \Omega_s can be split into groundstate \Omega_0 and a delta energy D\Omega
  /// \Omega_s=\Omega_0+D\Omega_s
  /// Z=boltz(\Omega_0)*SUM(boltz(D\Omega_s))  summing over all microstates
  ///
  /// For low temperatures we can approximate Z by truncating the sum after microstates that
  /// only involve point defects and no defects
  /// Z=boltz(\Omega_0)*SUM(boltz(D\Omega_s))  summing over all states with only point defects or no defects
  ///
  /// The free energy is
  /// Phi=-kB*T*ln(Z)
  /// Phi=-kB*T*(-\Omega_0/kBT+ln(SUM(boltz(D\Omega_s))    Sum is over point defects and no defects (in which case D\Omega_s == 0)
  /// Phi=(\Omega_0-kB*T(ln(SUM(boltz(D\Omega_s)))))/N
  ///
  double ChargeNeutralGrandCanonical::lte_grand_canonical_free_energy() const {
  // Zeyu: Everytime it will select 2 sites, and evaluate their energies when flip them, so it's a 2 flip instead of single-flip LTE
  // it should be called LTE2 
  // Actually, I tried to do LTE on NASICON however I found Phi_LTE is almost the same as the potential_energy since
  // Phi = E-TS-ux and S -> 0 when T->0; Phi_lowT -> E-ux = potential_energy
    const SiteExchanger &site_exch = m_site_swaps;
    const ConfigDoF &config_dof = configdof();
    ChargeNeutralGrandCanonicalEvent event = m_event;

    double tol = 1e-12;

    auto less = [&](const double & A, const double & B) {
      return A < B - tol;
    };

    std::map<double, unsigned long, decltype(less)> hist(less);

    // no defect case
    hist[0.0] = 1;

    // double sum_exp = 0.0;

    //Loop over sites that can change occupants
    int n_Na = 8;
    for(Index exch_ind_1 = 0; exch_ind_1 < site_exch.variable_sites().size(); exch_ind_1++) {
      for(Index exch_ind_2 = 0; exch_ind_2 < site_exch.variable_sites().size(); exch_ind_2++) {

        //Transform exchanger index to ConfigDoF index
        Index mutating_site_1 = site_exch.variable_sites()[exch_ind_1];
        int sublat_1 = site_exch.sublat()[exch_ind_1];
        int current_occupant_1 = config_dof.occ(mutating_site_1);

        Index mutating_site_2 = site_exch.variable_sites()[exch_ind_2];
        int sublat_2 = site_exch.sublat()[exch_ind_2];
        int current_occupant_2 = config_dof.occ(mutating_site_2);

        if ((((sublat_1 <= n_Na && sublat_2 > n_Na) || (sublat_1 > n_Na && sublat_2 <= n_Na)) && (current_occupant_1 == current_occupant_2))){
          //Loop over possible occupants for site that can change
          const auto &possible_1 = site_exch.possible_swap()[sublat_1][current_occupant_1];
          const auto &possible_2 = site_exch.possible_swap()[sublat_2][current_occupant_2];
          for(auto new_occ_it_1 = possible_1.begin(); new_occ_it_1 != possible_1.end(); ++new_occ_it_1) {
            for(auto new_occ_it_2 = possible_2.begin(); new_occ_it_2 != possible_2.end(); ++new_occ_it_2) {

              // Zeyu: creating pairs
              std::pair<Index,Index> mutating_sites (mutating_site_1,mutating_site_2);
              std::pair<Index,Index> sublats (sublat_1,sublat_2);
              std::pair<int,int> current_occupants (current_occupant_1,current_occupant_1);
              std::pair<int,int> new_occ_its (*new_occ_it_1,*new_occ_it_2);

              _update_deltas(event, mutating_sites, sublats, current_occupants, new_occ_its);

              //save the result
              double dpot_nrg = event.dEpot().first+event.dEpot().second;
              if(dpot_nrg < 0.0) {
                Log &err_log = default_err_log();
                err_log.error<Log::standard>("Calculating low temperature expansion");
                err_log << "  Defect lowered the potential energy. Your motif configuration "
                        << "is not the 0K ground state.\n" << std::endl;
                throw std::runtime_error("Error calculating low temperature expansion. Not in the ground state.");
              }


              auto it = hist.find(dpot_nrg);
              if(it == hist.end()) {
                hist[dpot_nrg] = 1;
              }
              else {
                it->second++;
              }
            }
          }
        }
      }
    }

    _log().results("Ground state and point defect potential energy details");
    _log() << "T: " << m_condition.temperature() << std::endl;
    _log() << "kT: " << 1.0 / m_condition.beta() << std::endl;
    _log() << "Beta: " << m_condition.beta() << std::endl << std::endl;

    _log() << std::setw(16) << "N/unitcell" << " "
           << std::setw(16) << "dPE" << " "
           << std::setw(24) << "N*exp(-dPE/kT)" << " "
           << std::setw(16) << "dphi" << " "
           << std::setw(16) << "phi" << std::endl;

    double tsum = 0.0;
    double phi = 0.0;
    double phi_prev;
    for(auto it = hist.rbegin(); it != hist.rend(); ++it) {
      phi_prev = phi;
      // Zeyu: remove double counting
      tsum += it->second * exp(-(it->first) * m_condition.beta())/2;
      phi = std::log(tsum) / m_condition.beta() / supercell().volume();

      if(almost_equal(it->first, 0.0, tol)) {
        _log() << std::setw(16) << "(gs)" << " ";
      }
      else {
        _log() << std::setw(16) << std::setprecision(8) << (1.0 * it->second) / supercell().volume() << " ";
      }
      _log() << std::setw(16) << std::setprecision(8) << it->first << " "
             << std::setw(24) << std::setprecision(8) << it->second *exp(-it->first * m_condition.beta()) << " "
             << std::setw(16) << std::setprecision(8) << phi - phi_prev << " "
             << std::setw(16) << std::setprecision(8) << potential_energy() - phi << std::endl;

    }

    _log() << "phi_LTE(1): " << std::setprecision(12) << potential_energy() - phi << std::endl << std::endl;

    return potential_energy() - phi;

  }

	/// \brief Write results to files
    void ChargeNeutralGrandCanonical::write_results(Index cond_index) const{
        CASM::write_results(settings(), *this, _log());
        write_conditions_json(settings(), *this, cond_index, _log());
        write_observations(settings(), *this, cond_index, _log());
        write_trajectory(settings(), *this, cond_index, _log());
        //write_pos_trajectory(settings(), *this, cond_index);
    }

  /// \brief Get potential energy
  ///
  /// - if(&config == &this->config()) { return potential_energy(); }, else
  ///   calculate potential_energy = formation_energy - comp_x.dot(param_chem_pot)
  double ChargeNeutralGrandCanonical::potential_energy(const Configuration &config) const {
    //if(&config == &this->config()) { return potential_energy(); }

    auto corr = correlations(config, _clexulator());
    double formation_energy = _eci() * corr.data();
    auto comp_x = primclex().composition_axes().param_composition(CASM::comp_n(config));
    return formation_energy - comp_x.dot(m_condition.param_chem_pot());
  }

  /// \brief Calculate delta correlations for an event
  void ChargeNeutralGrandCanonical::_set_dCorr(ChargeNeutralGrandCanonicalEvent &event,
                                  Index mutating_site,
                                  int sublat,
                                  int current_occupant,
                                  int new_occupant,
                                  bool use_deltas,
                                  bool all_correlations) const {

    // uses _clexulator(), nlist(), _configdof()

    // Point the Clexulator to the right neighborhood and right ConfigDoF
    _clexulator().set_config_occ(_configdof().occupation().begin());
    _clexulator().set_nlist(nlist().sites(nlist().unitcell_index(mutating_site)).data());

    if(use_deltas) {
      // Calculate the change in correlations due to this event
      if(all_correlations) {
        if (!event.is_swapped()) {
        _clexulator().calc_delta_point_corr(sublat,
                                            current_occupant,
                                            new_occupant,
                                            event.dCorr().first.data());
                                            }
        if (event.is_swapped()) {
        _clexulator().calc_delta_point_corr(sublat,
                                            current_occupant,
                                            new_occupant,
                                            event.dCorr().second.data());
                                            }
      }
      else {
        auto begin = _eci().index().data();
        auto end = begin + _eci().index().size();
        if (!event.is_swapped()){
        _clexulator().calc_restricted_delta_point_corr(sublat,
                                                       current_occupant,
                                                       new_occupant,
                                                       event.dCorr().first.data(),
                                                       begin,
                                                       end);
        }
        if (event.is_swapped()){
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
      Eigen::VectorXd before;
      Eigen::VectorXd after;
      if (!event.is_swapped()){
      before = Eigen::VectorXd::Zero(event.dCorr().first.size()); //Zeyu: is this correct????
      after = Eigen::VectorXd::Zero(event.dCorr().first.size());
      }
      if (event.is_swapped()){
      before = Eigen::VectorXd::Zero(event.dCorr().second.size());
      after = Eigen::VectorXd::Zero(event.dCorr().second.size());        
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
      if (!event.is_swapped()){
      event.dCorr().first = after - before;
      }
      if (event.is_swapped()){
      event.dCorr().second = after - before;  
      }

      // Unapply changes
      _configdof().occ(mutating_site) = current_occupant;
    }

    if(debug()) {
      if (!event.is_swapped()){
      _print_correlations(event.dCorr().first, "delta correlations", "dCorr", all_correlations);
      }
      if (event.is_swapped()){
      _print_correlations(event.dCorr().second, "delta correlations", "dCorr", all_correlations);
      }
    }
  }

  /// \brief Print correlations to _log()
  void ChargeNeutralGrandCanonical::_print_correlations(const Eigen::VectorXd &corr,
                                           std::string title,
                                           std::string colheader,
                                           bool all_correlations) const {

    _log().calculate(title);
    _log() << std::setw(12) << "i"
           << std::setw(16) << "ECI"
           << std::setw(16) << colheader
           << std::endl;

    for(int i = 0; i < corr.size(); ++i) {

      double eci = 0.0;
      bool calculated = true;
      Index index = find_index(_eci().index(), i);
      if(index != _eci().index().size()) {
        eci = _eci().value()[index];
      }
      if(!all_correlations && index == _eci().index().size()) {
        calculated = false;
      }

      _log() << std::setw(12) << i
             << std::setw(16) << std::setprecision(8) << eci;
      if(calculated) {
        _log() << std::setw(16) << std::setprecision(8) << corr[i];
      }
      else {
        _log() << std::setw(16) << "unknown";
      }
      _log() << std::endl;

    }
    _log() << std::endl;
  }

	/// This function needs to do all the math for energy and correlation deltas and store
	/// the results inside the containers hosted by event. 
    /// <- Zeyu: This is different from GrandCanonical.cc, not finished
    /// do this site by site and then calculate total dEpot and store in ChargeNeutralGrandCanonicalEvent
    /// and use it to for check()
	void ChargeNeutralGrandCanonical::_update_deltas(EventType &event, 
						std::pair<Index,Index> &mutating_sites,
						std::pair<Index,Index> &sublats,
						std::pair<int,int> &curr_occs,
						std::pair<int,int> &new_occs) const{
        // reset the flag
        event.set_is_swapped(false);

        // Site 1
        // ---- set OccMod --------------
        event.occupational_change().first.set(mutating_sites.first, sublats.first, new_occs.first);

        // ---- set dspecies --------------
        for(int i = 0; i < event.dN().first.size(); ++i) {
          event.set_dN(i, 0);
        }
        Index curr_species_1 = m_site_swaps.sublat_to_mol()[sublats.first][curr_occs.first];
        Index new_species_1 = m_site_swaps.sublat_to_mol()[sublats.first][new_occs.first];
        event.set_dN(curr_species_1, -1);
        event.set_dN(new_species_1, 1);

        // ---- set dcorr --------------
        _set_dCorr(event, mutating_sites.first, sublats.first, curr_occs.first, new_occs.first, m_use_deltas, m_all_correlations); // Zeyu: Shall we rewrite _set_dCorr?

        // ---- set dformation_energy --------------
        event.set_dEf(_eci() * event.dCorr().first.data());

        // ---- set dpotential_energy --------------
        double dEpot_1 = event.dEf().first - m_condition.exchange_chem_pot(new_species_1, curr_species_1);
        event.set_dEpot(dEpot_1);
        // back up site 1 occupation
        event.set_original_occ_first_swap(_configdof().occ(event.occupational_change().first.site_index()));
        // // Site 1 modification finished, update configuration ....
        _configdof().occ(event.occupational_change().first.site_index()) = event.occupational_change().first.to_value();
        // mark the changes of the first site
        event.set_is_swapped(true);
        // Site 2
        // ---- set OccMod --------------
        event.occupational_change().second.set(mutating_sites.second, sublats.second, new_occs.second);
        // ---- set dspecies --------------
        for(int i = 0; i < event.dN().first.size(); ++i) {
          event.set_dN(i, 0);
        }
        Index curr_species_2 = m_site_swaps.sublat_to_mol()[sublats.second][curr_occs.second];
        Index new_species_2 = m_site_swaps.sublat_to_mol()[sublats.second][new_occs.second];
        event.set_dN(curr_species_2, -1);
        event.set_dN(new_species_2, 1);
        // ---- set dcorr --------------
        _set_dCorr(event, mutating_sites.second, sublats.second, curr_occs.second, new_occs.second, m_use_deltas, m_all_correlations);
        // ---- set dformation_energy --------------
        event.set_dEf(_eci() * event.dCorr().second.data());
        // ---- set dpotential_energy --------------
        double dEpot_2 = event.dEf().second - m_condition.exchange_chem_pot(new_species_2, curr_species_2);
        event.set_dEpot(dEpot_2);
        // Calculate dEpot after two swaps
        event.set_dEpot_swapped_twice(dEpot_1+dEpot_2);
        // Zeyu: after get dEpot_swapped_twice, change configuration back to origin....
        _configdof().occ(event.occupational_change().first.site_index()) = event.original_occ_first_swap();
        event.set_is_swapped(false);
        
    }

  /// \brief Calculate properties given current conditions
  void ChargeNeutralGrandCanonical::_update_properties() {

    // initialize properties and store pointers to the data strucures
    _vector_properties()["corr"] = correlations_vec(_configdof(), supercell(), _clexulator());
    m_corr = &_vector_property("corr");

    _vector_properties()["comp_n"] = CASM::comp_n(_configdof(), supercell());
    m_comp_n = &_vector_property("comp_n");

    _scalar_properties()["formation_energy"] = _eci() * corr().data();
    m_formation_energy = &_scalar_property("formation_energy");

    _scalar_properties()["potential_energy"] = formation_energy() - primclex().composition_axes().param_composition(comp_n()).dot(m_condition.param_chem_pot());
    m_potential_energy = &_scalar_property("potential_energy");

    if(debug()) {

      _print_correlations(corr(), "correlations", "corr", m_all_correlations);

      auto origin = primclex().composition_axes().origin();
      auto exchange_chem_pot = m_condition.exchange_chem_pot();
      auto param_chem_pot = m_condition.param_chem_pot();
      auto comp_x = primclex().composition_axes().param_composition(comp_n());
      auto M = primclex().composition_axes().dmol_dparam();

      _log().custom("Calculate properties");
      _log() << "Semi-grand canonical ensemble: \n"
             << "  Thermodynamic potential (per unitcell), phi = -kT*ln(Z)/N \n"
             << "  Partition function, Z = sum_i exp(-N*potential_energy_i/kT) \n"
             << "  composition, comp_n = origin + M * comp_x \n"
             << "  parametric chemical potential, param_chem_pot = M.transpose() * chem_pot \n"
             << "  potential_energy (per unitcell) = formation_energy - param_chem_pot*comp_x \n\n"

             << "components: " << jsonParser(primclex().composition_axes().components()) << "\n"
             << "M:\n" << M << "\n"
             << "origin: " << origin.transpose() << "\n"
             << "comp_n: " << comp_n().transpose() << "\n"
             << "comp_x: " << comp_x.transpose() << "\n"
             << "param_chem_pot: " << param_chem_pot.transpose() << "\n"
             << "  param_chem_pot*comp_x: " << param_chem_pot.dot(comp_x)  << "\n"
             << "formation_energy: " << formation_energy() << "\n"
             << "  formation_energy - param_chem_pot*comp_x: " << formation_energy() - param_chem_pot.dot(comp_x) << "\n"
             << "potential_energy: " << potential_energy() << "\n" << std::endl;
    }

  }

  /// \brief Generate supercell filling ConfigDoF from default configuration
  ConfigDoF ChargeNeutralGrandCanonical::_default_motif() const {
    _log().set("DoF");
    _log() << "motif configname: default\n";
    _log() << "using configuration with default occupation...\n" << std::endl;
    return Configuration(_supercell(), jsonParser(), Array<int>(_supercell().num_sites(), 0)).configdof();
  }

  /// \brief Generate minimum potential energy ConfigDoF
  ///
  /// Raises exception if it doesn't tile the supercell
  std::pair<ConfigDoF, std::string> ChargeNeutralGrandCanonical::_auto_motif(const GrandCanonicalConditions &cond) const {

    _log().set("DoF");
    _log() << "motif configname: auto\n";
    _log() << "searching for minimum potential energy motif..." << std::endl;

    double tol = 1e-6;
    auto compare = [&](double A, double B) {
      return A < B - tol;
    };

    _log() << "using conditions: \n";
    _log() << cond << std::endl;

    std::multimap<double, const Configuration *, decltype(compare)> configmap(compare);
    for(auto it = primclex().config_begin(); it != primclex().config_end(); ++it) {
      configmap.insert(std::make_pair(_eci() * correlations(*it, _clexulator()) - cond.param_chem_pot().dot(CASM::comp(*it)), &(*it)));
    }

    const Configuration &min_config = *(configmap.begin()->second);
    double min_potential_energy = configmap.begin()->first;
    auto eq_range = configmap.equal_range(min_potential_energy);

    if(std::distance(eq_range.first, eq_range.second) > 1) {
      _log() << "Warning: Found degenerate ground states with potential energy: "
             << std::setprecision(8) << min_potential_energy << std::endl;
      for(auto it = eq_range.first; it != eq_range.second; ++it) {
        _log() << "  " << it->second->name() << std::endl;
      }
      _log() << "using: " << min_config.name() << "\n" << std::endl;
    }
    else {
      _log() << "using: " << min_config.name() << " with potential energy: "
             << std::setprecision(8) << min_potential_energy << "\n" << std::endl;
    }

    return std::make_pair(
             min_config.fill_supercell(_supercell(), primclex().get_prim().factor_group()).configdof(),
             min_config.name());
  }

  /// \brief Generate minimum potential energy ConfigDoF for this supercell
  std::pair<ConfigDoF, std::string> ChargeNeutralGrandCanonical::_restricted_auto_motif(const GrandCanonicalConditions &cond) const {

    _log().set("DoF");
    _log() << "motif configname: restricted_auto\n";
    _log() << "searching for minimum potential energy motif..." << std::endl;

    double tol = 1e-6;
    auto compare = [&](double A, double B) {
      return A < B - tol;
    };

    _log() << "using conditions: \n";
    _log() << cond << std::endl;

    std::multimap<double, const Configuration *, decltype(compare)> configmap(compare);
    for(auto it = primclex().config_begin(); it != primclex().config_end(); ++it) {
      configmap.insert(std::make_pair(_eci() * correlations(*it, _clexulator()) - cond.param_chem_pot().dot(CASM::comp(*it)), &(*it)));
    }

    // used to check if configurations can fill the monte carlo supercell
    const Lattice &scel_lat = supercell().get_real_super_lattice();
    const SymGroup &g = primclex().get_prim().factor_group();
    auto begin = g.begin();
    auto end = g.end();

    // save iterators pointing to configs that will fill the supercell
    std::vector<decltype(configmap)::const_iterator> allowed;

    // iterate through groups of degenerate configs
    auto next_it = configmap.begin();
    while(next_it != configmap.end()) {
      auto eq_range = configmap.equal_range(next_it->first);

      // save allowed configs
      for(auto it = eq_range.first; it != eq_range.second; ++it) {
        const Lattice &motif_lat = it->second->get_supercell().get_real_super_lattice();
        if(is_supercell(scel_lat, motif_lat, begin, end, TOL).first != end) {
          allowed.push_back(it);
        }
      }

      // if some found, break
      if(allowed.size()) {
        break;
      }

      // else, continue to next group
      next_it = eq_range.second;
    }

    if(!allowed.size()) {
      _log() << "Found no enumerated configurations that will fill the supercell\n";
      _log() << "using configuration with default occupation..." << std::endl;
      return std::make_pair(
               Configuration(_supercell(), jsonParser(), Array<int>(_supercell().num_sites(), 0)).configdof(),
               "default");
    }

    if(allowed.size() > 1) {
      _log() << "Warning: Found degenerate allowed configurations with potential energy: "
             << std::setprecision(8) << allowed[0]->first << std::endl;
      for(auto it = allowed.begin(); it != allowed.end(); ++it) {
        _log() << "  " << (*it)->second->name() << std::endl;
      }
      _log() << "using: " << allowed[0]->second->name() << "\n" << std::endl;
    }
    else {
      _log() << "using: " << allowed[0]->second->name() << " with potential energy: "
             << std::setprecision(8) << allowed[0]->first << "\n" << std::endl;
    }

    return std::make_pair(
             allowed[0]->second->fill_supercell(_supercell(), g).configdof(),
             allowed[0]->second->name());
  }

  /// \brief Generate supercell filling ConfigDoF from configuration
  ConfigDoF ChargeNeutralGrandCanonical::_configname_motif(const std::string &configname) const {

    _log().set("DoF");
    _log() << "motif configname: " << configname << "\n";
    _log() << "using configation: " << configname << "\n" << std::endl;

    const Configuration &config = primclex().configuration(configname);
    const SymGroup &g = primclex().get_prim().factor_group();
    return config.fill_supercell(_supercell(), g).configdof();
  }
}
