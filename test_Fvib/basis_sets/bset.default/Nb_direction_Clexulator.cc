#include <cstddef>
#include "casm/clex/Clexulator.hh"



/****** CLEXULATOR CLASS FOR PRIM ******
{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.268900000000, 0.513900000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.731100000000, 0.731100000000, 0.245000000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.268900000000, 1.000000000000, 0.513900000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 1.000000000000, 0.513900000000, 0.268900000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.268900000000, 0.513900000000, 0.000000000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.731100000000, 0.245000000000, 0.731100000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.245000000000, 0.731100000000, 0.731100000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.513900000000, 0.268900000000, 0.000000000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.513900000000, 1.000000000000, 0.268900000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.755000000000, 0.486100000000, 0.486100000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.486100000000, 0.486100000000, 0.755000000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 0.486100000000, 0.755000000000, 0.486100000000 ],
      "occupant_dof" : [ "Li" ]
    },
    {
      "coordinate" : [ 1.000000000000, 1.000000000000, 0.719420000000 ],
      "occupant_dof" : [ "Nb", "Ta" ]
    },
    {
      "coordinate" : [ 0.000000000000, 0.719420000000, 0.000000000000 ],
      "occupant_dof" : [ "Nb", "Ta" ]
    },
    {
      "coordinate" : [ 0.719420000000, 1.000000000000, 0.000000000000 ],
      "occupant_dof" : [ "Nb", "Ta" ]
    },
    {
      "coordinate" : [ 0.280580000000, 0.280580000000, 0.280580000000 ],
      "occupant_dof" : [ "Nb", "Ta" ]
    },
    {
      "coordinate" : [ 1.000000000000, 0.762100000000, 0.515100000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.237900000000, 0.237900000000, 0.753000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.762100000000, 1.000000000000, 0.515100000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.000000000000, 0.515100000000, 0.762100000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.762100000000, 0.515100000000, 0.000000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.237900000000, 0.753000000000, 0.237900000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.753000000000, 0.237900000000, 0.237900000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.515100000000, 0.762100000000, 0.000000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.515100000000, 1.000000000000, 0.762100000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.247000000000, 0.484900000000, 0.484900000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.484900000000, 0.484900000000, 0.247000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.484900000000, 0.247000000000, 0.484900000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 1.000000000000, 1.000000000000, 0.223800000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.000000000000, 0.223800000000, 0.000000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.223800000000, 0.000000000000, 0.000000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.776200000000, 0.776200000000, 0.776200000000 ],
      "occupant_dof" : [ "O" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "lattice_vectors" : [
    [ -4.221158000000, 4.221158000000, 4.221158000000 ],
    [ 4.221158000000, -4.221158000000, 4.221158000000 ],
    [ 4.221159000000, 4.221158000000, -4.221158000000 ]
  ],
  "title" : "Nb_direction"
}**/


/// \brief Returns a Clexulator_impl::Base* owning a Nb_direction_Clexulator
extern "C" CASM::Clexulator_impl::Base* make_Nb_direction_Clexulator();

namespace CASM {

  class Nb_direction_Clexulator : public Clexulator_impl::Base {

  public:

    Nb_direction_Clexulator();

    ~Nb_direction_Clexulator();

    /// \brief Clone the Nb_direction_Clexulator
    std::unique_ptr<Nb_direction_Clexulator> clone() const { 
      return std::unique_ptr<Nb_direction_Clexulator>(_clone()); 
    }

    /// \brief Calculate contribution to global correlations from one unit cell
    void calc_global_corr_contribution(double *corr_begin) const override;

    /// \brief Calculate contribution to select global correlations from one unit cell
    void calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;

    /// \brief Calculate point correlations about basis site 'b_index'
    void calc_point_corr(int b_index, double *corr_begin) const override;

    /// \brief Calculate select point correlations about basis site 'b_index'
    void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;

    /// \brief Calculate the change in point correlations due to changing an occupant
    void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const override;

    /// \brief Calculate the change in select point correlations due to changing an occupant
    void calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;


  private:

    /// \brief Clone the Clexulator
    virtual Nb_direction_Clexulator* _clone() const override {
      return new Nb_direction_Clexulator(*this);
    }

    // typedef for method pointers
    typedef double (Nb_direction_Clexulator::*BasisFuncPtr)() const;

    // typedef for method pointers
    typedef double (Nb_direction_Clexulator::*DeltaBasisFuncPtr)(int, int) const;

    // array of pointers to member functions for calculating basis functions
    BasisFuncPtr m_orbit_func_list[22];

    // array of pointers to member functions for calculating flower functions
    BasisFuncPtr m_flower_func_lists[32][22];

    // array of pointers to member functions for calculating DELTA flower functions
    DeltaBasisFuncPtr m_delta_func_lists[32][22];

    // Occupation Function tables for basis sites in asymmetric unit 1:
    //   - basis site 12:
    double m_occ_func_12_0[2];

    //   - basis site 13:
    double m_occ_func_13_0[2];

    //   - basis site 14:
    double m_occ_func_14_0[2];

    //   - basis site 15:
    double m_occ_func_15_0[2];

    // Occupation Function accessors for basis site 12:
    const double &occ_func_12_0(const int &nlist_ind)const{return m_occ_func_12_0[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}

    // Occupation Function accessors for basis site 13:
    const double &occ_func_13_0(const int &nlist_ind)const{return m_occ_func_13_0[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}

    // Occupation Function accessors for basis site 14:
    const double &occ_func_14_0(const int &nlist_ind)const{return m_occ_func_14_0[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}

    // Occupation Function accessors for basis site 15:
    const double &occ_func_15_0(const int &nlist_ind)const{return m_occ_func_15_0[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}

    //default functions for basis function evaluation 
    double zero_func() const{ return 0.0;};
    double zero_func(int,int) const{ return 0.0;};

    double eval_bfunc_0_0_0() const;

    double eval_bfunc_1_0_0() const;

    double site_eval_at_12_bfunc_1_0_0() const;

    double delta_site_eval_at_12_bfunc_1_0_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_1_0_0() const;

    double delta_site_eval_at_13_bfunc_1_0_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_1_0_0() const;

    double delta_site_eval_at_14_bfunc_1_0_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_1_0_0() const;

    double delta_site_eval_at_15_bfunc_1_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_0_0() const;

    double site_eval_at_12_bfunc_2_0_0() const;

    double delta_site_eval_at_12_bfunc_2_0_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_0_0() const;

    double delta_site_eval_at_13_bfunc_2_0_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_0_0() const;

    double delta_site_eval_at_14_bfunc_2_0_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_0_0() const;

    double delta_site_eval_at_15_bfunc_2_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_1_0() const;

    double site_eval_at_12_bfunc_2_1_0() const;

    double delta_site_eval_at_12_bfunc_2_1_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_1_0() const;

    double delta_site_eval_at_13_bfunc_2_1_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_1_0() const;

    double delta_site_eval_at_14_bfunc_2_1_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_1_0() const;

    double delta_site_eval_at_15_bfunc_2_1_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_2_0() const;

    double site_eval_at_12_bfunc_2_2_0() const;

    double delta_site_eval_at_12_bfunc_2_2_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_2_0() const;

    double delta_site_eval_at_13_bfunc_2_2_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_2_0() const;

    double delta_site_eval_at_14_bfunc_2_2_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_2_0() const;

    double delta_site_eval_at_15_bfunc_2_2_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_3_0() const;

    double site_eval_at_12_bfunc_2_3_0() const;

    double delta_site_eval_at_12_bfunc_2_3_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_3_0() const;

    double delta_site_eval_at_13_bfunc_2_3_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_3_0() const;

    double delta_site_eval_at_14_bfunc_2_3_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_3_0() const;

    double delta_site_eval_at_15_bfunc_2_3_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_4_0() const;

    double site_eval_at_12_bfunc_2_4_0() const;

    double delta_site_eval_at_12_bfunc_2_4_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_4_0() const;

    double delta_site_eval_at_13_bfunc_2_4_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_4_0() const;

    double delta_site_eval_at_14_bfunc_2_4_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_4_0() const;

    double delta_site_eval_at_15_bfunc_2_4_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_5_0() const;

    double site_eval_at_12_bfunc_2_5_0() const;

    double delta_site_eval_at_12_bfunc_2_5_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_5_0() const;

    double delta_site_eval_at_13_bfunc_2_5_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_5_0() const;

    double delta_site_eval_at_14_bfunc_2_5_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_5_0() const;

    double delta_site_eval_at_15_bfunc_2_5_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_6_0() const;

    double site_eval_at_12_bfunc_2_6_0() const;

    double delta_site_eval_at_12_bfunc_2_6_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_6_0() const;

    double delta_site_eval_at_13_bfunc_2_6_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_6_0() const;

    double delta_site_eval_at_14_bfunc_2_6_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_6_0() const;

    double delta_site_eval_at_15_bfunc_2_6_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_7_0() const;

    double site_eval_at_12_bfunc_2_7_0() const;

    double delta_site_eval_at_12_bfunc_2_7_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_7_0() const;

    double delta_site_eval_at_13_bfunc_2_7_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_7_0() const;

    double delta_site_eval_at_14_bfunc_2_7_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_7_0() const;

    double delta_site_eval_at_15_bfunc_2_7_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_8_0() const;

    double site_eval_at_12_bfunc_2_8_0() const;

    double delta_site_eval_at_12_bfunc_2_8_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_8_0() const;

    double delta_site_eval_at_13_bfunc_2_8_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_8_0() const;

    double delta_site_eval_at_14_bfunc_2_8_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_8_0() const;

    double delta_site_eval_at_15_bfunc_2_8_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_9_0() const;

    double site_eval_at_12_bfunc_2_9_0() const;

    double delta_site_eval_at_12_bfunc_2_9_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_9_0() const;

    double delta_site_eval_at_13_bfunc_2_9_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_9_0() const;

    double delta_site_eval_at_14_bfunc_2_9_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_9_0() const;

    double delta_site_eval_at_15_bfunc_2_9_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_10_0() const;

    double site_eval_at_12_bfunc_2_10_0() const;

    double delta_site_eval_at_12_bfunc_2_10_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_10_0() const;

    double delta_site_eval_at_13_bfunc_2_10_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_10_0() const;

    double delta_site_eval_at_14_bfunc_2_10_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_10_0() const;

    double delta_site_eval_at_15_bfunc_2_10_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_11_0() const;

    double site_eval_at_12_bfunc_2_11_0() const;

    double delta_site_eval_at_12_bfunc_2_11_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_11_0() const;

    double delta_site_eval_at_13_bfunc_2_11_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_11_0() const;

    double delta_site_eval_at_14_bfunc_2_11_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_11_0() const;

    double delta_site_eval_at_15_bfunc_2_11_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_12_0() const;

    double site_eval_at_12_bfunc_2_12_0() const;

    double delta_site_eval_at_12_bfunc_2_12_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_12_0() const;

    double delta_site_eval_at_13_bfunc_2_12_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_12_0() const;

    double delta_site_eval_at_14_bfunc_2_12_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_12_0() const;

    double delta_site_eval_at_15_bfunc_2_12_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_13_0() const;

    double site_eval_at_12_bfunc_2_13_0() const;

    double delta_site_eval_at_12_bfunc_2_13_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_2_13_0() const;

    double delta_site_eval_at_13_bfunc_2_13_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_2_13_0() const;

    double delta_site_eval_at_14_bfunc_2_13_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_2_13_0() const;

    double delta_site_eval_at_15_bfunc_2_13_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_0_0() const;

    double site_eval_at_12_bfunc_3_0_0() const;

    double delta_site_eval_at_12_bfunc_3_0_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_3_0_0() const;

    double delta_site_eval_at_13_bfunc_3_0_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_3_0_0() const;

    double delta_site_eval_at_14_bfunc_3_0_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_3_0_0() const;

    double delta_site_eval_at_15_bfunc_3_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_1_0() const;

    double site_eval_at_12_bfunc_3_1_0() const;

    double delta_site_eval_at_12_bfunc_3_1_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_3_1_0() const;

    double delta_site_eval_at_13_bfunc_3_1_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_3_1_0() const;

    double delta_site_eval_at_14_bfunc_3_1_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_3_1_0() const;

    double delta_site_eval_at_15_bfunc_3_1_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_2_0() const;

    double site_eval_at_12_bfunc_3_2_0() const;

    double delta_site_eval_at_12_bfunc_3_2_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_3_2_0() const;

    double delta_site_eval_at_13_bfunc_3_2_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_3_2_0() const;

    double delta_site_eval_at_14_bfunc_3_2_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_3_2_0() const;

    double delta_site_eval_at_15_bfunc_3_2_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_3_0() const;

    double site_eval_at_12_bfunc_3_3_0() const;

    double delta_site_eval_at_12_bfunc_3_3_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_3_3_0() const;

    double delta_site_eval_at_13_bfunc_3_3_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_3_3_0() const;

    double delta_site_eval_at_14_bfunc_3_3_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_3_3_0() const;

    double delta_site_eval_at_15_bfunc_3_3_0(int occ_i, int occ_f) const;

    double eval_bfunc_4_0_0() const;

    double site_eval_at_12_bfunc_4_0_0() const;

    double delta_site_eval_at_12_bfunc_4_0_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_4_0_0() const;

    double delta_site_eval_at_13_bfunc_4_0_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_4_0_0() const;

    double delta_site_eval_at_14_bfunc_4_0_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_4_0_0() const;

    double delta_site_eval_at_15_bfunc_4_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_4_1_0() const;

    double site_eval_at_12_bfunc_4_1_0() const;

    double delta_site_eval_at_12_bfunc_4_1_0(int occ_i, int occ_f) const;

    double site_eval_at_13_bfunc_4_1_0() const;

    double delta_site_eval_at_13_bfunc_4_1_0(int occ_i, int occ_f) const;

    double site_eval_at_14_bfunc_4_1_0() const;

    double delta_site_eval_at_14_bfunc_4_1_0(int occ_i, int occ_f) const;

    double site_eval_at_15_bfunc_4_1_0() const;

    double delta_site_eval_at_15_bfunc_4_1_0(int occ_i, int occ_f) const;


  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Nb_direction_Clexulator::Nb_direction_Clexulator() :
    Clexulator_impl::Base(59, 22) {
    m_occ_func_12_0[0] = -1.0000000000, m_occ_func_12_0[1] = 1.0000000000;

    m_occ_func_13_0[0] = -1.0000000000, m_occ_func_13_0[1] = 1.0000000000;

    m_occ_func_14_0[0] = -1.0000000000, m_occ_func_14_0[1] = 1.0000000000;

    m_occ_func_15_0[0] = -1.0000000000, m_occ_func_15_0[1] = 1.0000000000;

    m_orbit_func_list[0] = &Nb_direction_Clexulator::eval_bfunc_0_0_0;
    m_orbit_func_list[1] = &Nb_direction_Clexulator::eval_bfunc_1_0_0;
    m_orbit_func_list[2] = &Nb_direction_Clexulator::eval_bfunc_2_0_0;
    m_orbit_func_list[3] = &Nb_direction_Clexulator::eval_bfunc_2_1_0;
    m_orbit_func_list[4] = &Nb_direction_Clexulator::eval_bfunc_2_2_0;
    m_orbit_func_list[5] = &Nb_direction_Clexulator::eval_bfunc_2_3_0;
    m_orbit_func_list[6] = &Nb_direction_Clexulator::eval_bfunc_2_4_0;
    m_orbit_func_list[7] = &Nb_direction_Clexulator::eval_bfunc_2_5_0;
    m_orbit_func_list[8] = &Nb_direction_Clexulator::eval_bfunc_2_6_0;
    m_orbit_func_list[9] = &Nb_direction_Clexulator::eval_bfunc_2_7_0;
    m_orbit_func_list[10] = &Nb_direction_Clexulator::eval_bfunc_2_8_0;
    m_orbit_func_list[11] = &Nb_direction_Clexulator::eval_bfunc_2_9_0;
    m_orbit_func_list[12] = &Nb_direction_Clexulator::eval_bfunc_2_10_0;
    m_orbit_func_list[13] = &Nb_direction_Clexulator::eval_bfunc_2_11_0;
    m_orbit_func_list[14] = &Nb_direction_Clexulator::eval_bfunc_2_12_0;
    m_orbit_func_list[15] = &Nb_direction_Clexulator::eval_bfunc_2_13_0;
    m_orbit_func_list[16] = &Nb_direction_Clexulator::eval_bfunc_3_0_0;
    m_orbit_func_list[17] = &Nb_direction_Clexulator::eval_bfunc_3_1_0;
    m_orbit_func_list[18] = &Nb_direction_Clexulator::eval_bfunc_3_2_0;
    m_orbit_func_list[19] = &Nb_direction_Clexulator::eval_bfunc_3_3_0;
    m_orbit_func_list[20] = &Nb_direction_Clexulator::eval_bfunc_4_0_0;
    m_orbit_func_list[21] = &Nb_direction_Clexulator::eval_bfunc_4_1_0;


    m_flower_func_lists[0][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[0][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[1][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[1][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[2][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[2][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[3][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[3][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[4][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[4][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[5][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[5][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[6][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[6][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[7][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[7][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[8][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[8][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[9][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[9][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[10][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[10][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[11][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[11][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[12][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[12][1] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_1_0_0;
    m_flower_func_lists[12][2] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_0_0;
    m_flower_func_lists[12][3] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_1_0;
    m_flower_func_lists[12][4] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_2_0;
    m_flower_func_lists[12][5] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_3_0;
    m_flower_func_lists[12][6] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_4_0;
    m_flower_func_lists[12][7] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_5_0;
    m_flower_func_lists[12][8] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_6_0;
    m_flower_func_lists[12][9] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_7_0;
    m_flower_func_lists[12][10] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_8_0;
    m_flower_func_lists[12][11] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_9_0;
    m_flower_func_lists[12][12] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_10_0;
    m_flower_func_lists[12][13] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_11_0;
    m_flower_func_lists[12][14] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_12_0;
    m_flower_func_lists[12][15] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_2_13_0;
    m_flower_func_lists[12][16] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_3_0_0;
    m_flower_func_lists[12][17] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_3_1_0;
    m_flower_func_lists[12][18] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_3_2_0;
    m_flower_func_lists[12][19] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_3_3_0;
    m_flower_func_lists[12][20] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_4_0_0;
    m_flower_func_lists[12][21] = &Nb_direction_Clexulator::site_eval_at_12_bfunc_4_1_0;


    m_flower_func_lists[13][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[13][1] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_1_0_0;
    m_flower_func_lists[13][2] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_0_0;
    m_flower_func_lists[13][3] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_1_0;
    m_flower_func_lists[13][4] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_2_0;
    m_flower_func_lists[13][5] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_3_0;
    m_flower_func_lists[13][6] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_4_0;
    m_flower_func_lists[13][7] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_5_0;
    m_flower_func_lists[13][8] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_6_0;
    m_flower_func_lists[13][9] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_7_0;
    m_flower_func_lists[13][10] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_8_0;
    m_flower_func_lists[13][11] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_9_0;
    m_flower_func_lists[13][12] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_10_0;
    m_flower_func_lists[13][13] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_11_0;
    m_flower_func_lists[13][14] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_12_0;
    m_flower_func_lists[13][15] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_2_13_0;
    m_flower_func_lists[13][16] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_3_0_0;
    m_flower_func_lists[13][17] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_3_1_0;
    m_flower_func_lists[13][18] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_3_2_0;
    m_flower_func_lists[13][19] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_3_3_0;
    m_flower_func_lists[13][20] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_4_0_0;
    m_flower_func_lists[13][21] = &Nb_direction_Clexulator::site_eval_at_13_bfunc_4_1_0;


    m_flower_func_lists[14][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[14][1] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_1_0_0;
    m_flower_func_lists[14][2] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_0_0;
    m_flower_func_lists[14][3] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_1_0;
    m_flower_func_lists[14][4] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_2_0;
    m_flower_func_lists[14][5] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_3_0;
    m_flower_func_lists[14][6] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_4_0;
    m_flower_func_lists[14][7] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_5_0;
    m_flower_func_lists[14][8] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_6_0;
    m_flower_func_lists[14][9] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_7_0;
    m_flower_func_lists[14][10] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_8_0;
    m_flower_func_lists[14][11] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_9_0;
    m_flower_func_lists[14][12] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_10_0;
    m_flower_func_lists[14][13] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_11_0;
    m_flower_func_lists[14][14] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_12_0;
    m_flower_func_lists[14][15] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_2_13_0;
    m_flower_func_lists[14][16] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_3_0_0;
    m_flower_func_lists[14][17] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_3_1_0;
    m_flower_func_lists[14][18] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_3_2_0;
    m_flower_func_lists[14][19] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_3_3_0;
    m_flower_func_lists[14][20] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_4_0_0;
    m_flower_func_lists[14][21] = &Nb_direction_Clexulator::site_eval_at_14_bfunc_4_1_0;


    m_flower_func_lists[15][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[15][1] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_1_0_0;
    m_flower_func_lists[15][2] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_0_0;
    m_flower_func_lists[15][3] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_1_0;
    m_flower_func_lists[15][4] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_2_0;
    m_flower_func_lists[15][5] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_3_0;
    m_flower_func_lists[15][6] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_4_0;
    m_flower_func_lists[15][7] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_5_0;
    m_flower_func_lists[15][8] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_6_0;
    m_flower_func_lists[15][9] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_7_0;
    m_flower_func_lists[15][10] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_8_0;
    m_flower_func_lists[15][11] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_9_0;
    m_flower_func_lists[15][12] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_10_0;
    m_flower_func_lists[15][13] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_11_0;
    m_flower_func_lists[15][14] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_12_0;
    m_flower_func_lists[15][15] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_2_13_0;
    m_flower_func_lists[15][16] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_3_0_0;
    m_flower_func_lists[15][17] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_3_1_0;
    m_flower_func_lists[15][18] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_3_2_0;
    m_flower_func_lists[15][19] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_3_3_0;
    m_flower_func_lists[15][20] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_4_0_0;
    m_flower_func_lists[15][21] = &Nb_direction_Clexulator::site_eval_at_15_bfunc_4_1_0;


    m_flower_func_lists[16][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[16][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[17][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[17][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[18][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[18][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[19][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[19][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[20][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[20][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[21][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[21][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[22][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[22][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[23][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[23][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[24][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[24][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[25][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[25][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[26][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[26][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[27][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[27][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[28][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[28][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[29][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[29][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[30][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[30][21] = &Nb_direction_Clexulator::zero_func;


    m_flower_func_lists[31][0] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][1] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][2] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][3] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][4] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][5] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][6] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][7] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][8] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][9] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][10] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][11] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][12] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][13] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][14] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][15] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][16] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][17] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][18] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][19] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][20] = &Nb_direction_Clexulator::zero_func;
    m_flower_func_lists[31][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[0][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[0][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[1][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[1][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[2][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[2][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[3][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[3][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[4][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[4][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[5][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[5][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[6][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[6][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[7][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[7][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[8][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[8][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[9][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[9][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[10][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[10][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[11][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[11][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[12][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[12][1] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_1_0_0;
    m_delta_func_lists[12][2] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_0_0;
    m_delta_func_lists[12][3] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_1_0;
    m_delta_func_lists[12][4] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_2_0;
    m_delta_func_lists[12][5] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_3_0;
    m_delta_func_lists[12][6] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_4_0;
    m_delta_func_lists[12][7] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_5_0;
    m_delta_func_lists[12][8] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_6_0;
    m_delta_func_lists[12][9] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_7_0;
    m_delta_func_lists[12][10] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_8_0;
    m_delta_func_lists[12][11] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_9_0;
    m_delta_func_lists[12][12] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_10_0;
    m_delta_func_lists[12][13] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_11_0;
    m_delta_func_lists[12][14] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_12_0;
    m_delta_func_lists[12][15] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_13_0;
    m_delta_func_lists[12][16] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_0_0;
    m_delta_func_lists[12][17] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_1_0;
    m_delta_func_lists[12][18] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_2_0;
    m_delta_func_lists[12][19] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_3_0;
    m_delta_func_lists[12][20] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_4_0_0;
    m_delta_func_lists[12][21] = &Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_4_1_0;


    m_delta_func_lists[13][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[13][1] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_1_0_0;
    m_delta_func_lists[13][2] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_0_0;
    m_delta_func_lists[13][3] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_1_0;
    m_delta_func_lists[13][4] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_2_0;
    m_delta_func_lists[13][5] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_3_0;
    m_delta_func_lists[13][6] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_4_0;
    m_delta_func_lists[13][7] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_5_0;
    m_delta_func_lists[13][8] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_6_0;
    m_delta_func_lists[13][9] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_7_0;
    m_delta_func_lists[13][10] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_8_0;
    m_delta_func_lists[13][11] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_9_0;
    m_delta_func_lists[13][12] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_10_0;
    m_delta_func_lists[13][13] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_11_0;
    m_delta_func_lists[13][14] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_12_0;
    m_delta_func_lists[13][15] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_13_0;
    m_delta_func_lists[13][16] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_0_0;
    m_delta_func_lists[13][17] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_1_0;
    m_delta_func_lists[13][18] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_2_0;
    m_delta_func_lists[13][19] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_3_0;
    m_delta_func_lists[13][20] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_4_0_0;
    m_delta_func_lists[13][21] = &Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_4_1_0;


    m_delta_func_lists[14][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[14][1] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_1_0_0;
    m_delta_func_lists[14][2] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_0_0;
    m_delta_func_lists[14][3] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_1_0;
    m_delta_func_lists[14][4] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_2_0;
    m_delta_func_lists[14][5] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_3_0;
    m_delta_func_lists[14][6] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_4_0;
    m_delta_func_lists[14][7] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_5_0;
    m_delta_func_lists[14][8] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_6_0;
    m_delta_func_lists[14][9] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_7_0;
    m_delta_func_lists[14][10] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_8_0;
    m_delta_func_lists[14][11] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_9_0;
    m_delta_func_lists[14][12] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_10_0;
    m_delta_func_lists[14][13] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_11_0;
    m_delta_func_lists[14][14] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_12_0;
    m_delta_func_lists[14][15] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_13_0;
    m_delta_func_lists[14][16] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_0_0;
    m_delta_func_lists[14][17] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_1_0;
    m_delta_func_lists[14][18] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_2_0;
    m_delta_func_lists[14][19] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_3_0;
    m_delta_func_lists[14][20] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_4_0_0;
    m_delta_func_lists[14][21] = &Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_4_1_0;


    m_delta_func_lists[15][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[15][1] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_1_0_0;
    m_delta_func_lists[15][2] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_0_0;
    m_delta_func_lists[15][3] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_1_0;
    m_delta_func_lists[15][4] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_2_0;
    m_delta_func_lists[15][5] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_3_0;
    m_delta_func_lists[15][6] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_4_0;
    m_delta_func_lists[15][7] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_5_0;
    m_delta_func_lists[15][8] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_6_0;
    m_delta_func_lists[15][9] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_7_0;
    m_delta_func_lists[15][10] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_8_0;
    m_delta_func_lists[15][11] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_9_0;
    m_delta_func_lists[15][12] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_10_0;
    m_delta_func_lists[15][13] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_11_0;
    m_delta_func_lists[15][14] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_12_0;
    m_delta_func_lists[15][15] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_13_0;
    m_delta_func_lists[15][16] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_0_0;
    m_delta_func_lists[15][17] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_1_0;
    m_delta_func_lists[15][18] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_2_0;
    m_delta_func_lists[15][19] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_3_0;
    m_delta_func_lists[15][20] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_4_0_0;
    m_delta_func_lists[15][21] = &Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_4_1_0;


    m_delta_func_lists[16][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[16][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[17][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[17][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[18][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[18][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[19][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[19][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[20][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[20][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[21][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[21][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[22][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[22][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[23][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[23][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[24][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[24][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[25][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[25][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[26][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[26][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[27][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[27][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[28][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[28][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[29][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[29][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[30][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[30][21] = &Nb_direction_Clexulator::zero_func;


    m_delta_func_lists[31][0] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][1] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][2] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][3] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][4] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][5] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][6] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][7] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][8] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][9] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][10] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][11] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][12] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][13] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][14] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][15] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][16] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][17] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][18] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][19] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][20] = &Nb_direction_Clexulator::zero_func;
    m_delta_func_lists[31][21] = &Nb_direction_Clexulator::zero_func;


    m_weight_matrix.row(0) << 3, -1, -1;
    m_weight_matrix.row(1) << -1, 3, -1;
    m_weight_matrix.row(2) << -1, -1, 3;

    m_neighborhood = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -2, -2)},
      {UnitCellCoord(12, -2, -2, -1)},
      {UnitCellCoord(14, -2, -2, -1)},
      {UnitCellCoord(12, -2, -1, -2)},
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(13, -2, -1, -1)},
      {UnitCellCoord(14, -2, -1, -1)},
      {UnitCellCoord(15, -2, -1, -1)},
      {UnitCellCoord(12, -2, -1, 0)},
      {UnitCellCoord(14, -2, -1, 0)},
      {UnitCellCoord(12, -2, 0, -1)},
      {UnitCellCoord(15, -2, 0, -1)},
      {UnitCellCoord(12, -2, 0, 0)},
      {UnitCellCoord(14, -2, 0, 0)},
      {UnitCellCoord(12, -1, -2, -1)},
      {UnitCellCoord(13, -1, -2, -1)},
      {UnitCellCoord(14, -1, -2, -1)},
      {UnitCellCoord(15, -1, -2, -1)},
      {UnitCellCoord(12, -1, -2, 0)},
      {UnitCellCoord(14, -1, -2, 0)},
      {UnitCellCoord(12, -1, -1, -2)},
      {UnitCellCoord(13, -1, -1, -2)},
      {UnitCellCoord(14, -1, -1, -2)},
      {UnitCellCoord(15, -1, -1, -2)},
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(13, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(15, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(13, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(15, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 1)},
      {UnitCellCoord(12, -1, 0, -2)},
      {UnitCellCoord(15, -1, 0, -2)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(13, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(13, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, -1, 0, 1)},
      {UnitCellCoord(13, -1, 0, 1)},
      {UnitCellCoord(14, -1, 0, 1)},
      {UnitCellCoord(15, -1, 0, 1)},
      {UnitCellCoord(12, -1, 1, -1)},
      {UnitCellCoord(15, -1, 1, -1)},
      {UnitCellCoord(12, -1, 1, 0)},
      {UnitCellCoord(13, -1, 1, 0)},
      {UnitCellCoord(14, -1, 1, 0)},
      {UnitCellCoord(15, -1, 1, 0)},
      {UnitCellCoord(13, 0, -2, 0)},
      {UnitCellCoord(14, 0, -2, 0)},
      {UnitCellCoord(12, 0, -1, -1)},
      {UnitCellCoord(13, 0, -1, -1)},
      {UnitCellCoord(14, 0, -1, -1)},
      {UnitCellCoord(15, 0, -1, -1)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(15, 0, -1, 0)},
      {UnitCellCoord(12, 0, -1, 1)},
      {UnitCellCoord(13, 0, -1, 1)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(15, 0, -1, 1)},
      {UnitCellCoord(12, 0, 0, -2)},
      {UnitCellCoord(15, 0, 0, -2)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(14, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 0, 1)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 2)},
      {UnitCellCoord(12, 0, 1, -1)},
      {UnitCellCoord(13, 0, 1, -1)},
      {UnitCellCoord(14, 0, 1, -1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(13, 0, 1, 0)},
      {UnitCellCoord(14, 0, 1, 0)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(12, 0, 1, 1)},
      {UnitCellCoord(13, 0, 1, 1)},
      {UnitCellCoord(14, 0, 1, 1)},
      {UnitCellCoord(15, 0, 1, 1)},
      {UnitCellCoord(15, 0, 2, 0)},
      {UnitCellCoord(12, 1, -1, 0)},
      {UnitCellCoord(13, 1, -1, 0)},
      {UnitCellCoord(14, 1, -1, 0)},
      {UnitCellCoord(15, 1, -1, 0)},
      {UnitCellCoord(13, 1, -1, 1)},
      {UnitCellCoord(14, 1, -1, 1)},
      {UnitCellCoord(12, 1, 0, -1)},
      {UnitCellCoord(13, 1, 0, -1)},
      {UnitCellCoord(14, 1, 0, -1)},
      {UnitCellCoord(15, 1, 0, -1)},
      {UnitCellCoord(12, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(14, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(12, 1, 0, 1)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(15, 1, 0, 1)},
      {UnitCellCoord(13, 1, 0, 2)},
      {UnitCellCoord(14, 1, 0, 2)},
      {UnitCellCoord(13, 1, 1, -1)},
      {UnitCellCoord(15, 1, 1, -1)},
      {UnitCellCoord(12, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(14, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(12, 1, 1, 1)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(15, 1, 1, 1)},
      {UnitCellCoord(12, 1, 1, 2)},
      {UnitCellCoord(13, 1, 1, 2)},
      {UnitCellCoord(14, 1, 1, 2)},
      {UnitCellCoord(15, 1, 1, 2)},
      {UnitCellCoord(15, 1, 2, 0)},
      {UnitCellCoord(12, 1, 2, 1)},
      {UnitCellCoord(13, 1, 2, 1)},
      {UnitCellCoord(14, 1, 2, 1)},
      {UnitCellCoord(15, 1, 2, 1)},
      {UnitCellCoord(13, 2, 0, 0)},
      {UnitCellCoord(13, 2, 0, 1)},
      {UnitCellCoord(13, 2, 1, 0)},
      {UnitCellCoord(15, 2, 1, 0)},
      {UnitCellCoord(12, 2, 1, 1)},
      {UnitCellCoord(13, 2, 1, 1)},
      {UnitCellCoord(14, 2, 1, 1)},
      {UnitCellCoord(15, 2, 1, 1)},
      {UnitCellCoord(13, 2, 1, 2)},
      {UnitCellCoord(14, 2, 1, 2)},
      {UnitCellCoord(13, 2, 2, 1)},
      {UnitCellCoord(15, 2, 2, 1)},
      {UnitCellCoord(13, 2, 2, 2)},
      {UnitCellCoord(15, 2, 2, 2)}
    };


    m_orbit_neighborhood.resize(corr_size());
    m_orbit_neighborhood[0] = std::set<UnitCellCoord> {
    };

    m_orbit_neighborhood[1] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)}
    };

    m_orbit_neighborhood[2] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

    m_orbit_neighborhood[3] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)}
    };

    m_orbit_neighborhood[4] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(13, 2, 1, 1)}
    };

    m_orbit_neighborhood[5] = std::set<UnitCellCoord> {
      {UnitCellCoord(15, -1, -1, -1)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 0, 1)},
      {UnitCellCoord(13, 0, 1, 0)},
      {UnitCellCoord(14, 1, 0, 0)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

    m_orbit_neighborhood[6] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(13, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(13, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(15, 0, -1, 0)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(14, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 0, 1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(14, 0, 1, 0)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(12, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(12, 1, 1, 1)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)}
    };

    m_orbit_neighborhood[7] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -1, -2)},
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(14, -2, -1, -1)},
      {UnitCellCoord(12, -2, 0, -1)},
      {UnitCellCoord(12, -1, -2, -1)},
      {UnitCellCoord(14, -1, -2, 0)},
      {UnitCellCoord(12, -1, -1, -2)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(13, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 1)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 1)},
      {UnitCellCoord(15, -1, 1, 0)},
      {UnitCellCoord(12, 0, -1, -1)},
      {UnitCellCoord(13, 0, -1, -1)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(13, 0, -1, 1)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(12, 0, 1, -1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(14, 0, 1, 1)},
      {UnitCellCoord(15, 0, 1, 1)},
      {UnitCellCoord(13, 1, -1, 0)},
      {UnitCellCoord(13, 1, 0, -1)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(15, 1, 0, 1)},
      {UnitCellCoord(15, 1, 1, -1)},
      {UnitCellCoord(12, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 2)},
      {UnitCellCoord(15, 1, 2, 0)},
      {UnitCellCoord(15, 1, 2, 1)},
      {UnitCellCoord(13, 2, 0, 1)},
      {UnitCellCoord(13, 2, 1, 1)},
      {UnitCellCoord(15, 2, 1, 1)},
      {UnitCellCoord(13, 2, 1, 2)}
    };

    m_orbit_neighborhood[8] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(13, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(15, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(13, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, 0, -1, -1)},
      {UnitCellCoord(13, 0, -1, -1)},
      {UnitCellCoord(14, 0, -1, -1)},
      {UnitCellCoord(15, 0, -1, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 1, 1)},
      {UnitCellCoord(13, 0, 1, 1)},
      {UnitCellCoord(14, 0, 1, 1)},
      {UnitCellCoord(15, 0, 1, 1)},
      {UnitCellCoord(12, 1, 0, 1)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(15, 1, 0, 1)},
      {UnitCellCoord(12, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(14, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)}
    };

    m_orbit_neighborhood[9] = std::set<UnitCellCoord> {
      {UnitCellCoord(15, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 0, 1)},
      {UnitCellCoord(13, 0, 1, 0)},
      {UnitCellCoord(14, 1, 0, 0)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)}
    };

    m_orbit_neighborhood[10] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -2, -1)},
      {UnitCellCoord(14, -2, -1, 0)},
      {UnitCellCoord(14, -1, -2, -1)},
      {UnitCellCoord(12, -1, -1, -2)},
      {UnitCellCoord(13, -1, -1, -1)},
      {UnitCellCoord(12, -1, 0, -2)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 1)},
      {UnitCellCoord(15, -1, 1, -1)},
      {UnitCellCoord(12, -1, 1, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 1)},
      {UnitCellCoord(14, 0, 1, 0)},
      {UnitCellCoord(13, 1, -1, 0)},
      {UnitCellCoord(13, 1, -1, 1)},
      {UnitCellCoord(15, 1, 0, -1)},
      {UnitCellCoord(12, 1, 0, 0)},
      {UnitCellCoord(14, 1, 0, 2)},
      {UnitCellCoord(15, 1, 1, 1)},
      {UnitCellCoord(13, 1, 1, 2)},
      {UnitCellCoord(15, 1, 2, 1)},
      {UnitCellCoord(13, 2, 1, 0)},
      {UnitCellCoord(15, 2, 2, 1)}
    };

    m_orbit_neighborhood[11] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -2, -2)},
      {UnitCellCoord(14, -2, -2, -1)},
      {UnitCellCoord(14, -2, -1, 0)},
      {UnitCellCoord(14, -2, 0, 0)},
      {UnitCellCoord(13, -1, -2, -1)},
      {UnitCellCoord(12, -1, -1, -2)},
      {UnitCellCoord(12, -1, 0, -2)},
      {UnitCellCoord(14, -1, 0, 1)},
      {UnitCellCoord(12, -1, 1, -1)},
      {UnitCellCoord(14, -1, 1, 0)},
      {UnitCellCoord(13, 0, -2, 0)},
      {UnitCellCoord(12, 0, 0, -2)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 2)},
      {UnitCellCoord(15, 0, 2, 0)},
      {UnitCellCoord(13, 1, -1, 0)},
      {UnitCellCoord(13, 1, -1, 1)},
      {UnitCellCoord(12, 1, 0, -1)},
      {UnitCellCoord(13, 1, 0, 2)},
      {UnitCellCoord(15, 1, 1, 2)},
      {UnitCellCoord(15, 1, 2, 1)},
      {UnitCellCoord(13, 2, 0, 0)},
      {UnitCellCoord(15, 2, 1, 0)},
      {UnitCellCoord(15, 2, 2, 1)},
      {UnitCellCoord(15, 2, 2, 2)}
    };

    m_orbit_neighborhood[12] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -2, -2)},
      {UnitCellCoord(14, -2, -2, -1)},
      {UnitCellCoord(12, -2, -1, -2)},
      {UnitCellCoord(12, -2, -1, 0)},
      {UnitCellCoord(15, -2, 0, -1)},
      {UnitCellCoord(12, -2, 0, 0)},
      {UnitCellCoord(12, -1, -2, 0)},
      {UnitCellCoord(13, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, 1)},
      {UnitCellCoord(15, -1, 0, -2)},
      {UnitCellCoord(13, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, -1, 1, -1)},
      {UnitCellCoord(14, 0, -2, 0)},
      {UnitCellCoord(12, 0, -1, -1)},
      {UnitCellCoord(13, 0, -1, -1)},
      {UnitCellCoord(14, 0, -1, -1)},
      {UnitCellCoord(15, 0, -1, -1)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(15, 0, -1, 0)},
      {UnitCellCoord(15, 0, 0, -2)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(14, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(15, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 2)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(14, 0, 1, 0)},
      {UnitCellCoord(12, 0, 1, 1)},
      {UnitCellCoord(13, 0, 1, 1)},
      {UnitCellCoord(14, 0, 1, 1)},
      {UnitCellCoord(15, 0, 1, 1)},
      {UnitCellCoord(15, 0, 2, 0)},
      {UnitCellCoord(14, 1, -1, 1)},
      {UnitCellCoord(12, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 2)},
      {UnitCellCoord(13, 1, 1, -1)},
      {UnitCellCoord(12, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(15, 1, 2, 0)},
      {UnitCellCoord(13, 2, 0, 0)},
      {UnitCellCoord(13, 2, 0, 1)},
      {UnitCellCoord(15, 2, 1, 0)},
      {UnitCellCoord(14, 2, 1, 2)},
      {UnitCellCoord(13, 2, 2, 1)},
      {UnitCellCoord(13, 2, 2, 2)}
    };

    m_orbit_neighborhood[13] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -2, -1)},
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(15, -2, -1, -1)},
      {UnitCellCoord(12, -2, -1, 0)},
      {UnitCellCoord(14, -1, -2, -1)},
      {UnitCellCoord(15, -1, -1, -2)},
      {UnitCellCoord(15, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(15, -1, -1, 0)},
      {UnitCellCoord(15, -1, 0, -2)},
      {UnitCellCoord(13, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, -1, 0, 1)},
      {UnitCellCoord(15, -1, 1, -1)},
      {UnitCellCoord(12, -1, 1, 0)},
      {UnitCellCoord(14, 0, -1, -1)},
      {UnitCellCoord(15, 0, -1, -1)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(12, 0, -1, 1)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 0, 1)},
      {UnitCellCoord(13, 0, 1, -1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(13, 0, 1, 0)},
      {UnitCellCoord(12, 0, 1, 1)},
      {UnitCellCoord(13, 0, 1, 1)},
      {UnitCellCoord(14, 1, -1, 0)},
      {UnitCellCoord(14, 1, -1, 1)},
      {UnitCellCoord(15, 1, 0, -1)},
      {UnitCellCoord(14, 1, 0, 0)},
      {UnitCellCoord(12, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 2)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(14, 1, 1, 0)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(13, 1, 1, 2)},
      {UnitCellCoord(13, 1, 2, 1)},
      {UnitCellCoord(13, 2, 1, 0)},
      {UnitCellCoord(13, 2, 1, 1)},
      {UnitCellCoord(14, 2, 1, 1)},
      {UnitCellCoord(13, 2, 2, 1)}
    };

    m_orbit_neighborhood[14] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(13, -2, -1, -1)},
      {UnitCellCoord(12, -1, -2, -1)},
      {UnitCellCoord(14, -1, -2, -1)},
      {UnitCellCoord(13, -1, -1, -2)},
      {UnitCellCoord(14, -1, -1, -2)},
      {UnitCellCoord(13, -1, 0, 1)},
      {UnitCellCoord(15, -1, 0, 1)},
      {UnitCellCoord(12, -1, 1, 0)},
      {UnitCellCoord(15, -1, 1, 0)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(15, 0, -1, 1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 1, -1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 1, -1, 0)},
      {UnitCellCoord(15, 1, -1, 0)},
      {UnitCellCoord(13, 1, 0, -1)},
      {UnitCellCoord(15, 1, 0, -1)},
      {UnitCellCoord(13, 1, 1, 2)},
      {UnitCellCoord(14, 1, 1, 2)},
      {UnitCellCoord(12, 1, 2, 1)},
      {UnitCellCoord(14, 1, 2, 1)},
      {UnitCellCoord(12, 2, 1, 1)},
      {UnitCellCoord(13, 2, 1, 1)}
    };

    m_orbit_neighborhood[15] = std::set<UnitCellCoord> {
      {UnitCellCoord(14, -2, -1, -1)},
      {UnitCellCoord(15, -2, -1, -1)},
      {UnitCellCoord(13, -1, -2, -1)},
      {UnitCellCoord(15, -1, -2, -1)},
      {UnitCellCoord(12, -1, -1, -2)},
      {UnitCellCoord(15, -1, -1, -2)},
      {UnitCellCoord(12, -1, 0, 1)},
      {UnitCellCoord(14, -1, 0, 1)},
      {UnitCellCoord(13, -1, 1, 0)},
      {UnitCellCoord(14, -1, 1, 0)},
      {UnitCellCoord(12, 0, -1, 1)},
      {UnitCellCoord(13, 0, -1, 1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(12, 0, 1, -1)},
      {UnitCellCoord(13, 0, 1, -1)},
      {UnitCellCoord(13, 1, -1, 0)},
      {UnitCellCoord(14, 1, -1, 0)},
      {UnitCellCoord(12, 1, 0, -1)},
      {UnitCellCoord(14, 1, 0, -1)},
      {UnitCellCoord(12, 1, 1, 2)},
      {UnitCellCoord(15, 1, 1, 2)},
      {UnitCellCoord(13, 1, 2, 1)},
      {UnitCellCoord(15, 1, 2, 1)},
      {UnitCellCoord(14, 2, 1, 1)},
      {UnitCellCoord(15, 2, 1, 1)}
    };

    m_orbit_neighborhood[16] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

    m_orbit_neighborhood[17] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

    m_orbit_neighborhood[18] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(15, 1, 1, 1)},
      {UnitCellCoord(13, 2, 1, 1)}
    };

    m_orbit_neighborhood[19] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -2, -1, -1)},
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(15, -1, 0, 0)},
      {UnitCellCoord(12, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 1)},
      {UnitCellCoord(13, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 1)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, -1)},
      {UnitCellCoord(12, 0, 1, 0)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(15, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(14, 1, 1, 1)},
      {UnitCellCoord(13, 2, 1, 1)}
    };

    m_orbit_neighborhood[20] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

    m_orbit_neighborhood[21] = std::set<UnitCellCoord> {
      {UnitCellCoord(12, -1, -1, -1)},
      {UnitCellCoord(12, -1, -1, 0)},
      {UnitCellCoord(14, -1, -1, 0)},
      {UnitCellCoord(12, -1, 0, -1)},
      {UnitCellCoord(15, -1, 0, -1)},
      {UnitCellCoord(12, -1, 0, 0)},
      {UnitCellCoord(14, -1, 0, 0)},
      {UnitCellCoord(13, 0, -1, 0)},
      {UnitCellCoord(14, 0, -1, 0)},
      {UnitCellCoord(12, 0, 0, -1)},
      {UnitCellCoord(15, 0, 0, -1)},
      {UnitCellCoord(12, 0, 0, 0)},
      {UnitCellCoord(13, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 0)},
      {UnitCellCoord(15, 0, 0, 0)},
      {UnitCellCoord(14, 0, 0, 1)},
      {UnitCellCoord(15, 0, 1, 0)},
      {UnitCellCoord(13, 1, 0, 0)},
      {UnitCellCoord(13, 1, 0, 1)},
      {UnitCellCoord(14, 1, 0, 1)},
      {UnitCellCoord(13, 1, 1, 0)},
      {UnitCellCoord(15, 1, 1, 0)},
      {UnitCellCoord(13, 1, 1, 1)},
      {UnitCellCoord(15, 1, 1, 1)}
    };

  }

  Nb_direction_Clexulator::~Nb_direction_Clexulator(){
    //nothing here for now
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void Nb_direction_Clexulator::calc_global_corr_contribution(double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_orbit_func_list[i])();
    }
  }

  /// \brief Calculate contribution to select global correlations from one unit cell
  void Nb_direction_Clexulator::calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_orbit_func_list[*ind_list_begin])();
    }
  }

  /// \brief Calculate point correlations about basis site 'b_index'
  void Nb_direction_Clexulator::calc_point_corr(int b_index, double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_flower_func_lists[b_index][i])();
    }
  }

  /// \brief Calculate select point correlations about basis site 'b_index'
  void Nb_direction_Clexulator::calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_flower_func_lists[b_index][*ind_list_begin])();
    }
  }

  /// \brief Calculate the change in point correlations due to changing an occupant
  void Nb_direction_Clexulator::calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_delta_func_lists[b_index][i])(occ_i, occ_f);
    }
  }

  /// \brief Calculate the change in select point correlations due to changing an occupant
  void Nb_direction_Clexulator::calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_delta_func_lists[b_index][*ind_list_begin])(occ_i, occ_f);
    }
  }

  // Basis functions for empty cluster:
  double Nb_direction_Clexulator::eval_bfunc_0_0_0() const{
    return (1);
  }

  /**** Basis functions for orbit 1, 0****
#Points: 1
MaxLength: 0  MinLength: 0
               0.0000000    0.0000000    0.7194200 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_1_0_0() const{
    return ((occ_func_12_0(0)) + (occ_func_12_0(1)) + (occ_func_12_0(2)) + (occ_func_12_0(3)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_1_0_0() const{
    return ((occ_func_12_0(0)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_1_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((1))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_1_0_0() const{
    return ((occ_func_12_0(1)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_1_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((1))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_1_0_0() const{
    return ((occ_func_12_0(2)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_1_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((1))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_1_0_0() const{
    return ((occ_func_12_0(3)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_1_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((1))/4.0;
  }

  /**** Basis functions for orbit 2, 0****
#Points: 2
MaxLength: 3.3499115  MinLength: 3.3499115
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -0.0000000    1.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_0_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(22)) + (occ_func_12_0(1)*occ_func_14_0(10)) + (occ_func_12_0(0)*occ_func_14_0(53)) + (occ_func_12_0(2)*occ_func_14_0(59)) + (occ_func_12_0(0)*occ_func_14_0(35)) + (occ_func_12_0(1)*occ_func_14_0(27)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_0_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(22)) + (occ_func_12_0(0)*occ_func_14_0(53)) + (occ_func_12_0(0)*occ_func_14_0(35)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(22)) + (occ_func_14_0(53)) + (occ_func_14_0(35)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_0_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(10)) + (occ_func_12_0(40)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(27)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(10)) + (occ_func_12_0(40)) + (occ_func_14_0(27)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_0_0() const{
    return ((occ_func_12_0(16)*occ_func_14_0(2)) + (occ_func_12_0(29)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(59)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(16)) + (occ_func_12_0(29)) + (occ_func_14_0(59)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_0_0() const{
    return ((occ_func_12_0(38)*occ_func_14_0(3)) + (occ_func_12_0(4)*occ_func_14_0(3)) + (occ_func_12_0(13)*occ_func_14_0(3)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(38)) + (occ_func_12_0(4)) + (occ_func_12_0(13)))/6.0;
  }

  /**** Basis functions for orbit 2, 1****
#Points: 2
MaxLength: 4.9680021  MinLength: 4.9680021
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800    0.0000000    0.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)) + (occ_func_12_0(1)*occ_func_14_0(38)) + (occ_func_12_0(0)*occ_func_14_0(29)) + (occ_func_12_0(2)*occ_func_14_0(27)) + (occ_func_12_0(2)*occ_func_14_0(1)) + (occ_func_12_0(3)*occ_func_14_0(22)) + (occ_func_12_0(0)*occ_func_14_0(59)) + (occ_func_12_0(3)*occ_func_14_0(0)) + (occ_func_12_0(1)*occ_func_14_0(4)) + (occ_func_12_0(1)*occ_func_14_0(3)) + (occ_func_12_0(2)*occ_func_14_0(40)) + (occ_func_12_0(3)*occ_func_14_0(53)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)) + (occ_func_12_0(0)*occ_func_14_0(29)) + (occ_func_12_0(0)*occ_func_14_0(59)) + (occ_func_12_0(3)*occ_func_14_0(0)) + (occ_func_12_0(33)*occ_func_14_0(0)) + (occ_func_12_0(54)*occ_func_14_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(2)) + (occ_func_14_0(29)) + (occ_func_14_0(59)) + (occ_func_12_0(3)) + (occ_func_12_0(33)) + (occ_func_12_0(54)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_1_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(38)) + (occ_func_12_0(8)*occ_func_14_0(1)) + (occ_func_12_0(2)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(4)) + (occ_func_12_0(1)*occ_func_14_0(3)) + (occ_func_12_0(43)*occ_func_14_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(38)) + (occ_func_12_0(8)) + (occ_func_12_0(2)) + (occ_func_14_0(4)) + (occ_func_14_0(3)) + (occ_func_12_0(43)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)) + (occ_func_12_0(57)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(27)) + (occ_func_12_0(2)*occ_func_14_0(1)) + (occ_func_12_0(19)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(40)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(0)) + (occ_func_12_0(57)) + (occ_func_14_0(27)) + (occ_func_14_0(1)) + (occ_func_12_0(19)) + (occ_func_14_0(40)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_1_0() const{
    return ((occ_func_12_0(14)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(22)) + (occ_func_12_0(36)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(0)) + (occ_func_12_0(1)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(53)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(14)) + (occ_func_14_0(22)) + (occ_func_12_0(36)) + (occ_func_14_0(0)) + (occ_func_12_0(1)) + (occ_func_14_0(53)))/12.0;
  }

  /**** Basis functions for orbit 2, 2****
#Points: 2
MaxLength: 6.5191431  MinLength: 6.5191431
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -1.0000000    0.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_2_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)) + (occ_func_12_0(1)*occ_func_14_0(22)) + (occ_func_12_0(0)*occ_func_14_0(105)) + (occ_func_12_0(2)*occ_func_14_0(3)) + (occ_func_12_0(0)*occ_func_14_0(27)) + (occ_func_12_0(1)*occ_func_14_0(19)) + (occ_func_12_0(2)*occ_func_14_0(4)) + (occ_func_12_0(0)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(6)) + (occ_func_12_0(3)*occ_func_14_0(82)) + (occ_func_12_0(0)*occ_func_14_0(31)) + (occ_func_12_0(1)*occ_func_14_0(11)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_2_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)) + (occ_func_12_0(0)*occ_func_14_0(105)) + (occ_func_12_0(0)*occ_func_14_0(27)) + (occ_func_12_0(34)*occ_func_14_0(0)) + (occ_func_12_0(0)*occ_func_14_0(1)) + (occ_func_12_0(0)*occ_func_14_0(31)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(14)) + (occ_func_14_0(105)) + (occ_func_14_0(27)) + (occ_func_12_0(34)) + (occ_func_14_0(1)) + (occ_func_14_0(31)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_2_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(22)) + (occ_func_12_0(60)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(19)) + (occ_func_12_0(0)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(6)) + (occ_func_12_0(1)*occ_func_14_0(11)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(22)) + (occ_func_12_0(60)) + (occ_func_14_0(19)) + (occ_func_12_0(0)) + (occ_func_14_0(6)) + (occ_func_14_0(11)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_2_0() const{
    return ((occ_func_12_0(24)*occ_func_14_0(2)) + (occ_func_12_0(17)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(3)) + (occ_func_12_0(2)*occ_func_14_0(4)) + (occ_func_12_0(33)*occ_func_14_0(2)) + (occ_func_12_0(87)*occ_func_14_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(24)) + (occ_func_12_0(17)) + (occ_func_14_0(3)) + (occ_func_14_0(4)) + (occ_func_12_0(33)) + (occ_func_12_0(87)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_2_0() const{
    return ((occ_func_12_0(2)*occ_func_14_0(3)) + (occ_func_12_0(12)*occ_func_14_0(3)) + (occ_func_12_0(21)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(82)) + (occ_func_12_0(8)*occ_func_14_0(3)) + (occ_func_12_0(29)*occ_func_14_0(3)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(2)) + (occ_func_12_0(12)) + (occ_func_12_0(21)) + (occ_func_14_0(82)) + (occ_func_12_0(8)) + (occ_func_12_0(29)))/12.0;
  }

  /**** Basis functions for orbit 2, 3****
#Points: 2
MaxLength: 7.3112607  MinLength: 7.3112607
               0.0000000    0.0000000    0.7194200 Nb Ta
               0.0000000    0.0000000   -0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_3_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(16)) + (occ_func_12_0(1)*occ_func_12_0(13)) + (occ_func_12_0(2)*occ_func_12_0(10)) + (occ_func_12_0(3)*occ_func_12_0(35)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_3_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(16)) + (occ_func_12_0(20)*occ_func_12_0(0)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_12_0(16)) + (occ_func_12_0(20)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_3_0() const{
    return ((occ_func_12_0(1)*occ_func_12_0(13)) + (occ_func_12_0(25)*occ_func_12_0(1)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(13)) + (occ_func_12_0(25)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_3_0() const{
    return ((occ_func_12_0(2)*occ_func_12_0(10)) + (occ_func_12_0(30)*occ_func_12_0(2)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(10)) + (occ_func_12_0(30)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_3_0() const{
    return ((occ_func_12_0(3)*occ_func_12_0(35)) + (occ_func_12_0(7)*occ_func_12_0(3)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(35)) + (occ_func_12_0(7)))/4.0;
  }

  /**** Basis functions for orbit 2, 4****
#Points: 2
MaxLength: 7.3112607  MinLength: 7.3112607
               0.0000000    0.0000000    0.7194200 Nb Ta
              -1.0000000   -1.0000000   -0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_4_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(4)) + (occ_func_12_0(1)*occ_func_12_0(21)) + (occ_func_12_0(0)*occ_func_12_0(28)) + (occ_func_12_0(2)*occ_func_12_0(22)) + (occ_func_12_0(2)*occ_func_12_0(6)) + (occ_func_12_0(3)*occ_func_12_0(27)) + (occ_func_12_0(0)*occ_func_12_0(24)) + (occ_func_12_0(3)*occ_func_12_0(31)) + (occ_func_12_0(1)*occ_func_12_0(5)) + (occ_func_12_0(1)*occ_func_12_0(29)) + (occ_func_12_0(2)*occ_func_12_0(26)) + (occ_func_12_0(3)*occ_func_12_0(23)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_4_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(4)) + (occ_func_12_0(32)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(28)) + (occ_func_12_0(8)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(24)) + (occ_func_12_0(12)*occ_func_12_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_12_0(4)) + (occ_func_12_0(32)) + (occ_func_12_0(28)) + (occ_func_12_0(8)) + (occ_func_12_0(24)) + (occ_func_12_0(12)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_4_0() const{
    return ((occ_func_12_0(1)*occ_func_12_0(21)) + (occ_func_12_0(17)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(5)) + (occ_func_12_0(33)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(29)) + (occ_func_12_0(9)*occ_func_12_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(21)) + (occ_func_12_0(17)) + (occ_func_12_0(5)) + (occ_func_12_0(33)) + (occ_func_12_0(29)) + (occ_func_12_0(9)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_4_0() const{
    return ((occ_func_12_0(2)*occ_func_12_0(22)) + (occ_func_12_0(18)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(6)) + (occ_func_12_0(34)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(26)) + (occ_func_12_0(14)*occ_func_12_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(22)) + (occ_func_12_0(18)) + (occ_func_12_0(6)) + (occ_func_12_0(34)) + (occ_func_12_0(26)) + (occ_func_12_0(14)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_4_0() const{
    return ((occ_func_12_0(3)*occ_func_12_0(27)) + (occ_func_12_0(15)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(31)) + (occ_func_12_0(11)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(23)) + (occ_func_12_0(19)*occ_func_12_0(3)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(27)) + (occ_func_12_0(15)) + (occ_func_12_0(31)) + (occ_func_12_0(11)) + (occ_func_12_0(23)) + (occ_func_12_0(19)))/12.0;
  }

  /**** Basis functions for orbit 2, 5****
#Points: 2
MaxLength: 8.0421660  MinLength: 8.0421660
               0.0000000    0.0000000    0.7194200 Nb Ta
              -1.2805800   -1.0000000    0.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_5_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(38)) + (occ_func_12_0(1)*occ_func_14_0(74)) + (occ_func_12_0(0)*occ_func_14_0(189)) + (occ_func_12_0(2)*occ_func_14_0(35)) + (occ_func_12_0(2)*occ_func_14_0(45)) + (occ_func_12_0(3)*occ_func_14_0(10)) + (occ_func_12_0(0)*occ_func_14_0(103)) + (occ_func_12_0(3)*occ_func_14_0(44)) + (occ_func_12_0(1)*occ_func_14_0(112)) + (occ_func_12_0(1)*occ_func_14_0(59)) + (occ_func_12_0(2)*occ_func_14_0(84)) + (occ_func_12_0(3)*occ_func_14_0(81)) + (occ_func_12_0(2)*occ_func_14_0(68)) + (occ_func_12_0(0)*occ_func_14_0(13)) + (occ_func_12_0(1)*occ_func_14_0(62)) + (occ_func_12_0(3)*occ_func_14_0(134)) + (occ_func_12_0(0)*occ_func_14_0(50)) + (occ_func_12_0(0)*occ_func_14_0(107)) + (occ_func_12_0(3)*occ_func_14_0(40)) + (occ_func_12_0(1)*occ_func_14_0(51)) + (occ_func_12_0(1)*occ_func_14_0(16)) + (occ_func_12_0(2)*occ_func_14_0(53)) + (occ_func_12_0(3)*occ_func_14_0(89)) + (occ_func_12_0(2)*occ_func_14_0(183)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_5_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(38)) + (occ_func_12_0(0)*occ_func_14_0(189)) + (occ_func_12_0(0)*occ_func_14_0(103)) + (occ_func_12_0(51)*occ_func_14_0(0)) + (occ_func_12_0(197)*occ_func_14_0(0)) + (occ_func_12_0(82)*occ_func_14_0(0)) + (occ_func_12_0(98)*occ_func_14_0(0)) + (occ_func_12_0(0)*occ_func_14_0(13)) + (occ_func_12_0(0)*occ_func_14_0(50)) + (occ_func_12_0(0)*occ_func_14_0(107)) + (occ_func_12_0(55)*occ_func_14_0(0)) + (occ_func_12_0(21)*occ_func_14_0(0)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(38)) + (occ_func_14_0(189)) + (occ_func_14_0(103)) + (occ_func_12_0(51)) + (occ_func_12_0(197)) + (occ_func_12_0(82)) + (occ_func_12_0(98)) + (occ_func_14_0(13)) + (occ_func_14_0(50)) + (occ_func_14_0(107)) + (occ_func_12_0(55)) + (occ_func_12_0(21)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_5_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(74)) + (occ_func_12_0(120)*occ_func_14_0(1)) + (occ_func_12_0(50)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(112)) + (occ_func_12_0(1)*occ_func_14_0(59)) + (occ_func_12_0(87)*occ_func_14_0(1)) + (occ_func_12_0(24)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(62)) + (occ_func_12_0(1)*occ_func_14_0(51)) + (occ_func_12_0(1)*occ_func_14_0(16)) + (occ_func_12_0(42)*occ_func_14_0(1)) + (occ_func_12_0(79)*occ_func_14_0(1)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(74)) + (occ_func_12_0(120)) + (occ_func_12_0(50)) + (occ_func_14_0(112)) + (occ_func_14_0(59)) + (occ_func_12_0(87)) + (occ_func_12_0(24)) + (occ_func_14_0(62)) + (occ_func_14_0(51)) + (occ_func_14_0(16)) + (occ_func_12_0(42)) + (occ_func_12_0(79)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_5_0() const{
    return ((occ_func_12_0(56)*occ_func_14_0(2)) + (occ_func_12_0(93)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(35)) + (occ_func_12_0(2)*occ_func_14_0(45)) + (occ_func_12_0(31)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(84)) + (occ_func_12_0(2)*occ_func_14_0(68)) + (occ_func_12_0(105)*occ_func_14_0(2)) + (occ_func_12_0(179)*occ_func_14_0(2)) + (occ_func_12_0(44)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(53)) + (occ_func_12_0(2)*occ_func_14_0(183)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(56)) + (occ_func_12_0(93)) + (occ_func_14_0(35)) + (occ_func_14_0(45)) + (occ_func_12_0(31)) + (occ_func_14_0(84)) + (occ_func_14_0(68)) + (occ_func_12_0(105)) + (occ_func_12_0(179)) + (occ_func_12_0(44)) + (occ_func_14_0(53)) + (occ_func_14_0(183)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_5_0() const{
    return ((occ_func_12_0(6)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(10)) + (occ_func_12_0(64)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(44)) + (occ_func_12_0(37)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(81)) + (occ_func_12_0(3)*occ_func_14_0(134)) + (occ_func_12_0(60)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(40)) + (occ_func_12_0(45)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(89)) + (occ_func_12_0(130)*occ_func_14_0(3)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(6)) + (occ_func_14_0(10)) + (occ_func_12_0(64)) + (occ_func_14_0(44)) + (occ_func_12_0(37)) + (occ_func_14_0(81)) + (occ_func_14_0(134)) + (occ_func_12_0(60)) + (occ_func_14_0(40)) + (occ_func_12_0(45)) + (occ_func_14_0(89)) + (occ_func_12_0(130)))/24.0;
  }

  /**** Basis functions for orbit 2, 6****
#Points: 2
MaxLength: 8.4423170  MinLength: 8.4423170
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.0000000   -1.0000000   -0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_6_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(44)) + (occ_func_12_0(1)*occ_func_12_0(53)) + (occ_func_12_0(0)*occ_func_12_0(56)) + (occ_func_12_0(2)*occ_func_12_0(38)) + (occ_func_12_0(2)*occ_func_12_0(42)) + (occ_func_12_0(3)*occ_func_12_0(59)) + (occ_func_12_0(0)*occ_func_12_0(40)) + (occ_func_12_0(3)*occ_func_12_0(55)) + (occ_func_12_0(1)*occ_func_12_0(37)) + (occ_func_12_0(1)*occ_func_12_0(45)) + (occ_func_12_0(2)*occ_func_12_0(50)) + (occ_func_12_0(3)*occ_func_12_0(51)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_6_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(44)) + (occ_func_12_0(48)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(56)) + (occ_func_12_0(36)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(40)) + (occ_func_12_0(52)*occ_func_12_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_12_0(44)) + (occ_func_12_0(48)) + (occ_func_12_0(56)) + (occ_func_12_0(36)) + (occ_func_12_0(40)) + (occ_func_12_0(52)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_6_0() const{
    return ((occ_func_12_0(1)*occ_func_12_0(53)) + (occ_func_12_0(41)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(37)) + (occ_func_12_0(57)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(45)) + (occ_func_12_0(49)*occ_func_12_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(53)) + (occ_func_12_0(41)) + (occ_func_12_0(37)) + (occ_func_12_0(57)) + (occ_func_12_0(45)) + (occ_func_12_0(49)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_6_0() const{
    return ((occ_func_12_0(2)*occ_func_12_0(38)) + (occ_func_12_0(58)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(42)) + (occ_func_12_0(54)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(50)) + (occ_func_12_0(46)*occ_func_12_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(38)) + (occ_func_12_0(58)) + (occ_func_12_0(42)) + (occ_func_12_0(54)) + (occ_func_12_0(50)) + (occ_func_12_0(46)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_6_0() const{
    return ((occ_func_12_0(3)*occ_func_12_0(59)) + (occ_func_12_0(39)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(55)) + (occ_func_12_0(43)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(51)) + (occ_func_12_0(47)*occ_func_12_0(3)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(59)) + (occ_func_12_0(39)) + (occ_func_12_0(55)) + (occ_func_12_0(43)) + (occ_func_12_0(51)) + (occ_func_12_0(47)))/12.0;
  }

  /**** Basis functions for orbit 2, 7****
#Points: 2
MaxLength: 8.5893275  MinLength: 8.5893275
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.7194200   -0.7194200   -0.7194200 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_7_0() const{
    return ((occ_func_12_0(0)*occ_func_15_0(19)) + (occ_func_12_0(1)*occ_func_15_0(36)) + (occ_func_12_0(0)*occ_func_15_0(30)) + (occ_func_12_0(2)*occ_func_15_0(43)) + (occ_func_12_0(3)*occ_func_15_0(33)) + (occ_func_12_0(1)*occ_func_15_0(14)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_7_0() const{
    return ((occ_func_12_0(0)*occ_func_15_0(19)) + (occ_func_12_0(57)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_15_0(30)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_15_0(19)) + (occ_func_12_0(57)) + (occ_func_15_0(30)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_7_0() const{
    return ((occ_func_12_0(1)*occ_func_15_0(36)) + (occ_func_12_0(7)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_15_0(14)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_15_0(36)) + (occ_func_12_0(7)) + (occ_func_15_0(14)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_7_0() const{
    return ((occ_func_12_0(8)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_15_0(43)) + (occ_func_12_0(25)*occ_func_15_0(2)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(8)) + (occ_func_15_0(43)) + (occ_func_12_0(25)))/6.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_7_0() const{
    return ((occ_func_12_0(20)*occ_func_15_0(3)) + (occ_func_12_0(54)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_15_0(33)))/6.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(20)) + (occ_func_12_0(54)) + (occ_func_15_0(33)))/6.0;
  }

  /**** Basis functions for orbit 2, 8****
#Points: 2
MaxLength: 9.0826548  MinLength: 9.0826548
               0.0000000    0.0000000    0.7194200 Nb Ta
               0.0000000   -1.2805800    0.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_8_0() const{
    return ((occ_func_12_0(0)*occ_func_13_0(89)) + (occ_func_12_0(1)*occ_func_13_0(35)) + (occ_func_12_0(0)*occ_func_13_0(203)) + (occ_func_12_0(2)*occ_func_13_0(13)) + (occ_func_12_0(2)*occ_func_13_0(136)) + (occ_func_12_0(3)*occ_func_13_0(16)) + (occ_func_12_0(0)*occ_func_13_0(10)) + (occ_func_12_0(3)*occ_func_13_0(169)) + (occ_func_12_0(1)*occ_func_13_0(118)) + (occ_func_12_0(1)*occ_func_13_0(68)) + (occ_func_12_0(2)*occ_func_13_0(103)) + (occ_func_12_0(3)*occ_func_13_0(74)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_8_0() const{
    return ((occ_func_12_0(0)*occ_func_13_0(89)) + (occ_func_12_0(0)*occ_func_13_0(203)) + (occ_func_12_0(174)*occ_func_13_0(0)) + (occ_func_12_0(23)*occ_func_13_0(0)) + (occ_func_12_0(0)*occ_func_13_0(10)) + (occ_func_12_0(97)*occ_func_13_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_13_0(89)) + (occ_func_13_0(203)) + (occ_func_12_0(174)) + (occ_func_12_0(23)) + (occ_func_13_0(10)) + (occ_func_12_0(97)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_8_0() const{
    return ((occ_func_12_0(76)*occ_func_13_0(1)) + (occ_func_12_0(1)*occ_func_13_0(35)) + (occ_func_12_0(26)*occ_func_13_0(1)) + (occ_func_12_0(143)*occ_func_13_0(1)) + (occ_func_12_0(1)*occ_func_13_0(118)) + (occ_func_12_0(1)*occ_func_13_0(68)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(76)) + (occ_func_13_0(35)) + (occ_func_12_0(26)) + (occ_func_12_0(143)) + (occ_func_13_0(118)) + (occ_func_13_0(68)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_8_0() const{
    return ((occ_func_12_0(2)*occ_func_13_0(13)) + (occ_func_12_0(2)*occ_func_13_0(136)) + (occ_func_12_0(28)*occ_func_13_0(2)) + (occ_func_12_0(193)*occ_func_13_0(2)) + (occ_func_12_0(2)*occ_func_13_0(103)) + (occ_func_12_0(95)*occ_func_13_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_13_0(13)) + (occ_func_13_0(136)) + (occ_func_12_0(28)) + (occ_func_12_0(193)) + (occ_func_13_0(103)) + (occ_func_12_0(95)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_8_0() const{
    return ((occ_func_12_0(5)*occ_func_13_0(3)) + (occ_func_12_0(108)*occ_func_13_0(3)) + (occ_func_12_0(3)*occ_func_13_0(16)) + (occ_func_12_0(3)*occ_func_13_0(169)) + (occ_func_12_0(66)*occ_func_13_0(3)) + (occ_func_12_0(3)*occ_func_13_0(74)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(5)) + (occ_func_12_0(108)) + (occ_func_13_0(16)) + (occ_func_13_0(169)) + (occ_func_12_0(66)) + (occ_func_13_0(74)))/12.0;
  }

  /**** Basis functions for orbit 2, 9****
#Points: 2
MaxLength: 10.2309256  MinLength: 10.2309256
               0.0000000    0.0000000    0.7194200 Nb Ta
               1.2805800    1.2805800    2.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_9_0() const{
    return ((occ_func_12_0(0)*occ_func_15_0(235)) + (occ_func_12_0(1)*occ_func_15_0(136)) + (occ_func_12_0(0)*occ_func_15_0(74)) + (occ_func_12_0(2)*occ_func_15_0(216)) + (occ_func_12_0(2)*occ_func_15_0(203)) + (occ_func_12_0(3)*occ_func_15_0(213)) + (occ_func_12_0(0)*occ_func_15_0(169)) + (occ_func_12_0(3)*occ_func_15_0(118)) + (occ_func_12_0(1)*occ_func_15_0(103)) + (occ_func_12_0(1)*occ_func_15_0(210)) + (occ_func_12_0(2)*occ_func_15_0(89)) + (occ_func_12_0(3)*occ_func_15_0(68)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_9_0() const{
    return ((occ_func_12_0(0)*occ_func_15_0(235)) + (occ_func_12_0(173)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_15_0(74)) + (occ_func_12_0(222)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_15_0(169)) + (occ_func_12_0(99)*occ_func_15_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_15_0(235)) + (occ_func_12_0(173)) + (occ_func_15_0(74)) + (occ_func_12_0(222)) + (occ_func_15_0(169)) + (occ_func_12_0(99)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_9_0() const{
    return ((occ_func_12_0(1)*occ_func_15_0(136)) + (occ_func_12_0(227)*occ_func_15_0(1)) + (occ_func_12_0(140)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_15_0(103)) + (occ_func_12_0(1)*occ_func_15_0(210)) + (occ_func_12_0(78)*occ_func_15_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_15_0(136)) + (occ_func_12_0(227)) + (occ_func_12_0(140)) + (occ_func_15_0(103)) + (occ_func_15_0(210)) + (occ_func_12_0(78)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_9_0() const{
    return ((occ_func_12_0(92)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_15_0(216)) + (occ_func_12_0(2)*occ_func_15_0(203)) + (occ_func_12_0(195)*occ_func_15_0(2)) + (occ_func_12_0(229)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_15_0(89)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(92)) + (occ_func_15_0(216)) + (occ_func_15_0(203)) + (occ_func_12_0(195)) + (occ_func_12_0(229)) + (occ_func_15_0(89)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_9_0() const{
    return ((occ_func_12_0(204)*occ_func_15_0(3)) + (occ_func_12_0(110)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_15_0(213)) + (occ_func_12_0(3)*occ_func_15_0(118)) + (occ_func_12_0(65)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_15_0(68)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(204)) + (occ_func_12_0(110)) + (occ_func_15_0(213)) + (occ_func_15_0(118)) + (occ_func_12_0(65)) + (occ_func_15_0(68)))/12.0;
  }

  /**** Basis functions for orbit 2, 10****
#Points: 2
MaxLength: 10.6663929  MinLength: 10.6663929
               0.0000000    0.0000000    0.7194200 Nb Ta
              -1.2805800   -1.0000000   -1.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_10_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(6)) + (occ_func_12_0(1)*occ_func_14_0(134)) + (occ_func_12_0(0)*occ_func_14_0(229)) + (occ_func_12_0(2)*occ_func_14_0(51)) + (occ_func_12_0(2)*occ_func_14_0(5)) + (occ_func_12_0(3)*occ_func_14_0(50)) + (occ_func_12_0(0)*occ_func_14_0(183)) + (occ_func_12_0(3)*occ_func_14_0(28)) + (occ_func_12_0(1)*occ_func_14_0(204)) + (occ_func_12_0(1)*occ_func_14_0(31)) + (occ_func_12_0(2)*occ_func_14_0(140)) + (occ_func_12_0(3)*occ_func_14_0(173)) + (occ_func_12_0(2)*occ_func_14_0(112)) + (occ_func_12_0(0)*occ_func_14_0(45)) + (occ_func_12_0(1)*occ_func_14_0(110)) + (occ_func_12_0(3)*occ_func_14_0(222)) + (occ_func_12_0(0)*occ_func_14_0(26)) + (occ_func_12_0(0)*occ_func_14_0(195)) + (occ_func_12_0(3)*occ_func_14_0(24)) + (occ_func_12_0(1)*occ_func_14_0(23)) + (occ_func_12_0(1)*occ_func_14_0(44)) + (occ_func_12_0(2)*occ_func_14_0(21)) + (occ_func_12_0(3)*occ_func_14_0(189)) + (occ_func_12_0(2)*occ_func_14_0(227)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_10_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(6)) + (occ_func_12_0(0)*occ_func_14_0(229)) + (occ_func_12_0(0)*occ_func_14_0(183)) + (occ_func_12_0(11)*occ_func_14_0(0)) + (occ_func_12_0(233)*occ_func_14_0(0)) + (occ_func_12_0(170)*occ_func_14_0(0)) + (occ_func_12_0(198)*occ_func_14_0(0)) + (occ_func_12_0(0)*occ_func_14_0(45)) + (occ_func_12_0(0)*occ_func_14_0(26)) + (occ_func_12_0(0)*occ_func_14_0(195)) + (occ_func_12_0(15)*occ_func_14_0(0)) + (occ_func_12_0(49)*occ_func_14_0(0)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(6)) + (occ_func_14_0(229)) + (occ_func_14_0(183)) + (occ_func_12_0(11)) + (occ_func_12_0(233)) + (occ_func_12_0(170)) + (occ_func_12_0(198)) + (occ_func_14_0(45)) + (occ_func_14_0(26)) + (occ_func_14_0(195)) + (occ_func_12_0(15)) + (occ_func_12_0(49)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_10_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(134)) + (occ_func_12_0(208)*occ_func_14_0(1)) + (occ_func_12_0(34)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(204)) + (occ_func_12_0(1)*occ_func_14_0(31)) + (occ_func_12_0(139)*occ_func_14_0(1)) + (occ_func_12_0(48)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(110)) + (occ_func_12_0(1)*occ_func_14_0(23)) + (occ_func_12_0(1)*occ_func_14_0(44)) + (occ_func_12_0(18)*occ_func_14_0(1)) + (occ_func_12_0(123)*occ_func_14_0(1)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(134)) + (occ_func_12_0(208)) + (occ_func_12_0(34)) + (occ_func_14_0(204)) + (occ_func_14_0(31)) + (occ_func_12_0(139)) + (occ_func_12_0(48)) + (occ_func_14_0(110)) + (occ_func_14_0(23)) + (occ_func_14_0(44)) + (occ_func_12_0(18)) + (occ_func_12_0(123)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_10_0() const{
    return ((occ_func_12_0(32)*occ_func_14_0(2)) + (occ_func_12_0(177)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(51)) + (occ_func_12_0(2)*occ_func_14_0(5)) + (occ_func_12_0(47)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(140)) + (occ_func_12_0(2)*occ_func_14_0(112)) + (occ_func_12_0(201)*occ_func_14_0(2)) + (occ_func_12_0(219)*occ_func_14_0(2)) + (occ_func_12_0(12)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(21)) + (occ_func_12_0(2)*occ_func_14_0(227)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(32)) + (occ_func_12_0(177)) + (occ_func_14_0(51)) + (occ_func_14_0(5)) + (occ_func_12_0(47)) + (occ_func_14_0(140)) + (occ_func_14_0(112)) + (occ_func_12_0(201)) + (occ_func_12_0(219)) + (occ_func_12_0(12)) + (occ_func_14_0(21)) + (occ_func_14_0(227)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_10_0() const{
    return ((occ_func_12_0(46)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(50)) + (occ_func_12_0(128)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(28)) + (occ_func_12_0(9)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(173)) + (occ_func_12_0(3)*occ_func_14_0(222)) + (occ_func_12_0(116)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(24)) + (occ_func_12_0(17)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(189)) + (occ_func_12_0(214)*occ_func_14_0(3)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(46)) + (occ_func_14_0(50)) + (occ_func_12_0(128)) + (occ_func_14_0(28)) + (occ_func_12_0(9)) + (occ_func_14_0(173)) + (occ_func_14_0(222)) + (occ_func_12_0(116)) + (occ_func_14_0(24)) + (occ_func_12_0(17)) + (occ_func_14_0(189)) + (occ_func_12_0(214)))/24.0;
  }

  /**** Basis functions for orbit 2, 11****
#Points: 2
MaxLength: 11.2796759  MinLength: 11.2796759
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -1.0000000   -1.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_11_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(46)) + (occ_func_12_0(1)*occ_func_14_0(82)) + (occ_func_12_0(0)*occ_func_14_0(193)) + (occ_func_12_0(2)*occ_func_14_0(11)) + (occ_func_12_0(2)*occ_func_14_0(41)) + (occ_func_12_0(3)*occ_func_14_0(34)) + (occ_func_12_0(0)*occ_func_14_0(87)) + (occ_func_12_0(3)*occ_func_14_0(52)) + (occ_func_12_0(1)*occ_func_14_0(108)) + (occ_func_12_0(1)*occ_func_14_0(47)) + (occ_func_12_0(2)*occ_func_14_0(76)) + (occ_func_12_0(3)*occ_func_14_0(97)) + (occ_func_12_0(2)*occ_func_14_0(60)) + (occ_func_12_0(0)*occ_func_14_0(17)) + (occ_func_12_0(1)*occ_func_14_0(66)) + (occ_func_12_0(3)*occ_func_14_0(174)) + (occ_func_12_0(0)*occ_func_14_0(58)) + (occ_func_12_0(0)*occ_func_14_0(95)) + (occ_func_12_0(3)*occ_func_14_0(48)) + (occ_func_12_0(1)*occ_func_14_0(39)) + (occ_func_12_0(1)*occ_func_14_0(12)) + (occ_func_12_0(2)*occ_func_14_0(49)) + (occ_func_12_0(3)*occ_func_14_0(105)) + (occ_func_12_0(2)*occ_func_14_0(143)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_11_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(46)) + (occ_func_12_0(0)*occ_func_14_0(193)) + (occ_func_12_0(0)*occ_func_14_0(87)) + (occ_func_12_0(43)*occ_func_14_0(0)) + (occ_func_12_0(201)*occ_func_14_0(0)) + (occ_func_12_0(90)*occ_func_14_0(0)) + (occ_func_12_0(106)*occ_func_14_0(0)) + (occ_func_12_0(0)*occ_func_14_0(17)) + (occ_func_12_0(0)*occ_func_14_0(58)) + (occ_func_12_0(0)*occ_func_14_0(95)) + (occ_func_12_0(47)*occ_func_14_0(0)) + (occ_func_12_0(25)*occ_func_14_0(0)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(46)) + (occ_func_14_0(193)) + (occ_func_14_0(87)) + (occ_func_12_0(43)) + (occ_func_12_0(201)) + (occ_func_12_0(90)) + (occ_func_12_0(106)) + (occ_func_14_0(17)) + (occ_func_14_0(58)) + (occ_func_14_0(95)) + (occ_func_12_0(47)) + (occ_func_12_0(25)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_11_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(82)) + (occ_func_12_0(116)*occ_func_14_0(1)) + (occ_func_12_0(54)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(108)) + (occ_func_12_0(1)*occ_func_14_0(47)) + (occ_func_12_0(71)*occ_func_14_0(1)) + (occ_func_12_0(20)*occ_func_14_0(1)) + (occ_func_12_0(1)*occ_func_14_0(66)) + (occ_func_12_0(1)*occ_func_14_0(39)) + (occ_func_12_0(1)*occ_func_14_0(12)) + (occ_func_12_0(46)*occ_func_14_0(1)) + (occ_func_12_0(63)*occ_func_14_0(1)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(82)) + (occ_func_12_0(116)) + (occ_func_12_0(54)) + (occ_func_14_0(108)) + (occ_func_14_0(47)) + (occ_func_12_0(71)) + (occ_func_12_0(20)) + (occ_func_14_0(66)) + (occ_func_14_0(39)) + (occ_func_14_0(12)) + (occ_func_12_0(46)) + (occ_func_12_0(63)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_11_0() const{
    return ((occ_func_12_0(48)*occ_func_14_0(2)) + (occ_func_12_0(85)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(11)) + (occ_func_12_0(2)*occ_func_14_0(41)) + (occ_func_12_0(7)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(76)) + (occ_func_12_0(2)*occ_func_14_0(60)) + (occ_func_12_0(101)*occ_func_14_0(2)) + (occ_func_12_0(139)*occ_func_14_0(2)) + (occ_func_12_0(36)*occ_func_14_0(2)) + (occ_func_12_0(2)*occ_func_14_0(49)) + (occ_func_12_0(2)*occ_func_14_0(143)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(48)) + (occ_func_12_0(85)) + (occ_func_14_0(11)) + (occ_func_14_0(41)) + (occ_func_12_0(7)) + (occ_func_14_0(76)) + (occ_func_14_0(60)) + (occ_func_12_0(101)) + (occ_func_12_0(139)) + (occ_func_12_0(36)) + (occ_func_14_0(49)) + (occ_func_14_0(143)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_11_0() const{
    return ((occ_func_12_0(30)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(34)) + (occ_func_12_0(80)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(52)) + (occ_func_12_0(49)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(97)) + (occ_func_12_0(3)*occ_func_14_0(174)) + (occ_func_12_0(72)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(48)) + (occ_func_12_0(57)*occ_func_14_0(3)) + (occ_func_12_0(3)*occ_func_14_0(105)) + (occ_func_12_0(170)*occ_func_14_0(3)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(30)) + (occ_func_14_0(34)) + (occ_func_12_0(80)) + (occ_func_14_0(52)) + (occ_func_12_0(49)) + (occ_func_14_0(97)) + (occ_func_14_0(174)) + (occ_func_12_0(72)) + (occ_func_14_0(48)) + (occ_func_12_0(57)) + (occ_func_14_0(105)) + (occ_func_12_0(170)))/24.0;
  }

  /**** Basis functions for orbit 2, 12****
#Points: 2
MaxLength: 11.9392378  MinLength: 11.9392378
               0.0000000    0.0000000    0.7194200 Nb Ta
              -1.0000000    1.0000000    0.7194200 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_12_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(76)) + (occ_func_12_0(1)*occ_func_12_0(61)) + (occ_func_12_0(0)*occ_func_12_0(64)) + (occ_func_12_0(2)*occ_func_12_0(102)) + (occ_func_12_0(2)*occ_func_12_0(82)) + (occ_func_12_0(3)*occ_func_12_0(75)) + (occ_func_12_0(0)*occ_func_12_0(104)) + (occ_func_12_0(3)*occ_func_12_0(87)) + (occ_func_12_0(1)*occ_func_12_0(93)) + (occ_func_12_0(1)*occ_func_12_0(97)) + (occ_func_12_0(2)*occ_func_12_0(70)) + (occ_func_12_0(3)*occ_func_12_0(91)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_12_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(76)) + (occ_func_12_0(88)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(64)) + (occ_func_12_0(100)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(104)) + (occ_func_12_0(60)*occ_func_12_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_12_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_12_0(76)) + (occ_func_12_0(88)) + (occ_func_12_0(64)) + (occ_func_12_0(100)) + (occ_func_12_0(104)) + (occ_func_12_0(60)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_12_0() const{
    return ((occ_func_12_0(1)*occ_func_12_0(61)) + (occ_func_12_0(105)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(93)) + (occ_func_12_0(73)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(97)) + (occ_func_12_0(69)*occ_func_12_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_12_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(61)) + (occ_func_12_0(105)) + (occ_func_12_0(93)) + (occ_func_12_0(73)) + (occ_func_12_0(97)) + (occ_func_12_0(69)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_12_0() const{
    return ((occ_func_12_0(2)*occ_func_12_0(102)) + (occ_func_12_0(66)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(82)) + (occ_func_12_0(86)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(70)) + (occ_func_12_0(98)*occ_func_12_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_12_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(102)) + (occ_func_12_0(66)) + (occ_func_12_0(82)) + (occ_func_12_0(86)) + (occ_func_12_0(70)) + (occ_func_12_0(98)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_12_0() const{
    return ((occ_func_12_0(3)*occ_func_12_0(75)) + (occ_func_12_0(95)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(87)) + (occ_func_12_0(83)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(91)) + (occ_func_12_0(79)*occ_func_12_0(3)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_12_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(75)) + (occ_func_12_0(95)) + (occ_func_12_0(87)) + (occ_func_12_0(83)) + (occ_func_12_0(91)) + (occ_func_12_0(79)))/12.0;
  }

  /**** Basis functions for orbit 2, 13****
#Points: 2
MaxLength: 11.9392392  MinLength: 11.9392392
               0.0000000    0.0000000    0.7194200 Nb Ta
              -1.0000000   -1.0000000   -1.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_2_13_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(68)) + (occ_func_12_0(1)*occ_func_12_0(81)) + (occ_func_12_0(0)*occ_func_12_0(92)) + (occ_func_12_0(2)*occ_func_12_0(74)) + (occ_func_12_0(2)*occ_func_12_0(62)) + (occ_func_12_0(3)*occ_func_12_0(103)) + (occ_func_12_0(0)*occ_func_12_0(84)) + (occ_func_12_0(3)*occ_func_12_0(107)) + (occ_func_12_0(1)*occ_func_12_0(65)) + (occ_func_12_0(1)*occ_func_12_0(89)) + (occ_func_12_0(2)*occ_func_12_0(78)) + (occ_func_12_0(3)*occ_func_12_0(99)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_2_13_0() const{
    return ((occ_func_12_0(0)*occ_func_12_0(68)) + (occ_func_12_0(96)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(92)) + (occ_func_12_0(72)*occ_func_12_0(0)) + (occ_func_12_0(0)*occ_func_12_0(84)) + (occ_func_12_0(80)*occ_func_12_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_2_13_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_12_0(68)) + (occ_func_12_0(96)) + (occ_func_12_0(92)) + (occ_func_12_0(72)) + (occ_func_12_0(84)) + (occ_func_12_0(80)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_2_13_0() const{
    return ((occ_func_12_0(1)*occ_func_12_0(81)) + (occ_func_12_0(85)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(65)) + (occ_func_12_0(101)*occ_func_12_0(1)) + (occ_func_12_0(1)*occ_func_12_0(89)) + (occ_func_12_0(77)*occ_func_12_0(1)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_2_13_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(81)) + (occ_func_12_0(85)) + (occ_func_12_0(65)) + (occ_func_12_0(101)) + (occ_func_12_0(89)) + (occ_func_12_0(77)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_2_13_0() const{
    return ((occ_func_12_0(2)*occ_func_12_0(74)) + (occ_func_12_0(94)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(62)) + (occ_func_12_0(106)*occ_func_12_0(2)) + (occ_func_12_0(2)*occ_func_12_0(78)) + (occ_func_12_0(90)*occ_func_12_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_2_13_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(74)) + (occ_func_12_0(94)) + (occ_func_12_0(62)) + (occ_func_12_0(106)) + (occ_func_12_0(78)) + (occ_func_12_0(90)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_2_13_0() const{
    return ((occ_func_12_0(3)*occ_func_12_0(103)) + (occ_func_12_0(67)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(107)) + (occ_func_12_0(63)*occ_func_12_0(3)) + (occ_func_12_0(3)*occ_func_12_0(99)) + (occ_func_12_0(71)*occ_func_12_0(3)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_2_13_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(103)) + (occ_func_12_0(67)) + (occ_func_12_0(107)) + (occ_func_12_0(63)) + (occ_func_12_0(99)) + (occ_func_12_0(71)))/12.0;
  }

  /**** Basis functions for orbit 3, 0****
#Points: 3
MaxLength: 3.3499115  MinLength: 3.3499113
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -0.0000000    1.0000000 Nb Ta
               0.0000000   -0.2805800    1.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_3_0_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(22)*occ_func_13_0(53)) + (occ_func_12_0(1)*occ_func_14_0(10)*occ_func_13_0(27)) + (occ_func_12_0(0)*occ_func_14_0(53)*occ_func_13_0(35)) + (occ_func_12_0(3)*occ_func_14_0(38)*occ_func_13_0(4)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_3_0_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(22)*occ_func_13_0(53)) + (occ_func_12_0(0)*occ_func_14_0(53)*occ_func_13_0(35)) + (occ_func_12_0(35)*occ_func_14_0(22)*occ_func_13_0(0)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(22)*occ_func_13_0(53)) + (occ_func_14_0(53)*occ_func_13_0(35)) + (occ_func_12_0(35)*occ_func_14_0(22)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_3_0_0() const{
    return ((occ_func_12_0(40)*occ_func_14_0(10)*occ_func_13_0(1)) + (occ_func_12_0(1)*occ_func_14_0(10)*occ_func_13_0(27)) + (occ_func_12_0(40)*occ_func_14_0(1)*occ_func_13_0(27)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(40)*occ_func_14_0(10)) + (occ_func_14_0(10)*occ_func_13_0(27)) + (occ_func_12_0(40)*occ_func_13_0(27)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_3_0_0() const{
    return ((occ_func_12_0(16)*occ_func_14_0(2)*occ_func_13_0(29)) + (occ_func_12_0(29)*occ_func_14_0(2)*occ_func_13_0(59)) + (occ_func_12_0(59)*occ_func_14_0(2)*occ_func_13_0(16)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(16)*occ_func_13_0(29)) + (occ_func_12_0(29)*occ_func_13_0(59)) + (occ_func_12_0(59)*occ_func_13_0(16)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_3_0_0() const{
    return ((occ_func_12_0(13)*occ_func_14_0(38)*occ_func_13_0(3)) + (occ_func_12_0(4)*occ_func_14_0(13)*occ_func_13_0(3)) + (occ_func_12_0(3)*occ_func_14_0(38)*occ_func_13_0(4)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(13)*occ_func_14_0(38)) + (occ_func_12_0(4)*occ_func_14_0(13)) + (occ_func_14_0(38)*occ_func_13_0(4)))/4.0;
  }

  /**** Basis functions for orbit 3, 1****
#Points: 3
MaxLength: 4.9680025  MinLength: 3.3499113
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800    0.0000000    0.0000000 Nb Ta
               0.0000000   -0.2805800    0.0000000 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_3_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)) + (occ_func_12_0(0)*occ_func_14_0(29)*occ_func_13_0(59)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)) + (occ_func_12_0(2)*occ_func_14_0(1)*occ_func_13_0(40)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)) + (occ_func_12_0(0)*occ_func_14_0(59)*occ_func_13_0(2)) + (occ_func_12_0(3)*occ_func_14_0(0)*occ_func_13_0(53)) + (occ_func_12_0(1)*occ_func_14_0(4)*occ_func_13_0(38)) + (occ_func_12_0(1)*occ_func_14_0(3)*occ_func_13_0(4)) + (occ_func_12_0(2)*occ_func_14_0(40)*occ_func_13_0(27)) + (occ_func_12_0(3)*occ_func_14_0(53)*occ_func_13_0(22)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_3_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)) + (occ_func_12_0(0)*occ_func_14_0(29)*occ_func_13_0(59)) + (occ_func_12_0(54)*occ_func_14_0(53)*occ_func_13_0(0)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)) + (occ_func_12_0(0)*occ_func_14_0(59)*occ_func_13_0(2)) + (occ_func_12_0(3)*occ_func_14_0(0)*occ_func_13_0(53)) + (occ_func_12_0(33)*occ_func_14_0(0)*occ_func_13_0(22)) + (occ_func_12_0(33)*occ_func_14_0(35)*occ_func_13_0(0)) + (occ_func_12_0(54)*occ_func_14_0(0)*occ_func_13_0(35)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(2)*occ_func_13_0(29)) + (occ_func_14_0(29)*occ_func_13_0(59)) + (occ_func_12_0(54)*occ_func_14_0(53)) + (occ_func_12_0(3)*occ_func_14_0(22)) + (occ_func_14_0(59)*occ_func_13_0(2)) + (occ_func_12_0(3)*occ_func_13_0(53)) + (occ_func_12_0(33)*occ_func_13_0(22)) + (occ_func_12_0(33)*occ_func_14_0(35)) + (occ_func_12_0(54)*occ_func_13_0(35)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_3_1_0() const{
    return ((occ_func_12_0(8)*occ_func_14_0(10)*occ_func_13_0(1)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)) + (occ_func_12_0(8)*occ_func_14_0(1)*occ_func_13_0(27)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)) + (occ_func_12_0(2)*occ_func_14_0(1)*occ_func_13_0(40)) + (occ_func_12_0(43)*occ_func_14_0(40)*occ_func_13_0(1)) + (occ_func_12_0(1)*occ_func_14_0(4)*occ_func_13_0(38)) + (occ_func_12_0(1)*occ_func_14_0(3)*occ_func_13_0(4)) + (occ_func_12_0(43)*occ_func_14_0(1)*occ_func_13_0(10)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(8)*occ_func_14_0(10)) + (occ_func_14_0(38)*occ_func_13_0(3)) + (occ_func_12_0(8)*occ_func_13_0(27)) + (occ_func_12_0(2)*occ_func_14_0(27)) + (occ_func_12_0(2)*occ_func_13_0(40)) + (occ_func_12_0(43)*occ_func_14_0(40)) + (occ_func_14_0(4)*occ_func_13_0(38)) + (occ_func_14_0(3)*occ_func_13_0(4)) + (occ_func_12_0(43)*occ_func_13_0(10)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_3_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)) + (occ_func_12_0(57)*occ_func_14_0(2)*occ_func_13_0(59)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)) + (occ_func_12_0(2)*occ_func_14_0(1)*occ_func_13_0(40)) + (occ_func_12_0(19)*occ_func_14_0(2)*occ_func_13_0(16)) + (occ_func_12_0(0)*occ_func_14_0(59)*occ_func_13_0(2)) + (occ_func_12_0(57)*occ_func_14_0(16)*occ_func_13_0(2)) + (occ_func_12_0(2)*occ_func_14_0(40)*occ_func_13_0(27)) + (occ_func_12_0(19)*occ_func_14_0(29)*occ_func_13_0(2)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(0)*occ_func_13_0(29)) + (occ_func_12_0(57)*occ_func_13_0(59)) + (occ_func_14_0(27)*occ_func_13_0(1)) + (occ_func_14_0(1)*occ_func_13_0(40)) + (occ_func_12_0(19)*occ_func_13_0(16)) + (occ_func_12_0(0)*occ_func_14_0(59)) + (occ_func_12_0(57)*occ_func_14_0(16)) + (occ_func_14_0(40)*occ_func_13_0(27)) + (occ_func_12_0(19)*occ_func_14_0(29)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_3_1_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)) + (occ_func_12_0(36)*occ_func_14_0(13)*occ_func_13_0(3)) + (occ_func_12_0(14)*occ_func_14_0(3)*occ_func_13_0(13)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)) + (occ_func_12_0(36)*occ_func_14_0(3)*occ_func_13_0(38)) + (occ_func_12_0(3)*occ_func_14_0(0)*occ_func_13_0(53)) + (occ_func_12_0(1)*occ_func_14_0(3)*occ_func_13_0(4)) + (occ_func_12_0(14)*occ_func_14_0(4)*occ_func_13_0(3)) + (occ_func_12_0(3)*occ_func_14_0(53)*occ_func_13_0(22)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(1)*occ_func_14_0(38)) + (occ_func_12_0(36)*occ_func_14_0(13)) + (occ_func_12_0(14)*occ_func_13_0(13)) + (occ_func_14_0(22)*occ_func_13_0(0)) + (occ_func_12_0(36)*occ_func_13_0(38)) + (occ_func_14_0(0)*occ_func_13_0(53)) + (occ_func_12_0(1)*occ_func_13_0(4)) + (occ_func_12_0(14)*occ_func_14_0(4)) + (occ_func_14_0(53)*occ_func_13_0(22)))/12.0;
  }

  /**** Basis functions for orbit 3, 2****
#Points: 3
MaxLength: 6.5191431  MinLength: 3.3499113
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -1.0000000    0.0000000 Nb Ta
               0.2805800   -0.7194200    0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_3_2_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)*occ_func_15_0(31)) + (occ_func_12_0(1)*occ_func_14_0(22)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(105)*occ_func_15_0(34)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(4)) + (occ_func_12_0(2)*occ_func_14_0(17)*occ_func_15_0(87)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(29)) + (occ_func_12_0(0)*occ_func_14_0(27)*occ_func_15_0(1)) + (occ_func_12_0(3)*occ_func_14_0(12)*occ_func_15_0(82)) + (occ_func_12_0(1)*occ_func_14_0(60)*occ_func_15_0(11)) + (occ_func_12_0(1)*occ_func_14_0(19)*occ_func_15_0(6)) + (occ_func_12_0(2)*occ_func_14_0(24)*occ_func_15_0(33)) + (occ_func_12_0(3)*occ_func_14_0(21)*occ_func_15_0(8)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_3_2_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)*occ_func_15_0(31)) + (occ_func_12_0(1)*occ_func_14_0(22)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(105)*occ_func_15_0(34)) + (occ_func_12_0(34)*occ_func_14_0(35)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(27)*occ_func_15_0(1)) + (occ_func_12_0(27)*occ_func_14_0(0)*occ_func_15_0(22)) + (occ_func_12_0(105)*occ_func_14_0(0)*occ_func_15_0(35)) + (occ_func_12_0(14)*occ_func_14_0(0)*occ_func_15_0(53)) + (occ_func_12_0(31)*occ_func_14_0(53)*occ_func_15_0(0)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(14)*occ_func_15_0(31)) + (occ_func_12_0(1)*occ_func_14_0(22)) + (occ_func_14_0(105)*occ_func_15_0(34)) + (occ_func_12_0(34)*occ_func_14_0(35)) + (occ_func_14_0(27)*occ_func_15_0(1)) + (occ_func_12_0(27)*occ_func_15_0(22)) + (occ_func_12_0(105)*occ_func_15_0(35)) + (occ_func_12_0(14)*occ_func_15_0(53)) + (occ_func_12_0(31)*occ_func_14_0(53)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_3_2_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(22)*occ_func_15_0(0)) + (occ_func_12_0(60)*occ_func_14_0(1)*occ_func_15_0(10)) + (occ_func_12_0(22)*occ_func_14_0(1)*occ_func_15_0(27)) + (occ_func_12_0(11)*occ_func_14_0(10)*occ_func_15_0(1)) + (occ_func_12_0(0)*occ_func_14_0(27)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_14_0(60)*occ_func_15_0(11)) + (occ_func_12_0(1)*occ_func_14_0(19)*occ_func_15_0(6)) + (occ_func_12_0(6)*occ_func_14_0(40)*occ_func_15_0(1)) + (occ_func_12_0(19)*occ_func_14_0(1)*occ_func_15_0(40)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(22)*occ_func_15_0(0)) + (occ_func_12_0(60)*occ_func_15_0(10)) + (occ_func_12_0(22)*occ_func_15_0(27)) + (occ_func_12_0(11)*occ_func_14_0(10)) + (occ_func_12_0(0)*occ_func_14_0(27)) + (occ_func_14_0(60)*occ_func_15_0(11)) + (occ_func_14_0(19)*occ_func_15_0(6)) + (occ_func_12_0(6)*occ_func_14_0(40)) + (occ_func_12_0(19)*occ_func_15_0(40)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_3_2_0() const{
    return ((occ_func_12_0(24)*occ_func_14_0(2)*occ_func_15_0(59)) + (occ_func_12_0(17)*occ_func_14_0(2)*occ_func_15_0(16)) + (occ_func_12_0(4)*occ_func_14_0(29)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(4)) + (occ_func_12_0(2)*occ_func_14_0(17)*occ_func_15_0(87)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(29)) + (occ_func_12_0(87)*occ_func_14_0(16)*occ_func_15_0(2)) + (occ_func_12_0(33)*occ_func_14_0(59)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(24)*occ_func_15_0(33)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(24)*occ_func_15_0(59)) + (occ_func_12_0(17)*occ_func_15_0(16)) + (occ_func_12_0(4)*occ_func_14_0(29)) + (occ_func_14_0(3)*occ_func_15_0(4)) + (occ_func_14_0(17)*occ_func_15_0(87)) + (occ_func_12_0(3)*occ_func_15_0(29)) + (occ_func_12_0(87)*occ_func_14_0(16)) + (occ_func_12_0(33)*occ_func_14_0(59)) + (occ_func_14_0(24)*occ_func_15_0(33)))/12.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_3_2_0() const{
    return ((occ_func_12_0(8)*occ_func_14_0(38)*occ_func_15_0(3)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(4)) + (occ_func_12_0(82)*occ_func_14_0(13)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(29)) + (occ_func_12_0(12)*occ_func_14_0(3)*occ_func_15_0(13)) + (occ_func_12_0(3)*occ_func_14_0(12)*occ_func_15_0(82)) + (occ_func_12_0(29)*occ_func_14_0(4)*occ_func_15_0(3)) + (occ_func_12_0(21)*occ_func_14_0(3)*occ_func_15_0(38)) + (occ_func_12_0(3)*occ_func_14_0(21)*occ_func_15_0(8)))/12.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(8)*occ_func_14_0(38)) + (occ_func_12_0(2)*occ_func_15_0(4)) + (occ_func_12_0(82)*occ_func_14_0(13)) + (occ_func_14_0(2)*occ_func_15_0(29)) + (occ_func_12_0(12)*occ_func_15_0(13)) + (occ_func_14_0(12)*occ_func_15_0(82)) + (occ_func_12_0(29)*occ_func_14_0(4)) + (occ_func_12_0(21)*occ_func_15_0(38)) + (occ_func_14_0(21)*occ_func_15_0(8)))/12.0;
  }

  /**** Basis functions for orbit 3, 3****
#Points: 3
MaxLength: 6.5191431  MinLength: 4.9680020
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -1.0000000    0.0000000 Nb Ta
              -0.7194200   -0.7194200    0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_3_3_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)*occ_func_15_0(3)) + (occ_func_12_0(1)*occ_func_14_0(22)*occ_func_15_0(8)) + (occ_func_12_0(0)*occ_func_14_0(105)*occ_func_15_0(54)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(2)*occ_func_14_0(17)*occ_func_15_0(19)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(0)*occ_func_14_0(27)*occ_func_15_0(33)) + (occ_func_12_0(3)*occ_func_14_0(12)*occ_func_15_0(14)) + (occ_func_12_0(1)*occ_func_14_0(60)*occ_func_15_0(43)) + (occ_func_12_0(1)*occ_func_14_0(19)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(24)*occ_func_15_0(57)) + (occ_func_12_0(3)*occ_func_14_0(21)*occ_func_15_0(36)) + (occ_func_12_0(2)*occ_func_14_0(4)*occ_func_15_0(19)) + (occ_func_12_0(0)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(1)*occ_func_14_0(6)*occ_func_15_0(43)) + (occ_func_12_0(3)*occ_func_14_0(82)*occ_func_15_0(36)) + (occ_func_12_0(0)*occ_func_14_0(34)*occ_func_15_0(33)) + (occ_func_12_0(0)*occ_func_14_0(31)*occ_func_15_0(54)) + (occ_func_12_0(3)*occ_func_14_0(8)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_14_0(11)*occ_func_15_0(8)) + (occ_func_12_0(1)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(33)*occ_func_15_0(0)) + (occ_func_12_0(3)*occ_func_14_0(29)*occ_func_15_0(14)) + (occ_func_12_0(2)*occ_func_14_0(87)*occ_func_15_0(57)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_3_3_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)*occ_func_15_0(3)) + (occ_func_12_0(29)*occ_func_14_0(54)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(105)*occ_func_15_0(54)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(27)*occ_func_15_0(33)) + (occ_func_12_0(27)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(105)*occ_func_14_0(0)*occ_func_15_0(59)) + (occ_func_12_0(14)*occ_func_14_0(0)*occ_func_15_0(29)) + (occ_func_12_0(59)*occ_func_14_0(33)*occ_func_15_0(0)) + (occ_func_12_0(34)*occ_func_14_0(0)*occ_func_15_0(59)) + (occ_func_12_0(0)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(59)*occ_func_14_0(54)*occ_func_15_0(0)) + (occ_func_12_0(0)*occ_func_14_0(34)*occ_func_15_0(33)) + (occ_func_12_0(0)*occ_func_14_0(31)*occ_func_15_0(54)) + (occ_func_12_0(31)*occ_func_14_0(0)*occ_func_15_0(29)) + (occ_func_12_0(29)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(1)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(33)*occ_func_15_0(0)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_3_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(14)*occ_func_15_0(3)) + (occ_func_12_0(29)*occ_func_14_0(54)) + (occ_func_14_0(105)*occ_func_15_0(54)) + (occ_func_12_0(2)*occ_func_14_0(3)) + (occ_func_14_0(27)*occ_func_15_0(33)) + (occ_func_12_0(27)*occ_func_15_0(2)) + (occ_func_12_0(105)*occ_func_15_0(59)) + (occ_func_12_0(14)*occ_func_15_0(29)) + (occ_func_12_0(59)*occ_func_14_0(33)) + (occ_func_12_0(34)*occ_func_15_0(59)) + (occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(59)*occ_func_14_0(54)) + (occ_func_14_0(34)*occ_func_15_0(33)) + (occ_func_14_0(31)*occ_func_15_0(54)) + (occ_func_12_0(31)*occ_func_15_0(29)) + (occ_func_12_0(29)*occ_func_14_0(3)) + (occ_func_12_0(1)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(33)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_3_3_0() const{
    return ((occ_func_12_0(1)*occ_func_14_0(22)*occ_func_15_0(8)) + (occ_func_12_0(60)*occ_func_14_0(1)*occ_func_15_0(38)) + (occ_func_12_0(22)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(4)*occ_func_14_0(43)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_14_0(60)*occ_func_15_0(43)) + (occ_func_12_0(1)*occ_func_14_0(19)*occ_func_15_0(2)) + (occ_func_12_0(38)*occ_func_14_0(8)*occ_func_15_0(1)) + (occ_func_12_0(19)*occ_func_14_0(1)*occ_func_15_0(4)) + (occ_func_12_0(0)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(1)*occ_func_14_0(6)*occ_func_15_0(43)) + (occ_func_12_0(4)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(3)*occ_func_14_0(8)*occ_func_15_0(1)) + (occ_func_12_0(1)*occ_func_14_0(11)*occ_func_15_0(8)) + (occ_func_12_0(1)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(6)*occ_func_14_0(1)*occ_func_15_0(4)) + (occ_func_12_0(11)*occ_func_14_0(1)*occ_func_15_0(38)) + (occ_func_12_0(38)*occ_func_14_0(43)*occ_func_15_0(1)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_3_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_14_0(22)*occ_func_15_0(8)) + (occ_func_12_0(60)*occ_func_15_0(38)) + (occ_func_12_0(22)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_14_0(2)) + (occ_func_12_0(4)*occ_func_14_0(43)) + (occ_func_14_0(60)*occ_func_15_0(43)) + (occ_func_14_0(19)*occ_func_15_0(2)) + (occ_func_12_0(38)*occ_func_14_0(8)) + (occ_func_12_0(19)*occ_func_15_0(4)) + (occ_func_12_0(0)*occ_func_15_0(3)) + (occ_func_14_0(6)*occ_func_15_0(43)) + (occ_func_12_0(4)*occ_func_14_0(2)) + (occ_func_12_0(3)*occ_func_14_0(8)) + (occ_func_14_0(11)*occ_func_15_0(8)) + (occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(6)*occ_func_15_0(4)) + (occ_func_12_0(11)*occ_func_15_0(38)) + (occ_func_12_0(38)*occ_func_14_0(43)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_3_3_0() const{
    return ((occ_func_12_0(24)*occ_func_14_0(2)*occ_func_15_0(27)) + (occ_func_12_0(17)*occ_func_14_0(2)*occ_func_15_0(40)) + (occ_func_12_0(40)*occ_func_14_0(57)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(2)*occ_func_14_0(17)*occ_func_15_0(19)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(27)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(1)*occ_func_14_0(19)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(24)*occ_func_15_0(57)) + (occ_func_12_0(2)*occ_func_14_0(4)*occ_func_15_0(19)) + (occ_func_12_0(33)*occ_func_14_0(2)*occ_func_15_0(27)) + (occ_func_12_0(87)*occ_func_14_0(2)*occ_func_15_0(40)) + (occ_func_12_0(4)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(40)*occ_func_14_0(19)*occ_func_15_0(2)) + (occ_func_12_0(1)*occ_func_14_0(0)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(33)*occ_func_15_0(0)) + (occ_func_12_0(27)*occ_func_14_0(57)*occ_func_15_0(2)) + (occ_func_12_0(2)*occ_func_14_0(87)*occ_func_15_0(57)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_3_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(24)*occ_func_15_0(27)) + (occ_func_12_0(17)*occ_func_15_0(40)) + (occ_func_12_0(40)*occ_func_14_0(57)) + (occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_14_0(17)*occ_func_15_0(19)) + (occ_func_12_0(3)*occ_func_15_0(1)) + (occ_func_12_0(27)*occ_func_14_0(0)) + (occ_func_12_0(1)*occ_func_14_0(19)) + (occ_func_14_0(24)*occ_func_15_0(57)) + (occ_func_14_0(4)*occ_func_15_0(19)) + (occ_func_12_0(33)*occ_func_15_0(27)) + (occ_func_12_0(87)*occ_func_15_0(40)) + (occ_func_12_0(4)*occ_func_15_0(1)) + (occ_func_12_0(40)*occ_func_14_0(19)) + (occ_func_12_0(1)*occ_func_14_0(0)) + (occ_func_14_0(33)*occ_func_15_0(0)) + (occ_func_12_0(27)*occ_func_14_0(57)) + (occ_func_14_0(87)*occ_func_15_0(57)))/24.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_3_3_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(14)*occ_func_15_0(3)) + (occ_func_12_0(2)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(22)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(12)*occ_func_14_0(3)*occ_func_15_0(53)) + (occ_func_12_0(3)*occ_func_14_0(12)*occ_func_15_0(14)) + (occ_func_12_0(53)*occ_func_14_0(36)*occ_func_15_0(3)) + (occ_func_12_0(21)*occ_func_14_0(3)*occ_func_15_0(22)) + (occ_func_12_0(3)*occ_func_14_0(21)*occ_func_15_0(36)) + (occ_func_12_0(22)*occ_func_14_0(36)*occ_func_15_0(3)) + (occ_func_12_0(0)*occ_func_14_0(1)*occ_func_15_0(3)) + (occ_func_12_0(53)*occ_func_14_0(14)*occ_func_15_0(3)) + (occ_func_12_0(3)*occ_func_14_0(82)*occ_func_15_0(36)) + (occ_func_12_0(8)*occ_func_14_0(3)*occ_func_15_0(22)) + (occ_func_12_0(3)*occ_func_14_0(8)*occ_func_15_0(1)) + (occ_func_12_0(29)*occ_func_14_0(3)*occ_func_15_0(0)) + (occ_func_12_0(3)*occ_func_14_0(29)*occ_func_15_0(14)) + (occ_func_12_0(82)*occ_func_14_0(3)*occ_func_15_0(53)))/24.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_3_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(0)*occ_func_14_0(14)) + (occ_func_12_0(2)*occ_func_15_0(0)) + (occ_func_12_0(22)*occ_func_14_0(1)) + (occ_func_14_0(2)*occ_func_15_0(1)) + (occ_func_12_0(12)*occ_func_15_0(53)) + (occ_func_14_0(12)*occ_func_15_0(14)) + (occ_func_12_0(53)*occ_func_14_0(36)) + (occ_func_12_0(21)*occ_func_15_0(22)) + (occ_func_14_0(21)*occ_func_15_0(36)) + (occ_func_12_0(22)*occ_func_14_0(36)) + (occ_func_12_0(0)*occ_func_14_0(1)) + (occ_func_12_0(53)*occ_func_14_0(14)) + (occ_func_14_0(82)*occ_func_15_0(36)) + (occ_func_12_0(8)*occ_func_15_0(22)) + (occ_func_14_0(8)*occ_func_15_0(1)) + (occ_func_12_0(29)*occ_func_15_0(0)) + (occ_func_14_0(29)*occ_func_15_0(14)) + (occ_func_12_0(82)*occ_func_15_0(53)))/24.0;
  }

  /**** Basis functions for orbit 4, 0****
#Points: 4
MaxLength: 3.3499117  MinLength: 3.3499113
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800   -0.0000000    1.0000000 Nb Ta
               0.0000000   -0.2805800    1.0000000 Nb Ta
               0.2805800    0.2805800    1.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_4_0_0() const{
    return (occ_func_12_0(0)*occ_func_14_0(22)*occ_func_13_0(53)*occ_func_15_0(35));
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_4_0_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(22)*occ_func_13_0(53)*occ_func_15_0(35)))/1.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_4_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(22)*occ_func_13_0(53)*occ_func_15_0(35)))/1.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_4_0_0() const{
    return ((occ_func_12_0(40)*occ_func_14_0(10)*occ_func_13_0(1)*occ_func_15_0(27)))/1.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_4_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(40)*occ_func_14_0(10)*occ_func_15_0(27)))/1.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_4_0_0() const{
    return ((occ_func_12_0(16)*occ_func_14_0(2)*occ_func_13_0(29)*occ_func_15_0(59)))/1.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_4_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(16)*occ_func_13_0(29)*occ_func_15_0(59)))/1.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_4_0_0() const{
    return ((occ_func_12_0(4)*occ_func_14_0(38)*occ_func_13_0(13)*occ_func_15_0(3)))/1.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_4_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(4)*occ_func_14_0(38)*occ_func_13_0(13)))/1.0;
  }

  /**** Basis functions for orbit 4, 1****
#Points: 4
MaxLength: 4.9680025  MinLength: 3.3499113
               0.0000000    0.0000000    0.7194200 Nb Ta
              -0.2805800    0.0000000    0.0000000 Nb Ta
               0.0000000   -0.2805800    0.0000000 Nb Ta
               0.2805800    0.2805800    0.2805800 Nb Ta
****/
  double Nb_direction_Clexulator::eval_bfunc_4_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)*occ_func_15_0(59)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)*occ_func_15_0(4)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)*occ_func_15_0(40)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)*occ_func_15_0(53)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_12_bfunc_4_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)*occ_func_15_0(59)) + (occ_func_12_0(33)*occ_func_14_0(22)*occ_func_13_0(35)*occ_func_15_0(0)) + (occ_func_12_0(54)*occ_func_14_0(35)*occ_func_13_0(53)*occ_func_15_0(0)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)*occ_func_15_0(53)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_12_bfunc_4_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_12_0[occ_f] - m_occ_func_12_0[occ_i])*((occ_func_14_0(2)*occ_func_13_0(29)*occ_func_15_0(59)) + (occ_func_12_0(33)*occ_func_14_0(22)*occ_func_13_0(35)) + (occ_func_12_0(54)*occ_func_14_0(35)*occ_func_13_0(53)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_15_0(53)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_13_bfunc_4_1_0() const{
    return ((occ_func_12_0(8)*occ_func_14_0(10)*occ_func_13_0(1)*occ_func_15_0(27)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)*occ_func_15_0(4)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)*occ_func_15_0(40)) + (occ_func_12_0(43)*occ_func_14_0(10)*occ_func_13_0(40)*occ_func_15_0(1)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_13_bfunc_4_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_13_0[occ_f] - m_occ_func_13_0[occ_i])*((occ_func_12_0(8)*occ_func_14_0(10)*occ_func_15_0(27)) + (occ_func_14_0(38)*occ_func_13_0(3)*occ_func_15_0(4)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_15_0(40)) + (occ_func_12_0(43)*occ_func_14_0(10)*occ_func_13_0(40)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_14_bfunc_4_1_0() const{
    return ((occ_func_12_0(0)*occ_func_14_0(2)*occ_func_13_0(29)*occ_func_15_0(59)) + (occ_func_12_0(57)*occ_func_14_0(2)*occ_func_13_0(59)*occ_func_15_0(16)) + (occ_func_12_0(2)*occ_func_14_0(27)*occ_func_13_0(1)*occ_func_15_0(40)) + (occ_func_12_0(19)*occ_func_14_0(2)*occ_func_13_0(16)*occ_func_15_0(29)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_14_bfunc_4_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_14_0[occ_f] - m_occ_func_14_0[occ_i])*((occ_func_12_0(0)*occ_func_13_0(29)*occ_func_15_0(59)) + (occ_func_12_0(57)*occ_func_13_0(59)*occ_func_15_0(16)) + (occ_func_14_0(27)*occ_func_13_0(1)*occ_func_15_0(40)) + (occ_func_12_0(19)*occ_func_13_0(16)*occ_func_15_0(29)))/4.0;
  }

  double Nb_direction_Clexulator::site_eval_at_15_bfunc_4_1_0() const{
    return ((occ_func_12_0(36)*occ_func_14_0(38)*occ_func_13_0(13)*occ_func_15_0(3)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_13_0(3)*occ_func_15_0(4)) + (occ_func_12_0(14)*occ_func_14_0(3)*occ_func_13_0(13)*occ_func_15_0(4)) + (occ_func_12_0(3)*occ_func_14_0(22)*occ_func_13_0(0)*occ_func_15_0(53)))/4.0;
  }

  double Nb_direction_Clexulator::delta_site_eval_at_15_bfunc_4_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_15_0[occ_f] - m_occ_func_15_0[occ_i])*((occ_func_12_0(36)*occ_func_14_0(38)*occ_func_13_0(13)) + (occ_func_12_0(1)*occ_func_14_0(38)*occ_func_15_0(4)) + (occ_func_12_0(14)*occ_func_13_0(13)*occ_func_15_0(4)) + (occ_func_14_0(22)*occ_func_13_0(0)*occ_func_15_0(53)))/4.0;
  }

}


extern "C" {
  /// \brief Returns a Clexulator_impl::Base* owning a Nb_direction_Clexulator
  CASM::Clexulator_impl::Base* make_Nb_direction_Clexulator() {
    return new CASM::Nb_direction_Clexulator();
  }

}

