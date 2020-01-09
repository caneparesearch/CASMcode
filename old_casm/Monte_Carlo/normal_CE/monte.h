#include<iostream>
#include<fstream>
#include<istream>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<stdarg.h> //Added by John
#include<limits.h>
#include<fnmatch.h>
#include<sys/types.h>
#include<unistd.h>
#include<dirent.h>
#include<iomanip.h>
#include<time.h>
#include "Array.h"

using namespace std;

//NOTE: This is licensed software. See LICENSE.txt for license terms.
//written by Anton Van der Ven, John Thomas, Qingchuan Xu, and Jishnu Bhattacharya
//please see CASMdocumentation.pdf for a tutorial for using the code.
//Version as of May 26, 2010

/*Changes since clusters10.0.h
  ====================
  x  Average correlations for MC simulation
  x  Fix int_to_string and add double_to_string
  x  Generalized Susceptibility for multiple sublattices
  x  Symmetry classification (currently only output for point_group)
  x  Tensor Class
  x  read_mc_input routine
*/
//**************************************************************
//**************************************************************
double tol = 1.0e-3;
double kb=0.00008617;
////////////////////////////////////////////////////////////////////////////////
//swoboda
char basis_type;
////////////////////////////////////////////////////////////////////////////////
void coord_trans_mat(double lat[3][3], double FtoC[3][3], double CtoF[3][3]);


class vec;
class sym_op;
class tensor;
class specie;
class atompos;
class cluster;
class orbit;
class multiplet;
class structure;
class concentration;
class arrangement;
class superstructure;
class configurations;
class facet;
class hull;
class chempot;
class trajectory;
class fluctuation;
class hop;
class mc_index;
class Monte_Carlo;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class vec{
 public:
  bool frac_on,cart_on;
  double fcoord[3],ccoord[3];
  double length;

  vec() {frac_on=false; cart_on=false;}

  double calc_dist();
  vec apply_sym(sym_op op);
  void print_frac(ostream &stream);
  void print_cart(ostream &stream);
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class sym_op{
 public:

  bool frac_on,cart_on;
  double fsym_mat[3][3],csym_mat[3][3];
  double ftau[3],ctau[3];
  double lat[3][3];
  double FtoC[3][3],CtoF[3][3];
  vec eigenvec; //Added by John to hold rotation axis for rotation or mirror plane normal
  short int sym_type;  //Added by John, Not calculated=-2, Inversion=-1, Identity=0, Rotation=1, Mirror=2, Rotoinversion=3
  short int op_angle; //Added by John to hold rotation angle for rotation op


  sym_op(){frac_on=false; cart_on=false;}

  void get_sym_type(); //Added by John, populates eigenvec and sym_type;
  void print_fsym_mat(ostream &stream);
  void print_csym_mat(ostream &stream);
  void get_trans_mat(){ coord_trans_mat(lat,FtoC,CtoF); }
  void get_csym_mat();
  void get_fsym_mat();
  void update();
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//Edited by John
//generalized for any tensor, though needs to be cleaned up.  Some constructors may be redundant.
// ?*?*? May decide to replace unsized arrays with vectors for stability.

class tensor{
 public:
  int rank, size;
  int *dim, *mult;
  double *K;
  tensor();  //Default constructor
  tensor(tensor const& ttens);  //Copy constructor
  tensor& operator=(const tensor& ttens);  //Assignment operator
  ~tensor(){delete [] dim; delete [] mult; delete [] K;};  //Destructor
  tensor(int trank, ...);  //Constrctor taking rank, and size along each dimension.
  tensor(int trank, int *tdim); //Constructor taking rank, and array of sizes along the dimensions
  double get_elem(int ind, ...);  //get element at set of indeces, separated by commas
  double get_elem(int *inds);  //get element at indeces contained in array *inds
  void set_elem(double new_elem, ...);   //set element to new_elem at indeces, separated by commas
  void set_elem(double new_elem, int *inds);  //set elem to new_elem at indeces specified by array *inds
  tensor apply_sym(sym_op op);            // applies op.csym_mat^T* K * op.csym_mat
  void print(ostream &stream);

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class specie{
 public:
  //char name[2];  // commented by jishnu
  string name; // jishnu
  int spin;
  double mass;
  double magmom;  // jishnu
  double U;	  // jishnu  // this is for LDA+U calculations
  double J;	  // jishnu  // this is for LDA+U calculations
	
  specie(){name =""; spin =0; mass =0.0; magmom =0.0; U =0.0;J = 0.0;}
  void print(ostream &stream);
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class atompos{
 public:
  //Cluster expansion related variables
  int bit;
  specie occ;
  vector<specie> compon;
  double fcoord[3],ccoord[3];     //atom coordinates
  double dfcoord[3],dccoord[3];   //difference from ideal atom position
  double delta;
  ////////////////////////////////////////////////////////////////////////////////
  //added for occupation basis by Ben Swoboda
  vector<int> p_vec;              //occupation basis vector
  vector<int> spin_vec;           //vector of spins to be used as basis_vec if flagged
  vector<int> basis_vec;          //vector that stores the values of the basis in use
  char basis_flag;                //charachter 0,1,2,etc. that indicates which basis to use
                                  //nothing or 0=spin-basis, 1=occ-basis
  ////////////////////////////////////////////////////////////////////////////////

  //Cluster expansion related functions
  atompos();
  atompos apply_sym(sym_op op);
  void get_cart(double FtoC[3][3]);
  void get_frac(double CtoF[3][3]);
  void readf(istream &stream);
  void readc(istream &stream);
  void print(ostream &stream);
  void assign_spin();              //the first compon[0] has the highest spin
  int get_spin(string name);      //returns the spin of name[2]

  //Monte Carlo related variables
  int shift[4];                    //gives coordinates of unit cell and basis relative to prim
  vector< vector<int> > flip;      //for each compon, this gives the spins of the other components
  vector<double> mu;               //contains the chemical potentials for each specie with i<compon.size()-1, while mu[compon.size()-1]=0

  //Monte Carlo related functions
  void assemble_flip();
  int iflip(int spin);             //given a spin this gives the index in flip and dmu for that specie

};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class cluster{
 public:
  double min_leng,max_leng; 
  vector<atompos> point;
  vector<sym_op> clust_group;
  //double clustmat[3][3];   // could be a force constant matrix, or a measure of strain for that cluster
  // we want to generalize this to a tensor object

  cluster(){min_leng=0; max_leng=0;} 
  cluster apply_sym(sym_op op);
  void get_dimensions();
  void get_cart(double FtoC[3][3]);
  void get_frac(double CtoF[3][3]);
  void readf(istream &stream, int np);
  void readc(istream &stream, int np);
  void print(ostream &stream);
  void write_clust_group(ostream &stream);
  void determine_site_attributes(structure prim);
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class orbit{
 public:
  vector<cluster> equiv;
  double eci;
  int stability;  // added by jishnu to determine environment stability in monte carlo
  //vector<sym_op> orb_group;  // this will be a set of matrices that link all equiv to the first one


  orbit(){eci=0.0;stability=0;};
  void get_cart(double FtoC[3][3]);
  void get_frac(double CtoF[3][3]);
  void readf(istream &stream, int np, int mult);
  void readc(istream &stream, int np, int mult);
  void print(ostream &stream);
  void determine_site_attributes(structure prim);
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class multiplet{
 public:
  vector< vector<orbit> > orb;

  vector<int> size;
  vector<int> order;
  vector< vector<int> > index; 
  vector< vector<int> > subcluster;

  void readf(istream &stream);
  void readc(istream &stream);
  void print(ostream &stream);
  void get_cart(double FtoC[3][3]);
  void get_frac(double CtoF[3][3]);
  void sort(int np);
  void get_index();
  void get_hierarchy();
  void print_hierarchy(ostream &stream);
  void read_eci(istream &stream);
  void determine_site_attributes(structure prim);
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class structure{
 public:

  //lattice related variables

  char title[200];
  double scale;
  double lat[3][3];                   // cartesian coordinates of the lattice vectors (rows)
  double ilat[3][3];                  // ideal cartesian coordinates of the lattice vectors (rows) - either unrelaxed, or slat[][]*prim.lat[][]
  double slat[3][3];                  // supercell coordinates in terms of a primitive lattice (could be identity matrix)
  double strain[3][3];                //
  double latparam[3],latangle[3];
  int permut[3];
  double ilatparam[3],ilatangle[3];
  int ipermut[3];
  double FtoC[3][3],CtoF[3][3];       // Fractional to Cartesian trans mat and vice versa
  double PtoS[3][3],StoP[3][3];       // Primitive to Supercell trans mat and vice versa
  vector<sym_op> point_group;
  vector<vec> prim_grid;

  //basis related variables

  bool frac_on,cart_on;
  vector<int> num_each_specie;
  vector<specie> compon;
  vector<atompos> atom;
  vector<sym_op> factor_group;

  //reciprocal lattice variables

  double recip_lat[3][3];             // cartesian coordinates of the reciprocal lattice
  double recip_latparam[3],recip_latangle[3];
  int recip_permut[3];



  structure();

  //lattice related routines

  void get_trans_mat(){ coord_trans_mat(lat,FtoC,CtoF); coord_trans_mat(slat,StoP,PtoS);}
  void get_latparam();
  void get_ideal_latparam();
  void calc_point_group();
  void update_lat();
  void generate_3d_supercells(vector<structure> &supercell, int max_vol);
  void generate_2d_supercells(vector<structure> &supercell, int max_vol, int excluded_axis);
  void generate_3d_reduced_cell();
  void generate_2d_reduced_cell(int excluded_axis);
  void generate_slat(structure prim);                     //generates slat from lat
  void generate_lat(structure prim);                     //generates slat from lat   // added by jishnu
  void generate_slat(structure prim, double rescale);     //generates slat from lat after rescaling with rescale
  void generate_ideal_slat(structure prim);             //generates slat from ilat
  void generate_ideal_slat(structure prim, double rescale);
  void calc_strain();
  void generate_prim_grid();

  void read_lat_poscar(istream &stream);
  void write_lat_poscar(ostream &stream);
  void write_point_group();


  //lattice+basis related routines

  void calc_fractional();
  void calc_cartesian();
  void bring_in_cell();
  void calc_factor_group();
  void expand_prim_basis(structure prim);
  void map_on_expanded_prim_basis(structure prim);
  void map_on_expanded_prim_basis(structure prim, arrangement &conf);
  void idealize();
  void expand_prim_clust(multiplet basiplet, multiplet &super_basiplet);
  void collect_components();
  void collect_relax(string dir_name);      // reads in POS and CONTCAR and fills lat, rlat, atom etc.
  void update_struc();


  ////////////////////////////////////////////////////////////////////////////////
  //added by Ben Swoboda
  void read_struc_prim(istream &stream);
  ////////////////////////////////////////////////////////////////////////////////
  void read_species();     // added by jishnu
  void read_struc_poscar(istream &stream);
  void write_struc_poscar(ostream &stream);
  void write_struc_xyz(ostream &stream);
  void write_struc_xyz(ostream &stream, concentration out_conc);
  void write_factor_group();

  //reciprocal lattice related routines

  void calc_recip_lat();
  void get_recip_latparam();

};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class concentration{
 public:
  vector< vector<specie> > compon;
  vector< vector<double> > occup;
  vector< vector<double> > mu;

  void collect_components(structure &prim);
  void calc_concentration(structure &struc);
  void print_concentration(ostream &stream);
  void print_concentration_without_names(ostream &stream);
  void print_names(ostream &stream);
  void get_occup(istream &stream);  // added by jishnu
    
  //Monte Carlo related functions
  void set_zero();
  void increment(concentration conc);
  void normalize(int n);

};




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class arrangement{
 public:
  int ns,nc;                       // supercell index and configuration index
  vector<int> bit;
  concentration conc;
  vector<double> correlations;
  string name;
  double energy,fenergy,cefenergy,fpfenergy;  // first principles energy, formation energy, cluster expanded formation energy // added by jishnu
  int fp,ce,te; // if fp=0, fenergy != fpfenergy and fp=1, fenergy = fpfenergy and same is true for ce and te // te stands for total energy  // added by jishnu
  double delE_from_facet; // added by jishnu // this is energy differnce from the hull
  //double norm_dist_from_facet;  // added by jishnu // this is normal distance from the hull
  double weight,reduction;
  vector<double> coordinate;   // contains concentration and energy of the arrangement
  // vector<double> coordinate_CE;  //contains concentration and CE_enenrgy of the arrangement // added by jishnu
  bool calculated,make,got_fenergy,got_cefenergy;   // got_fenergy is true when formation energy calculation is done for that arranegment // added by jishnu
  int relax_step;  //added by jishnu // to get the no of relaxation steps in the final vasp calculation
	
  arrangement(){calculated = false; make = false; got_fenergy = false; got_cefenergy =false; weight = 0.0; fp =0; ce=0; te =0; reduction=0;}
  void assemble_coordinate_fenergy();          //assemble coordinate vector using first-principles formation energy // added by jishnu
  // void assemble_coordinate_CE();          //assemble coordinate vector using Cluster Expansion formation energy // added by jishnu
  void print_bit(ostream &stream);
  void print_correlations(ostream &stream);
  void print_coordinate(ostream &stream);
  void print_coordinate_ternary(ostream &stream);  // added by jishnu
  void get_bit(istream &stream);  // added by jishnu
  void update_ce(); // to be made to convert fenergy = cefenergy and update coordinate, flag to indicate what is in fenergy // added by jishnu
  void update_fp(); // same but fenergy=fpfenergy   // added by jishnu
  void update_te(); // same but fenergy = energy   // added by jishnu
  //need routine that copy fenergy into fpfenergy or cefenergy 
  void print_in_energy_file(ostream &stream); // added by jishnu
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class superstructure{
 public:
  structure struc,curr_struc;
  vector<arrangement> conf;
  double kmesh[3];
  int nodes,ppn,walltime;  // added by jishnu
  string queue,parent_directory;  // added by jishnu
  vector< vector< vector< vector< int > > > > corr_to_atom_vec;  //added by John - associates basis functions with various curr_struc.atom[i] 
  vector< vector< vector< int > > > basis_to_bit_vec;  //added by John - associates various curr_struc.atom[i].bit with appropriate spin/occupation basis values


  void decorate_superstructure(arrangement conf);
  void determine_kpoint_grid(double kpoint_dens);
  void print(string dir_name, string file_name);
  void print_incar(string dir_name);    // jishnu
  void print_potcar(string dir_name);
  void print_kpoint(string dir_name);
  void print_yihaw(string dir_name);   // added by jishnu // do not use this routine
  void read_yihaw_input();  // added by jishnu // do not use this one
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class facet{
 public:
  vector<arrangement> corner;
  vector < double > normal_vec; // added by jishnu
  double offset; // added by jishnu
  vector<double> mu;
	
  void get_norm_vec(); // finds the coefficients of the facet plane equation ax+by+cz+d = 0; // added by jishnu
  bool find_endiff(arrangement arr, double &delE_from_facet); // sees whether that facet contains an arrangement and if yes, finds out the energy on facet // added by jishnu
  void get_mu(); // added by jishnu  
};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class hull{
 public:
  vector<arrangement> point;
  //vector<structure> struc;
  vector<facet> face;
  //svector< vector<int> > point_to_face;
	
  void sort_conc();
  // bool below_hull(arrangement conf);
	
  //void assemble_coordinate_fenergy();
  void write_hull();  // added by jishnu
  void write_clex_hull();  // added by jishnu	
  void clear_arrays();  // added by jishnu	

  //routine to calculate the cluster expanded energy for each hull point
  //routines that update all coordinates with either the FP energy or the CE energy
	
  //for each point on the hull keep track of all facets that contain that point
  //the chemical potentials that stabilize the facets are the chemical potential bounds for the point
  //have something to make a chemical potential stability map (phase diagram)
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class configurations{
 public:
  structure prim;
  multiplet basiplet;
  vector<superstructure> superstruc;
  vector<arrangement> reference;       // contains the reference states to calculate formation energies
  hull chull;
	
  void generate_configurations(vector<structure> supercells);
  void generate_configurations_fast(vector<structure> supercells); //Added by John
  void generate_vasp_input_directories();
  void print_con();
  void print_con_old();    // added by jishnu
  void read_con();    // added by jishnu
  void read_corr();    // added by jishnu
  void print_corr();
  void read_energy_and_corr(); //added by jishnu
  void print_corr_old();    // added by jishnu
  void print_make_dirs();    // added by jishnu
  void read_make_dirs();    // added by jishnu
  void collect_reference();
  void collect_energies();
  void collect_energies_fast(); //Added by John
  void calculate_formation_energy();
  // void find_CE_formation_energies();  // jishnu
  // void get_CE_hull(); // added by jishnu // this is not the actual hull but the recalculation of FP-hull with CE.
  // void write_below_hull(); // added by jishnu // find distance from CE-hull,write it in below hull,label them as fitting/calculated/non-fitting but calculated etc
  //	// also do the mapbelowhull part in write_below_hull
  void assemble_hull();  // modified and rewritten by jishnu
  void CEfenergy_analysis();  // added by jishnu
  void get_delE_from_hull();  // added by jishnu // you can use this to weight according to distance from hull // right now not being used
  void get_delE_from_hull_w_clexen();  // added by jishnu 
  void print_eci_inputfiles();
  void print_eci_inputfiles_old();  // added by jishnu
  void assemble_coordinate_fenergy();
  void cluster_expanded_energy();
  void reconstruct_from_read_files();   //added by jishnu
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class chempot{
 public:
  vector< vector<double> > m;
  vector< vector<specie> > compon;

  void initialize(concentration conc);      //sets up the vector structure of mu that is compatible with the concentration
  void initialize(vector<vector< specie > > init_compon);
  void set(facet face);                     //set the values of mu to be equal to those that stabilize the given facet
  void increment(chempot dmu);
  void print(ostream &stream);
  void print_compon(ostream &stream);

};



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class trajectory{
 public:
  vector< vector<double> > Rx;
  vector< vector<double> > Ry;
  vector< vector<double> > Rz;
  vector< vector<double> > R2;
  vector< vector<specie> > elements;
  vector< vector<int> > spin;


  void initialize(concentration conc);      //sets up the vector structure of mu that is compatible with the concentration
  void set_zero();
  void increment(trajectory R);
  void normalize(double D);
  void normalize(concentration conc);
  void print(ostream &stream);
  void print_elements(ostream &stream);
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class fluctuation{
 public:
  vector< vector<double> > f;
  vector< vector<specie> > compon;


  void initialize(concentration conc);
  void set_zero();
  void evaluate(concentration conc);
  void evaluate(trajectory R);
  void increment(fluctuation FF);
  void decrement(fluctuation FF);
  void normalize(double n);
  void normalize(double n, concentration conc);
  void print(ostream &stream);
  void print_elements(ostream &stream);
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//hop class collects possible hops for a particular basis site

class hop{
 public:
  int b;                                   // index of the basis
  int vac_spin_init;                       // spin of the vacancy on the initial site
  vector<int> vac_spin;                    // spin of the vacancy on the final site
  atompos initial;                         // initial site of the hop
  vector<cluster> endpoints;               // collects all clusters of possible hops (e.g. all equivalent nearest neighbors)
  vector<vec> jump_vec;
  vector<double> jump_leng;
  vector< mc_index > final;                // shifts of the final states of the hops
  vector< mc_index > activated;            // shifts of the activated states of the hop cluster
  vector< vector< mc_index > > reach;      // list of update sites for each endpoints cluster

  void get_reach(vector<multiplet> montiplet, bool clear, vector<atompos> basis);  // takes the montiplet and determines the reach for each endpoints cluster
  void print_hop_info(ostream &stream);

  // need to allow for repeated application of get_reach in case we have several cluster expansions simultaneously
  // have a boolean flag that clears the reach list, or enlarges it
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class mc_index{
 public:
  string name;
  int l;
  int shift[4];
  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  int num_specie;
  char basis_flag;
  ////////////////////////////////////////////////////////////////////////////////

  void print(ostream &stream);
  void print_name(ostream &stream);
  void print_shift(ostream &stream);
  void create_name(int i);
};




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Monte_Carlo{
 public:
  vector<atompos> basis;
  concentration conc,num_atoms,sublat_conc,num_hops;   // work concentration variable that has the right dimensions
  vector<int> basis_to_conc;            // links each basis site to a concentration unit
  vector< vector<int> > conc_to_basis;  // links each concentration unit to all of its basis sites
		
  vector<int> basis_to_sublat;
  vector< vector<int> > sublat_to_basis;
		
  structure prim;
  structure Monte_Carlo_cell;
  multiplet basiplet;
  vector<multiplet> montiplet;
  vector< vector <hop> > jumps;    // for each basis site, there is a vector of hop objects
  vector<double> AVcorr;		   // Vector containing average correlations
  int di,dj,dk,db,ind1,ind2,ind3;
  int nuc;                         // number of unit cells
  int si,sj,sk;                    // shift variables used in index to insure a positive modulo
  int corr_flag;                   // Flag to indicate whether to calculate average correlations
  int *mcL;                        // unrolled Monte Carlo cell in a linear array
  int *ltoi,*ltoj,*ltok,*ltob;     // given the index of a site in the unrolled Monte Carlo cell mcL, it gives back i,j,k, or b
  double *Rx,*Ry,*Rz;                 // the trajectory of atom at site l - these have same length as mcL
	
  int *s;                          // unrolled shift factors, sequentially for each basis site, in a linear array
  int *nums;                       // number of sites for each cluster
  int *nump;                       // number of products for each eci
  int *mult;                       // multiplicity for each eci
  int *startend;                   // start and end indices for each basis
  double *eci;                     // eci vector
  double eci_empty;
  double *prob;                    //contains the probabilities for all hop events
  double *cprob;                   //contains the cumulative hop probabilities
  double tot_prob;
  int *ptol;                    // p is the index of prob and l is the index of mcL
  int *ptoj;                    // p is the index of prob and j is the index of a jump type for site l in ptol
  int *ptoh;                    // p is the index of prob and h is the index of the hop for site l in ptol
  int *ltosp;                   // l is the index of mcL sp=start p index for hops of site l
  int *arrayi, *arrayj, *arrayk;
	
  int nmcL,neci,ns,nse,np;           // need explicit sizes for all these arrays
  int idum;                        // random number integer used as seed
	
  concentration AVconc,AVnum_atoms,AVsublat_conc;
  fluctuation AVSusc,Susc;
  double AVenergy,heatcap,flipfreq;
	
  fluctuation AVkinL,kinL;
  trajectory Dtrace,corrfac;
  trajectory AVDtrace,AVcorrfac;
  trajectory R;
  double hop_leng;
	
  //Thermfac = inverse of Susc matrix (for each distinct sublattice) (inverse taken after multiplying Susc by kT)
	
  Monte_Carlo(structure in_prim, structure in_struc, multiplet in_basiplet, int idim, int jdim, int kdim);
  ~Monte_Carlo(){delete [] mcL; delete [] eci; delete [] s; delete [] nums; delete [] nump; delete [] mult; delete [] startend;
    delete [] arrayi; delete [] arrayj; delete[] arrayk;
  };
	
  void collect_basis();
  void collect_sublat();                 //maps basis sites onto crystallographically distinct sublattice sites
  void assemble_conc_basis_links();      //assembles the basis_to_conc and conc_to_basis vectors
  void update_mu(chempot mu);
	
	
  inline int index(int i, int j, int k, int b){
    return (b)*ind1+mdi(i,di)*ind2+mdj(j,dj)*ind3+mdk(k,dk);
  };
  inline int mdi(int i, int di) {return arrayi[i>=0 ? i : i+di];};  // -di <= i <= 2di - 1, e.g. di = 12, i is [-12, 23]
  inline int mdj(int j, int dj) {return arrayj[j>=0 ? j : j+dj];};
  inline int mdk(int k, int dk) {return arrayk[k>=0 ? k : k+dk];};
  inline void invert_index(int &i, int &j, int &k, int &b, int l);
  void generate_eci_arrays();
  void write_point_energy(ofstream &out);
  void write_point_corr(ofstream &out); //Added by John
  void write_normalized_point_energy(ofstream &out);	
  void write_monte_h(string class_file);
  void write_monte_xyz(ostream &stream);
  void write_monte_poscar(ostream &stream);
	
	
  double pointenergy(int i, int j, int k, int b);
  void  pointcorr(int i, int j, int k, int b);
  double normalized_pointenergy(int i, int j, int k, int b);
  void calc_energy(double &energy);
  void calc_concentration();
  void calc_num_atoms();
  void calc_sublat_concentration();
  void update_num_hops(int l, int ll, int b, int bb);
  double calc_grand_canonical_energy(chempot mu);
	
	
	
	
  //conventional Monte Carlo routines
  bool compatible(structure init_struc, int &ndi, int &ndj, int &ndk);
  void initialize(structure init_struc);
  void initialize(concentration conc);
  ////////////////////////////////////////////////////////////////////////////////
  //added by Ben Swoboda
  void initialize_1vac(concentration conc);  //1 vacancy in monte cell
  void initialize_1_specie(double in_conc) ; //Added by Aziz
  ////////////////////////////////////////////////////////////////////////////////
  void grand_canonical(double beta, chempot mu, int n_pass, int n_equil_pass);
  void canonical_1_species(double beta, int Temp, int n_pass, int n_equil_pass, int n_pass_output); //Added by Aziz
  double lte(double beta, chempot mu);
  //  void canonical(double beta, double n_pass, double n_equil_pass);
  //  void n_fold_grand_canonical(double beta, chempot mu, double n_pass, double n_equil_pass);
	
	
  //kinetic Monte Carlo routines
  void initialize_kmc();
  void extend_reach();
  void get_hop_prob(int i, int j, int k, int b, double beta);
  double calc_barrier(int i, int j, int k, int b, int ii, int jj, int kk, int bb, int l, int ll, int ht, int h);
  void initialize_prob(double beta);
  int pick_hop();
  void update_prob(int i, int j, int k, int b, int ht, int h, double beta);
  void kinetic(double beta, double n_pass, double n_equil_pass);
	
  void collect_R();
	
  void output_Monte_Carlo_cell();
	
};




//********************************************************************
//Routines

double determinant(double mat[3][3]);
void inverse(double mat[3][3], double invmat[3][3]);
void matrix_mult(double mat1[3][3], double mat2[3][3], double mat3[3][3]);
void get_perp(double vec1[3], double vec2[3]); //Added by John
void get_perp(double vec1[3], double vec2[3], double vec3[3]); //Added by John
bool normalize(double vec1[3], double length); //Added by John -- returns false if vector is null
void lat_dimension(double lat[3][3], double radius, int dim[3]);
bool compare(double mat1[3][3], double mat2[3][3]);
bool compare(double vec1[3], double vec2[3]);
bool compare(double vec1[3], double vec2[3], int trans[3]);
bool compare(vector<double> vec1, vector<double> vec2);
bool compare(char name1[2], char name2[2]);
bool compare(specie compon1, specie compon2);
bool compare(vector<specie> compon1, vector<specie> compon2);
bool compare(atompos &atom1, atompos &atom2);
bool compare(atompos atom1, atompos atom2, int trans[3]);
bool compare(cluster &clust1, cluster &clust2);
////////////////////////////////////////////////////////////////////////////////
//added by anton
bool compare(orbit orb1, orbit orb2);
////////////////////////////////////////////////////////////////////////////////
bool compare(concentration conc1, concentration conc2);
bool compare(mc_index m1, mc_index m2);
bool new_mc_index(vector<mc_index> v1, mc_index m2);
bool is_integer(double vec[3]);
bool is_integer(double mat[3][3]);
void within(double fcoord[3]);
void within(atompos &atom);
void within(cluster &clust);
void within(cluster &clust, int n);
////////////////////////////////////////////////////////////////////////////////
//added by Ben Swoboda
void within(structure &struc);
////////////////////////////////////////////////////////////////////////////////
void latticeparam(double lat[3][3], double latparam[3], double latangle[3]);
void latticeparam(double lat[3][3], double latparam[3], double latangle[3], int permut[3]);
void conv_AtoB(double AtoB[3][3], double Acoord[3], double Bcoord[3]);
double distance(atompos atom1, atompos atom2);
bool update_bit(vector<int> max_bit, vector<int> &bit, int &last);
void get_equiv(orbit &orb, vector<sym_op> &op);
bool new_clust(cluster clust, orbit &orb);
bool new_clust(cluster clust, vector<orbit> &orbvec);
void get_loc_equiv(orbit &orb, vector<sym_op> &op);
bool new_loc_clust(cluster clust, orbit orb);
bool new_loc_clust(cluster clust, vector<orbit> torbvec);
void calc_correlations(structure struc, multiplet super_basiplet, arrangement &conf);
void get_super_basis_vec(structure &superstruc, vector < vector < vector < int > > > &super_basis_vec); //Added by John
void get_corr_vector(structure &struc, multiplet &super_basiplet, vector < vector < vector < vector < int > > > > &corr_vector); //Added by John
bool new_conf(arrangement &conf,superstructure &superstruc);
bool new_conf(arrangement &conf,vector<superstructure> &superstruc);
void get_shift(atompos &atom, vector<atompos> basis);
void double_to_string(double n, string &a, int dec_places); //Added by John
void int_to_string(int i, string &a, int base);
void generate_ext_clust(structure struc, int min_num_compon, int max_num_points,vector<double> max_radius, multiplet &clustiplet);
void generate_ext_basis(structure struc, multiplet clustiplet, multiplet &basiplet);
void generate_ext_monteclust(vector<atompos> basis, multiplet basiplet, vector<multiplet> &montiplet);
////////////////////////////////////////////////////////////////////////////////
//added by anton - filters a multiplet for clusters containing just one activated site (with occupation basis = 1)
void filter_activated_clust(multiplet clustiplet);
void merge_multiplets(multiplet clustiplet1, multiplet clustiplet2, multiplet &clustiplet3); 
void write_clust(multiplet clustiplet, string out_file);
void write_fclust(multiplet clustiplet, string out_file);
////////////////////////////////////////////////////////////////////////////////

bool scandirectory(string dirname, string filename);
bool read_oszicar(string dirname, double& e0);
bool read_oszicar(string dirname, double& e0, int &count);   // added by jishnu
bool read_mc_input(string cond_file, int &n_pass, int &n_equil_pass, int &nx, int &ny, int &nz, chempot &muinit, chempot &mu_min, chempot &mu_max, vector<chempot> &muinc, double $Tinit, double &Tmin, double &Tmax, double &Tinc, int &xyz_step, int &corr_flag, int &temp_chem);
double ran0(int &idum);

////////////////////////////////////////////////////////////////////////////////
//added by Ben Swoboda to utilize multiple basis
void get_clust_func(atompos atom1, atompos atom2, double &clust_func);
void get_basis_vectors(atompos &atom);
////////////////////////////////////////////////////////////////////////////////
void read_junk(istream &stream);  // added by jishnu
//int hullfinder_bi(double); // added by jishnu // this generates binary hull and cutts off the high energy structres // made member functions of cofigurations
// int hullfinder_ter(double); // added by jishnu // this generates ternary hull and cutts off the high energy structres // made member functions of cofigurations

//************************************************************

double vec::calc_dist(){
  if(cart_on==false){
    cout << "no cartesian coordinates for vec \n";
    exit(1);
  }

  length=0;
  for(int i=0; i<3; i++)
    length=length+ccoord[i]*ccoord[i];
  length=sqrt(length);

  return length;
}


//************************************************************

vec vec::apply_sym(sym_op op){
  int i,j,k;
  vec tlat;

  if(op.frac_on == false || op.cart_on == false)op.update();


  for(i=0; i<3; i++){
    tlat.fcoord[i]=op.ftau[i];
    tlat.ccoord[i]=op.ctau[i];
    for(j=0; j<3; j++){
      tlat.fcoord[i]=tlat.fcoord[i]+op.fsym_mat[i][j]*fcoord[j];
      tlat.ccoord[i]=tlat.ccoord[i]+op.csym_mat[i][j]*ccoord[j];
    }
  }
  return tlat;
}


//************************************************************

void vec::print_frac(ostream &stream){

  for(int i=0; i<3; i++){
    stream.precision(5);stream.width(10);stream.setf(ios::showpoint);
    stream << fcoord[i] << " ";
  }
  stream << "\n";
}

//************************************************************

void vec::print_cart(ostream &stream){

  for(int i=0; i<3; i++){
    stream.precision(5);stream.width(10);stream.setf(ios::showpoint);
    stream << ccoord[i] << " ";
  }
  stream << "\n";
}

//************************************************************
//Edited by John
tensor tensor::apply_sym(sym_op op){
  // calculates the op.csym_mat^T * K * op.csym_mat

  tensor ttensor(rank, dim);
  bool D_flag=true;
  for(int i=0; i<rank; i++) D_flag=D_flag&&(dim[i]==3);
  if(!D_flag){
    cout << "ERROR:  Attempting apply 3D transformation to tensor with improper dimensionality.\n";
    return ttensor;
  }
  int *elem_array=new int[rank];
  int *sub_array=new int[rank];
  bool elem_cont=true;
  for(int i=0; i<rank; i++)
    elem_array[i]=0;
  while(elem_cont){
    double telem=0.0;
    for(int i=0; i<rank; i++)
      sub_array[i]=0;
    bool sub_cont=true;
    while(sub_cont){
      double ttelem=1.0;
      for(int i=0; i<rank; i++){
        ttelem*=op.csym_mat[elem_array[i]][sub_array[i]];
      }
      ttelem*=get_elem(sub_array);
      telem+=ttelem;
      for(int i=0; i<rank; i++){
	sub_array[i]+=1;
	if(sub_array[i]==dim[i])
          sub_array[i]=0;
	else break;
        if(i==rank-1) sub_cont=false;
      }
    }
    ttensor.set_elem(telem, elem_array);

    for(int i=0; i<rank; i++){
      elem_array[i]=elem_array[i]+1;
      if(elem_array[i]>=dim[i])
        elem_array[i]=0;
      else break;
      if(i==rank-1) elem_cont=false;
    }
  }
  return ttensor;
}
//\Edited by John    

//************************************************************

//Added by John
//Copy constructor
tensor::tensor(tensor const& ttens){
  rank=ttens.rank;
  size=ttens.size;
  dim=new int[rank];
  mult=new int[rank];
  K=new double[size];
  for(int i=0; i<rank; i++){
    dim[i]=ttens.dim[i];
    mult[i]=ttens.mult[i];
  }
  for(int j=0; j<size; j++)
    K[j]=ttens.K[j];
}

//************************************************************
  
//Added by John
//Default constructor
tensor::tensor(){
  dim=0;
  K=0;
  mult=0;
}
  
//************************************************************

//Added by John
tensor::tensor(int trank, ...){
  size=1;
  rank=trank;
  va_list argPtr;
  va_start( argPtr, trank );
  dim=new int[rank];
  mult=new int[rank];
  for(int i=0; i<rank; i++){
    dim[i]=va_arg( argPtr, int );
    mult[i]=size;
    size*=dim[i];
  }
  va_end(argPtr);
  K=new double[size];
  for(int i=0; i<size; i++)
    K[i]=0.0;
}

//\Added by John

//************************************************************

//Added by John
tensor::tensor(int trank, int *tdim){
  size=1;
  rank=trank;
  dim=new int[rank];
  mult=new int[rank];
  for(int i=0; i<rank; i++){
    dim[i]=tdim[i];
    mult[i]=size;
    size*=dim[i];
  }
  K=new double[size];
  for(int i=0; i<size; i++)
    K[i]=0.0;
}

//\Added by John

//************************************************************

//Added by John
//Assignment operator
tensor& tensor::operator=(const tensor& ttens){
  rank = ttens.rank;
  size = ttens.size;
  delete dim;
  delete mult;
  delete K;
  dim = new int[rank];
  mult = new int[rank];
  K = new double[size];
  for(int i=0; i<rank; i++){
    dim[i]=ttens.dim[i];
    mult[i]=ttens.mult[i];
  }
  for(int i=0; i<size; i++)
    K[i]=ttens.K[i];

  return *this;
}

//************************************************************

//Added by John
  double tensor::get_elem(int ind, ...){
    va_list argPtr;
    va_start(argPtr, ind);
    if(!(ind>-1 && ind<dim[0])){
      cout << "WARNING:  Attempted to acess tensor element out of bounds.";
      return 0.0;
    }
    int ctr=ind*mult[0];

    for(int i=1; i<rank; i++){
      int tind=va_arg( argPtr, int );
      if(tind>-1 && tind <dim[i])
	ctr+=tind*mult[i];
      else{
	cout << "WARNING:  Attempted to acess tensor element out of bounds.";
	return 0.0;
      }
    }
    va_end(argPtr);
    return K[ctr];
  }
//\Added by John

//************************************************************

//Added by John
double tensor::get_elem(int *inds){
  int ctr=0;
  for(int i=0; i<rank; i++){
    if(inds[i]>-1 && inds[i] <dim[i])
      ctr+=inds[i]*mult[i];
    else{
      cout << "WARNING:  Attempted to acess tensor element out of bounds.";
      return 0.0;
    }
  }
  return K[ctr];
}
//\Added by John
  
//************************************************************
//Added by John         
void tensor::set_elem(double new_elem, ...){        
  va_list argPtr;       
  va_start(argPtr, new_elem);                         
  int ctr=0;            
  for(int i=0; i<rank; i++){                          
    int tind=va_arg(argPtr, int);                     
    if(tind>-1 && tind <dim[i])                       
      ctr+=tind*mult[i];                              
    else{               
      cout << "WARNING:  Attempted to acess tensor element out of bounds.";         
      return;           
    }                   
  }                     
  va_end(argPtr);       
  K[ctr]=new_elem;      
  return;               
}
//\Added by John

//************************************************************      
//Added by John               
void tensor::set_elem(double new_elem, int *inds){
  int ctr=0;
  //  cout << "Inside set_elem, new_elem =" << new_elem << ";  rank =" << rank << "\n";                   
  for(int i=0; i<rank; i++){

    if(inds[i]>=0 && inds[i] <dim[i]){
      ctr=ctr+inds[i]*mult[i];
      //      cout << "ctr=" << ctr << " and inds[" << i << "]=" << inds[i] << "\n";                      
    }
    else{
      cout << "WARNING:  Attempted to acess tensor element out of bounds.";
      return;
    }
  }
  //  cout << "\n";           
  K[ctr]=new_elem;
  return;
}
//\Added by John 

//************************************************************                                                                                  

//Added by John                                                                                                                                 
void tensor::print(ostream &stream){

  for(int i=0; i<size; i++){
    stream << "   " << K[i];
    for(int j=0; j<rank; j++){
      if(!((i+1)%mult[j])&&mult[j]!=1)
        stream << "\n";
    }
  }
}
//\Addition                                                                                                                                     


//************************************************************
//Added by John
void sym_op::get_sym_type(){
  
  int i, j;
  double vec_sum, vec_mag;
  if(!cart_on){
    update();
  }

  double det=0.0;
  double trace=0.0;
  double tmat[3][3];

  //Copy csym_mat to temporary
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      tmat[i][j]=csym_mat[i][j];
    }
  }
  //Get Determinant and trace
  det=determinant(tmat);
  
  for(i=0; i<3; i++){
    trace+=tmat[i][i];
  }


  if(abs(trace-3.0)<tol){ //Sym_op is identity
    sym_type=0;
    return;
  }

  if(abs(trace+3.0)<tol){ //Sym_op is inversion
    sym_type=-1;
    return;
  }
  

  if(det<0 && abs(trace-1.0)<tol){ //operation is mirror
    //The trace criterion can be shown by noting that a mirror
    //is a 180 degree rotation composed with an inversion

    sym_type=2;

    //Mirror planes have eigenvalues 1,1,-1; the Eigenvectors form an
    //orthonormal basis into which any member of R^3 can be decomposed.
    //For mirror operation S and test vector v, we take the vector
    //w=v-S*v to be the eigenvector with eigenvalue -1.  We must test as many
    //as 3 cases to ensure that test vector v is not coplanar with mirror
    double vec1[3], vec2[3];
    vec1[0]=1;  vec1[1]=1;  vec1[2]=1;
    conv_AtoB(tmat, vec1, vec2);
    vec_sum=0.0;
    for(i=0; i<3; i++){
      vec2[i]-=vec1[i];
      if(abs(vec_sum)<tol && abs(vec2[i])>tol)
	vec_sum=vec2[i]/abs(vec2[i]);
      eigenvec.ccoord[i]=vec2[i];
    }
    if(normalize(eigenvec.ccoord, vec_sum)){
      conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
      eigenvec.frac_on=true;
      eigenvec.cart_on=true;
      return;
    }
    
    vec1[2]=-2;
    conv_AtoB(tmat, vec1, vec2);
    vec_sum=0.0;
    for(i=0; i<3; i++){
      vec2[i]-=vec1[i];
      if(abs(vec_sum)<tol && abs(vec2[i])>tol)
	vec_sum=vec2[i]/abs(vec2[i]);
      eigenvec.ccoord[i]=vec2[i];
    }
    if(normalize(eigenvec.ccoord, vec_sum)){
      conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
      eigenvec.frac_on=true;
      eigenvec.cart_on=true;
      return;
    }
    
    eigenvec.ccoord[0]=1.0/sqrt(2);
    eigenvec.ccoord[1]=-1.0/sqrt(2);
    eigenvec.ccoord[2]=0.0;
    conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
    eigenvec.frac_on=true;
    eigenvec.cart_on=true;
    return;
  }//\End Mirror Plane Conditions

  else { //operation is either rotation or rotoinversion
    if(det<(-tol)){ //operation is rotoinversion
      trace*=-1;
      det*=-1;
      sym_type=3;
      //compose rotoinversion with inversion so that angle and axis may be extracted
      //in same was as proper rotation
      for(i=0; i<3; i++){
	for(j=0; j<3; j++){
	  tmat[i][j]*=-1;
	}
      }
    }
    else sym_type=1; //operation is rotation
    
    if(abs(trace+1)<tol){ //rotation is 180 degrees, which requires special care, since rotation matrix becomes symmetric
      //180 rotation can be decomposed into two orthogonal mirror planes, so we use similar method as above, but finding +1 eigenvector
      double vec1[3], vec2[3];
      op_angle=180;
      vec1[0]=1;  vec1[1]=1;  vec1[2]=1;
      conv_AtoB(tmat, vec1, vec2);
      double vec_sum=0.0;
      for(i=0; i<3; i++){
	vec2[i]+=vec1[i];
	if(abs(vec_sum)<tol && abs(vec2[i])>tol)
	  vec_sum=vec2[i]/abs(vec2[i]);
	eigenvec.ccoord[i]=vec2[i];
      }
      if(normalize(eigenvec.ccoord, vec_sum)){
	conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
	eigenvec.frac_on=true;
	eigenvec.cart_on=true;
	return;
      }

      vec1[2]=-2;
      conv_AtoB(tmat, vec1, vec2);
      vec_sum=0.0;
      for(i=0; i<3; i++){
	vec2[i]+=vec1[i];
	if(abs(vec_sum)<tol && abs(vec2[i])>tol)
	  vec_sum=vec2[i]/abs(vec2[i]);
	eigenvec.ccoord[i]=vec2[i];
      }
      if(normalize(eigenvec.ccoord, vec_sum)){
	conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
	eigenvec.frac_on=true;
	eigenvec.cart_on=true;
	return;
      }

      eigenvec.ccoord[0]=1.0/sqrt(2);
      eigenvec.ccoord[1]=-1.0/sqrt(2);
      eigenvec.ccoord[2]=0.0;
      conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
      eigenvec.frac_on=true;
      eigenvec.cart_on=true;
      return;
    }
    
    // Following only evaluates if we have non-180 proper rotation
    // Method uses inversion of axis-angle interpretation of a rotation matrix R
    // With axis v=(x,y,z) and angle ÅŒÅ∏, with ||v||=1
    //  c = cos(ÅŒÅ∏); s = sin(ÅŒÅ∏); C = 1-c
    //      [ x*xC+c   xyC-zs   zxC+ys ]
    //  R = [ xyC+zs   y*yC+c   yzC-xs ]
    //      [ zxC-ys   yzC+xs   z*zC+c ]
    double tangle;
    vec_sum=0.0;
    vec_mag=0.0;
    for(i=0; i<3; i++){
      eigenvec.ccoord[i]=tmat[(i+2)%3][(i+1)%3]-tmat[(i+1)%3][(i+2)%3];
      vec_mag += eigenvec.ccoord[i]*eigenvec.ccoord[i];
      if(abs(vec_sum)<tol && abs(eigenvec.ccoord[i])>tol)
	vec_sum=eigenvec.ccoord[i]/abs(eigenvec.ccoord[i]);
    }
    vec_mag=sqrt(vec_mag);
    tangle=round((180.0/3.141592654)*atan2(vec_mag,trace-1));
    op_angle=int(tangle);
    normalize(eigenvec.ccoord,vec_sum);
    if(vec_sum<0){
      op_angle=360-op_angle;
    }

    conv_AtoB(CtoF, eigenvec.ccoord, eigenvec.fcoord);
    eigenvec.frac_on=true;
    eigenvec.cart_on=true;
    return;
  }
}
//\End Addition

//************************************************************

void sym_op::print_fsym_mat(ostream &stream){
  //Added by John
  
  if(sym_type==-1) stream << "Inversion Operation: \n";
  if(!sym_type) stream << "Identity Operation: \n";
  if(sym_type==1){
    stream << op_angle << " degree Rotation (or screw) Operation about axis: ";
    eigenvec.print_frac(stream);
  }
  if(sym_type==2){
    stream << "Mirror (or glide) Operation with plane normal: ";
    eigenvec.print_frac(stream);
  }
  if(sym_type==3){
    stream << op_angle << " degree Rotoinversion (or screw) Operation about axis: ";
    eigenvec.print_frac(stream);
  }
  //\End Addition
  stream << "        symmetry operation matrix                  shift \n";
  for(int i=0; i<3; i++){
    stream.precision(9);stream.width(15);stream.setf(ios::showpoint);
    for(int j=0; j<3; j++) stream << fsym_mat[i][j] << "  ";
    stream << "       " << ftau[i] << "\n";
  }
}

//************************************************************

void sym_op::print_csym_mat(ostream &stream){
  //Added by John
  if(sym_type==-1) stream << "Inversion Operation: \n";
  if(!sym_type) stream << "Identity Operation: \n";
  if(sym_type==1){
    stream << op_angle << " degree Rotation (or screw) Operation about axis: ";
    eigenvec.print_cart(stream);
  }
  if(sym_type==2){
    stream << "Mirror (or glide) Operation with plane normal: ";
    eigenvec.print_cart(stream);
  }
  if(sym_type==3){
    stream << op_angle << " degree Rotoinversion (or screw) Operation about axis: ";
    eigenvec.print_cart(stream);
  }
  //\End Addition
  stream << "        symmetry operation matrix                  shift \n";
  for(int i=0; i<3; i++){
    stream.precision(9);stream.width(15);stream.setf(ios::showpoint);
    for(int j=0; j<3; j++) stream << csym_mat[i][j] << "  ";
    stream << "       " << ctau[i] << "\n";
  }
}


//************************************************************

void sym_op::get_csym_mat(){
  int i,j,k;
  double temp[3][3];


  if(cart_on == true) return;

  if(cart_on == false) {
    if(frac_on == false){
      cout << "No sym_op initialized - cannot get_csym_mat\n";
      return;
    }

    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	temp[i][j]=0;
	for(k=0; k<3; k++) temp[i][j]=temp[i][j]+fsym_mat[i][k]*CtoF[k][j];
      }
    }

    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	csym_mat[i][j]=0.0;
	for(k=0; k<3; k++) csym_mat[i][j]=csym_mat[i][j]+FtoC[i][k]*temp[k][j];
      }
    }

    for(i=0; i<3; i++){
      ctau[i]=0.0;
      for(j=0; j<3; j++) ctau[i]=ctau[i]+FtoC[i][j]*ftau[j];
    }
  }
  cart_on=true;
  return;
}



//************************************************************

void sym_op::get_fsym_mat(){
  int i,j,k;
  double temp[3][3];

  if(frac_on == true) return;
  if(frac_on == false) {
    if(cart_on == false){
      cout << "No sym_op initialized - cannot get_fsym_mat\n";
      return;
    }

    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	temp[i][j]=0;
	for(k=0; k<3; k++) temp[i][j]=temp[i][j]+csym_mat[i][k]*FtoC[k][j];
      }
    }

    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	fsym_mat[i][j]=0.0;
	for(k=0; k<3; k++) fsym_mat[i][j]=fsym_mat[i][j]+CtoF[i][k]*temp[k][j];
      }
    }

    for(i=0; i<3; i++){
      ftau[i]=0.0;
      for(j=0; j<3; j++) ftau[i]=ftau[i]+CtoF[i][j]*ctau[j];
    }
  }
  frac_on=true;
  return;
}


//************************************************************

void sym_op::update(){
  get_trans_mat();
  get_csym_mat();
  get_fsym_mat();
  return;
}


//************************************************************
structure::structure(){
  int i,j;
  for(i=0; i<200; i++) title[i]=0;
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      lat[i][j]=0.0;
      slat[i][j]=0.0;
      FtoC[i][j]=0.0;
      CtoF[i][j]=0.0;
    }
    slat[i][i]=1.0;
  }

  frac_on=false;
  cart_on=false;

}


//************************************************************
void structure::get_latparam(){

  latticeparam(lat,latparam,latangle,permut);

  return;

}


//************************************************************
void structure::get_ideal_latparam(){

  latticeparam(ilat,ilatparam,ilatangle,ipermut);

  return;

}


//************************************************************
//read a POSCAR like file and collect all the structure variable

void structure::read_lat_poscar(istream &stream){
  int i,j;
  char ch;
  char buff[200];
  char sp[200];


  stream.getline(title,199);
  stream >> scale;
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      stream >> lat[i][j];
    }
  }

  //normalize the scale to 1.0 and adjust lattice vectors accordingly

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)lat[i][j]=scale*lat[i][j];
  scale=1.0;

  return;

}



//************************************************************
//write the structure to a file in POSCAR like format

void structure::write_lat_poscar(ostream &stream) {
  int i,j;

  stream << title <<"\n";

  stream.precision(7);stream.width(12);stream.setf(ios::showpoint);
  stream << scale <<"\n";

  for(int i=0; i<3; i++){
    stream << "  ";
    for(int j=0; j<3; j++){

      stream.precision(9);stream.width(15);stream.setf(ios::showpoint);
      stream << lat[i][j] << " ";

    }
    stream << "\n";
  }

  return;

}




//************************************************************

void structure::calc_point_group(){
  int i,j,k;
  int dim[3];
  double radius;
  vec temp;
  vector<vec> gridlat;
  double tlat[3][3],tlatparam[3],tlatangle[3];
  double tcsym_mat[3][3];
  sym_op tpoint_group;

  get_latparam();
  int num_point_group=0;

  //make a lattice parallelepiped that encompasses a sphere with radius = 2.1*largest latparam

  radius=1.5*latparam[permut[2]];

  lat_dimension(lat,radius,dim);

  for(i=0; i<3; i++){
    if(dim[i] > 3) dim[i]=3;
  }

  cout << "inside calc_point group \n";
  cout << "dimension = ";
  for(i=0; i<3; i++)cout << dim[i] << " ";
  cout << "\n";

  for(i=-dim[0]; i<=dim[0]; i++){
    for(j=-dim[1]; j<=dim[1]; j++){
      for(k=-dim[2]; k<=dim[2]; k++){
	if(!(i==0 && j==0 && k==0)){
	  temp.ccoord[0]=i*lat[0][0]+j*lat[1][0]+k*lat[2][0];
	  temp.ccoord[1]=i*lat[0][1]+j*lat[1][1]+k*lat[2][1];
	  temp.ccoord[2]=i*lat[0][2]+j*lat[1][2]+k*lat[2][2];
	  temp.cart_on=true;
	  temp.calc_dist();
	}

	//keep only the lattice points within the sphere with radius

	if(temp.length < radius){
	  gridlat.push_back(temp);
	}
      }
    }
  }


  cout << "made the grid \n";
  cout << "number of sites in the grid = " << gridlat.size() << "\n";

  //for each set of three lattice points within the sphere see which one has the
  //same sets of lengths and angles as the original lattice unit cell vectors.

  for(i=0; i<gridlat.size(); i++){
    for(j=0; j<gridlat.size(); j++){
      for(k=0; k<gridlat.size(); k++){

	if(i!=j && i!=k && j!=k){
	  for(int ii=0; ii<3; ii++){
	    tlat[0][ii]=gridlat[i].ccoord[ii];
	    tlat[1][ii]=gridlat[j].ccoord[ii];
	    tlat[2][ii]=gridlat[k].ccoord[ii];
	  }
	  latticeparam(tlat,tlatparam,tlatangle);

	  //compare the tlat... and lat... to see if they are the same lattice
	  // that is do the lattice vectors have the same lengths and the same angles

	  if(abs(latparam[0]-tlatparam[0]) < tol*latparam[permut[0]] &&
	     abs(latparam[1]-tlatparam[1]) < tol*latparam[permut[0]] &&
	     abs(latparam[2]-tlatparam[2]) < tol*latparam[permut[0]] &&
	     abs(latangle[0]-tlatangle[0]) < tol &&
	     abs(latangle[1]-tlatangle[1]) < tol &&
	     abs(latangle[2]-tlatangle[2]) < tol ){

	    // get the matrix that relates the two lattice vectors


	    for(int ii=0; ii<3; ii++){
	      for(int jj=0; jj<3; jj++){
		tcsym_mat[ii][jj]=0.0;
		for(int kk=0; kk<3; kk++)
		  tcsym_mat[ii][jj]=tcsym_mat[ii][jj]+tlat[kk][ii]*CtoF[kk][jj];
	      }
	    }

	    // check whether this symmetry operation is new or not

	    int ll=0;
	    for(int ii=0; ii<num_point_group; ii++)
	      if(compare(tcsym_mat,point_group[ii].csym_mat))break;
	      else ll++;

	    // if the symmetry operation is new, add it to the pointgroup array
	    // and update all info about the sym_op object

	    if(num_point_group == 0 || ll == num_point_group){
	      for(int jj=0; jj<3; jj++){
		tpoint_group.frac_on=false;
		tpoint_group.cart_on=false;
		tpoint_group.ctau[jj]=0.0;
		for(int kk=0; kk<3; kk++){
		  tpoint_group.csym_mat[jj][kk]=tcsym_mat[jj][kk];
		  tpoint_group.lat[jj][kk]=lat[jj][kk];
		}
	      }

	      tpoint_group.cart_on=true;
	      tpoint_group.update();
 	      tpoint_group.get_sym_type(); // Added by John
	      point_group.push_back(tpoint_group);
	      num_point_group++;
	    }


	  }
	}
      }
    }
  }
  cout << "finished finding all point group operartions \n";
}


//************************************************************
void structure::update_lat(){
  get_trans_mat();
  get_latparam();
  calc_point_group();
  return;
}


//************************************************************
void structure::write_point_group(){
  int pg;

  ofstream out("point_group");
  if(!out){
    cout << "Cannot open file.\n";
    return;
  }

  cout << " number of point group ops " << point_group.size() << "\n";

  for(pg=0; pg<point_group.size(); pg++){
    out << "point group operation " << pg << " \n";
    point_group[pg].print_fsym_mat(out);
    out << "\n";
  }

  out.close();
}


//************************************************************
void structure::generate_3d_supercells(vector<structure> &supercell, int max_vol){
  int vol,pg,i,j,k;
  int tslat[3][3];


  //algorithm relayed to me by Tim Mueller
  //make upper triangular matrix where the product of the diagonal elements equals the volume
  //then for the elements above the diagonal choose all values less than the diagonal element
  //for each lattice obtained this way, apply point group symmetry operations
  //see if the transformed superlattice can be written as a linear combination of superlattices already found
  //if not add it to the list


  for(vol = 1; vol <= max_vol; vol++){
    vector<structure> tsupercell;

    //generate all tslat[][] matrices that are upper diagonal where the product of the
    //diagonal equals the current volume vol and where the upper diagonal elements take
    //all values less than the diagonal below it

    //initialize the superlattice vectors to zero
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)tslat[i][j]=0;

    for(tslat[0][0]=1; tslat[0][0]<=vol; tslat[0][0]++){
      if(vol%tslat[0][0] == 0){
	for(tslat[1][1]=1; tslat[1][1]<=vol/tslat[0][0]; tslat[1][1]++){
	  if((vol/tslat[0][0])%tslat[1][1] == 0){
	    tslat[2][2]=(vol/tslat[0][0])/tslat[1][1];
	    for(tslat[0][1]=0; tslat[0][1]<tslat[1][1]; tslat[0][1]++){
	      for(tslat[0][2]=0; tslat[0][2]<tslat[2][2]; tslat[0][2]++){
		for(tslat[1][2]=0; tslat[1][2]<tslat[2][2]; tslat[1][2]++){

		  //copy the superlattice vectors into lattice_point objects
		  //and get their cartesian coordinates

		  vec lat_point[3];

		  for(i=0; i<3; i++){
		    for(j=0; j<3; j++) lat_point[i].fcoord[j]=tslat[i][j];
		    lat_point[i].frac_on=true;
		    conv_AtoB(FtoC, lat_point[i].fcoord, lat_point[i].ccoord);
		    lat_point[i].cart_on=true;
		  }

		  //if no supercells have been added for this volume, add it
		  //else perform all point group operations to the supercell and
		  //see if it is a linear combination of already found supercells with the same volume

		  if(tsupercell.size() == 0){
		    structure tsup_lat;
		    strcpy(tsup_lat.title,"supercell of ");
		    int leng=strlen(tsup_lat.title);
		    for(i=0; title[i]!=0 && i<199-leng; i++)tsup_lat.title[i+leng]=title[i];
		    tsup_lat.scale = scale;
		    for(i=0; i<3; i++){
		      for(j=0; j<3; j++){
			tsup_lat.lat[i][j]= lat_point[i].ccoord[j];
			tsup_lat.slat[i][j]= lat_point[i].fcoord[j];
		      }
		    }

		    tsup_lat.get_trans_mat();
		    tsupercell.push_back(tsup_lat);


		  }
		  else{
		    //apply each point group to this superlattice stored in lat_point

		    int num_miss1=0;
		    for(pg=0; pg<point_group.size(); pg++){
		      vec tlat_point[3];

		      for(i=0; i<3; i++) tlat_point[i]=lat_point[i].apply_sym(point_group[pg]);

		      //see if tlat_point[] can be expressed as a linear combination of any
		      //superlattice with volume vol already found

		      int num_miss2=0;
		      for(int ts=0; ts<tsupercell.size(); ts++){
			double lin_com[3][3];    //contains the coefficients relating the two lattices
			for(i=0; i<3; i++){
			  for(j=0; j<3; j++){
			    lin_com[i][j]=0.0;
			    for(k=0; k<3; k++)
			      lin_com[i][j]=lin_com[i][j]+tsupercell[ts].CtoF[i][k]*tlat_point[j].ccoord[k];
			  }
			}

			//check whether lin_com[][] are strictly integer
			//if so, the transformed superlattice is a linear combination of a previous one

			if(!is_integer(lin_com))num_miss2++;
		      }
		      if(num_miss2 == tsupercell.size())num_miss1++;
		    }
		    if(num_miss1 == point_group.size()){
		      structure tsup_lat;
		      strcpy(tsup_lat.title,"supercell of ");
		      int leng=strlen(tsup_lat.title);
		      for(i=0; title[i]!=0 && i<199-leng; i++)tsup_lat.title[i+leng]=title[i];
		      tsup_lat.scale = scale;
		      for(i=0; i<3; i++){
			for(j=0; j<3; j++){
			  tsup_lat.lat[i][j]= lat_point[i].ccoord[j];
			  tsup_lat.slat[i][j]= lat_point[i].fcoord[j];
			}
		      }
		      tsup_lat.get_trans_mat();
		      tsupercell.push_back(tsup_lat);

		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    for(i=0; i<tsupercell.size(); i++)supercell.push_back(tsupercell[i]);
  }
  return;
}



//************************************************************
void structure::generate_2d_supercells(vector<structure> &supercell, int max_vol, int excluded_axis){
  int vol,pg,i,j,k;
  int tslat[3][3];
  int tslat2d[2][2];


  //algorithm relayed to me by Tim Mueller
  //make upper triangular matrix where the product of the diagonal elements equals the volume
  //then for the elements above the diagonal choose all values less than the diagonal element
  //for each lattice obtained this way, apply point group symmetry operations
  //see if the transformed superlattice can be written as a linear combination of superlattices already found
  //if not add it to the list


  for(vol = 1; vol <= max_vol; vol++){
    vector<structure> tsupercell;


    //generate all tslat[][] matrices that are upper diagonal where the product of the
    //diagonal equals the current volume vol and where the upper diagonal elements take
    //all values less than the diagonal below it

    //initialize the superlattice vectors to zero
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)tslat[i][j]=0;

    for(tslat2d[0][0]=1; tslat2d[0][0]<=vol; tslat2d[0][0]++){
      if(vol%tslat2d[0][0] == 0){
	tslat2d[1][1]=(vol/tslat2d[0][0]);
	tslat2d[1][0]=0;
	for(tslat2d[0][1]=0; tslat2d[0][1]<tslat2d[1][1]; tslat2d[0][1]++){

	  if(excluded_axis=0){
	    tslat[0][0]=1; tslat[0][1]=0; tslat[0][2]=0;
	    tslat[1][0]=0; tslat[1][1]=tslat2d[0][0]; tslat[1][2]=tslat2d[0][1];
	    tslat[2][0]=0; tslat[2][1]=tslat2d[1][0]; tslat[2][2]=tslat2d[1][1];
	  }
	  if(excluded_axis=1){
	    tslat[0][0]=tslat2d[0][0]; tslat[0][1]=0; tslat[0][2]=tslat2d[0][1];
	    tslat[1][0]=0; tslat[1][1]=1; tslat[1][2]=0;
	    tslat[2][0]=tslat2d[1][0]; tslat[2][1]=0; tslat[2][2]=tslat2d[1][1];
	  }
	  if(excluded_axis=2){
	    tslat[0][0]=tslat2d[0][0]; tslat[0][1]=tslat2d[0][1]; tslat[0][2]=0;
	    tslat[1][0]=tslat2d[1][0]; tslat[1][1]=tslat2d[1][1]; tslat[1][2]=0;
	    tslat[2][0]=0; tslat[2][1]=0; tslat[2][2]=1;
	  }

	  //copy the superlattice vectors into lattice_point objects
	  //and get their cartesian coordinates

	  vec lat_point[3];

	  for(i=0; i<3; i++){
	    for(j=0; j<3; j++) lat_point[i].fcoord[j]=tslat[i][j];
	    lat_point[i].frac_on=true;
	    conv_AtoB(FtoC, lat_point[i].fcoord, lat_point[i].ccoord);
	    lat_point[i].cart_on=true;
	  }

	  //if no supercells have been added for this volume, add it
	  //else perform all point group operations to the supercell and
	  //see if it is a linear combination of already found supercells with the same volume


	  if(tsupercell.size() == 0){
	    structure tsup_lat;
	    tsup_lat.scale = scale;
	    for(i=0; i<3; i++){
	      for(j=0; j<3; j++){
		tsup_lat.lat[i][j]= lat_point[i].ccoord[j];
		tsup_lat.slat[i][j]= lat_point[i].fcoord[j];
	      }
	    }
	    tsup_lat.get_trans_mat();
	    tsupercell.push_back(tsup_lat);
	  }
	  else{
	    //apply each point group to this superlattice stored in lat_point

	    int num_miss1=0;
	    for(pg=0; pg<point_group.size(); pg++){
	      vec tlat_point[3];

	      for(i=0; i<3; i++) tlat_point[i]=lat_point[i].apply_sym(point_group[pg]);

	      //see if tlat_point[] can be expressed as a linear combination of any
	      //superlattice with volume vol already found

	      int num_miss2=0;
	      for(int ts=0; ts<tsupercell.size(); ts++){
		double lin_com[3][3];    //contains the coefficients relating the two lattices
		for(i=0; i<3; i++){
		  for(j=0; j<3; j++){
		    lin_com[i][j]=0.0;
		    for(k=0; k<3; k++)
		      lin_com[i][j]=lin_com[i][j]+tsupercell[ts].CtoF[i][k]*tlat_point[j].ccoord[k];
		  }
		}

		//check whether lin_com[][] are strictly integer
		//if so, the transformed superlattice is a linear combination of a previous one

		if(!is_integer(lin_com))num_miss2++;
	      }
	      if(num_miss2 == tsupercell.size())num_miss1++;
	    }
	    if(num_miss1 == point_group.size()){
	      structure tsup_lat;
	      tsup_lat.scale = scale;
	      for(i=0; i<3; i++){
		for(j=0; j<3; j++){
		  tsup_lat.lat[i][j]= lat_point[i].ccoord[j];
		  tsup_lat.slat[i][j]= lat_point[i].fcoord[j];
		}
	      }
	      tsup_lat.get_trans_mat();
	      tsupercell.push_back(tsup_lat);

	    }
	  }
	}
      }
    }
    for(i=0; i<tsupercell.size(); i++)supercell.push_back(tsupercell[i]);
  }
  return;
}
//************************************************************

void structure::generate_lat(structure prim){   /// added by jishnu

  //lat[][]=slat[][]*prim.lat[][]    in matrix form
  //this routine determines lat[][] given slat and prim.lat
  //then it rounds all elements of lat[][] to the nearest integer

  //

  
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      lat[i][j]=0.0;
      for(int k=0; k<3; k++){
	lat[i][j]=lat[i][j]+slat[i][k]*prim.lat[k][j]*prim.scale;
      }
    }
  }
  scale=1.0;
  
}    // end of s/r



//************************************************************

void structure::generate_slat(structure prim){

  //lat[][]=slat[][]*prim.lat[][]    in matrix form
  //this routine determines slat[][]
  //then it rounds all elements of slat[][] to the nearest integer



  double inv_lat[3][3];
  inverse(prim.lat,inv_lat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      slat[i][j]=0.0;
      for(int k=0; k<3; k++){
	slat[i][j]=slat[i][j]+lat[i][k]*inv_lat[k][j];
      }
    }
  }

  //round the elements of slat[][] to the closest integer

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(slat[i][j]-floor(slat[i][j]) < 0.5) slat[i][j]=floor(slat[i][j]);
      else slat[i][j]=ceil(slat[i][j]);
    }
  }

}


//************************************************************

void structure::generate_slat(structure prim, double rescale){

  //lat[][]=slat[][]*prim.lat[][]    in matrix form
  //this routine determines slat[][]
  //then it rounds all elements of slat[][] to the nearest integer


  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      lat[i][j]=rescale*lat[i][j];
    }

  double inv_lat[3][3];
  inverse(prim.lat,inv_lat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      slat[i][j]=0.0;
      for(int k=0; k<3; k++){
	slat[i][j]=slat[i][j]+lat[i][k]*inv_lat[k][j];
      }
    }
  }

  //round the elements of slat[][] to the closest integer

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(slat[i][j]-floor(slat[i][j]) < 0.5) slat[i][j]=floor(slat[i][j]);
      else slat[i][j]=ceil(slat[i][j]);
    }
  }


}


//************************************************************

void structure::generate_ideal_slat(structure prim){

  //ilat[][]=slat[][]*prim.lat[][]    in matrix form
  //this routine determines slat[][]
  //then it rounds all elements of slat[][] to the nearest integer

  double inv_lat[3][3];
  inverse(prim.lat,inv_lat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      slat[i][j]=0.0;
      for(int k=0; k<3; k++){
	slat[i][j]=slat[i][j]+ilat[i][k]*inv_lat[k][j];
      }
    }
  }

  //round the elements of slat[][] to the closest integer

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(slat[i][j]-floor(slat[i][j]) < 0.5) slat[i][j]=floor(slat[i][j]);
      else slat[i][j]=ceil(slat[i][j]);
    }
  }


}


//************************************************************

void structure::generate_ideal_slat(structure prim, double rescale){

  //ilat[][]=slat[][]*prim.lat[][]    in matrix form
  //this routine determines slat[][]
  //then it rounds all elements of slat[][] to the nearest integer

  double tilat[3][3];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      tilat[i][j]=rescale*ilat[i][j];
    }

  double inv_lat[3][3];
  inverse(prim.lat,inv_lat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      slat[i][j]=0.0;
      for(int k=0; k<3; k++){
	slat[i][j]=slat[i][j]+tilat[i][k]*inv_lat[k][j];
      }
    }
  }

  //round the elements of slat[][] to the closest integer

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(slat[i][j]-floor(slat[i][j]) < 0.5) slat[i][j]=floor(slat[i][j]);
      else slat[i][j]=ceil(slat[i][j]);
    }
  }


}

//************************************************************

void structure::calc_strain(){

  //get the matrix the relates lat[][]^transpose=deform[][]*ilat[][]^transpose
  //get the symmetric part of deform as 1/2(deform[][]+deform[][]^transpose)

  double tilat[3][3],tlat[3][3],inv_tilat[3][3],deform[3][3];

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      tlat[i][j]=lat[j][i];
      tilat[i][j]=ilat[j][i];
    }
  }

  inverse(tilat,inv_tilat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      deform[i][j]=0.0;
      for(int k=0; k<3; k++)
	deform[i][j]=deform[i][j]+tlat[i][k]*inv_tilat[k][j];
    }
    //    deform[i][i]=deform[i][i]-1.0;
  }

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) strain[i][j]=0.5*(deform[i][j]+deform[j][i]);


}

//************************************************************

void structure::generate_prim_grid(){
  int i,j,k;

  prim_grid.clear();

  //create a mesh of primitive lattice points that encompasses the supercell

  //first determine the extrema of all corners of the supercell projected onto the
  //different axes of the primitive cell

  int min[3],max[3];
  int corner[8][3];

  //generate the corners of the supercell

  for(j=0; j<3; j++) corner[0][j]=0;
  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      if(slat[i][j] <= 0.0) corner[i+1][j]=int(floor(slat[i][j]));
      if(slat[i][j] > 0.0) corner[i+1][j]=int(ceil(slat[i][j]));
    }

  //add up pairs of lattice vectors
  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      if(slat[(i+1)%3][j]+slat[(i+2)%3][j] <= 0.0) corner[i+4][j]=int(floor(slat[(i+1)%3][j]+slat[(i+2)%3][j]));
      if(slat[(i+1)%3][j]+slat[(i+2)%3][j] > 0.0) corner[i+4][j]=int(ceil(slat[(i+1)%3][j]+slat[(i+2)%3][j]));
    }

  //add up all three lattice vectors
  for(j=0; j<3; j++){
    if(slat[0][j]+slat[1][j]+slat[2][j] <= 0.0) corner[7][j]=int(floor(slat[0][j]+slat[1][j]+slat[2][j]));
    if(slat[0][j]+slat[1][j]+slat[2][j] > 0.0) corner[7][j]=int(ceil(slat[0][j]+slat[1][j]+slat[2][j]));
  }


  //get the extrema of the coordinates projected on the primitive

  for(j=0; j<3; j++){
    min[j]=corner[0][j];
    max[j]=corner[0][j];
  }

  for(i=1; i<8; i++){
    for(j=0; j<3; j++){
      if(min[j] > corner[i][j]) min[j]=corner[i][j];
      if(max[j] < corner[i][j]) max[j]=corner[i][j];
    }
  }


  //generate a grid of primitive lattice sites that encompasses the supercell
  //keep only those primitive lattice sites that reside within the supercell

  for(i=min[0]; i <= max[0]; i++){
    for(j=min[1]; j <= max[1]; j++){
      for(k=min[2]; k <= max[2]; k++){
	double ptemp[3],stemp[3],ctemp[3];

	ptemp[0]=i;
	ptemp[1]=j;
	ptemp[2]=k;

	conv_AtoB(PtoS,ptemp,stemp);

	int ll=0;
	for(int ii=0; ii<3; ii++)
	  if(stemp[ii] >=0.0-tol && stemp[ii] < 1.0-tol) ll++;

	if(ll == 3){
	  vec tlat_point;
	  conv_AtoB(FtoC,stemp,ctemp);
	  for(int jj=0; jj<3; jj++){
	    tlat_point.fcoord[jj]=stemp[jj];
	    tlat_point.ccoord[jj]=ctemp[jj];
	  }
	  tlat_point.frac_on=true;
	  tlat_point.cart_on=true;

	  prim_grid.push_back(tlat_point);
	}
      }

    }
  }

}


//************************************************************
//algorithm taken from B. Z. Yanchitsky, A. N. Timoshevskii,
//Computer Physics Communications vol 139 (2001) 235-242


void structure::generate_3d_reduced_cell(){
  int i,j,k;
  double rlat[3][3];
  double rslat[3][3];
  double leng[3],angle[3];

  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      rlat[i][j]=lat[i][j];
      rslat[i][j]=slat[i][j];
    }

  double vol=determinant(lat);
  latticeparam(lat,leng,angle);

  //get the largest lattice parameter
  //for all diagonals, find that smallest one that reduces the cell
  //use that diagonal to reduce the cell
  //if among all the diagonals none are smaller than the existing lattice vectors, the cell is reduced
  //then put it into Niggli unique form

  bool small_diag=true;
  while(small_diag){

    int replace;
    double min_leng=leng[0];
    double min_diag[3],min_sdiag[3];

    for(i=0; i<3; i++) if(min_leng < leng[i])min_leng=leng[i];

    small_diag=false;
    for(int i0=-1; i0<=1; i0++){
      for(int i1=-1; i1<=1; i1++){
	for(int i2=-1; i2<=1; i2++){
	  if(!(i0==0 && i1==0) && !(i0==0 && i2==0) && !(i1==0 && i2==0)){

	    double diag[3];
	    double sdiag[3];
	    double tleng=0.0;

	    for(j=0; j<3; j++){
	      diag[j]=i0*rlat[0][j]+i1*rlat[1][j]+i2*rlat[2][j];
	      tleng=tleng+diag[j]*diag[j];
	      sdiag[j]=i0*rslat[0][j]+i1*rslat[1][j]+i2*rslat[2][j];
	    }
	    tleng=sqrt(tleng);
	    if(tleng < min_leng){

	      for(i=0; i<3; i++){

		if(tleng < leng[i]-tol){
		  double tlat[3][3],tslat[3][3];
		  for(j=0; j<3; j++)
		    for(k=0; k<3; k++){
		      tlat[j][k]=rlat[j][k];
		      tslat[j][k]=rslat[j][k];
		    }
		  for(k=0; k<3; k++){
		    tlat[i][k]=diag[k];
		    tslat[i][k]=sdiag[k];
		  }

		  if(abs(determinant(tlat)-vol) < tol){
		    min_leng=tleng;
		    replace=i;
		    for(k=0; k<3; k++){
		      min_diag[k]=diag[k];
		      min_sdiag[k]=sdiag[k];
		    }
		    small_diag=true;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(small_diag){
      for(k=0; k<3; k++){
	rlat[replace][k]=min_diag[k];
	rslat[replace][k]=min_sdiag[k];
      }
      latticeparam(rlat,leng,angle);
    }
  }

  //rearrange cell so a < b < c

  int rpermut[3];
  latticeparam(rlat,leng,angle,rpermut);
  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      lat[i][j]=rlat[rpermut[i]][j];
      slat[i][j]=rslat[rpermut[i]][j];
    }

  //make sure that the angles are either all obtuse or all acute.
  bool found=false;
  for(int i0=1; i0>=-1; i0--){
    for(int i1=1; i1>=-1; i1--){
      for(int i2=1; i2>=-1; i2--){
	if(!found && (i0!=0) && (i1!=0) && (i2!=0)){
	  for(j=0; j<3; j++){
	    rlat[0][j]=i0*lat[0][j];
	    rslat[0][j]=i0*slat[0][j];
	    rlat[1][j]=i1*lat[1][j];
	    rslat[1][j]=i1*slat[1][j];
	    rlat[2][j]=i2*lat[2][j];
	    rslat[2][j]=i2*slat[2][j];
	  }
	  if(determinant(rlat) > 0){
	    latticeparam(rlat,leng,angle);
	    if(abs(angle[0]-90) < tol && abs(angle[1]-90) < tol && abs(angle[2]-90) < tol) found=true;
	    if(angle[0]-90 > -tol && angle[1]-90 > -tol && angle[2]-90 > -tol) found=true;
	    if(angle[0]-90 < tol && angle[1]-90 < tol && angle[2]-90 < tol) found=true;
	  }
	}
      }
    }
  }

  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      lat[i][j]=rlat[i][j];
      slat[i][j]=rslat[i][j];
    }

  //put the reduced cell into Niggli form

}


//************************************************************
//algorithm taken from B. Z. Yanchitsky, A. N. Timoshevskii,
//Computer Physics Communications vol 139 (2001) 235-242


void structure::generate_2d_reduced_cell(int excluded_axis){
  int i,j,k;
  double rlat[3][3];
  double rslat[3][3];
  double leng[3],angle[3];

  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      rlat[i][j]=lat[i][j];
      rslat[i][j]=slat[i][j];
    }

  double vol=determinant(lat);
  latticeparam(lat,leng,angle);

  //get the largest lattice parameter
  //for all diagonals, find that smallest one that reduces the cell
  //use that diagonal to reduce the cell
  //if among all the diagonals none are smaller than the existing lattice vectors, the cell is reduced
  //then put it into Niggli unique form

  bool small_diag=true;
  while(small_diag){

    int replace;
    double min_leng=0;
    double min_diag[3],min_sdiag[3];

    for(i=0; i<3; i++)
      if(i != excluded_axis)
	if(min_leng < leng[i])min_leng=leng[i];

    small_diag=false;
    for(int ii0=-1; ii0<=1; ii0++){
      for(int ii1=-1; ii1<=1; ii1++){
	for(int ii2=-1; ii2<=1; ii2++){
	  int i0=ii0;
	  int i1=ii1;
	  int i2=ii2;
	  if(excluded_axis == 0) i0=0;
	  if(excluded_axis == 1) i1=0;
	  if(excluded_axis == 2) i2=0;

	  if(!(i0==0 && i1==0) && !(i0==0 && i2==0) && !(i1==0 && i2==0)){

	    double diag[3];
	    double sdiag[3];
	    double tleng=0.0;

	    for(j=0; j<3; j++){
	      diag[j]=i0*rlat[0][j]+i1*rlat[1][j]+i2*rlat[2][j];
	      tleng=tleng+diag[j]*diag[j];
	      sdiag[j]=i0*rslat[0][j]+i1*rslat[1][j]+i2*rslat[2][j];
	    }
	    tleng=sqrt(tleng);
	    if(tleng < min_leng){

	      for(i=0; i<3; i++){
		if(i != excluded_axis){
		  if(tleng < leng[i]-tol){
		    double tlat[3][3],tslat[3][3];
		    for(j=0; j<3; j++)
		      for(k=0; k<3; k++){
			tlat[j][k]=rlat[j][k];
			tslat[j][k]=rslat[j][k];
		      }
		    for(k=0; k<3; k++){
	              tlat[i][k]=diag[k];
		      tslat[i][k]=sdiag[k];
		    }

		    if(abs(determinant(tlat)-vol) < tol){
		      min_leng=tleng;
		      replace=i;
		      for(k=0; k<3; k++){
			min_diag[k]=diag[k];
			min_sdiag[k]=sdiag[k];
		      }
		      small_diag=true;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if(small_diag){
      for(k=0; k<3; k++){
	rlat[replace][k]=min_diag[k];
	rslat[replace][k]=min_sdiag[k];
      }
      latticeparam(rlat,leng,angle);
    }
  }

  for(i=0; i<3; i++)
    for(j=0; j<3; j++){
      lat[i][j]=rlat[i][j];
      slat[i][j]=rslat[i][j];
    }

}



//************************************************************

void specie::print(ostream &stream){
  //for(int i=0; i<2; i++) stream << name[i];   // commented by jishnu
  stream << name;   // jishnu
  return;
}



//************************************************************

atompos::atompos(){

  bit=0;
  specie tcompon;
  compon.push_back(tcompon);

  for(int i=0; i<3; i++){
    fcoord[i]=0.0;
    ccoord[i]=0.0;
    dfcoord[i]=0.0;
    dccoord[i]=0.0;
  }
}


//************************************************************

atompos atompos::apply_sym(sym_op op){
  int i,j,k;
  atompos tatom;

  if(op.frac_on == false && op.cart_on == false){
    cout << "no coordinates available to perform a symmetry operation on\n";
    return tatom;
  }

  if(op.frac_on == false || op.cart_on == false)op.update();

  tatom.bit=bit;
  tatom.occ=occ;

  // added by anton
  for(i=0; i< compon.size(); i++) tatom.compon.push_back(compon[i]);
  for(i=0; i< p_vec.size(); i++) tatom.p_vec.push_back(p_vec[i]);
  for(i=0; i< spin_vec.size(); i++) tatom.spin_vec.push_back(spin_vec[i]);  
  for(i=0; i< basis_vec.size(); i++) tatom.basis_vec.push_back(basis_vec[i]);
  tatom.basis_flag=basis_flag;

  tatom.compon.clear();
  for(i=0; i<compon.size(); i++) tatom.compon.push_back(compon[i]);

  for(i=0; i<3; i++){
    tatom.fcoord[i]=op.ftau[i];
    tatom.ccoord[i]=op.ctau[i];
    for(j=0; j<3; j++){
      tatom.fcoord[i]=tatom.fcoord[i]+op.fsym_mat[i][j]*fcoord[j];
      tatom.ccoord[i]=tatom.ccoord[i]+op.csym_mat[i][j]*ccoord[j];
    }
  }
  return tatom;
}


//************************************************************

void atompos::get_cart(double FtoC[3][3]){
  for(int i=0; i<3; i++){
    ccoord[i]=0.0;
    for(int j=0; j<3; j++){
      ccoord[i]=ccoord[i]+FtoC[i][j]*fcoord[j];
    }
  }
}



//************************************************************

void atompos::get_frac(double CtoF[3][3]){
  for(int i=0; i<3; i++){
    fcoord[i]=0.0;
    for(int j=0; j<3; j++){
      fcoord[i]=fcoord[i]+CtoF[i][j]*ccoord[j];
    }
  }
}



//************************************************************

void atompos::readf(istream &stream){

  for(int i=0; i<3; i++) stream >> fcoord[i];

  bit=0;   // is set to zero, unless we read a different number after the specie list

  char buff[200];
  char sp[200];
  stream.getline(buff,199);

  int on=0;
  int off=1;
  int count=0;
  int ii;
  compon.clear();   //clear the blanck specie in the component vector
  for(ii=0; buff[ii]!=0; ii++){
    if(buff[ii] != ' ' && buff[ii] != '\t'){
      if(off == 1){
	on=1;
	off=0;
	count=0;
	for(int jj=0; jj< sizeof(sp); jj++) sp[jj]=' ';
	sp[count]=buff[ii];
      }
      else{
	count++;
	sp[count]=buff[ii];

      }
    }
    else{
      if(on == 1){
	on=0;
	off=1;
	if(isdigit(sp[0])){
	  bit=atoi(sp);
	}
	else{
	  specie tcompon;
	  //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
	  int jj=0;
	  do{
	    tcompon.name.push_back(sp[jj]);
	    jj++;
	  }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
	  compon.push_back(tcompon);
	}
      }
    }
  }
  if(buff[ii] == 0 && on == 1){
    on=0;
    off=1;
    if(isdigit(sp[0])){
      bit=atoi(sp);
    }
    else{
      specie tcompon;
      //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
      int jj=0;
      do{
	tcompon.name.push_back(sp[jj]);
	jj++;
      }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
      compon.push_back(tcompon);
    }
  }

}


//************************************************************

void atompos::readc(istream &stream){
  for(int i=0; i<3; i++) stream >> ccoord[i];

  bit=0;   // is set to zero, unless we read a different number after the specie list

  char buff[200];
  char sp[200];
  stream.getline(buff,199);

  int on=0;
  int off=1;
  int count=0;
  int ii;
  compon.clear();   //clear the blanck specie in the component vector
  for(ii=0; buff[ii]!=0; ii++){
    if(buff[ii] != ' ' && buff[ii] != '\t'){
      if(off == 1){
	on=1;
	off=0;
	count=0;
	for(int jj=0; jj< sizeof(sp); jj++) sp[jj]=' ';
	sp[count]=buff[ii];
      }
      else{
	count++;
	sp[count]=buff[ii];

      }
    }
    else{
      if(on == 1){
	on=0;
	off=1;
	if(isdigit(sp[0])){
	  bit=atoi(sp);
	}
	else{
	  specie tcompon;
	  //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
	  int jj=0;
	  do{
	    tcompon.name.push_back(sp[jj]);
	    jj++;
	  }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
	  compon.push_back(tcompon);
	}
      }
    }
  }
  if(buff[ii] == 0 && on == 1){
    on=0;
    off=1;
    if(isdigit(sp[0])){
      bit=atoi(sp);
    }
    else{
      specie tcompon;
      //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
      int jj=0;
      do{
	tcompon.name.push_back(sp[jj]);
	jj++;
      }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
      compon.push_back(tcompon);
    }
  }

}



//************************************************************

void atompos::print(ostream &stream){
  for(int i=0; i<3; i++){
    stream.precision(5);stream.width(10);stream.setf(ios::showpoint);
    stream << fcoord[i] << "  ";
  }
  for(int i=0; i<compon.size(); i++){
    //for(int j=0; j<2; j++) stream << compon[i].name[j]; // commented by jishnu
    stream << compon[i].name;   // jishnu
    stream << "  ";
  }
  stream<< bit << " ";
  stream << "\n";
}


//************************************************************

void atompos::assign_spin(){

  int remain=compon.size()%2;
  int bound=(compon.size()-remain)/2;

  int j=bound;
  for(int i=0; i<compon.size(); i++){
    if(j==0 && remain == 0) j--;
    compon[i].spin=j;
    j--;
  }

  return;
}


//************************************************************

void atompos::assemble_flip(){
  for(int i=0; i<compon.size(); i++){
    vector<int> tflip;
    for(int j=1; j<compon.size(); j++){
      tflip.push_back(compon[(i+j)%compon.size()].spin);
    }
    flip.push_back(tflip);
  }
}


//************************************************************

int atompos::iflip(int spin){
  int i;
  for(i=0; i<compon.size(); i++){
    if(compon[i].spin == spin) return i;
  }
  cout << "spin = " << spin << "  is not among those for this sublattice \n";
}

//************************************************************

int atompos::get_spin(string name){
  for(int i=0; i<compon.size(); i++){
    if(compon[i].name.compare(name) == 0) return compon[i].spin;
  }
  //cout << name[0] << name[1] << " is not among the list of components\n";   // commented by jishnu
  cout << name << " is not among the list of components\n";

}




//************************************************************

void structure::bring_in_cell(){
  for(int na=0; na<atom.size(); na++)within(atom[na].fcoord);
  calc_cartesian();
  return;
}



//************************************************************

void structure::calc_factor_group(){
  int pg,i,j,k,n,m,num_suc_maps;
  atompos hatom;
  double shift[3],temp[3];
  sym_op tfactor_group;



  //all symmetry operations are done within the fractional coordinate system
  //since translations back into the unit cell are straightforward

  //for each point group operation, apply it to the crystal and store the transformed
  //coordinates in tatom[]


  for(pg=0; pg<point_group.size(); pg++){
    vector<atompos> tatom;
    for(i=0; i<atom.size(); i++){
      hatom=atom[i].apply_sym(point_group[pg]);
      tatom.push_back(hatom);
    }


    //consider all internal shifts that move an atom of the original structure (e.g. the first
    //atom) to an atom of the transformed structure. Then see if that translation maps the
    //transformed crystal onto the original crystal.


    for(i=0; i<atom.size(); i++){
      if(compare(atom[0].compon, atom[i].compon)){
	for(j=0; j<3; j++) shift[j]=atom[0].fcoord[j]-tatom[i].fcoord[j];

	num_suc_maps=0;
	for(n=0; n<atom.size(); n++){
	  for(m=0; m<atom.size(); m++){
	    if(compare(atom[n].compon, atom[m].compon)){
	      for(j=0; j<3; j++) temp[j]=atom[n].fcoord[j]-tatom[m].fcoord[j]-shift[j];
	      within(temp);

	      k=0;
	      for(j=0; j<3; j++)
		if(abs(temp[j]) < 0.00005 ) k++;
	      if(k==3)num_suc_maps++;
	    }
	  }
	}

	if(num_suc_maps == atom.size()){
	  within(shift);

	  //check whether the symmetry operation already exists in the factorgroup array

	  int ll=0;
	  for(int ii=0; ii<factor_group.size(); ii++)
	    if(compare(point_group[pg].fsym_mat,factor_group[ii].fsym_mat)
	       && compare(shift,factor_group[ii].ftau) )break;
	    else ll++;

	  // if the symmetry operation is new, add it to the factorgroup array
	  // and update all info about the sym_op object

	  if(factor_group.size() == 0 || ll == factor_group.size()){
	    for(int jj=0; jj<3; jj++){
	      tfactor_group.frac_on=false;
	      tfactor_group.cart_on=false;
	      tfactor_group.ftau[jj]=shift[jj];
	      for(int kk=0; kk<3; kk++){
		tfactor_group.fsym_mat[jj][kk]=point_group[pg].fsym_mat[jj][kk];
		tfactor_group.lat[jj][kk]=lat[jj][kk];
	      }
	    }
	    tfactor_group.frac_on=true;
	    tfactor_group.update();
	    tfactor_group.get_sym_type(); // Added by John
	    factor_group.push_back(tfactor_group);
	  }
	}
      }
    }
    tatom.clear();
  }

  return;
}




//************************************************************

void structure::expand_prim_basis(structure prim){
  int i,j,k;


  //add all the prim_grid lattice points to the basis within prim
  //transform the coordinates using PtoS to get the fractional coordinates within
  //current lattice frame

  get_trans_mat();
  generate_prim_grid();

  num_each_specie.clear();
  int ii=0;
  for(i=0; i<prim.num_each_specie.size(); i++){
    int tnum_each_specie=0;


    for(j=0; j<prim.num_each_specie[i]; j++){

      for(k=0; k<prim_grid.size(); k++){
	atompos hatom;
	hatom=prim.atom[ii];
	double temp[3];
	for(int jj=0; jj<3; jj++){
	  temp[jj]=0;
	  for(int kk=0; kk<3; kk++) temp[jj]=temp[jj]+PtoS[jj][kk]*hatom.fcoord[kk];
	}
	for(int kk=0; kk<3; kk++) hatom.fcoord[kk]=temp[kk]+prim_grid[k].fcoord[kk];
	within(hatom);
	hatom.assign_spin();
	atom.push_back(hatom);
	tnum_each_specie++;
      }
      ii++;
    }
    num_each_specie.push_back(tnum_each_specie);
  }
  frac_on=1;
  calc_cartesian();
  //  cout << "calc_cartesian inside expand_prim_cell \n";
  return;
}


//************************************************************

void structure::map_on_expanded_prim_basis(structure prim){

  //takes the current structure and maps the atom positions on the ideal positions

  //prim is expanded according to slat[][]: the expanded structure is ideal_struc
  //this means that the ideal coordinates (within ideal_struc) are those as if simply a volume
  //      relaxation had been done, freezing the internal coordinates

  //the atoms of the current structure are then mapped on to those of ideal_struc
  //the deviation from the ideal positions are stored in dfcoord[] and dccoord[] of atom


  structure ideal_struc;

  //copy title, scale, lat and slat into ideal_struc

  for(int i=0; i<200; i++)ideal_struc.title[i]=title[i];
  ideal_struc.scale=scale;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ideal_struc.lat[i][j]=lat[i][j];
      ideal_struc.slat[i][j]=slat[i][j];
    }
  }

  //put ideal lattice parameters in ilat
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ilat[i][j]=0.0;
      for(int k=0; k<3; k++){
	ilat[i][j]=ilat[i][j]+slat[i][k]*prim.lat[k][j];
      }
    }
  }


  // expand the primitive cell within the slat of ideal_struc
  // expand_prim_basis will convert to cartesian in the lat[][] coordinate system (which is relaxed)

  ideal_struc.expand_prim_basis(prim);

  //if you want distances (cartesian coordinates) in the unrelaxed coordinate system
  //  coord_trans_mat(ilat,ideal_struc.FtoC,ideal_struc.CtoF);
  //  ideal_struc.cart_on=0;
  //  ideal_struc.calc_cartesian();


  // make sure that the current structure has cartesian coordinates (i.e. the one with relaxed coordinates)

  get_trans_mat();
  //  coord_trans_mat(ilat,FtoC,CtoF);      //to have cartesian coordinates and therefore distances in the unrelaxed coordinate system
  calc_cartesian();

  //find the positions in ideal_struc that are closest to those in the current structure
  //   if the position is found, increase atom[].bit from 0 to 1
  //make sure that the first component of structure coincides with one of the components in the component list of ideal_struc
  //if so update occ of ideal_struc and put delta in atom of ideal_struc
  //if some ideal_struc positions are not matched - check if there are vacancies in the component list if so, make the occ = Va

  //transcribe the atom objects of ideal_struc into those of the current structure

  for(int na=0; na<atom.size(); na++){
    if(atom[na].compon.size() == 0){
      cout << "the structure file has no atomic labels next to the coordinates\n";
      cout << "exiting \n";
      exit(1);
    }
    int min_ind=-1;
    int comp_ind=-1;
    int min_i,min_j,min_k;
    double min_dist=1.0e10;

    //look for the closest ideal position
    for(int ina=0; ina<ideal_struc.atom.size(); ina++){

      //first make sure that atom[na].compon[0] belongs to one of ideal_struc.atom[ina].compon
      for(int c=0; c<ideal_struc.atom[ina].compon.size(); c++){

	//if(compare(atom[na].compon[0],ideal_struc.atom[ina].compon[c])){

	if(compare(atom[na].occ,ideal_struc.atom[ina].compon[c])){

	  //need to translate atom to all edges of the unit cell and get the minimal distance
	  //to ideal_struc.atom[ina]
	  double tmin_dist=1.0e10;
	  int tmin_i,tmin_j,tmin_k;
	  for(int i=-1; i<=1; i++){
	    for(int j=-1; j<=1; j++){
	      for(int k=-1; k<=1; k++){
		atompos hatom=atom[na];
		hatom.fcoord[0]=atom[na].fcoord[0]+i;
		hatom.fcoord[1]=atom[na].fcoord[1]+j;
		hatom.fcoord[2]=atom[na].fcoord[2]+k;
		conv_AtoB(FtoC,hatom.fcoord,hatom.ccoord);
		double dist=0.0;
		for(int ii=0; ii<3; ii++){
		  dist=dist+(hatom.ccoord[ii]-ideal_struc.atom[ina].ccoord[ii])*(hatom.ccoord[ii]-ideal_struc.atom[ina].ccoord[ii]);
		}
		dist=sqrt(dist);

		if(dist<tmin_dist){
		  tmin_dist=dist;
		  tmin_i=i;
		  tmin_j=j;
		  tmin_k=k;
		}
	      }
	    }
	  }
	  if(tmin_dist<min_dist){
	    min_dist=tmin_dist;
	    min_ind=ina;
	    comp_ind=c;
	    min_i=tmin_i;
	    min_j=tmin_j;
	    min_k=tmin_k;
	  }
	  break;
	}
      }

    }
    if(min_ind >= 0 && ideal_struc.atom[min_ind].bit==0){
      ideal_struc.atom[min_ind].assign_spin();
      ideal_struc.atom[min_ind].occ=ideal_struc.atom[min_ind].compon[comp_ind];
      atom[na].fcoord[0]=atom[na].fcoord[0]+min_i;
      atom[na].fcoord[1]=atom[na].fcoord[1]+min_j;
      atom[na].fcoord[2]=atom[na].fcoord[2]+min_k;
      conv_AtoB(FtoC,atom[na].fcoord,atom[na].ccoord);

      for(int i=0; i<3; i++){
	ideal_struc.atom[min_ind].dfcoord[i]=atom[na].fcoord[i]-ideal_struc.atom[min_ind].fcoord[i];
	ideal_struc.atom[min_ind].dccoord[i]=atom[na].ccoord[i]-ideal_struc.atom[min_ind].ccoord[i];
      }
      ideal_struc.atom[min_ind].delta=min_dist;
      ideal_struc.atom[min_ind].bit=1;    // this means this sites has already been taken
    }
    else {
      if(min_ind == -1){
	cout << "it was not possible to map atom\n";
	atom[na].print(cout);
	cout << "onto the ideal structure\n";
      }
      if(ideal_struc.atom[min_ind].bit == 1){
	cout << "it was not possible to map atom \n";
	atom[na].print(cout);
	cout << "onto the closest ideal position\n";
	ideal_struc.atom[min_ind].print(cout);
	cout << "since this ideal atom position has already been claimed\n";
      }
    }

  }

  // run through ideal_struc to find sites that have not been claimed yet
  // check if these sites can hold vacancies, if so, occ becomes the vacancy
  // if not, there is a problem

  for(int ina=0; ina<ideal_struc.atom.size(); ina++){
    if(ideal_struc.atom[ina].bit == 0){
      for(int c=0; c< ideal_struc.atom[ina].compon.size(); c++){
	if(ideal_struc.atom[ina].compon[c].name[0] == 'V' && ideal_struc.atom[ina].compon[c].name[1] == 'a'){
	  ideal_struc.atom[ina].occ=ideal_struc.atom[ina].compon[c];
	  ideal_struc.atom[ina].delta=0.0;
	  ideal_struc.atom[ina].bit=1;
	  break;
	}
      }
      if(ideal_struc.atom[ina].bit == 0){
	cout << "was not able to map any atom onto this position\n";
	ideal_struc.atom[ina].print(cout);
      }
    }
  }

  //copy the ideal_struc atom positions into the current structure

  num_each_specie.clear();
  atom.clear();


  //for now, we just list the total number of atoms in num_each_specie (later we may
  //modify this part to do the same as is done when structures are generated from scratch

  num_each_specie.push_back(ideal_struc.atom.size());
  for(int ina=0; ina<ideal_struc.atom.size(); ina++){
    atom.push_back(ideal_struc.atom[ina]);
    //assign to atom.bit the index of the component at that site and assign spins to each component

  }

  //modified by anton
  for(int i=0; i<atom.size(); i++) get_basis_vectors(atom[i]);

  return;

}


//************************************************************

void structure::map_on_expanded_prim_basis(structure prim, arrangement &conf){

  //takes the current structure and maps the atom positions on the ideal positions

  //prim is expanded according to slat[][]: the expanded structure is ideal_struc
  //this means that the ideal coordinates (within ideal_struc) are those as if simply a volume
  //      relaxation had been done, freezing the internal coordinates

  //the atoms of the current structure are then mapped on to those of ideal_struc
  //the deviation from the ideal positions are stored in dfcoord[] and dccoord[] of atom


  structure ideal_struc;

  //copy title, scale, lat and slat into ideal_struc

  for(int i=0; i<200; i++)ideal_struc.title[i]=title[i];
  ideal_struc.scale=scale;
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ideal_struc.lat[i][j]=lat[i][j];
      ideal_struc.slat[i][j]=slat[i][j];
    }
  }

  // expand the primitive cell within the slat of ideal_struc
  // expand_prim_basis will convert to cartesian in the lat[][] coordinate system (which is relaxed)

  ideal_struc.expand_prim_basis(prim);

  //initialize the bit of each atom to be -1 (means has not been mapped onto yet)

  for(int na=0; na<ideal_struc.atom.size(); na++) ideal_struc.atom[na].bit=-1;

  //if you want distances (cartesian coordinates) in the unrelaxed coordinate system
  //  coord_trans_mat(ilat,ideal_struc.FtoC,ideal_struc.CtoF);
  //  ideal_struc.cart_on=0;
  //  ideal_struc.calc_cartesian();


  // make sure that the current structure has cartesian coordinates (i.e. the one with relaxed coordinates)

  get_trans_mat();
  //  coord_trans_mat(ilat,FtoC,CtoF);      //to have cartesian coordinates and therefore distances in the unrelaxed coordinate system
  calc_cartesian();

  //find the positions in ideal_struc that are closest to those in the current structure
  //   if the position is found, change atom[].bit from -1 to the index of the component at that site
  //make sure that the first component of structure coincides with one of the components in the component list of ideal_struc
  //if so update occ of ideal_struc and put delta in atom of ideal_struc
  //if some ideal_struc positions are not matched - check if there are vacancies in the component list if so, make the occ = Va

  //transcribe the atom objects of ideal_struc into those of the current structure

  for(int na=0; na<atom.size(); na++){
    if(atom[na].compon.size() == 0){
      cout << "the structure file has no atomic labels next to the coordinates\n";
      cout << "exiting \n";
      exit(1);
    }
    int min_ind=-1;
    int comp_ind=-1;
    int min_i,min_j,min_k;
    double min_dist=1.0e10;

    //look for the closest ideal position
    for(int ina=0; ina<ideal_struc.atom.size(); ina++){

      //first make sure that atom[na].compon[0] belongs to one of ideal_struc.atom[ina].compon
      for(int c=0; c<ideal_struc.atom[ina].compon.size(); c++){


	if(compare(atom[na].compon[0],ideal_struc.atom[ina].compon[c])){

	  //need to translate atom to all edges of the unit cell and get the minimal distance
	  //to ideal_struc.atom[ina]
	  double tmin_dist=1.0e10;
	  int tmin_i,tmin_j,tmin_k;
	  for(int i=-1; i<=1; i++){
	    for(int j=-1; j<=1; j++){
	      for(int k=-1; k<=1; k++){
		atompos hatom=atom[na];
		hatom.fcoord[0]=atom[na].fcoord[0]+i;
		hatom.fcoord[1]=atom[na].fcoord[1]+j;
		hatom.fcoord[2]=atom[na].fcoord[2]+k;
		conv_AtoB(FtoC,hatom.fcoord,hatom.ccoord);
		double dist=0.0;
		for(int ii=0; ii<3; ii++){
		  dist=dist+(hatom.ccoord[ii]-ideal_struc.atom[ina].ccoord[ii])*(hatom.ccoord[ii]-ideal_struc.atom[ina].ccoord[ii]);
		}
		dist=sqrt(dist);

		if(dist<tmin_dist){
		  tmin_dist=dist;
		  tmin_i=i;
		  tmin_j=j;
		  tmin_k=k;
		}
	      }
	    }
	  }
	  if(tmin_dist<min_dist){
	    min_dist=tmin_dist;
	    min_ind=ina;
	    comp_ind=c;
	    min_i=tmin_i;
	    min_j=tmin_j;
	    min_k=tmin_k;
	  }
	  break;
	}
      }

    }
    if(min_ind >= 0 && ideal_struc.atom[min_ind].bit==-1){
      ideal_struc.atom[min_ind].assign_spin();
      ideal_struc.atom[min_ind].occ=ideal_struc.atom[min_ind].compon[comp_ind];
      ideal_struc.atom[min_ind].bit=comp_ind;   // also indicates that this sites has already been taken
      atom[na].fcoord[0]=atom[na].fcoord[0]+min_i;
      atom[na].fcoord[1]=atom[na].fcoord[1]+min_j;
      atom[na].fcoord[2]=atom[na].fcoord[2]+min_k;
      conv_AtoB(FtoC,atom[na].fcoord,atom[na].ccoord);

      for(int i=0; i<3; i++){
	ideal_struc.atom[min_ind].dfcoord[i]=atom[na].fcoord[i]-ideal_struc.atom[min_ind].fcoord[i];
	ideal_struc.atom[min_ind].dccoord[i]=atom[na].ccoord[i]-ideal_struc.atom[min_ind].ccoord[i];
      }
      ideal_struc.atom[min_ind].delta=min_dist;
    }
    else {
      if(min_ind == -1){
	cout << "it was not possible to map atom\n";
	atom[na].print(cout);
	cout << "onto the ideal structure\n";
      }
      if(ideal_struc.atom[min_ind].bit != -1){
	cout << "it was not possible to map atom \n";
	atom[na].print(cout);
	cout << "onto the closest ideal position\n";
	ideal_struc.atom[min_ind].print(cout);
	cout << "since this ideal atom position has already been claimed\n";
      }
    }

  }

  // run through ideal_struc to find sites that have not been claimed yet
  // check if these sites can hold vacancies, if so, occ becomes the vacancy
  // if not, there is a problem

  for(int ina=0; ina<ideal_struc.atom.size(); ina++){
    if(ideal_struc.atom[ina].bit == -1){
      for(int c=0; c< ideal_struc.atom[ina].compon.size(); c++){
	if(ideal_struc.atom[ina].compon[c].name[0] == 'V' && ideal_struc.atom[ina].compon[c].name[1] == 'a'){
	  ideal_struc.atom[ina].occ=ideal_struc.atom[ina].compon[c];
	  ideal_struc.atom[ina].delta=0.0;
	  ideal_struc.atom[ina].bit=c;
	  break;
	}
      }
      if(ideal_struc.atom[ina].bit == -1){
	cout << "was not able to map any atom onto this position\n";
	ideal_struc.atom[ina].print(cout);
      }
    }
  }

  //copy the ideal_struc atom positions into the current structure

  num_each_specie.clear();
  atom.clear();
  conf.bit.clear();

  //for now, we just list the total number of atoms in num_each_specie (later we may
  //modify this part to do the same as is done when structures are generated from scratch

  num_each_specie.push_back(ideal_struc.atom.size());
  for(int ina=0; ina<ideal_struc.atom.size(); ina++){
    atom.push_back(ideal_struc.atom[ina]);
    conf.bit.push_back(ideal_struc.atom[ina].bit);
  }

  //modified by anton
  for(int i=0; i<atom.size(); i++) get_basis_vectors(atom[i]);

  return;

}



//************************************************************

void structure::idealize(){
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) ilat[i][j]=lat[i][j];

  get_latparam();
  get_trans_mat();
  calc_fractional();
  calc_cartesian();
}


//************************************************************

void structure::expand_prim_clust(multiplet basiplet, multiplet &super_basiplet){
  int nm,no,nc,np,ng;
  int i,j,k;

  get_trans_mat();

  generate_prim_grid();

  //transform the cluster coordinates (in prim system) to the current lattice coordinate system
  //add prim_grid points to each cluster and store the new clusters in super_basis

  super_basiplet.orb.push_back(basiplet.orb[0]);


  for(nm=1; nm<basiplet.orb.size(); nm++){
    vector<orbit> torbvec;
    for(no=0; no<basiplet.orb[nm].size(); no++){
      orbit torb1,torb2;
      for(nc=0; nc<basiplet.orb[nm][no].equiv.size(); nc++){
	cluster tclust;
	for(np=0; np<basiplet.orb[nm][no].equiv[nc].point.size(); np++){
	  atompos tpoint=basiplet.orb[nm][no].equiv[nc].point[np];
	  for(i=0; i<3; i++){
	    tpoint.fcoord[i]=0.0;
	    for(j=0; j<3; j++)
	      tpoint.fcoord[i]=tpoint.fcoord[i]+PtoS[i][j]*basiplet.orb[nm][no].equiv[nc].point[np].fcoord[j];
	  }
	  tclust.point.push_back(tpoint);
	}
	torb2.equiv.push_back(tclust);
      }
      for(ng=0; ng<prim_grid.size(); ng++){
	for(nc=0; nc<torb2.equiv.size(); nc++){
	  cluster tclust;
	  for(np=0; np<torb2.equiv[nc].point.size(); np++){
	    atompos tpoint=torb2.equiv[nc].point[np];
	    for(i=0; i<3; i++){
	      tpoint.fcoord[i]=tpoint.fcoord[i]+prim_grid[ng].fcoord[i];
	    }
	    within(tpoint);
	    tclust.point.push_back(tpoint);
	  }
	  torb1.equiv.push_back(tclust);
	}
      }
      torbvec.push_back(torb1);
    }
    super_basiplet.orb.push_back(torbvec);
  }
  return;
}



//************************************************************

void structure::write_factor_group(){
  int fg;

  ofstream out("factor_group");
  if(!out){
    cout << "Cannot open file.\n";
    return;
  }

  out << " number of factor group ops " << factor_group.size() << "\n";
  out << " *** Please Note:  Depending on translation vectors, rotations and mirrors may correspond to screw axes and glide planes.\n";
  for(fg=0; fg<factor_group.size(); fg++){
    out << "factor group operation " << fg << " \n";
    factor_group[fg].print_fsym_mat(out);
    out << "\n";
  }

  out.close();
}


//************************************************************
//read a POSCAR like file and collect all the structure variables
//modified to read PRIM file and determine which basis to use
//modified by Ben Swoboda

void structure::read_struc_prim(istream &stream){
  int i,j;
  char ch;
  char buff[200];
  char sp[200];
  atompos hatom;
	
  stream.getline(title,199);
  stream >> scale;
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      stream >> lat[i][j];
    }
  }
	
  //normalize the scale to 1.0 and adjust lattice vectors accordingly
	
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)lat[i][j]=scale*lat[i][j];
  scale=1.0;
	
	
  stream.getline(buff,199);
  stream.getline(buff,199);
	
  //Figure out how many species
	
  int on=0;
  int off=1;
  int count=0;
	
  for(i=0; buff[i]!=0; i++){
    if(buff[i] != ' ' && buff[i] != '\t'){
      if(off == 1){
	on=1;
	off=0;
	count=0;
	for(int ii=0; ii< sizeof(sp); ii++) sp[ii]=' ';
	sp[count]=buff[i];
      }
      else{
	count++;
	sp[count]=buff[i];
      }
    }
    else{
      if(on == 1){
	on=0;
	off=1;
				
	int ll=atoi(sp);
	num_each_specie.push_back(ll);
      }
    }
  }
  if(buff[i] == 0 && on == 1){
    on=0;
    off=1;
    int ll=atoi(sp);
    num_each_specie.push_back(ll);
  }
	
	
	
  int num_atoms=0;
  for(i=0; i<num_each_specie.size(); i++)num_atoms +=num_each_specie[i];
	
  // fractional coordinates or cartesian
	
  stream.getline(buff,199);
  i=0;
  while(buff[i] == ' ')i++;
	
  //first check for selective dynamics
	
  if(buff[i] == 'S' || buff[i] == 's') {
    stream.getline(buff,199);
    i=0;
    while(buff[i] == ' ')i++;
  }
	
  frac_on=false; cart_on=false;
	
  if(buff[i] == 'D' || buff[i] == 'd'){
    frac_on=true;
    cart_on=false;
  }
  else
    if(buff[i] == 'C' || buff[i] == 'c'){
      frac_on=false;
      cart_on=true;
    }
    else{
      cout << "ERROR in input\n";
      cout << "If not using selective dynamics the 7th line should be Direct or Cartesian.\n";
      cout << "Otherwise the 8th line should be Direct or Cartesian \n";
      exit(1);
    }
	
	
  //read the coordinates
	
  if(atom.size() != 0 ){
    cout << "the structure is going to be overwritten";
    atom.clear();
  }
	

	
  // The following part written by jishnu to take care of the spaces in prim file (spaces do not matter any more)	
  for(i=0; i<num_atoms; i++){   				
    for(j=0; j<3; j++){			
      if(frac_on == true) stream >> hatom.fcoord[j];
      else stream >> hatom.ccoord[j];
    }	
    // cout << "atom #" << i << "coordinates are : " <<  hatom.fcoord[0] << "   " <<  hatom.fcoord[1] << "   "<<  hatom.fcoord[2] << "\n";
    // determine the number and name of species that can occupy this atomposition
    hatom.compon.clear();   //clear the blank specie in the component vector
    do{
      if(!(((stream.peek()>='A')&&(stream.peek()<='Z'))||((stream.peek()>='a')&&(stream.peek()<='z')))) stream.get(ch);
      if(((stream.peek()>='A')&&(stream.peek()<='Z'))||((stream.peek()>='a')&&(stream.peek()<='z'))){
	specie tcompon;
	while(!(stream.peek()==' ')) {
	  stream.get(ch);
	  tcompon.name.push_back(ch);
	  if(stream.peek()=='\n') break;
	}
	hatom.compon.push_back(tcompon);
      }
      if((stream.peek()>='0')&&(stream.peek()<='9')){
	if((stream.peek()>='0')&&(stream.peek()<='1')) stream >> hatom.basis_flag;
	else cout << "Check the prim file; there is an invalid basis_flag \n";			
      }
    }while(!(stream.peek()=='\n'));       
    hatom.occ=hatom.compon[0];
    atom.push_back(hatom);		       
  }	
	
}


//************************************************************
//read a POSCAR like file and collect all the structure variables

void structure::read_struc_poscar(istream &stream){
  int i,j;
  char ch;
  char buff[200];
  char sp[200];
  atompos hatom;

  stream.getline(title,199);
  stream >> scale;
  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      stream >> lat[i][j];
    }
  }

  //normalize the scale to 1.0 and adjust lattice vectors accordingly

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)lat[i][j]=scale*lat[i][j];
  scale=1.0;


  stream.getline(buff,199);
  stream.getline(buff,199);

  //Figure out how many species

  int on=0;
  int off=1;
  int count=0;

  for(i=0; buff[i]!=0; i++){
    if(buff[i] != ' ' && buff[i] != '\t'){
      if(off == 1){
	on=1;
	off=0;
	count=0;
	for(int ii=0; ii< sizeof(sp); ii++) sp[ii]=' ';
	sp[count]=buff[i];
      }
      else{
	count++;
	sp[count]=buff[i];
      }
    }
    else{
      if(on == 1){
	on=0;
	off=1;

	int ll=atoi(sp);
	num_each_specie.push_back(ll);
      }
    }
  }
  if(buff[i] == 0 && on == 1){
    on=0;
    off=1;
    int ll=atoi(sp);
    num_each_specie.push_back(ll);
  }



  int num_atoms=0;
  for(i=0; i<num_each_specie.size(); i++)num_atoms +=num_each_specie[i];

  // fractional coordinates or cartesian

  stream.getline(buff,199);
  i=0;
  while(buff[i] == ' ')i++;

  //first check for selective dynamics

  if(buff[i] == 'S' || buff[i] == 's') {
    stream.getline(buff,199);
    i=0;
    while(buff[i] == ' ')i++;
  }

  frac_on=false; cart_on=false;

  if(buff[i] == 'D' || buff[i] == 'd'){
    frac_on=true;
    cart_on=false;
  }
  else
    if(buff[i] == 'C' || buff[i] == 'c'){
      frac_on=false;
      cart_on=true;
    }
    else{
      cout << "ERROR in input\n";
      cout << "If not using selective dynamics the 7th line should be Direct or Cartesian.\n";
      cout << "Otherwise the 8th line should be Direct or Cartesian \n";
      exit(1);
    }


  //read the coordinates

  if(atom.size() != 0 ){
    cout << "the structure is going to be overwritten";
    atom.clear();
  }

  for(i=0; i<num_atoms; i++){
    for(j=0; j<3; j++)
      if(frac_on == true)stream >> hatom.fcoord[j];
      else stream >> hatom.ccoord[j];

    // determine the number and name of species that can occupy this atomposition

    stream.getline(buff,199);

    int on=0;
    int off=1;
    int count=0;
    int ii;
    hatom.compon.clear();   //clear the blanck specie in the component vector
    for(ii=0; buff[ii]!=0; ii++){
      if(buff[ii] != ' ' && buff[ii] != '\t'){
	if(off == 1){
	  on=1;
	  off=0;
	  count=0;
	  for(int jj=0; jj< sizeof(sp); jj++) sp[jj]=' ';
	  sp[count]=buff[ii];
	}
	else{
	  count++;
	  sp[count]=buff[ii];
	}
      }
      else{
	if(on == 1){
	  on=0;
	  off=1;
	  specie tcompon;
	  //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
	  int jj=0;
	  do{
	    tcompon.name.push_back(sp[jj]);
	    jj++;
	  }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
	  hatom.compon.push_back(tcompon);
	}
      }
    }
    if(buff[ii] == 0 && on == 1){
      on=0;
      off=1;
      specie tcompon;
      //for(int jj=0; jj<2; jj++) tcompon.name[jj]=sp[jj];
      int jj=0;
      do{
	tcompon.name.push_back(sp[jj]);
	jj++;
      }while(!(sp[jj]==' ') && !(sp[jj]=='\n'));
      hatom.compon.push_back(tcompon);
    }
    //assign the first component to the occupation slot
    hatom.occ=hatom.compon[0];
    atom.push_back(hatom);
  }

}


//************************************************************
//************************************************************
//write the structure to a file in POSCAR like format

void structure::write_struc_poscar(ostream &stream) {
  int i,j;


  stream << title <<"\n";

  stream.precision(7);stream.width(12);stream.setf(ios::showpoint);
  stream << scale <<"\n";

  for(int i=0; i<3; i++){
    stream << "  ";
    for(int j=0; j<3; j++){

      stream.precision(9);stream.width(15);stream.setf(ios::showpoint);
      stream << lat[i][j] << "  ";

    }
    stream << "\n";
  }
  for(i=0; i < num_each_specie.size(); i++) stream << " " << num_each_specie[i] ;
  stream << "\n";

  stream << "Direct\n";

  for(i=0; i<atom.size(); i++){
    for(j=0; j<3; j++){
      stream.precision(9);stream.width(15);stream.setf(ios::showpoint);
      stream << atom[i].fcoord[j] << "  ";
    }
    //    for(int ii=0; ii<atom[i].compon.size(); ii++){
    stream << atom[i].occ.name;   // jishnu
    //      stream << "  ";
    //    }
    stream << "\n";
  }
}


//************************************************************
void structure::write_struc_xyz(ostream &stream){
  calc_cartesian();
  int tot_num_atoms=0;
  for(int i=0; i<atom.size(); i++){
    if(!(atom[i].occ.name[0] == 'V' && atom[i].occ.name[1] == 'a'))
      tot_num_atoms++;
  }
  
  stream << tot_num_atoms << "\n";
  stream << title << "\n";
  for(int i=0; i<atom.size(); i++){
    if(!(atom[i].occ.name[0] == 'V' && atom[i].occ.name[1] == 'a')){
      atom[i].occ.print(stream);
      stream << "  ";
      for(int j=0; j< 3; j++) stream << atom[i].ccoord[j] << "  ";
      stream << "\n";
    }
  }
  return;

}

//************************************************************
void structure::write_struc_xyz(ostream &stream, concentration out_conc){
  calc_cartesian();
  int tot_num_atoms=0;
  for(int i=0; i<atom.size(); i++){
    if(!(atom[i].occ.name[0] == 'V' && atom[i].occ.name[1] == 'a'))
      tot_num_atoms++;
  }
  
  stream << tot_num_atoms << "\n";
  stream << "Configuration with concentrations: ";
  out_conc.print_concentration(stream);
  stream << "\n";
  for(int i=0; i<atom.size(); i++){
    if(!(atom[i].occ.name[0] == 'V' && atom[i].occ.name[1] == 'a')){
      atom[i].occ.print(stream);
      stream << "  ";
      for(int j=0; j< 3; j++) stream << atom[i].ccoord[j] << "  ";
      stream << "\n";
    }
  }
  return;

}

//************************************************************

void structure::calc_cartesian(){
  int i,j,k;

  if(cart_on == true) return;
  if(cart_on == false) {
    if(frac_on == false){
      cout << "No structure initialized - cannot calc_cartesian\n";
      return;
    }
    for(i=0; i<atom.size(); i++) conv_AtoB(FtoC, atom[i].fcoord, atom[i].ccoord);

    cart_on=true;
  }
}



//************************************************************

void structure::calc_fractional(){
  int i,j,k;

  if(frac_on == true) return;
  if(frac_on == false) {
    if(cart_on == false){
      cout << "No structure initialized - cannot calc_fractional\n";
      return;
    }
    for(i=0; i<atom.size(); i++) conv_AtoB(CtoF, atom[i].ccoord, atom[i].fcoord);

    frac_on=true;
  }
}


//************************************************************

void structure::collect_components(){

  if(atom.size() == 0){
    cout << "cannot collect_components since no atoms in structure \n";
    return;
  }

  //find the first atom that has at least one component
  int i=0;
  while(atom[i].compon.size() < 1 && i<atom.size()){
    i++;
  }
  if(i==atom.size()){
    cout << "no atoms with at least one component\n";
    cout << "not collecting components\n";
    return;
  }


  compon.clear();
  compon.push_back(atom[i].compon[0]);
  for(int i=0; i<atom.size(); i++){
    if(atom[i].compon.size() >= 1){
      for(int j=0; j<atom[i].compon.size(); j++){
        int l=0;
        for(int k=0; k<compon.size(); k++)
          if(!compare(atom[i].compon[j],compon[k]))l++;
        if(l==compon.size())compon.push_back(atom[i].compon[j]);
      }
    }
  }

  return;
}


//************************************************************

void structure::collect_relax(string dir_name){

  //reads a POS and a CONTCAR, collects the info from both and places it
  //all in the current structure object
  //then prints some relevant info in a RELAX file

  //-puts the relaxed cell vectors in rlat (from CONTCAR), puts the original cell vectors in lat (from POS)
  //-puts the unrelaxed coordinates from POS in atompos together with the atom labels
  //-get the difference between the unrelaxed coordinates and the relaxed coordinates (from CONTCAR) and
  //place in dfcoord, dccoord and get the distance delta
  //-after printing out the info about the relaxations in RELAX, replace the atomic coordinates with
  //relaxed coordinates from CONTCAR, set dfcoord, dccoord and delta all to zero.


  //create a string for the POS filename and the CONTCAR filename
  //define structure objects for POS and for CONTCAR

  structure pos;

  string pos_filename=dir_name;
  pos_filename.append("/POS");
  ifstream in_pos;
  in_pos.open(pos_filename.c_str());
  if(!in_pos){
    cout << "cannot open file " << pos_filename << "\n";
    return;
  }

  pos.read_struc_poscar(in_pos);

  in_pos.close();

  structure contcar;

  string contcar_filename=dir_name;
  contcar_filename.append("/CONTCAR");
  ifstream in_contcar;
  in_contcar.open(contcar_filename.c_str());
  if(!in_contcar){
    cout << "cannot open file " << contcar_filename << "\n";
    return;
  }

  contcar.read_struc_poscar(in_contcar);

  in_contcar.close();

  if(pos.atom.size() != contcar.atom.size()){
    cout << "POS and CONTCAR in " << dir_name << "\n";
    cout << "are incompatible \n";
    cout << "quitting collect_relax() \n";
    return;
  }


  //collect the data from the two structures and place it in the current structure

  for(int i=0; i<200; i++) title[i]=pos.title[i];
  scale=pos.scale;

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      ilat[i][j]=pos.lat[i][j];
      lat[i][j]=contcar.lat[i][j];
    }
  }

  atom.clear();

  for(int na=0; na<pos.atom.size(); na++) atom.push_back(pos.atom[na]);
  frac_on=1;

  get_trans_mat();
  calc_cartesian();


  //go through contcar and translate the atoms so they are closest to the pos atom positions
  //we get cartesian coordinates of contcar by transforming the fractional into the
  //relaxed coordinate system

  for(int na=0; na<atom.size(); na++){
    double min_dist=1.0e10;
    int min_i,min_j,min_k;
    for(int i=-1; i<=1; i++){
      for(int j=-1; j<=1; j++){
	for(int k=-1; k<=1; k++){
	  atompos hatom=contcar.atom[na];
	  hatom.fcoord[0]=contcar.atom[na].fcoord[0]+i;
	  hatom.fcoord[1]=contcar.atom[na].fcoord[1]+j;
	  hatom.fcoord[2]=contcar.atom[na].fcoord[2]+k;
	  conv_AtoB(FtoC,hatom.fcoord,hatom.ccoord);
	  double dist=0.0;
	  for(int ii=0; ii<3; ii++)
	    dist=dist+(hatom.ccoord[ii]-atom[na].ccoord[ii])*(hatom.ccoord[ii]-atom[na].ccoord[ii]);
	  dist=sqrt(dist);

	  if(dist<min_dist){
	    min_dist=dist;
	    min_i=i;
	    min_j=j;
	    min_k=k;
	  }
	}
      }
    }
    atompos hatom=contcar.atom[na];
    hatom.fcoord[0]=contcar.atom[na].fcoord[0]+min_i;
    hatom.fcoord[1]=contcar.atom[na].fcoord[1]+min_j;
    hatom.fcoord[2]=contcar.atom[na].fcoord[2]+min_k;
    conv_AtoB(FtoC,hatom.fcoord,hatom.ccoord);

    for(int i=0; i<3; i++){
      atom[na].dfcoord[i]=hatom.fcoord[i]-atom[na].fcoord[i];
      atom[na].dccoord[i]=hatom.ccoord[i]-atom[na].ccoord[i];
    }
    atom[na].delta=min_dist;
  }


  //calculate the strain matrix and the lattice parameters

  calc_strain();
  get_latparam();
  get_ideal_latparam();


  //write some of this info into a RELAX file in the directory containing POS and CONTCAR
  string relax_file=dir_name;
  relax_file.append("/RELAX");
  ofstream out_relax;
  out_relax.open(relax_file.c_str());
  if(!out_relax){
    cout << "cannot open file " << relax_file << "\n";
    return;
  }

  out_relax << "CONFIGURATION = " << dir_name << "\n";
  out_relax << "CHANGE IN VOLUME = " << determinant(strain) << "\n";
  out_relax << "\n";
  out_relax << "Original lattice parameters \n";
  for(int ii=0; ii<3; ii++) out_relax << ilatparam[ii] << "  ";
  for(int ii=0; ii<3; ii++) out_relax << ilatangle[ii] << "  ";
  out_relax << "\n";
  out_relax << "\n";

  out_relax << "Lattice parameters after relaxation \n";
  for(int ii=0; ii<3; ii++) out_relax << latparam[ii] << "  ";
  for(int ii=0; ii<3; ii++) out_relax << latangle[ii] << "  ";
  out_relax << "\n";
  out_relax << "\n";

  out_relax << "STRAIN MATRIX \n";
  for(int ii=0; ii<3; ii++){
    for(int jj=0; jj<3; jj++){
      out_relax.precision(8);out_relax.width(15);out_relax.setf(ios::showpoint);
      out_relax << strain[ii][jj] << "  ";
    }
    out_relax << "\n";
  }
  out_relax << "\n";

  out_relax << " Atomic relaxations (in the reference system of the relaxed unit cell) \n";
  int count=0;
  for(int na=0; na<atom.size(); na++){
    if(!(atom[na].occ.name[0] == 'V' && atom[na].occ.name[1] == 'a')){
      count++;
      out_relax << "Atom " << count << " relaxation distance = " << atom[na].delta << " Angstrom\n";
    }
  }

  out_relax.close();


  //now place the relaxed coordinates in the structure and set dfcoord=0, dccoord=0, and delta=0.

  for(int na=0; na<atom.size(); na++){
    for(int j=0; j<3; j++){
      atom[na].fcoord[j]=contcar.atom[na].fcoord[j];
      atom[na].ccoord[j]=contcar.atom[na].ccoord[j];
      atom[na].dfcoord[j]=0.0;
      atom[na].dccoord[j]=0.0;
    }
    atom[na].delta=0.0;
  }


  return;

}


//************************************************************

void structure::update_struc(){
  update_lat();
  calc_fractional();
  calc_cartesian();
  calc_factor_group();
  calc_recip_lat();
  get_recip_latparam();

  //test print out lattice parameters and angles
  cout << "lattice parameters a b c \n";
  for(int i=0; i<3; i++)cout << latparam[i] << " ";
  cout << "\n";
  cout << "lattice angle alpha beta gamma\n";
  for(int i=0; i<3; i++)cout << latangle[i] << " ";
  cout << "\n";
  cout << "lattice parameters in order of descending length \n";
  for(int i=0; i<3; i++)cout << latparam[permut[i]] << " ";
  cout << "\n";

}


//************************************************************

void structure::get_recip_latparam(){

  latticeparam(recip_lat, recip_latparam, recip_latangle, recip_permut);

}



//************************************************************

void concentration::collect_components(structure &prim){

  if(prim.atom.size() == 0){
    cout << "cannot collect_components since no atoms in structure \n";
    return;
  }

  for(int i=0; i<prim.atom.size(); i++){
    prim.atom[i].assign_spin();
  }

  //find the first atom that has a min_num_components=2
  int i=0;
  while(prim.atom[i].compon.size() < 2 && i<prim.atom.size()){
    i++;
  }
  if(i==prim.atom.size()){
    cout << "no atoms with more than or equal to 2 components\n";
    cout << "not collecting components\n";
    return;
  }

  compon.clear();
  compon.push_back(prim.atom[i].compon);

  //fill the occup vector with zeros
  occup.clear();
  mu.clear();
  vector<double> toccup;
  for(int k=0; k<prim.atom[i].compon.size(); k++){
    toccup.push_back(0);
  }
  occup.push_back(toccup);
  mu.push_back(toccup);

  for(int j=i+1; j<prim.atom.size(); j++){
    if(prim.atom[j].compon.size() >= 2){
      int k=0;
      for(k=0; k<compon.size(); k++)
	if(compare(prim.atom[j].compon,compon[k])) break;
      if(k==compon.size()){
	compon.push_back(prim.atom[j].compon);
	vector<double> toccup;
	for(int l=0; l<prim.atom[j].compon.size(); l++){
	  toccup.push_back(0);
	}
	occup.push_back(toccup);
	mu.push_back(toccup);
      }
    }
  }

  return;
}



//************************************************************

void concentration::calc_concentration(structure &struc){

  occup.clear();
  for(int i=0; i<compon.size(); i++){
    vector<double> toccup;
    double total=0.0;
    double correction=0.0;
    for(int j=0; j<compon[i].size(); j++){
      double conc=0.0;
      for(int k=0; k<struc.atom.size(); k++){
	//modified by anton
	if(compare(struc.atom[k].compon,compon[i]) && compare(struc.atom[k].occ,compon[i][j])){
	  if(struc.atom[k].basis_flag != '1'){
	    conc=conc+1.0;
	    total=total+1.0;
	  }
	  else{ //Uncomment followin lines to neglect vacancies when calculating concentrations in the occupation basis
	    //if(!(compon[i][j].name[0] == 'V' && compon[i][j].name[1] == 'a')){ //Commented by John.  
	    conc=conc+1.0;
	    total=total+1.0; //Added by John
	    //  correction=correction-1.0; //Commented by John.  
	    //} //Commented by John.  
	  }
	}
      }
      toccup.push_back(conc);
    }
    if(total > tol)
      for(int j=0; j<compon[i].size(); j++){
	//modified by anton - not fool proof
	if(compon[i][j].name[0] == 'V' && compon[i][j].name[1] == 'a') toccup[j]=toccup[j]+correction;
	toccup[j]=toccup[j]/total;
      }
    occup.push_back(toccup);
  }
}



//************************************************************

void concentration::print_concentration(ostream &stream){
  for(int i=0; i<compon.size(); i++){
    for(int j=0; j<compon[i].size()-1; j++){
      compon[i][j].print(stream);
      stream << "=" << occup[i][j] << "  ";
    }
  }
}


//************************************************************

void concentration::print_concentration_without_names(ostream &stream){
  for(int i=0; i<compon.size(); i++){
    for(int j=0; j<compon[i].size()-1; j++){
      stream << occup[i][j] << "  ";
    }
  }
}

//************************************************************

void concentration::print_names(ostream &stream){
  for(int i=0; i<compon.size(); i++){
    for(int j=0; j<compon[i].size()-1; j++){
      compon[i][j].print(stream);
      stream << "  ";
    }
  }
}

//************************************************************

void concentration::get_occup(istream &stream){  // added by jishnu
  for(int i=0; i<compon.size(); i++){
    double sum=0.0;
    for(int j=0; j<compon[i].size()-1; j++){			
      stream >> occup[i][j];
      sum=sum+occup[i][j];
    }
    occup[i][compon[i].size()-1]=1.0-sum;
  }
}

//************************************************************

void arrangement::get_bit(istream &stream){  // added by jishnu
  int bit1;
  char ch;
  while(!(stream.peek()=='\n')){			
    stream.get(ch);
    if(ch==' ') bit1=110;
    if(ch=='0') bit1=0;
    if(ch=='1') bit1=1;
    if(ch=='2') bit1=2;
    if(ch=='3') bit1=3;
    if(ch=='4') bit1=4;
    if(ch=='5') bit1=5;
    if(ch=='6') bit1=6;
    if(ch=='7') bit1=7;
    if(ch=='8') bit1=8;
    if(ch=='9') bit1=9;			
    if(bit1!=110){
      bit.push_back(bit1);
    }			
  };			
			
}

//************************************************************
//************************************************************

void arrangement::update_ce(){  // added by jishnu
  if (ce==1) return;
  fenergy = cefenergy;
  assemble_coordinate_fenergy();
  ce = 1;
  fp = 0;
  te = 0;
  return;
}
//************************************************************

void arrangement::update_fp(){  // added by jishnu
  if (fp==1) return;
  fenergy = fpfenergy;
  assemble_coordinate_fenergy();
  ce = 0;
  fp = 1;
  te = 0;
  return;
}//************************************************************

void arrangement::update_te(){  // added by jishnu
  if (te==1) return;
  fenergy = energy;
  assemble_coordinate_fenergy();
  ce = 0;
  fp = 0;
  te = 1;
  return;
}
//************************************************************

void concentration::set_zero(){
  for(int i=0; i<occup.size(); i++){
    for(int j=0; j<occup[i].size(); j++){
      occup[i][j]=0.0;
    }
  }
}

//************************************************************

void concentration::normalize(int n){
  for(int i=0; i<occup.size(); i++){
    for(int j=0; j<occup[i].size(); j++){
      occup[i][j]=occup[i][j]/n;
    }
  }

}




//************************************************************

void concentration::increment(concentration conc){
  if(conc.occup.size() != occup.size()){
    cout << "in increment, concentrations are not compatible\n";
    return;
  }
  for(int i=0; i<conc.occup.size(); i++){
    if(conc.occup[i].size() != occup[i].size()){
      cout << "in increment, concentrations are not compatible\n";
    }
    for(int j=0; j<occup[i].size(); j++){
      occup[i][j]=occup[i][j]+conc.occup[i][j];
    }
  }
}




//************************************************************


void cluster::get_dimensions(){
  double dist,diff;
  min_leng=1.0e20;
  max_leng=0.0;
  for(int i=0; i<point.size(); i++)
    for(int j=0; j<point.size(); j++)
      if(i != j){
	dist=0;
	for(int k=0; k<3; k++){
	  diff=(point[i].ccoord[k]-point[j].ccoord[k]);
	  dist=dist+diff*diff;
	}
	dist=sqrt(dist);
	if(dist < min_leng)min_leng=dist;
	if(dist > max_leng)max_leng=dist;
      }
}


//************************************************************

cluster cluster::apply_sym(sym_op op){
  cluster tclust;
  tclust.min_leng=min_leng;
  tclust.max_leng=max_leng;
  for(int i=0; i<point.size(); i++){
    atompos tatom;
    tatom=point[i].apply_sym(op);
    tclust.point.push_back(tatom);
  }
  tclust.clust_group.clear();
  return tclust;
}


//************************************************************

void cluster::get_cart(double FtoC[3][3]){
  int np;
  for(np=0; np<point.size(); np++)
    conv_AtoB(FtoC,point[np].fcoord,point[np].ccoord);
}


//************************************************************

void cluster::get_frac(double CtoF[3][3]){
  int np;
  for(np=0; np<point.size(); np++)
    conv_AtoB(CtoF,point[np].ccoord,point[np].fcoord);
}


//************************************************************

void cluster::readf(istream &stream, int np){
  //clear out the points before reading new ones
  point.clear();
  for(int i=0; i<np; i++){
    atompos tatom;
    tatom.readf(stream);
    point.push_back(tatom);
  }

}


//************************************************************

void cluster::readc(istream &stream, int np){
  //clear out the points before reading new ones
  point.clear();
  for(int i=0; i<np; i++){
    atompos tatom;
    tatom.readc(stream);
    point.push_back(tatom);
  }

}


//************************************************************

void cluster::print(ostream &stream){
  for(int i=0; i<point.size(); i++){
    point[i].print(stream);
  }
}


//************************************************************

void cluster::write_clust_group(ostream &stream){
  int cg;

  stream << "cluster group for cluster \n";
  print(stream);
  stream << "\n";
  stream << " number of cluster group ops " << clust_group.size() << "\n";

  for(cg=0; cg<clust_group.size(); cg++){
    stream << "cluster group operation " << cg << " \n";
    clust_group[cg].print_fsym_mat(stream);
    stream << "\n";
  }
}


//************************************************************

void cluster::determine_site_attributes(structure prim){
  for(int i=0; i<point.size(); i++){
    //determine which prim site this point maps onto
    //when there is a match - copy all the attributes from prim onto that site

    for(int j=0; j < prim.atom.size(); j++){
      int trans[3];
      if(compare(point[i],prim.atom[j],trans)){
	//copy attributes from prim.atom[j] onto point[i]
	point[i]=prim.atom[j];
	for(int k=0; k<3; k++) point[i].fcoord[k]=point[i].fcoord[k]+trans[k];
      }
    }

  }

}





//************************************************************

void orbit::readf(istream &stream, int np, int mult){
  //clear out the equivalent clusters before reading new ones
  equiv.clear();
  for(int nm=0; nm<mult; nm++){
    char buff[200];
    stream.getline(buff,199);
    cluster tclust;
    tclust.readf(stream,np);
    equiv.push_back(tclust);
  }

}


//************************************************************

void orbit::readc(istream &stream, int np, int mult){
  //clear out the equivalent clusters before reading new ones
  equiv.clear();
  for(int nm=0; nm<mult; nm++){
    char buff[200];
    stream.getline(buff,199);
    cluster tclust;
    tclust.readc(stream,np);
    equiv.push_back(tclust);
  }


}


//************************************************************

void orbit::print(ostream &stream){
  for(int i=0; i<equiv.size(); i++){
    stream << "equivalent cluster " << i+1 << "\n";
    equiv[i].print(stream);
  }
}


//************************************************************

void orbit::get_cart(double FtoC[3][3]){
  for(int ne=0; ne < equiv.size(); ne++){
    equiv[ne].get_cart(FtoC);
    equiv[ne].get_dimensions();
  }
  return;
}



//************************************************************

void orbit::get_frac(double CtoF[3][3]){
  for(int ne=0; ne < equiv.size(); ne++){
    equiv[ne].get_frac(CtoF);
  }
  return;
}


//************************************************************

void orbit::determine_site_attributes(structure prim){
  for(int ne=0; ne < equiv.size(); ne++){
    equiv[ne].determine_site_attributes(prim);
  }
}



//************************************************************

void multiplet::readf(istream &stream){
  char buff[200];
  char bull;

  //clear out the orbit before reading clusters
  orb.clear();
  size.clear();
  order.clear();
  index.clear();

  vector<orbit> torbvec;

  //make the empty cluster and put it in orb
  {
    cluster tclust;
    tclust.max_leng=0;
    tclust.min_leng=0;
    orbit torb;
    torb.equiv.push_back(tclust);
    torbvec.push_back(torb);
    //the first index i.e. 0 is always for the empty cluster
    //    int i=0;
    size.push_back(0);
    order.push_back(0);
  }


  //read the orbits and collect them torbvec

  int max_np=0;
  int nc;
  stream >> nc;
  stream.get(bull);
  for(int i=1; i<= nc; i++){
    stream.getline(buff,199);

    int dummy,np,mult;
    stream >> dummy;
    stream >> np;
    stream >> mult;
    stream >> dummy;
    stream.getline(buff,199);
    if(np > max_np) max_np=np;
    orbit torb;
    torb.readf(stream,np,mult);
    size.push_back(np);
    order.push_back(0);
    torbvec.push_back(torb);
  }



  //then for all cluster sizes less than or equal to max_np, collect all orbits of the same size
  //we also keep track of the indexing so we remember the order in which the clusters were input
  //(necessary for eci match up)

  for(int np=0; np<=max_np; np++){
    
    vector<orbit> orbvec;
    vector<int> tindex;
    for(int i=0; i<=nc; i++){
      if(size[i] == np){
	orbvec.push_back(torbvec[i]);
	tindex.push_back(i);
	order[i]=tindex.size()-1;
      }
    }
    orb.push_back(orbvec);
    index.push_back(tindex);
  }

}



//************************************************************

void multiplet::readc(istream &stream){

}



//************************************************************

void multiplet::print(ostream &stream){
  for(int i=0; i<orb.size(); i++)
    for(int j=0; j<orb[i].size(); j++)
      orb[i][j].print(stream);
}


//************************************************************

void multiplet::sort(int np){

  for(int i=0; i<orb[np].size(); i++){
    for(int j=i+1; j<orb[np].size(); j++){
      if(orb[np][i].equiv[0].max_leng > orb[np][j].equiv[0].max_leng){
	orbit torb=orb[np][j];
	orb[np][j]=orb[np][i];
	orb[np][i]=torb;
      }
    }
  }
}


//************************************************************

void multiplet::get_index(){
  int count=0;

  size.clear();
  order.clear();
  index.clear();

  //first the emtpy cluster
  size.push_back(0);
  order.push_back(0);
  vector <int> tindex;
  tindex.push_back(count);
  index.push_back(tindex);

  for(int np=1; np<orb.size(); np++){
    vector <int> tindex;
    for(int nc=0; nc<orb[np].size(); nc++){
      count++;
      tindex.push_back(count);
      size.push_back(np);
      order.push_back(nc);
    }
    index.push_back(tindex);
  }
  
  return;

}


//************************************************************

void multiplet::get_hierarchy(){

  int count=0;
 
  size.clear();
  order.clear();
  index.clear();
  subcluster.clear();

  //first the emtpy cluster
  size.push_back(0);
  order.push_back(0);
  vector <int> tindex;
  tindex.push_back(count);
  index.push_back(tindex);

  for(int np=1; np<orb.size(); np++){
    vector <int> tindex;
    for(int nc=0; nc<orb[np].size(); nc++){
      count++;
      tindex.push_back(count);
      size.push_back(np);
      order.push_back(nc);
    }
    index.push_back(tindex);
  }


  // make the subcluster table for the empty cluster
  {
    vector<int>temptysubcluster;
    subcluster.push_back(temptysubcluster);
  }

  // make theh subcluster tables for the non-empty clusters
  for(int np=1; np<orb.size(); np++){

    for(int nc=0; nc<orb[np].size(); nc++){
      vector<int> tsubcluster;

      //enumerate all subclusters of the cluster orbit with index (np,nc)
      //find which cluster this subcluster is equivalent to
      //record the result in tsubcluster

      for(int snp=1; snp<np; snp++){
	vector<int> sc;            // contains the indices of the subcluster
	for(int i=0; i<snp; i++) sc.push_back(i);

	while(sc[0]<=(np-snp)){
	  while(sc[snp-1]<=(np-1)){

	    //BLOCK WHERE SUB CLUSTER IS FOUND AMONG THE LIST

	    cluster tclust;
	    for(int i=0; i<snp; i++)
	      tclust.point.push_back(orb[np][nc].equiv[0].point[sc[i]]);
	    within(tclust);



	    //compare the subclusters with all clusters of the same size

	    for(int i=0; i<orb[snp].size(); i++){
	      if(!new_clust(tclust,orb[snp][i])){
		// snp,i is a subcluster
		// check whether it already exists among the subcluster list

		int j;
		for(j=0; j<tsubcluster.size(); j++){
		  if(index[snp][i] == tsubcluster[j])break;
		}
		if(j== tsubcluster.size()) tsubcluster.push_back(index[snp][i]);

	      }

	    }

	    //END OF BLOCK TO DETERMINE SUBCLUSTERS

	    sc[snp-1]++;
	  }
	  int j=snp-2;
	  if(j>-1){
	    while(sc[j] == (np-(snp-j)) && j>-1) j--;

	    if(j>-1){
	      sc[j]++;
	      for(int k=j+1; k<snp; k++) sc[k]=sc[k-1]+1;
	    }
	    else break;

	  }
	}
      }
      // extra for pair clusters (and possibly point clusters)

      if(np == 2){
	if(nc>0){
	  int k=nc-1;
	  while(k>=0 && abs(orb[np][k].equiv[0].max_leng-orb[np][nc].equiv[0].max_leng) < 0.0001)k--;
	  if(k>=0 && orb[np][nc].equiv[0].max_leng-orb[np][k].equiv[0].max_leng >= 0.0001){
	    tsubcluster.push_back(index[np][k]);
	  }
	}

      }
		
      subcluster.push_back(tsubcluster);
    }
  }

}


//************************************************************

void multiplet::print_hierarchy(ostream &stream){
  stream << "label    weight    mult    size    length    heirarchy \n";
  for(int i=0; i<subcluster.size(); i++){
    stream << i << "   " << "   0    " << orb[size[i]][order[i]].equiv.size() << "   ";
    stream << orb[size[i]][order[i]].equiv[0].point.size() << "   ";
    stream << orb[size[i]][order[i]].equiv[0].max_leng << "   ";
    stream << subcluster[i].size() << "   ";
    for(int j=0; j<subcluster[i].size(); j++) stream << subcluster[i][j] << "   ";
    stream << "\n";
  }

}


//************************************************************

void multiplet::read_eci(istream &stream){
  char buff[200];
  double eci1,eci2;
  int index;
  get_hierarchy();
  for(int i=0; i<7; i++)stream.getline(buff,199);
  while(!stream.eof()){
    stream >> eci1;
    stream >> eci2;
    stream >> index;
    if(index >=size.size()){
      cout << "WARNING:  eci.out has cluster indeces larger than total number of clusters.  Please check for source of incompatibility.  Exiting...\n";
      exit(1);
    }
    stream.getline(buff,199);
    //cout << eci1 << "  " << eci2 << "  " << index << "\n";
    orb[size[index]][order[index]].eci=eci2;
  }

}


//************************************************************

void multiplet::get_cart(double FtoC[3][3]){
  for(int np=1; np < orb.size(); np++){
    for(int ne=0; ne < orb[np].size(); ne++){
      orb[np][ne].get_cart(FtoC);
    }
  }
  return;
}





//************************************************************

void multiplet::get_frac(double CtoF[3][3]){
  for(int np=1; np < orb.size(); np++){
    for(int ne=0; ne < orb[np].size(); ne++){
      orb[np][ne].get_frac(CtoF);
    }
  }
  return;
}



//************************************************************

void multiplet::determine_site_attributes(structure prim){
  for(int np=0; np<orb.size(); np++){
    for(int no=0; no<orb[np].size(); no++){
      orb[np][no].determine_site_attributes(prim);
    }
  }

}


//************************************************************

void arrangement::assemble_coordinate_fenergy(){
  coordinate.clear();
  for(int i=0; i<conc.occup.size(); i++){
    for(int j=0; j<conc.occup[i].size()-1; j++){
      coordinate.push_back(conc.occup[i][j]);
    }
  }
  coordinate.push_back(fenergy);

}


//************************************************************
//************************************************************
/*
  void arrangement::assemble_coordinate_CE(){
  coordinate_CE.clear();
  for(int i=0; i<conc.occup.size(); i++){
  for(int j=0; j<conc.occup[i].size()-1; j++){
  coordinate_CE.push_back(conc.occup[i][j]);
  }
  }
  coordinate_CE.push_back(cefenergy);

  }*/


//************************************************************
void arrangement::print_bit(ostream &stream){
  for(int i=0; i<bit.size(); i++) stream << bit[i] << " ";
  stream << "\n";
}



//************************************************************

void arrangement::print_correlations(ostream &stream){
  for(int i=0; i<correlations.size(); i++){
    stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
    stream << correlations[i] << " ";
  }
  stream << "\n";

}


//************************************************************

void arrangement::print_coordinate(ostream &stream){
	
	
  for(int i=0; i<coordinate.size(); i++){
    stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
    stream << coordinate[i] << " ";
  }
  stream << name;
  stream << "\n";

}


//************************************************************
//************************************************************

void arrangement::print_in_energy_file(ostream &stream){    //added by jishnu
  stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
  stream << fenergy << "   ";
  stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
  stream << weight << "   ";
  for(int i=0; i<coordinate.size()-1; i++){
    stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
    stream << coordinate[i] << " ";
  }
  stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
  stream << delE_from_facet << "   ";
  stream.precision(8);stream.width(12);stream.setf(ios::showpoint);
  stream << name;   
  stream << "\n";

}


//************************************************************

void superstructure::decorate_superstructure(arrangement conf){

  if(struc.atom.size() != conf.bit.size()){
    cout << "inside decorate_superstructure and conf.bit.size() is \n";
    cout << "not equal to the number of atoms in the structure \n";
    exit(1);
  }


  for(int i=0; i<struc.atom.size(); i++){
    struc.atom[i].occ=struc.atom[i].compon[conf.bit[i]];
  }


  //collect the components within the structure
  struc.collect_components();
  curr_struc=struc;
  curr_struc.num_each_specie.clear();


  int curr_ind=0;
  for(int nc=0; nc<struc.compon.size(); nc++){
    int num=0;
    for(int i=0; i<struc.atom.size(); i++){
      if(compare(struc.compon[nc],struc.atom[i].occ)){
	curr_struc.atom[curr_ind]=struc.atom[i];
	curr_ind++;
	num++;
      }
    }
    curr_struc.num_each_specie.push_back(num);
  }


}


//************************************************************

void superstructure::determine_kpoint_grid(double kpoint_dens){

  struc.calc_recip_lat();
  struc.get_recip_latparam();
  double recip_vol=determinant(struc.recip_lat);

  double mleng=struc.recip_latparam[struc.recip_permut[0]];

  int i=0;
  do{
    i++;
    double delta=mleng/i;
    for(int j=0; j<3; j++)
      kmesh[j]=int(ceil(struc.recip_latparam[j]/delta));
  }while((kmesh[0]*kmesh[1]*kmesh[2])/recip_vol < kpoint_dens && i < 99);

  if(i > 99){
    cout << "k-point grid is unusually high \n";
  }


}


//************************************************************

void superstructure::print_incar(string dir_name){
	
	
  if(!scandirectory(dir_name,"INCAR")){	
		
    string outfile = dir_name;    
    outfile.append("/INCAR");
        
    string infile = "INCAR";		
		
    ifstream in;
    in.open(infile.c_str());
		
    if(!in){
      cout << "cannot open parent INCAR file \n";
      cout << "no INCAR can be created for the configurations.\n";
      return;
    }
		
    ofstream out(outfile.c_str());
		
    if(!out){
      cout << "no INCAR created for the configuration " <<  dir_name <<".\n";
      return;
    }
		
    string line;
    while(getline(in,line)){			
      // adding magnetic moments
      bool spin_pol=false;
      string check=line.substr(0,5);
      string check2=line.substr(0,6);
      if (  (check == "ISPIN") && (check2 != "ISPIND") ){
	for(int i=0;i<line.size();i++){
	  if(line[i]=='1') {spin_pol=false;break;}
	  if(line[i]=='2') {spin_pol=true;break;}
	}
      }
      if(spin_pol) {
	out << "MAGMOM = ";
	for(int i=0; i<curr_struc.atom.size(); i++){					
	  if(!(curr_struc.atom[i].occ.name[0] == 'V' && curr_struc.atom[i].occ.name[1] == 'a')){
	    out << curr_struc.atom[i].occ.magmom;  
	    out << " ";
	  }
	}				
	out << "\n";				
      }  // end of adding magmoms
      // writing L's
      bool Lline=false;
      check=line.substr(0,5);			
      if (check == "LDAUL") {	
	Lline = true;
	out << "LDAUL = ";
	for(int i=0; i < curr_struc.num_each_specie.size(); i++){
	  if(!(curr_struc.compon[i].name[0] == 'V' && curr_struc.compon[i].name[1] == 'a') && curr_struc.num_each_specie[i] != 0)
	    out << "  2";
	}				
	out << "\n";				
      }  // end of writing L's
      // writing U's
      bool Uline=false;
      check=line.substr(0,5);			
      if (check == "LDAUU") {	
	Uline = true;
	out << "LDAUU = ";
	for(int i=0; i < curr_struc.num_each_specie.size(); i++){
	  if(!(curr_struc.compon[i].name[0] == 'V' && curr_struc.compon[i].name[1] == 'a') && curr_struc.num_each_specie[i] != 0)
	    out << " " << curr_struc.compon[i].U;
	}				
	out << "\n";				
      }  // end of writing U's
      // writing J's
      bool Jline=false;
      check=line.substr(0,5);			
      if (check == "LDAUJ") {	
	Jline = true;
	out << "LDAUJ = ";
	for(int i=0; i < curr_struc.num_each_specie.size(); i++){
	  if(!(curr_struc.compon[i].name[0] == 'V' && curr_struc.compon[i].name[1] == 'a') && curr_struc.num_each_specie[i] != 0)
	    out << " " << curr_struc.compon[i].J;
	}				
	out << "\n";				
      }  // end of writing J's
			
      if(!Lline && !Uline && !Jline) {
	out << line;
	out << "\n";
      }
			
			
    }	
		
    in.close();		
    out.close();			
		
    return;
  }
  return;
	
}


//************************************************************
//************************************************************

void superstructure::print(string dir_name, string file_name){

  if(!scandirectory(dir_name,file_name)){

    string file=dir_name;
    file.append("/");
    file.append(file_name);

    ofstream out(file.c_str());

    if(!out){
      cout << "cannot open " << file << "\n";
      return;
    }

    out << curr_struc.title <<"\n";

    out.precision(7);out.width(12);out.setf(ios::showpoint);
    out << curr_struc.scale <<"\n";

    for(int i=0; i<3; i++){
      out << "  ";
      for(int j=0; j<3; j++){

	out.precision(9);out.width(15);out.setf(ios::showpoint);
	out << curr_struc.lat[i][j] << " ";

      }
      out << "\n";
    }
    for(int i=0; i < curr_struc.num_each_specie.size(); i++){
      if(!(curr_struc.compon[i].name[0] == 'V' && curr_struc.compon[i].name[1] == 'a') && curr_struc.num_each_specie[i] != 0)
	out << " " << curr_struc.num_each_specie[i] ;
    }
    out << "\n";

    out << "Direct\n";

    for(int i=0; i<curr_struc.atom.size(); i++){
      if(!(curr_struc.atom[i].occ.name[0] == 'V' && curr_struc.atom[i].occ.name[1] == 'a')){
	for(int j=0; j<3; j++){
	  out.precision(9);out.width(15);out.setf(ios::showpoint);
	  out << curr_struc.atom[i].fcoord[j] << " ";
	}
	out << curr_struc.atom[i].occ.name;  // jishnu
	out << "\n";
      }
    }
    out.close();
    return;
  }

}


//************************************************************

void superstructure::print_potcar(string dir_name){

  string potcar_file=dir_name;
  potcar_file.append("/POTCAR");

  //  ofstream out;
  //  out.open(potcar_file.c_str());
  //  if(!out){
  //    cout << "cannot open POTCAR file.\n";
  //    return;
  //  }
  //  out.close();

  //look at the first atom
  string last_element;

  string command = "cat ";

  if(!(curr_struc.atom[0].occ.name[0] == 'V' && curr_struc.atom[0].occ.name[1] == 'a')){
    string element=curr_struc.atom[0].occ.name;  // jishnu
    string potcar="POTCAR_";
    potcar.append(element);
    command.append(potcar);

    last_element=element;
  }

  for(int i=1; i<curr_struc.atom.size(); i++){
    if(!(curr_struc.atom[i].occ.name[0] == 'V' && curr_struc.atom[i].occ.name[1] == 'a')){
      string element=curr_struc.atom[i].occ.name;  // jishnu
      if(element != last_element){
	string potcar="POTCAR_";
	potcar.append(element);
	command.append(" ");
	command.append(potcar);

	last_element=element;
      }

    }
  }

  command.append(" > ");
  command.append(potcar_file);
  //    cout << command << "\n";

  int s=system(command.c_str());
  if(s == -1){cout << "was unable to perform system command\n";}


  return;
}





//************************************************************

void superstructure::print_kpoint(string dir_name){


  if(!scandirectory(dir_name,"KPOINTS")){

    string file_name=dir_name;
    file_name.append("/KPOINTS");

    ofstream out;
    out.open(file_name.c_str());
    if(!out){
      cout << "cannot open " << file_name << "\n";
      return;
    }
    out << dir_name << "\n";
    out << " 0 \n";
    out << "Gamma point shift\n";
    for(int j=0; j<3; j++)out << " " << kmesh[j];
    out << "\n";
    out << " 0 0 0 \n";
    out.close();
  }

}


//************************************************************
//************************************************************
//read the configuration folder names to check whether a particular one is already present   //added by jishnu
bool read_vasp_list_file(string name) {

  string s;
  ifstream rf;	      
  rf.open("vasp_list_file",ios::out);
  do{
    rf>>s;
    if (s==name) {rf.close();return true;}					  
  }while (!rf.eof());	
  return false;
}
//************************************************************
//************************************************************
//write the configuration folder names to a file so that automatic submitvasp can work   //added by jishnu
void write_vasp_list_file(string name) {

  bool ifthere;
  if(!scandirectory(".","vasp_list_file")){
    string file_name="vasp_list_file";
    ofstream out;
    out.open(file_name.c_str());
    if(!out){
      cout << "cannot open " << file_name << "\n";
      return;
    }
    out.close();
  }

  ifthere=read_vasp_list_file(name);

  if(!ifthere) {
    ofstream out;
    out.open("vasp_list_file",ios::app);
    out << "D "<< name <<"\n";
    out.close();
  }  // end of if(!ifthere)
  
  return;

}
//************************************************************
//************************************************************
//write the structure to a file in yihaw like format   //added by jishnu

void superstructure::print_yihaw(string dir_name) {


  if(!scandirectory(dir_name,"yihaw")){

    string file_name=dir_name;
    file_name.append("/yihaw");

    ofstream out;
    out.open(file_name.c_str());
    if(!out){
      cout << "cannot open " << file_name << "\n";
      return;
    }
         
    out << "#!/bin/sh" <<"\n";
    out << "#PBS -S /bin/sh" <<"\n";
    out << "#PBS -N vasp" <<"\n";
    out << "#PBS -l nodes="<< nodes <<":ppn="<< ppn <<",walltime=" << walltime <<":00:00"<<"\n";
    out << "#PBS -q "<< queue <<"\n";
    out << "\n";
    out << "#PBS -o " << parent_directory << "/" <<dir_name<<"\n";
    out << "#PBS -e " << parent_directory << "/" <<dir_name<<"\n";
    out << "#PBS -joe" <<"\n";
    out << "#PBS -V" <<"\n";
    out << "#" <<"\n";
    out << "\n";
    //out << "echo ""<< "l ran on:"<<"""<<"\n";
    out << "cat $PBS_NODEFILE" <<"\n";
    out << "#" <<"\n";
    out << "# Change to your execution directory." <<"\n";
    out << "cd " << parent_directory << "/" <<dir_name<<"\n";
    out << "#" <<"\n";
    out << "\n";
    out << "lamboot" <<"\n";
    out << "\n";
    out << "#" <<"\n";
    out << "# Use mpirun to run with "<< ppn << " cpu for "<< walltime <<" hours" <<"\n";
    out << "\n";
    out << "mpirun  -np "<< ppn << " vasp" <<"\n";
    out << "\n";
    out << "lamhalt" <<"\n";
    out << "#" <<"\n";
  
    out.close();
  
  }
}



//************************************************************
//************************************************************
//read the yihaw_input file   //added by jishnu

void superstructure::read_yihaw_input() {

  if(!scandirectory(".","yihaw_input")){
    cout << "No yihaw_input file to open  \n";
  }

  ifstream readfrom;
  int n;
  char ch;
  readfrom.open("yihaw_input",ios::out);
  n=0;
  do{
    readfrom.get(ch);
    if(ch=='=') n=n+1;
  }while(n<1);
  readfrom>>nodes;
  n=0;
  do{
    readfrom.get(ch);
    if(ch=='=') n=n+1;
  }while(n<1);
  readfrom>>ppn;
  n=0;
  do{
    readfrom.get(ch);
    if(ch=='=') n=n+1;
  }while(n<1);
  readfrom>>walltime;
  n=0;
  do{
    readfrom.get(ch);
    if(ch=='=') n=n+1;
  }while(n<1);
  readfrom>>queue;
  n=0;
  do{
    readfrom.get(ch);
    if(ch=='=') n=n+1;
  }while(n<1);
  readfrom>>parent_directory;		

  readfrom.close(); 
    
}

//************************************************************

void structure::calc_recip_lat(){
  double pi=3.141592654;
  double vol=determinant(lat);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) recip_lat[i][j]=2.0*pi*(lat[(i+1)%3][(j+1)%3]*lat[(i+2)%3][(j+2)%3]-
						   lat[(i+1)%3][(j+2)%3]*lat[(i+2)%3][(j+1)%3])/vol;
  }
}



//************************************************************
////////////////////////////////////////////////////////////////////////////////
//added by Ben Swoboda
//takes atompos as input and determines the basis vectors for occupation and spin methods
void get_basis_vectors(atompos &atom){

  int tspin;
  atom.spin_vec.clear();
  atom.p_vec.clear();

  tspin=1;
  //modified by Anton (...compon.size()-1 ...)
  for(int i=0; i < atom.compon.size()-1; i++){
    tspin=tspin*atom.occ.spin;
    atom.spin_vec.push_back(tspin);
  }


  //  atom.p_vec.push_back(1);
  //modified by Anton (i=0 and compon.size()-1 instead of i=1 and compon.size() )
  for(int i=0; i<atom.compon.size()-1; i++){
    if(compare(atom.occ,atom.compon[i]))atom.p_vec.push_back(1);
    else atom.p_vec.push_back(0);
  }

  return;

}
////////////////////////////////////////////////////////////////////////////////
//************************************************************

void configurations::generate_configurations(vector<structure> supercells){
  int ns,i;

  for(ns=0; ns<supercells.size(); ns++){
    superstructure tsuperstruc;
    multiplet super_basiplet;
    //used to be copy_lattice below
    tsuperstruc.struc=supercells[ns];
    tsuperstruc.struc.expand_prim_basis(prim);
    tsuperstruc.struc.expand_prim_clust(basiplet,super_basiplet);

    //generate the different bit combinations

    int last=0;

    while(last == 0){
      tsuperstruc.struc.atom[0].bit++;
      for(i=0; i<(tsuperstruc.struc.atom.size()-1); i++){
	if(tsuperstruc.struc.atom[i].bit !=0 &&
	   tsuperstruc.struc.atom[i].bit%(tsuperstruc.struc.atom[i].compon.size()) == 0){
	  tsuperstruc.struc.atom[i+1].bit++;
	  tsuperstruc.struc.atom[i].bit=0;
	}
      }
      if(tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit !=0 &&
	 tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit%(tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].compon.size()) == 0){
	last=last+1;
	tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit=0;
      }

      //for each atom position, assign its spin and specie for this configuration of bits

      arrangement tconf;

      for(i=0; i<tsuperstruc.struc.atom.size(); i++){
	tsuperstruc.struc.atom[i].occ=tsuperstruc.struc.atom[i].compon[tsuperstruc.struc.atom[i].bit];
	////////////////////////////////////////////////////////////////////////////////
	//modified by Ben Swoboda
	//assign the spin to the atompos object for use in basis
	get_basis_vectors(tsuperstruc.struc.atom[i]);


	//cout  << "p-vec[" << i << "]\n";
	//for(int x=0; x<tsuperstruc.struc.atom[i].p_vec.size(); x++){
	//    cout << tsuperstruc.struc.atom[i].p_vec[x] << "\t";
	//}
	//cout << "\n";
	//cout << "bit[" << i << "]: " << tsuperstruc.struc.atom[i].bit << "\tspin: " << tsuperstruc.struc.atom[i].occ.spin << "\tspecie: "
	//<< tsuperstruc.struc.atom[i].occ.name << "\n";
	////////////////////////////////////////////////////////////////////////////////
	tconf.bit.push_back(tsuperstruc.struc.atom[i].bit);
      }
      calc_correlations(tsuperstruc.struc,super_basiplet,tconf);


      //calculate the concentration over the sites with more than min_num_components = 2
      tconf.conc.collect_components(prim);
      tconf.conc.calc_concentration(tsuperstruc.struc);


      if(new_conf(tconf,superstruc) && new_conf(tconf,tsuperstruc)){

	//give this newly found configuration its name
	tconf.name="con";
	string scel_num;
	string period=".";
	int_to_string(ns,scel_num,10);
	tconf.name.append(scel_num);
	tconf.name.append(period);
	string conf_num;
	int_to_string(tsuperstruc.conf.size(),conf_num,10);
	tconf.name.append(conf_num);

	//record the indices of this configuration
	tconf.ns=ns;
	tconf.nc=tsuperstruc.conf.size();

	//add the new configuration to the list for this superstructure
	tsuperstruc.conf.push_back(tconf);

      }

    }
    superstruc.push_back(tsuperstruc);
  }

  return;
}

//************************************************************

void configurations::generate_configurations_fast(vector<structure> supercells){
  int ns,i,j,k;
  double tcorr, tclust_func;

  for(ns=0; ns<supercells.size(); ns++){
    superstructure tsuperstruc;
    multiplet super_basiplet;
    //used to be copy_lattice below
    tsuperstruc.struc=supercells[ns];
    tsuperstruc.struc.expand_prim_basis(prim);
    tsuperstruc.struc.expand_prim_clust(basiplet,super_basiplet);

    //get cluster function and basis info for this supercell
    get_corr_vector(tsuperstruc.struc, super_basiplet, tsuperstruc.corr_to_atom_vec);
    get_super_basis_vec(tsuperstruc.struc, tsuperstruc.basis_to_bit_vec);

    //generate the different bit combinations

    int last=0;

    while(last == 0){
      tsuperstruc.struc.atom[0].bit++;
      for(i=0; i<(tsuperstruc.struc.atom.size()-1); i++){
	if(tsuperstruc.struc.atom[i].bit !=0 &&
	   tsuperstruc.struc.atom[i].bit%(tsuperstruc.struc.atom[i].compon.size()) == 0){
	  tsuperstruc.struc.atom[i+1].bit++;
	  tsuperstruc.struc.atom[i].bit=0;
	}
      }
      if(tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit !=0 &&
	 tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit%(tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].compon.size()) == 0){
	last=last+1;
	tsuperstruc.struc.atom[tsuperstruc.struc.atom.size()-1].bit=0;
      }

      //for each atom position, assign its spin and specie for this configuration of bits

      arrangement tconf;

      for(i=0; i<tsuperstruc.struc.atom.size(); i++){
	tsuperstruc.struc.atom[i].occ=tsuperstruc.struc.atom[i].compon[tsuperstruc.struc.atom[i].bit];
	tconf.bit.push_back(tsuperstruc.struc.atom[i].bit);
      }
      tconf.correlations.push_back(1.0);
      for(i=0; i<tsuperstruc.corr_to_atom_vec.size(); i++){
	tcorr=0.0;
	for(j=0; j<tsuperstruc.corr_to_atom_vec[i].size(); j++){
	  tclust_func=1.0;
	  for(k=0; k<tsuperstruc.corr_to_atom_vec[i][j].size(); k++){
	    tclust_func*=tsuperstruc.basis_to_bit_vec[tsuperstruc.corr_to_atom_vec[i][j][k][0]][tconf.bit[tsuperstruc.corr_to_atom_vec[i][j][k][0]]][tsuperstruc.corr_to_atom_vec[i][j][k][1]];
	  }
	  tcorr+=tclust_func;
	}
	tconf.correlations.push_back(tcorr/tsuperstruc.corr_to_atom_vec[i].size());
      }




      if(new_conf(tconf,superstruc) && new_conf(tconf,tsuperstruc)){

	//calculate the concentration over the sites with more than min_num_components = 2
	tconf.conc.collect_components(prim);
	tconf.conc.calc_concentration(tsuperstruc.struc);

	//give this newly found configuration its name
	tconf.name="con";
	string scel_num;
	string period=".";
	int_to_string(ns,scel_num,10);
	tconf.name.append(scel_num);
	tconf.name.append(period);
	string conf_num;
	int_to_string(tsuperstruc.conf.size(),conf_num,10);
	tconf.name.append(conf_num);

	//record the indices of this configuration
	tconf.ns=ns;
	tconf.nc=tsuperstruc.conf.size();

	//add the new configuration to the list for this superstructure
	tsuperstruc.conf.push_back(tconf);

      }

    }
    superstruc.push_back(tsuperstruc);
  }

  return;
}



//************************************************************

void configurations::generate_vasp_input_directories(){

  //first determine the kpoints density

  double kpt_dens;
  {
    int kmesh[3];

    // read in the kpoints for the primitive cell
    string kpoint_file="KPOINTS";
    ifstream in_kpoints;
    in_kpoints.open(kpoint_file.c_str());
    if(!in_kpoints){
      cout << "cannot open " << kpoint_file << "\n";
      return;
    }
    char buff[200];
    for(int i=0; i<3; i++) in_kpoints.getline(buff,199);
    for(int j=0; j<3; j++)in_kpoints >> kmesh[j];
    in_kpoints.close();

    double recip_vol=determinant(prim.recip_lat);
    kpt_dens=(kmesh[0]*kmesh[1]*kmesh[2])/recip_vol;
  }


  //create vasp files for all the configurations

  for(int sc=0; sc<superstruc.size(); sc++){

    //make the kpoint-mesh for this superstructure

    superstruc[sc].determine_kpoint_grid(kpt_dens);

    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){

      if(!scandirectory(".",superstruc[sc].conf[nc].name) && superstruc[sc].conf[nc].make){
        write_vasp_list_file(superstruc[sc].conf[nc].name);  // added by jishnu for auto submission of vasp runs
        string command="mkdir ";
        command.append(superstruc[sc].conf[nc].name);
        int s=system(command.c_str());
        if(s == -1){cout << "was unable to perform system command\n";}
      }

      if(scandirectory(".",superstruc[sc].conf[nc].name)){

	// string command="cp INCAR ";    // This part is from the time when we just copied INCAR not write it
	// command.append(superstruc[sc].conf[nc].name);
	// int s=system(command.c_str());
	// if(s == -1){cout << "was unable to perform system command\n";}

	superstruc[sc].decorate_superstructure(superstruc[sc].conf[nc]);

	superstruc[sc].print(superstruc[sc].conf[nc].name,"POSCAR");

	superstruc[sc].print(superstruc[sc].conf[nc].name,"POS");

	superstruc[sc].print_incar(superstruc[sc].conf[nc].name);      // added by jishnu  // this is to explicitly write the INCAR based upon the parent INCAR in parent directory
	
	superstruc[sc].print_potcar(superstruc[sc].conf[nc].name);

	superstruc[sc].print_kpoint(superstruc[sc].conf[nc].name);
	
	//superstruc[sc].read_yihaw_input();   // added by jishnu
	
	//superstruc[sc].print_yihaw(superstruc[sc].conf[nc].name);   // added by jishnu
      }
    }
  }

}


//************************************************************

void configurations::print_con_old(){
  int i,j,k;

  ofstream out;
  out.open("CON");
  if(!out){
    cout << "cannot open CON \n";
    return;
  }

  out << "Structures generated within supercells \n";
  out << "\n";
  out << superstruc.size() << "  supercells considered\n";
  for(int sc=0; sc<superstruc.size(); sc++){
    out << "\n";
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	out.precision(5);out.width(5);
	out << superstruc[sc].struc.slat[i][j] << " ";
      }
      out << "  ";
    }
    out << "\n";
    out << superstruc[sc].conf.size() << " configurations in this supercell\n";
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      superstruc[sc].conf[nc].conc.print_concentration(out);
      superstruc[sc].conf[nc].print_bit(out);
    }
  }

}

//************************************************************

void configurations::read_energy_and_corr(){    // added by jishnu
	
  //first open a log file for mismatching names
  string log_file="energy_names.log";
  ofstream enname;
  enname.open(log_file.c_str());
  if(!enname){
    cout << "cannot open " << log_file << "\n";
    return;
  }
	
	
	
  char buff[300];
  ifstream enin;
  ifstream corrin;
  if(!scandirectory(".","energy")) {
    cout << "No energy file fto read from \n";
    exit(1);
  }
  else {
    enin.open("energy");
  }
  if(!scandirectory(".","corr.in")) {
    cout << "No corr.in file fto read from \n";
    exit(1);
  }
  else {
    corrin.open("corr.in");
  }	
  // collect info from energy file
  enin.getline(buff,299);	
  double fe,wt,co,dfh;
  vector <double> coo;
  string nm;vector <double> vecfe,vecwt,vecdfh;
  vector< vector<double> > veccoo;
  vector<string> vecnm;	
  while(!enin.eof()) {
    enin >> fe >> wt;
    vecfe.push_back(fe);
    vecwt.push_back(wt);	
    for (int i=0; i<superstruc[0].conf[0].coordinate.size()-1; i++) {
      enin >> co;
      coo.push_back(co);
    }
    veccoo.push_back(coo);
    coo.clear();
    enin >> dfh;
    vecdfh.push_back(dfh);
    enin >> nm;
    if(nm.size() > 1) vecnm.push_back(nm);
    nm.erase();
    enin.getline(buff,299);
  }
  //collect info from corr.in file
  int neci,nconf;
  corrin >> neci;
  corrin.getline(buff,299);	
  corrin >> nconf;
  corrin.getline(buff,299);
  corrin.getline(buff,299);
  double elem;
  vector <double> cor;
  vector < vector<double> > veccor;	
  for(int i=0;i<nconf;i++) {
    for (int j=0; j<neci; j++) {
      corrin >> elem;
      cor.push_back(elem);
    }		
    veccor.push_back(cor);				
    cor.clear();		
  }
  // check if the energy and corr.in files are compatible to each other
  cout << "vecnm.size() =" << vecnm.size() << "\n";
  cout << "veccor.size() =" << veccor.size() << "\n";
  if(vecnm.size() != veccor.size()) {
    cout << " The size of the corr.in file is not same as the size in energy file ; Please check what's wrong. \n";
  }
	
  //put these info in superstruc	 -- new approach
  for(int sc=0; sc<superstruc.size(); sc++){				
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      int index = -1;			
      for(int nv=0;nv<vecnm.size();nv++) {
	if(superstruc[sc].conf[nc].name == vecnm[nv]) index = nv;
	if(index != -1) break;
      }
      if(index != -1){
	superstruc[sc].conf[nc].calculated = true;
	superstruc[sc].conf[nc].fenergy = vecfe[index];
	superstruc[sc].conf[nc].fpfenergy = vecfe[index];
	superstruc[sc].conf[nc].weight = vecwt[index];
	superstruc[sc].conf[nc].correlations.clear();
	for(int i=0; i<veccor[index].size(); i++){
	  superstruc[sc].conf[nc].correlations.push_back(veccor[index][i]);
	}	
				
	int i=0;
	for(int j=0; j<superstruc[sc].conf[nc].conc.occup.size(); j++){
	  double sum=0.0;
	  for(int k=0; k<superstruc[sc].conf[nc].conc.occup[j].size()-1; k++){      					
	    superstruc[sc].conf[nc].conc.occup[j][k] = veccoo[index][i];
	    sum=sum+superstruc[sc].conf[nc].conc.occup[j][k];
	    i++;
	  }
	  superstruc[sc].conf[nc].conc.occup[j][superstruc[sc].conf[nc].conc.occup[j].size()-1]=1.0-sum;
	}
	superstruc[sc].conf[nc].assemble_coordinate_fenergy();				
				
      }
    }
  }
	
  string file_name="custom_structures";
  ifstream in_dir;
  in_dir.open(file_name.c_str());
	
  if(!in_dir){
    cout << "cannot open file " << file_name << "\n";
    return;
  }	
  while(!in_dir.eof()){
    string dir_name;
    in_dir >> dir_name;
    if(dir_name.size() > 0){	
      int index = -1;			
      for(int nv=0;nv<vecnm.size();nv++) {
	if(dir_name == vecnm[nv]) index = nv;
	if(index != -1) break;
      }
      if(index != -1){
	superstructure tsup;
	arrangement tarr;				
	tarr.calculated = true;
	tarr.name = vecnm[index];
	tarr.fenergy = vecfe[index];
	tarr.fpfenergy = vecfe[index];
	tarr.weight = vecwt[index];
	tarr.delE_from_facet = vecdfh[index];
	tarr.correlations.clear();
	for(int i=0; i<veccor[index].size(); i++){
	  tarr.correlations.push_back(veccor[index][i]);
	}
				
	for(int j=0; j<superstruc[0].conf[0].conc.compon.size(); j++){  // set the concentration.compon structure same as the genrated structure
	  vector <specie> dummy;
	  for(int k=0; k<superstruc[0].conf[0].conc.compon[j].size(); k++){
	    dummy.push_back(superstruc[0].conf[0].conc.compon[j][k]);						
	  }
	  tarr.conc.compon.push_back(dummy);
	  dummy.clear();  					    
	}
	for(int j=0; j<superstruc[0].conf[0].conc.occup.size(); j++){  // set the concentration.occup structure of the custom_structure same as the generated ones
	  vector <double> dummy;
	  for(int k=0; k<superstruc[0].conf[0].conc.occup[j].size(); k++){
	    dummy.push_back(0.0);
	  }
	  tarr.conc.occup.push_back(dummy);
	  dummy.clear();  					    
	}
  				
	int i=0;
	for(int j=0; j<tarr.conc.occup.size(); j++){
	  double sum=0.0;
	  for(int k=0; k<tarr.conc.occup[j].size()-1; k++){
	    tarr.conc.occup[j][k] = veccoo[index][i];
	    sum=sum+tarr.conc.occup[j][k];
	    i++;
	  }
	  tarr.conc.occup[j][tarr.conc.occup[j].size()-1]=1.0-sum;
	}	
	tarr.assemble_coordinate_fenergy();			
	tsup.conf.push_back(tarr);
	superstruc.push_back(tsup);	
      }
    }
  }
	
	
  in_dir.close();
  enname.close();
	
}
//************************************************************

void configurations::print_corr_old(){
  int i,j,k;
  int num_basis=0;
  int num_conf=0;

  ofstream out;
  out.open("CON.CORR");
  if(!out){
    cout << "cannot open CON.CORR \n";
    return;
  }

  for(i=0; i<basiplet.orb.size(); i++)
    num_basis=num_basis+basiplet.orb[i].size();

  for(i=0; i<superstruc.size(); i++)
    num_conf=num_conf+superstruc[i].conf.size();

  out << num_basis << "  number of basis function correlations\n";
  out << num_conf << "  number of configurations\n";
  out << "correlations \n";
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++)
      superstruc[sc].conf[nc].print_correlations(out);
  }
}




//************************************************************

void configurations::print_con(){
  int i,j,k;

  ofstream out;
  out.open("configuration");
  if(!out){
    cout << "cannot open configuration \n";
    return;
  }

  out << "Structures generated within supercells \n";
  out << "\n";
  out << superstruc.size() << "  supercells considered\n";
  for(int sc=0; sc<superstruc.size(); sc++){
    out << "\n";
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	out.precision(5);out.width(5);
	out << superstruc[sc].struc.slat[i][j] << " ";
      }
      out << "  ";
    }
    out << "\n";
    out << superstruc[sc].conf.size() << " configurations in this supercell\n";
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      out << superstruc[sc].conf[nc].name << "  ";  // added by jishnu
      superstruc[sc].conf[nc].conc.print_concentration_without_names(out);
      superstruc[sc].conf[nc].print_bit(out);
    }
  }

}

//************************************************************

void configurations::read_con(){  // added by jishnu (this s/r)
 
  string junk;
  double value;
	
   
  //----------------
  //test if prim sturcture is there or not   
  if(prim.atom.size()==0){
    cout << "No prim structure read \n";
    return;
  }
  //----------------

  ifstream in;
  in.open("configuration");
  if(!in){
    cout << "cannot open configuration \n";
    return;
  }

  read_junk(in);
  int superstruc_size;
  in >> superstruc_size;   
  read_junk(in);
  for(int sc=0;sc<superstruc_size;sc++){
    superstructure tsuperstructure;	
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	in >> tsuperstructure.struc.slat[i][j];
      }
    }
    tsuperstructure.struc.generate_lat(prim);
    tsuperstructure.struc.expand_prim_basis(prim);
				
    int conf_size;
    in >> conf_size;	
    read_junk(in);
    for (int nc=0;nc<conf_size;nc++){
      arrangement tarrangement;				
      tarrangement.conc.collect_components(prim);
      in >> tarrangement.name;	
      tarrangement.conc.get_occup(in);				
      tarrangement.get_bit(in);
      tsuperstructure.conf.push_back(tarrangement);						  
    }			
    superstruc.push_back(tsuperstructure);
  }	
  
}   // end of s/r read_con()

//************************************************************

void configurations::print_corr(){
  int i,j,k;
  int num_basis=0;
  int num_conf=0;

  ofstream out;
  out.open("configuration.corr");
  if(!out){
    cout << "cannot open configuration.corr \n";
    return;
  }

  for(i=0; i<basiplet.orb.size(); i++)
    num_basis=num_basis+basiplet.orb[i].size();

  for(i=0; i<superstruc.size(); i++)
    num_conf=num_conf+superstruc[i].conf.size();

  out << num_basis << "  number of basis function correlations\n";
  out << num_conf << "  number of configurations\n";
  out << "correlations \n";
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      out << superstruc[sc].conf[nc].name;  // added by jishnu
      superstruc[sc].conf[nc].print_correlations(out);
    }
  }
}  // end of s/r
//************************************************************

void configurations::read_corr(){  // added by jishnu (this s/r) // this must be called after read_con

  string junk;
  int corr_size;
  
  ifstream in;
  in.open("configuration.corr");
  if(!in){
    cout << "cannot open configuration.corr \n";
    return;
  }
  
  in >> corr_size;
  for(int i=0;i<3;i++){
    read_junk(in);
  } 
  double value;
  for(int sc=0; sc<superstruc.size(); sc++){	
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){	
      in >> junk;  // the name of the configuration is aleady read in from configurations file , so no need to do that again or overwrite
      for(int i=0;i<corr_size;i++){						
	in >> value;
	superstruc[sc].conf[nc].correlations.push_back(value);
      }
    }
  }
    
}  // end of the s/r

//************************************************************
void configurations::reconstruct_from_read_files(){   // added by jishnu
	
  read_con();
  read_corr();

}   // end of the s/r
//************************************************************

void hull::clear_arrays(){   // added by jishnu
  point.clear();
  face.clear();		
}
//************************************************************

void configurations::print_make_dirs(){   // modified by jishnu
  if(!scandirectory(".","make_dirs")){
    ofstream out("make_dirs");
    out << "#    name      make      concentrations  \n";
    
    for(int ns=0; ns<superstruc.size(); ns++){
      for(int nc=0; nc<superstruc[ns].conf.size(); nc++){
	out << superstruc[ns].conf[nc].name;
	if(superstruc[ns].conf[nc].make) out << "  1     ";
	else out << "  0     ";
	superstruc[ns].conf[nc].assemble_coordinate_fenergy(); 
	for(int i=0; i<superstruc[ns].conf[nc].coordinate.size()-1; i++){
	  out << superstruc[ns].conf[nc].coordinate[i] << "  ";
	}
	out << " \n";	
      }
    }
    out.close();
  }


}


//************************************************************

void configurations::read_make_dirs(){   // modified by jishnu

  ifstream in("make_dirs");
  if(!in){
    cout << "cannot open file make_dirs\n";
    return;
  }
  char buff[300];
  in.getline(buff,299);
  
  while(!in.eof()){
    string struc_name;
    // double weight;
    int make;    
    in >> struc_name;
    if(struc_name.size() > 0){
      // in >> weight;
      in >> make;
      in.getline(buff,299);

      //among all the structures in the superstructure vector, find that with the same name

      bool match=false;
      if(make){
	for(int ns=0; ns<superstruc.size(); ns++){
	  for(int nc=0; nc<superstruc[ns].conf.size(); nc++){
	    if(struc_name.compare(superstruc[ns].conf[nc].name) == 0){
	      match=true;
	      //superstruc[ns].conf[nc].weight=weight;
	      //if(make == 1) superstruc[ns].conf[nc].make=true;
	      //else superstruc[ns].conf[nc].make=false;
	      superstruc[ns].conf[nc].make=true;
	      break;
	    }
	  }
	  if(match)break;
	}
      }
    }


  }
  in.close();

  //remove the make_dirs file so that the most recent weights can be written

  int s=system("rm make_dirs");
  if(s == -1){
    cout << "was unable to remove make_dirs \n";
  }

}


//************************************************************

void configurations::collect_reference(){

  // reference states are those structures where either all concentration variables are zero or
  // at most one concentration variable is 1.0

  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){

      // go through the concentration object for each arrangement conf[] and search for
      // those in which all concentrations are either zero, or all are zero except for one

      int zero=0; 
      int one=0; 
      int total=0; 
      for(int k=0; k<superstruc[sc].conf[nc].conc.occup.size(); k++){
	for(int l=0; l<superstruc[sc].conf[nc].conc.occup[k].size()-1; l++){
	  total++; 
	  if(abs(superstruc[sc].conf[nc].conc.occup[k][l]) < tol)zero++; 
	  if(abs(superstruc[sc].conf[nc].conc.occup[k][l]-1.0) < tol)one++; 
	}
      }
      if(zero == total){
	reference.push_back(superstruc[sc].conf[nc]); 
      }
      if(zero == total-1 && one == 1){
	// this part added by jishnu
	int count=0; 
	for(int er=0; er<reference.size(); er++){		  
	  if(!compare(superstruc[sc].conf[nc].correlations,reference[er].correlations)) count++; 
	  else{
	    if(superstruc[sc].conf[nc].calculated && reference[er].calculated) {
	      double difference=superstruc[sc].conf[nc].energy-reference[er].energy; 
	      if(difference >=tol) {cout << "CAUTION!! The reference structures are same but the energies are significantly different and we are not adding the new reference\n"; } 
	    }
	    else if(superstruc[sc].conf[nc].calculated && !reference[er].calculated) reference[er]=superstruc[sc].conf[nc]; 
	  }  
	}	
	if(count == reference.size())   reference.push_back(superstruc[sc].conf[nc]); 
	// end of added by jishnu
      }
    }
  }

  if(!scandirectory(".","reference")){
    ofstream out("reference"); 

    out << "concentration and energy of the reference states\n"; 
    out << "\n"; 

    for(int i=0; i<reference.size(); i++){
      out << "reference " << i << "  "; 
      out <<  reference[i].name << "\n"; // added by jishnu
      for(int ii=0; ii<reference[i].conc.compon.size(); ii++){
	out << "sublattice " << ii << "\n"; 
	for(int j=0; j<reference[i].conc.compon[ii].size(); j++){
	  for(int k=0; k<2; k++)
	    out << reference[i].conc.compon[ii][j].name[k]; 
	  out << "  "; 
	}
	out << "\n"; 

	for(int j=0; j<reference[i].conc.compon[ii].size(); j++)
	  out << reference[i].conc.occup[ii][j] << "  "; 
	out << "\n"; 

	out << reference[i].energy << "  energy \n"; 
      }
      out << "\n"; 

    }
    out.close(); 
  }

  ifstream in("reference"); 

  char buff[200]; 
  in.getline(buff,199); 
  
  for(int i=0; i<reference.size(); i++){
    if(in.eof()){
      cout << "reference file is not compatible with current system\n"; 
      cout << "using conventional reference \n"; 
      return; 
    }
    in.getline(buff,199); 
    in.getline(buff,199); 
    for(int ii=0; ii<reference[i].conc.compon.size(); ii++){
      if(in.eof()){
	cout << "reference file is not compatible with current system\n"; 
	cout << "using conventional reference \n"; 
	return; 
      }
      in.getline(buff,199); 
      in.getline(buff,199); 
      for(int j=0; j<reference[i].conc.compon[ii].size(); j++){
	in >> reference[i].conc.occup[ii][j]; 
      }
      in.getline(buff,199); 
      in >> reference[i].energy; 		
      in.getline(buff,199); 
    }
  }
}


//************************************************************

void configurations::collect_energies(){
	
  //first open a log file with problematic relaxations
  string log_file="relax.log";
  ofstream relaxlog;
  relaxlog.open(log_file.c_str());
  if(!relaxlog){
    cout << "cannot open " << log_file << "\n";
    return;
  }
	
	
  for(int sc=0; sc<superstruc.size(); sc++){
    multiplet super_basiplet;
    superstruc[sc].struc.expand_prim_clust(basiplet,super_basiplet);
		
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(scandirectory(".",superstruc[sc].conf[nc].name)){
	if(scandirectory(superstruc[sc].conf[nc].name, "OSZICAR")){
					
	  double energy;
	  int relax_step;    // added by jishnu
					
	  //extract the energy from the OSZICAR file in directory= superstruc[sc].conf[nc].name
					
	  // if(read_oszicar(superstruc[sc].conf[nc].name, energy)){
	  if(read_oszicar(superstruc[sc].conf[nc].name, energy, relax_step)){   // changed by jishnu // used overloaded function to count the no of relax step
						
	    //normalize the energy by the number of primitive unit cells
						
	    double vol=determinant(superstruc[sc].struc.slat);
	    if(abs(vol) > tol){
	      superstruc[sc].conf[nc].energy= energy/abs(vol);
	      superstruc[sc].conf[nc].calculated = true;
	      superstruc[sc].conf[nc].relax_step = relax_step;    // added by jishnu
							
	      //check first if there is a POS file, if not, create it for this configuration				
	      structure relaxed;
	      relaxed.collect_relax(superstruc[sc].conf[nc].name);
							
	      double rescale=pow(1.0/determinant(relaxed.strain),1.0/3.0);
							
	      relaxed.generate_slat(prim,rescale);
							
	      if(!compare(relaxed.slat,superstruc[sc].struc.slat)){
		relaxlog << superstruc[sc].conf[nc].name << " the relaxed cell is not the same supercell as the original supercell\n";
	      }
	      else{
		arrangement relaxed_conf;
		relaxed.map_on_expanded_prim_basis(prim,relaxed_conf);
		calc_correlations(relaxed, super_basiplet, relaxed_conf);
		if(!compare(relaxed_conf.correlations,superstruc[sc].conf[nc].correlations)){
		  relaxlog << superstruc[sc].conf[nc].name << " relaxed to a new structure\n";
									
		  //----- added by jishnu ----------
		  for(int ncc=0; ncc<superstruc[sc].conf.size(); ncc++){
		    if(compare(relaxed_conf.correlations,superstruc[sc].conf[ncc].correlations)){
		      relaxlog << relaxed_conf.name << " and is a duplicate of " << superstruc[sc].conf[ncc].name << "\n";
		      break;
		    }
		  }
		  //----- added by jishnu ----------			
									
		  superstruc[sc].conf[nc].weight=0.0;
		  for(int i=0; i<superstruc[sc].conf[nc].correlations.size(); i++){
		    superstruc[sc].conf[nc].correlations[i]=relaxed_conf.correlations[i];
		  }
		}
		else{
		  superstruc[sc].conf[nc].weight=1.0;
		}
	      }
	    }
	  }
	}
      }
			
    }
  }
	
  //Read in the names of directories of manually made configurations, collect them and add them to the
  //configs object
	
  // -read in the POSCAR and the CONTCAR, determine slat, and compare to all supercells already found
  // - also compare the slat from the POSCAR with that from the CONTCAR
  // -determine the arrangement and add it to that structure (if it is not already included)
  // -if the arrangement already exists, keep a log of overlapping structures
	
	
	
  string file_name="custom_structures";
  ifstream in_dir;
  in_dir.open(file_name.c_str());
	
  if(!in_dir){
    cout << "cannot open file " << file_name << "\n";
    return;
  }
	
	
  while(!in_dir.eof()){
    string dir_name;
    in_dir >> dir_name;
    if(dir_name.size() > 0){
			
      cout << "WORKING ON " << dir_name << "\n";
			
      //first check whether the original POS can be mapped onto a supercell of PRIM
			
			
      if(scandirectory(dir_name,"POS")){
				
	structure prerelaxed;
				
	string pos_file=dir_name;
	pos_file.append("/POS");
	ifstream in_pos;
	in_pos.open(pos_file.c_str());
	if(!in_pos){
	  cout << "cannot open file " << in_pos << "\n";
	  return;
	}
				
				
	prerelaxed.read_struc_poscar(in_pos);
				
	prerelaxed.generate_slat(prim);
				
	prerelaxed.idealize();
	arrangement prerelaxed_conf;
	prerelaxed_conf.name=dir_name;
	prerelaxed.map_on_expanded_prim_basis(prim, prerelaxed_conf);
				
	prerelaxed_conf.conc.collect_components(prim);
	prerelaxed_conf.conc.calc_concentration(prerelaxed);
				
	//	prerelaxed.map_on_expanded_prim_basis(prim, prerelaxed_conf);
				
				
	multiplet super_basiplet;
				
	prerelaxed.expand_prim_clust(basiplet,super_basiplet);
				
	calc_correlations(prerelaxed, super_basiplet, prerelaxed_conf);
				
	//read the energy from the OSZICAR file if it exists
				
	if(scandirectory(dir_name,"OSZICAR")){
	  double energy;
	  int relax_step; // added by jishnu
	  // if(read_oszicar(dir_name,energy)){
	  if(read_oszicar(dir_name, energy, relax_step)){   // changed by jishnu // used overloaded function to count the no of relax step
	    double vol=determinant(prerelaxed.slat);
	    if(abs(vol) > tol){
	      prerelaxed_conf.energy=energy/abs(vol);
	      prerelaxed_conf.calculated = true;
	      prerelaxed_conf.relax_step = relax_step;   // added by jishnu
	    }
	  }
	}
				
	//if there is a CONTCAR - read it and compare slat and the configuration
				
	if(scandirectory(dir_name,"CONTCAR")){
	  structure relaxed;
	  relaxed.collect_relax(dir_name);
	  double rescale=pow(1.0/determinant(relaxed.strain),1.0/3.0);
	  relaxed.generate_slat(prim,rescale);
					
					
	  if(!compare(relaxed.slat,prerelaxed.slat)){
	    relaxlog << prerelaxed_conf.name << " the relaxed cell is not the same supercell as the original supercell\n";
	  }
	  else{
	    arrangement relaxed_conf;
	    relaxed.map_on_expanded_prim_basis(prim,relaxed_conf);
	    calc_correlations(relaxed, super_basiplet, relaxed_conf);
	    if(!compare(relaxed_conf.correlations,prerelaxed_conf.correlations)){
	      relaxlog << prerelaxed_conf.name << " relaxed to a new structure\n";
	      for(int i=0; i< prerelaxed_conf.correlations.size(); i++){
		prerelaxed_conf.correlations[i]=relaxed_conf.correlations[i];
	      }
	    }
	  }
	}
				
	//go through all current superstructures from configs and see if this supercell is already there
	//if not, add the super cell
	//else see if the conf is already there
	//if not add it, other wise make a note in the log file
	//collect the energy
				
	int non_match_sc=0;
	for(int sc=0; sc<superstruc.size(); sc++){
	  if(compare(prerelaxed.slat,superstruc[sc].struc.slat)){
						
	    //see if the configuration already exists
	    int non_match_nc=0;
	    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
	      if(compare(prerelaxed_conf.correlations,superstruc[sc].conf[nc].correlations)){
		relaxlog << prerelaxed_conf.name << " is a duplicate of " << superstruc[sc].conf[nc].name << "\n";
		//Edited by jishnu
		if((superstruc[sc].conf[nc].calculated)&& (prerelaxed_conf.calculated)){
		  if(abs(superstruc[sc].conf[nc].energy-prerelaxed_conf.energy) > tol) {
		    relaxlog << "CAUTION!!!!  The energy of custom structure " << prerelaxed_conf.name <<" is diffferent from its duplicate structure (" << superstruc[sc].conf[nc].name 
			     << ") by " << prerelaxed_conf.energy-superstruc[sc].conf[nc].energy << " eV.\n"; 
		    cout << "CAUTION!!!!  The energy of custom structure " << prerelaxed_conf.name <<" is diffferent from its duplicate structure (" << superstruc[sc].conf[nc].name 
			 << ") by " << prerelaxed_conf.energy-superstruc[sc].conf[nc].energy << " eV.\n"; 

		  }
		}
		if((!superstruc[sc].conf[nc].calculated) && (prerelaxed_conf.calculated)){
		  superstruc[sc].conf[nc].name=prerelaxed_conf.name;  // The custom structure replaces generated structure 
		  superstruc[sc].conf[nc].energy=prerelaxed_conf.energy;   
		  superstruc[sc].conf[nc].calculated=true;										
		}															
		break;
	      }
	      else{
		non_match_nc++;
	      }
	    }
	    if(non_match_nc == superstruc[sc].conf.size()){
	      cout << "THIS STRUCTURE IS UNIQUE - we are including it \n";
	      cout << "It's name is " << prerelaxed_conf.name << " \n";
	      cout << "\n";
	      prerelaxed_conf.ns=sc;
	      prerelaxed_conf.nc=superstruc[sc].conf.size();
	      superstruc[sc].conf.push_back(prerelaxed_conf);
	      break;
	    }
	  }
	  else{
	    non_match_sc++;
	  }
	}
	if(non_match_sc == superstruc.size()){
	  cout << "THIS STRUCTURE IS UNIQUE - we are including it \n";
	  cout << "It's name is " << prerelaxed_conf.name << " \n";
	  cout << "\n";
	  superstructure tsuperstruc;
	  tsuperstruc.struc=prerelaxed;
	  prerelaxed_conf.ns=superstruc.size();
	  prerelaxed_conf.nc=0;
	  superstruc.push_back(tsuperstruc);
	  superstruc[non_match_sc].conf.push_back(prerelaxed_conf);
	}
				
      }
			
    }
  }
	
	
  in_dir.close();
  relaxlog.close();
	
	
  //read in the weights for all these structures
	
	
	
}

//************************************************************

void configurations::collect_energies_fast(){
	
  //first open a log file with problematic relaxations
  string log_file="relax.log";
  ofstream relaxlog;
  relaxlog.open(log_file.c_str());
  if(!relaxlog){
    cout << "cannot open " << log_file << "\n";
    return;
  }
	
  double tclust_func, tcorr;

  for(int sc=0; sc<superstruc.size(); sc++){
		
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(scandirectory(".",superstruc[sc].conf[nc].name)){
	if(scandirectory(superstruc[sc].conf[nc].name, "OSZICAR")){
	  if(!(superstruc[sc].corr_to_atom_vec.size() && superstruc[sc].corr_to_atom_vec.size())){
	    //If cluster function and basis information do not exist for this supercell, populate the vectors
	    multiplet super_basiplet;
	    superstruc[sc].struc.expand_prim_clust(basiplet,super_basiplet);
	    get_corr_vector(superstruc[sc].struc, super_basiplet, superstruc[sc].corr_to_atom_vec);
	    get_super_basis_vec(superstruc[sc].struc, superstruc[sc].basis_to_bit_vec);
	  }

	  double energy;
	  int relax_step;    // added by jishnu
	  
	  //extract the energy from the OSZICAR file in directory= superstruc[sc].conf[nc].name
	  // if(read_oszicar(superstruc[sc].conf[nc].name, energy)){
	  if(read_oszicar(superstruc[sc].conf[nc].name, energy, relax_step)){   // changed by jishnu // used overloaded function to count the no of relax step
						
	    //normalize the energy by the number of primitive unit cells
	    
	    double vol=determinant(superstruc[sc].struc.slat);
	    if(abs(vol) > tol){
	      superstruc[sc].conf[nc].energy= energy/abs(vol);
	      superstruc[sc].conf[nc].calculated = true;
	      superstruc[sc].conf[nc].relax_step = relax_step;    // added by jishnu
	      
	      //check first if there is a POS file, if not, create it for this configuration				
	      structure relaxed;
	      relaxed.collect_relax(superstruc[sc].conf[nc].name);
	      
	      double rescale=pow(1.0/determinant(relaxed.strain),1.0/3.0);
	      
	      relaxed.generate_slat(prim,rescale);
	      
	      if(!compare(relaxed.slat,superstruc[sc].struc.slat)){
		relaxlog << superstruc[sc].conf[nc].name << " the relaxed cell is not the same supercell as the original supercell\n";
	      }
	      else{
		arrangement relaxed_conf;
		relaxed.map_on_expanded_prim_basis(prim,relaxed_conf);
		
		//Calculate correlations
		relaxed_conf.correlations.push_back(1.0);
		int atom_ind, bit_ind;
		for(int i=0; i<superstruc[sc].corr_to_atom_vec.size(); i++){
		  tcorr=0.0;
		  for(int j=0; j<superstruc[sc].corr_to_atom_vec[i].size(); j++){
		    tclust_func=1.0;
		    for(int k=0; k<superstruc[sc].corr_to_atom_vec[i][j].size(); k++){
		      atom_ind=superstruc[sc].corr_to_atom_vec[i][j][k][0];
		      bit_ind=superstruc[sc].corr_to_atom_vec[i][j][k][1];
		      tclust_func*=superstruc[sc].basis_to_bit_vec[atom_ind][relaxed_conf.bit[atom_ind]][bit_ind];
		    }
		    tcorr+=tclust_func;
		  }
		  relaxed_conf.correlations.push_back(tcorr/superstruc[sc].corr_to_atom_vec[i].size());
		}

		if(!compare(relaxed_conf.correlations,superstruc[sc].conf[nc].correlations)){
		  relaxlog << superstruc[sc].conf[nc].name << " relaxed to a new structure\n";
									
		  //----- added by jishnu ----------
		  for(int ncc=0; ncc<superstruc[sc].conf.size(); ncc++){
		    if(compare(relaxed_conf.correlations,superstruc[sc].conf[ncc].correlations)){
		      relaxlog << relaxed_conf.name << " and is a duplicate of " << superstruc[sc].conf[ncc].name << "\n";
		      break;
		    }
		  }
		  //----- added by jishnu ----------			
									
		  superstruc[sc].conf[nc].weight=0.0;
		  for(int i=0; i<superstruc[sc].conf[nc].correlations.size(); i++){
		    superstruc[sc].conf[nc].correlations[i]=relaxed_conf.correlations[i];
		  }
		}
		else{
		  superstruc[sc].conf[nc].weight=1.0;
		}
	      }
	    }
	  }
	}
      }
			
    }
  }
	
  //Read in the names of directories of manually made configurations, collect them and add them to the
  //configs object
	
  // -read in the POSCAR and the CONTCAR, determine slat, and compare to all supercells already found
  // - also compare the slat from the POSCAR with that from the CONTCAR
  // -determine the arrangement and add it to that structure (if it is not already included)
  // -if the arrangement already exists, keep a log of overlapping structures
	
	
	
  string file_name="custom_structures";
  ifstream in_dir;
  in_dir.open(file_name.c_str());
	
  if(!in_dir){
    cout << "cannot open file " << file_name << "\n";
    return;
  }
	
	
  while(!in_dir.eof()){
    string dir_name;
    in_dir >> dir_name;
    if(dir_name.size() > 0){
			
      cout << "WORKING ON " << dir_name << "\n";
			
      //first check whether the original POS can be mapped onto a supercell of PRIM
			
			
      if(scandirectory(dir_name,"POS")){
				
	structure prerelaxed;
				
	string pos_file=dir_name;
	pos_file.append("/POS");
	ifstream in_pos;
	in_pos.open(pos_file.c_str());
	if(!in_pos){
	  cout << "cannot open file " << in_pos << "\n";
	  return;
	}
				
				
	prerelaxed.read_struc_poscar(in_pos);
	prerelaxed.generate_slat(prim);
	prerelaxed.idealize();
	arrangement prerelaxed_conf;
	prerelaxed_conf.name=dir_name;
	prerelaxed.map_on_expanded_prim_basis(prim, prerelaxed_conf);
				
	prerelaxed_conf.conc.collect_components(prim);
	prerelaxed_conf.conc.calc_concentration(prerelaxed);
				

	//Edited by John - find supercell number first, ensure existence of corr_to_atom_vec and basis_to_bit_vec
	int sc_ind=-1;
	for(int sc=0; sc<superstruc.size(); sc++){
	  if(compare(prerelaxed.slat,superstruc[sc].struc.slat)){
	    sc_ind=sc;
	    break;
	  }
	}
	if(sc_ind==-1){
	  // cout << "New supercell encountered...";
	  sc_ind=superstruc.size();

	  superstructure tsuperstruc;
	  tsuperstruc.struc=prerelaxed;
	  superstruc.push_back(tsuperstruc);
	  // cout << "Added. \n";
	}

	int atom_ind, bit_ind;

	if(!(superstruc[sc_ind].corr_to_atom_vec.size() && superstruc[sc_ind].basis_to_bit_vec.size())){
	  // cout << "Cluster calculation vectors not present for this supercell... ";
	  multiplet super_basiplet;
	  superstruc[sc_ind].struc.expand_prim_clust(basiplet,super_basiplet);
	  get_corr_vector(superstruc[sc_ind].struc, super_basiplet, superstruc[sc_ind].corr_to_atom_vec);
	  get_super_basis_vec(superstruc[sc_ind].struc, superstruc[sc_ind].basis_to_bit_vec);
	  // cout << "Added.\n"
	}

	//Calculate correlations of POS
	prerelaxed_conf.correlations.push_back(1.0);
	for(int i=0; i<superstruc[sc_ind].corr_to_atom_vec.size(); i++){
	  tcorr=0.0;
	  for(int j=0; j<superstruc[sc_ind].corr_to_atom_vec[i].size(); j++){
	    tclust_func=1.0;
	    for(int k=0; k<superstruc[sc_ind].corr_to_atom_vec[i][j].size(); k++){
	      atom_ind=superstruc[sc_ind].corr_to_atom_vec[i][j][k][0];
	      bit_ind=superstruc[sc_ind].corr_to_atom_vec[i][j][k][1];
	      tclust_func*=superstruc[sc_ind].basis_to_bit_vec[atom_ind][prerelaxed_conf.bit[atom_ind]][bit_ind];
	    }
	    tcorr+=tclust_func;
	  }
	  prerelaxed_conf.correlations.push_back(tcorr/superstruc[sc_ind].corr_to_atom_vec[i].size());
	}
	  
	//\End edit by John
      
	//read the energy from the OSZICAR file if it exists
				
	if(scandirectory(dir_name,"OSZICAR")){
	  double energy;
	  int relax_step; // added by jishnu
	  // if(read_oszicar(dir_name,energy)){
	  if(read_oszicar(dir_name, energy, relax_step)){   // changed by jishnu // used overloaded function to count the no of relax step
	    double vol=determinant(prerelaxed.slat);
	    if(abs(vol) > tol){
	      prerelaxed_conf.energy=energy/abs(vol);
	      prerelaxed_conf.calculated = true;
	      prerelaxed_conf.relax_step = relax_step;   // added by jishnu
	    }
	  }
	}
				
	//if there is a CONTCAR - read it and compare slat and the configuration
				
	if(scandirectory(dir_name,"CONTCAR")){
	  structure relaxed;
	  relaxed.collect_relax(dir_name);
	  double rescale=pow(1.0/determinant(relaxed.strain),1.0/3.0);
	  relaxed.generate_slat(prim,rescale);
					
					
	  if(!compare(relaxed.slat,prerelaxed.slat)){
	    relaxlog << prerelaxed_conf.name << " the relaxed cell is not the same supercell as the original supercell\n";
	  }
	  else{
	    arrangement relaxed_conf;
	    relaxed.map_on_expanded_prim_basis(prim,relaxed_conf);

	    //Calculate correlations of CONTCAR
	    relaxed_conf.correlations.push_back(1.0);
	    for(int i=0; i<superstruc[sc_ind].corr_to_atom_vec.size(); i++){
	      tcorr=0.0;
	      for(int j=0; j<superstruc[sc_ind].corr_to_atom_vec[i].size(); j++){
		tclust_func=1.0;
		for(int k=0; k<superstruc[sc_ind].corr_to_atom_vec[i][j].size(); k++){
		  atom_ind=superstruc[sc_ind].corr_to_atom_vec[i][j][k][0];
		  bit_ind=superstruc[sc_ind].corr_to_atom_vec[i][j][k][1];
		  tclust_func*=superstruc[sc_ind].basis_to_bit_vec[atom_ind][relaxed_conf.bit[atom_ind]][bit_ind];
		}
		tcorr+=tclust_func;
	      }
	      relaxed_conf.correlations.push_back(tcorr/superstruc[sc_ind].corr_to_atom_vec[i].size());
	    }


	    if(!compare(relaxed_conf.correlations,prerelaxed_conf.correlations)){
	      relaxlog << prerelaxed_conf.name << " relaxed to a new structure\n";
	      for(int i=0; i< prerelaxed_conf.correlations.size(); i++){
		prerelaxed_conf.correlations[i]=relaxed_conf.correlations[i];
	      }
	    }

	  }
	}
				
	//go through all current superstructures and configs and see if this
	//configuration is already there.  If not, add the configuration.
	//other wise make a note in the log file
	//collect the energy
				
	//see if the configuration already exists
	bool new_flag=true;
	for(int ns=0; ns<superstruc.size(); ns++){
	  for(int nc=0; nc<superstruc[ns].conf.size(); nc++){
	    if(compare(prerelaxed_conf.correlations,superstruc[ns].conf[nc].correlations)){
	      new_flag=false;
	      relaxlog << prerelaxed_conf.name << " is a duplicate of " << superstruc[ns].conf[nc].name << "\n";
	      //Edited by jishnu
	      if((superstruc[ns].conf[nc].calculated)&& (prerelaxed_conf.calculated)){
		if(abs(superstruc[ns].conf[nc].energy-prerelaxed_conf.energy) > tol) {
		  relaxlog << "CAUTION!!!!  The energy of custom structure (" << prerelaxed_conf.name 
			   <<") is very diffferent from that of the generated structure(" 
			   << superstruc[ns].conf[nc].name << ")\n";
		  relaxlog << "   Energy of " << prerelaxed_conf.name << " differs by " 
			   << prerelaxed_conf.energy-superstruc[ns].conf[nc].energy
			   << " eV.\n";
		  cout << "CAUTION!!!!  The energy custom structure (" << prerelaxed_conf.name 
		       <<") is very diffferent from that of the generated structure(" 
		       << superstruc[ns].conf[nc].name << ")\n";
		  cout << "   Energy of " << prerelaxed_conf.name << " differs by " 
		       << prerelaxed_conf.energy-superstruc[ns].conf[nc].energy
		       << " eV.\n";
		}
	      }
	      if((!superstruc[ns].conf[nc].calculated) && (prerelaxed_conf.calculated)){
		superstruc[ns].conf[nc].name=prerelaxed_conf.name;  // The custom structure replaces generated structure 
		superstruc[ns].conf[nc].energy=prerelaxed_conf.energy;   
		superstruc[ns].conf[nc].calculated=true;										
	      }															
	      // break;  //Commented by John.  could cause issues if more than two structures are identical.
	    }
	  }
	}
	if(new_flag){
	  cout << "THIS STRUCTURE IS UNIQUE - we are including it \n";
	  cout << "Its name is " << prerelaxed_conf.name << " \n";
	  cout << "\n";
	  prerelaxed_conf.ns=sc_ind;
	  prerelaxed_conf.nc=superstruc[sc_ind].conf.size();
	  superstruc[sc_ind].conf.push_back(prerelaxed_conf);
	}
      }
			
    }
  }
	
	
  in_dir.close();
  relaxlog.close();
	
	
  //read in the weights for all these structures
	
	
	
}


//************************************************************
//************************************************************

void configurations::calculate_formation_energy(){  // modified by jishnu
	
  facet ref_plane;
  for(int i=0;i<reference.size();i++)	ref_plane.corner.push_back(reference[i]);	
  for(int i=0;i<reference.size();i++) ref_plane.corner[i].update_te();   // very IMPORTANT statement	
  ref_plane.get_norm_vec();
	
  for(int sc=0; sc<superstruc.size();  sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){				
	superstruc[sc].conf[nc].update_te();	// again very IMPORTANT				
	superstruc[sc].conf[nc].got_fenergy = ref_plane.find_endiff(superstruc[sc].conf[nc],superstruc[sc].conf[nc].fpfenergy); // here got_fenergy does not mean anything, whether it is true or false, the formation energy value is always accepted.
	superstruc[sc].conf[nc].got_fenergy = true; // here we are forcing it to be true
	superstruc[sc].conf[nc].update_fp();
      }
    }
  }
}

//************************************************************
void configurations::assemble_hull(){    // Modified by jishnu

  // of all calculated configurations, keep only the lowest energy one for each concentration
  // copy these into an array
  // feed the array to the hull finder
  // keep track of indexing

  vector<arrangement> tconf;
  for(int ns=0; ns<superstruc.size(); ns++){
    for(int nc=0; nc<superstruc[ns].conf.size(); nc++){
      if(superstruc[ns].conf[nc].calculated || ((superstruc[ns].conf[nc].ce == 1) && (superstruc[ns].conf[nc].fp == 0) && (superstruc[ns].conf[nc].te == 0))){                // first condition works for the fp hull and second condition works for the CE hull // jishnu

	//compare this config with the tconf already collected
	//if the concentration is already present, keep the one with the lowest energy
	//otherwise add this point to the list

	int i;
	for(i=0; i<tconf.size(); i++){
	  if(compare(superstruc[ns].conf[nc].conc,tconf[i].conc)){
	    if(superstruc[ns].conf[nc].fenergy < tconf[i].fenergy){
	      tconf[i]=superstruc[ns].conf[nc];
	      tconf[i].assemble_coordinate_fenergy();
	    }
	    break;
	  }
	}
	if(i == tconf.size()){
	  tconf.push_back(superstruc[ns].conf[nc]);
	  tconf[i].assemble_coordinate_fenergy();
	}
      }
    }
  }
  
  if(tconf.size() == 0){
    cout << "No configurations available to determine convex hull\n";
    cout << "quitting assemble_hull() \n";
    return;
  }
  
  
  //determine the number of independent concentration variables
  int dim=0;
  for(int i=0; i<tconf[0].conc.occup.size(); i++){
    for(int j=0; j<tconf[0].conc.occup[i].size()-1; j++){
      dim++;
    }
  }


  if(dim > 2){
    cout << "At this point we can only determine the convex hull for at most a ternary system\n";
    cout << "quitting assemble_hull() \n";
    return;
  }
  
  if (dim == 1) {   // i.e. if the system is binary
  
    double *matrix = new double [tconf.size()*3];  // matrix stores label, conc and FP energy
    for (int i=0;i<tconf.size();i++){
      matrix[i*3+0] = i;  // keeping track of the unique structure details
      matrix[i*3+1] = tconf[i].coordinate[0];
      matrix[i*3+2] = tconf[i].coordinate[1];
		
    }
		
    Array mat(tconf.size(),3);
    mat.setArray(matrix);
    Array gs(3,3);  // first index is no of rows,which can be anything, but second colum is no of columns which must be 3 for binary.
    Array edge(3,2); //again u can put any no instead of 3, but keep 2 fixed.
    gs = mat.hullnd(2,edge); // write the hull points to gs Array, the 1st column is label, 2nd concentration, 3rd E		
				
    gs.assort_a(2);  // this sorts the gs array according to the i th column (starting from 1) // here it is concentration column
	
    int nr = gs.num_row();
    double *gs_all = new double [nr*3];
    for(int i=0;i<nr;i++){
      gs_all[i*3+0] = gs.elem(i+1,1);
      gs_all[i*3+1] = gs.elem(i+1,2);
      gs_all[i*3+2] = gs.elem(i+1,3);
    }  //gs_all is equivalent to matrix (before finding the hull), it contains all the hull points 	(label, concentration and energy in order)			
	
    // finding the left and right end points of the binary hull
    double *left_end = new double [3];
    double *right_end = new double [3];
    for(int i=0;i<3;i++)	left_end[i]=gs_all[0+i];	
    for(int i=0;i<3;i++)	right_end[i]=gs_all[0+i];
										
    for(int i=1;i<nr;i++){
      if (gs_all[i*3+1] < left_end[1]) for(int j=0;j<3;j++)	{left_end[j]=gs_all[i*3+j];}
		
      if (gs_all[i*3+1] > right_end[1]) for(int j=0;j<3;j++)	{right_end[j]=gs_all[i*3+j];}
				
    } // end of finding the ends
	
    // finding the striaght line joining the end points	
    double slope = (right_end[2]-left_end[2])/(right_end[1]-left_end[1]);
    // double intercept = right_end[1] - slope*right_end[0]; // no need to calculate this
    // end of finding the straight line	
	
    // finding the lower half of the hull
    for(int j=0;j<tconf.size();j++){  // first put the left end into the pool of hull points
      if(j == left_end[0]) chull.point.push_back(tconf[j]);					
    }
		
    for (int i=0;i<nr;i++){
      if((gs_all[3*i+1]-left_end[1]) != 0.0){
	double slope_p=(gs_all[3*i+2]-left_end[2])/(gs_all[3*i+1]-left_end[1]);
	if(slope_p <= slope){
	  for(int j=0;j<tconf.size();j++){
	    if(j == gs_all[3*i+0]) chull.point.push_back(tconf[j]);					
	  }
	}
      }			
    }
    // end of finding and saving the lower half of the hull
    // save the facet info
    for (int i=0;i<(chull.point.size()-1);i++){   // "-1" is there because there will be 7 facets(or edges in case of binary) for 8 points
      facet tfacet;
      tfacet.corner.push_back(chull.point[i]);  // the points are already sorted according to the concentration
      tfacet.corner.push_back(chull.point[i+1]);
      chull.face.push_back(tfacet);
    }																															
			
  }  // end of dim == 1 loop
  
  if (dim == 2) {   // i.e. if the system is ternary  
  
    double *matrix = new double [tconf.size()*4];  // matrix stores conc1, conc2, FP energy, and label
    for (int i=0;i<tconf.size();i++){
      matrix[i*4+0] = i;  // keeping track of the unique structure details
      matrix[i*4+1] = tconf[i].coordinate[0];
      matrix[i*4+2] = tconf[i].coordinate[1];
      matrix[i*4+3] = tconf[i].coordinate[2];
		
    }
			
    Array mat(tconf.size(),4);
    mat.setArray(matrix);
    mat.assort_a2(2,3);    //assort data by concentration in ascending order
    Array gs(3,4);  // first index is no of rows,which can be anything, but second colum is no of columns which must be 4 for ternary.
    Array fa(3,3); //again u can put any no instead of first 3, but keep second 3 fixed.
    gs = mat.half_hull3(fa); // write the hull points to gs Array, the 1st column is label, 2nd and 3rd concentration, 4th E // this has already eliminated the upper portion of the hull
    int nr = gs.num_row();
    int nf = fa.num_row();  // number of facets
    //cout << "no of facets = " << nf << "\n";
	
    for (int i=0;i<nr;i++){
      for (int j=0;j<tconf.size();j++){
	if(gs.elem(i+1,1) == j) chull.point.push_back(tconf[j]);						
      }
    }	
		
    // save the facet info	
    for (int i=0;i<nf;i++){
      facet tfacet;
      for(int j=0;j<(dim+1);j++){
	for(int k=0;k<nr;k++){
	  if(fa.elem(i+1,j+1) == gs.elem(k+1,1)) tfacet.corner.push_back(tconf[gs.elem(k+1,1)]);
	}
      }
      chull.face.push_back(tfacet);
    }
	
  }  // end of dim == 2 loop 
  
  
 
}     // end of assemble_hull
//************************************************************
//************************************************************
void configurations::get_delE_from_hull(){     // added by jishnu
	
  for (int i=0;i<chull.face.size();i++){
    chull.face[i].get_norm_vec();
  }	
	
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){			
	for (int i=0;i<chull.face.size();i++){
	  if(chull.face[i].find_endiff(superstruc[sc].conf[nc],superstruc[sc].conf[nc].delE_from_facet)) {  // dont need to do anything becuase the correct one is the last one and that is the only delE_from_facet that will be saved.
	    break;   // there is no need to continue on i loop once u find the right facet // not only "no need" but also, the correct value of delE_from_facet will be overwritten
	  }
	}	
      }
    }
  }
	
	
} // end of get_delE_from_hull 
//************************************************************
//************************************************************
void configurations::get_delE_from_hull_w_clexen(){     // added by jishnu
	
  for (int i=0;i<chull.face.size();i++){
    chull.face[i].get_norm_vec();
  }	
	
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].got_cefenergy){
	int i;
	for (i=0;i<chull.face.size();i++){
	  if(chull.face[i].find_endiff(superstruc[sc].conf[nc],superstruc[sc].conf[nc].delE_from_facet)) {  // dont need to do anything becuase the correct one is the last one and that is the only delE_from_facet that will be saved.
	    break;   // there is no need to continue on i loop once u find the right facet // not only "no need" but also, the correct value of delE_from_facet will be overwritten
	  }					
	}		
      }
    }
  }
	
	
} // end of get_delE_from_clex_hull 
//************************************************************

//************************************************************
void hull::write_clex_hull(){   // added by jishnu
  
  string hull_clex_file="hull.clex";
  ofstream hullclex;
  hullclex.open(hull_clex_file.c_str());
  if(!hullclex){
    cout << "cannot open hull.clex file.\n";
    return;
  }
  hullclex << "# formation_energy        meaningless          concentrations            meaningless            name \n";
  for(int i=0; i<point.size(); i++){
    hullclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << point[i].coordinate[point[i].coordinate.size()-1];
    hullclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;
    for(int j=0;j<(point[i].coordinate.size()-1);j++){
      hullclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  point[i].coordinate[j];
    }
    hullclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;
    hullclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0) << point[i].name << "\n";
  }
  
  hullclex.close();
  
  
  string facet_clex_file="facet.clex";
  ofstream facettclex;
  facettclex.open(facet_clex_file.c_str());
  if(!facettclex){
    cout << "cannot open facet.clex file.\n";
    return;
  }
  
  for(int i=0;i<face.size();i++){
    for(int j=0;j<face[i].corner.size();j++){
      facettclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0)<< face[i].corner[j].name;
    }
    facettclex << "\n";
  }
  
  facettclex.close();
  
  
}
//************************************************************
//************************************************************

void configurations::cluster_expanded_energy(){    // modified by jishnu
  // cout << "inside cluster_expanded_energy \n";
  //basiplet.get_index();  //commented by jishnu
  string filename_eciout = "eci.out";  // added by jishnu
  ifstream eciout; // added by jishnu
  eciout.open(filename_eciout.c_str()); // added by jishnu  
  basiplet.read_eci(eciout);  // added by jishnu  
  
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      //test that the number of correlations matches the number of clusters in basiplet
      if(basiplet.size.size() != superstruc[sc].conf[nc].correlations.size()){
	cout << "Cannot calculate cluster expanded energy since the correlations and clusters are not compatible\n";
	return;
      }
      superstruc[sc].conf[nc].cefenergy=0.0;
      for(int i=0; i<basiplet.size.size(); i++){
	int s=basiplet.size[i];
	int o=basiplet.order[i];
	superstruc[sc].conf[nc].cefenergy=superstruc[sc].conf[nc].cefenergy+basiplet.orb[s][o].equiv.size()*basiplet.orb[s][o].eci*superstruc[sc].conf[nc].correlations[i];
      }
      superstruc[sc].conf[nc].got_cefenergy = true;
    }
  }      
  return;
}

//************************************************************
// ***********************************************************
void configurations::CEfenergy_analysis(){  // added by jishnu

  cluster_expanded_energy();
  int dim = chull.face[0].corner.size()-1;
  
  int ssc,nnc;	
  for(int i=0;i<chull.point.size();i++){  
    bool namematch = false;
    for(int sc=0; sc<superstruc.size(); sc++){
      for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
	if(chull.point[i].name == superstruc[sc].conf[nc].name) {namematch = true;ssc = sc; nnc = nc; break;}
      }
      if(namematch == true) break;
    }    		
    chull.point[i].fpfenergy = superstruc[ssc].conf[nnc].fpfenergy;
    chull.point[i].cefenergy = superstruc[ssc].conf[nnc].cefenergy;     
    chull.point[i].update_ce();     
  }
  
  for(int i=0;i<chull.face.size(); i++){
    for(int j=0;j<chull.face[i].corner.size();j++){ 
      bool namematch = false;
      for(int sc=0; sc<superstruc.size(); sc++){
	for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
	  if(chull.face[i].corner[j].name == superstruc[sc].conf[nc].name) {namematch = true;ssc = sc; nnc = nc; break;}
	}
	if(namematch == true) break;
      }       
    
      chull.face[i].corner[j].fpfenergy = superstruc[ssc].conf[nnc].fpfenergy;
      chull.face[i].corner[j].cefenergy = superstruc[ssc].conf[nnc].cefenergy;
      chull.face[i].corner[j].update_ce();
    }
  }

  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      superstruc[sc].conf[nc].update_ce();
    }
  }
	
  // ------------- this writes the vasp hull with clex energy--------
  string FP_hull_clex_en_file="FPhull.clex";
  ofstream FPhullclexen;
  FPhullclexen.open(FP_hull_clex_en_file.c_str());
  if(!FPhullclexen){
    cout << "cannot open FPhull.clex file.\n";
    return;
  }
  FPhullclexen << "# clex_form_en        FP_form_en         concentrations             meaningless           name \n";
  for(int i=0; i<chull.point.size(); i++){
    FPhullclexen <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << chull.point[i].coordinate[chull.point[i].coordinate.size()-1];
    FPhullclexen <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << chull.point[i].fpfenergy;
    for(int j=0;j<(chull.point[i].coordinate.size()-1);j++){
      FPhullclexen <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  chull.point[i].coordinate[j];
    }
    FPhullclexen <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;
    FPhullclexen <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0) << chull.point[i].name << "\n";
  }
	
  FPhullclexen.close();
	
  //----------------------------------------------------------------	

  get_delE_from_hull_w_clexen();
	
  string energyclex = "energy.clex";
  ofstream enclex;
  enclex.open(energyclex.c_str());
  if(!enclex){
    cout << "cannot open energy.clex file \n";
    return;
  }

  string belowhull_file = "below.hull";
  ofstream belowhull;
  belowhull.open(belowhull_file.c_str());
  if(!belowhull){
    cout << "cannot open below.hull file \n";
    return;
  }

  enclex << "formation energy           calculated/not          concentrations            dist_from_hull            name \n";
  
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){ 
      enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].cefenergy;
      if(superstruc[sc].conf[nc].calculated)enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 1;   
      else enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;
      for(int j=0;j<dim;j++){
	enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  superstruc[sc].conf[nc].coordinate[j];
      }
      enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].delE_from_facet;
      enclex <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].name;
      enclex << "\n";
    }
  }
  enclex.close();
  
  belowhull << "formation energy           calculated/not          concentrations            dist_from_hull            name \n";
  int number_belowhull = 0;    
      
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].delE_from_facet < (-tol)){
	if (superstruc[sc].conf[nc].calculated) {  
	  number_belowhull++;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].cefenergy;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 1;   
	  for(int j=0;j<dim;j++){
	    belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  superstruc[sc].conf[nc].coordinate[j];
	  }
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].delE_from_facet;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].name;
	  belowhull << "\n";
        }
      }   
    }
  }
  
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].delE_from_facet < (-tol)){
	if (!superstruc[sc].conf[nc].calculated) {  
	  number_belowhull++;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].cefenergy;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;   
	  for(int j=0;j<dim;j++){
	    belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  superstruc[sc].conf[nc].coordinate[j];
	  }
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].delE_from_facet;
	  belowhull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << superstruc[sc].conf[nc].name;
	  belowhull << "\n";
        }
      }   
    }
  }
  
  belowhull << "Total No of below Hull points = "  << number_belowhull << "\n";
  belowhull.close();

  chull.clear_arrays(); 
  assemble_hull();
  chull.write_clex_hull();

}
//************************************************************

void configurations::print_eci_inputfiles_old(){
	
  //first determine how many calculated structures there are
  int num_calc=0;
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){
	num_calc++;
      }
    }
  }
	
  //determine how many basis clusters there are
  basiplet.get_hierarchy();
	
  //now print out the corr.in and ener.in files
  string corr_file="corr.in";
  string ener_master_file="ener_master.in";
  string ener_file="ener.in";
  string energyold_file="energy_old";
	
  ofstream corr;
  corr.open(corr_file.c_str());
  if(!corr){
    cout << "cannot open corr.in file.\n";
    return;
  }
	
	
  ofstream ener;
  ener.open(ener_file.c_str());
  if(!ener){
    cout << "cannot open ener.in file.\n";
    return;
  }
	
  ofstream ener_master;
  ener_master.open(ener_master_file.c_str());
  if(!ener_master){
    cout << "cannot open ener_master.in file.\n";
    return;
  }
	
  ofstream energyold;
  energyold.open(energyold_file.c_str());
  if(!energyold){
    cout << "cannot open energy_old file.\n";
    return;
  }
	
	
  corr << basiplet.size.size() << " # number of clusters\n";
  corr << num_calc << " # number of configurations\n";
  corr << "clusters \n";
	
  ener_master << "       exact_ener    weight   structure name  \n";
	
  ener << "       exact_ener    weight   structure name  \n";
	
  energyold << "#   form_energy			weight			concentrations			dist_from_hull          name\n";
	
  num_calc=0;
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){
	num_calc++;
	superstruc[sc].conf[nc].print_correlations(corr);
	ener_master << num_calc << "  " << superstruc[sc].conf[nc].fenergy << "   " << superstruc[sc].conf[nc].weight 
		    << "    " << superstruc[sc].conf[nc].name << "\n";
	ener << num_calc << "  " << superstruc[sc].conf[nc].fenergy << "   " << superstruc[sc].conf[nc].weight 
	     << "    " << superstruc[sc].conf[nc].name << "\n";
	superstruc[sc].conf[nc].print_in_energy_file(energyold); 				
      }
    }
  }
	
  corr.close();
  ener_master.close();
  ener.close();
  energyold.close();
	
	
}




//************************************************************

void configurations::print_eci_inputfiles(){   //changed by jishnu (ener.in and energy files are in energy now) (corr.in is in energy.corr now)
	
  //first determine how many calculated structures there are
  int num_calc=0;
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){
	num_calc++;
      }
    }
  }
	
  //determine how many basis clusters there are
  basiplet.get_hierarchy();
	
  //now print out the energy.corr and enery files
  string corr_file="corr.in";
  string energy_file="energy";
	
  ofstream corr;
  corr.open(corr_file.c_str());
  if(!corr){
    cout << "cannot open corr.in file.\n";
    return;
  }  
	
  ofstream energy;
  energy.open(energy_file.c_str());
  if(!energy){
    cout << "cannot open energy file.\n";
    return;
  }
	
	
  corr << basiplet.size.size() << " # number of clusters\n";
  corr << num_calc << " # number of configurations\n";
  corr << "clusters \n";
	
  energy << "# formation energy   weight        concentrations         dist_from_hull        name \n";
  get_delE_from_hull();
  num_calc=0;
  for(int sc=0; sc<superstruc.size(); sc++){
    for(int nc=0; nc<superstruc[sc].conf.size(); nc++){
      if(superstruc[sc].conf[nc].calculated){
	num_calc++;
	// corr << superstruc[sc].conf[nc].name;
	superstruc[sc].conf[nc].print_correlations(corr);
	//energy << num_calc << "      ";
	//energy << superstruc[sc].conf[nc].relax_step << "        ";
	//energy << superstruc[sc].conf[nc].weight << "       ";
	superstruc[sc].conf[nc].print_in_energy_file(energy);  	
      }
    }
  }
	
  corr.close();
  energy.close();
	
}



//************************************************************

void configurations::assemble_coordinate_fenergy(){
  for(int ns=0; ns<superstruc.size(); ns++){
    for(int nc=0; nc<superstruc[ns].conf.size(); nc++){
      superstruc[ns].conf[nc].assemble_coordinate_fenergy();
    }
  }
}




//************************************************************



//************************************************************
void facet::get_norm_vec(){  // added by jishnu  // finding a-b-c of ax+by+cz+d=0; // d is the offset
	
  normal_vec.clear();
  int num_row = corner.size() -1 ;
  int num_col = corner.size() -1 ;
  double vec_mag = 0.0;
  offset = 0.0;
  for(int i=0; i<corner.size(); i++){
    Array det_mat(num_row,num_col);
    double *matrix = new double [num_row*num_col];
    int mat_ind =0;
		
    // begin filling matrix
    for(int j=1;j<corner.size();j++){
      for(int k=0; k<corner[j].coordinate.size();k++){
	if(k!=i){
	  matrix[mat_ind]=corner[j].coordinate[k]-corner[0].coordinate[k];
	  mat_ind++;
	}
      }
    }
    det_mat.setArray(matrix);
    // matrix filled; find determinant
		
    double det_val=det_mat.det();
    det_val*=pow(-1.0,i+2.0);
    normal_vec.push_back(det_val);
    offset-=corner[0].coordinate[i]*det_val;
    vec_mag+=det_val*det_val;
  }
	
  vec_mag=pow(vec_mag,0.5);
  for(int i=0;i<normal_vec.size();i++){
    normal_vec[i]/=vec_mag;
  }
	
  offset/=vec_mag;
	
	
} // end of get_norm_vec


//************************************************************
//************************************************************
bool facet::find_endiff(arrangement arr, double &delE_from_facet) { // to find which facet contains which structure and to find the corresponding energy on the facet // added by jishnu
	
  vector <double> phase_frac;
  int dim_whole = corner.size();
  int dim = dim_whole - 1;
  vector <double> trow;
  vector <vector <double> > coord;
  vector <double> con_vec;
	
  for (int i=0;i<dim;i++){
    for (int j=0;j<dim;j++){
      trow.push_back(corner[i].coordinate[j]-corner[dim].coordinate[j]);
    }
    coord.push_back(trow);
    trow.clear();
  }
	
  for (int i=0;i<dim;i++){
    con_vec.push_back(arr.coordinate[i]-corner[dim].coordinate[i]);
  }
	
  double *ccoord = new double [dim*dim];	
  for (int i=0;i<dim;i++){
    for (int j=0;j<dim;j++){
      ccoord[i*dim+j] = coord[i][j];				
    }
  }	
	
  double *ccon_vec = new double [dim];
  for (int i=0;i<dim;i++){
    ccon_vec[i] = con_vec[i];
  }		
	
  Array cccoord(dim,dim);
  cccoord.setArray(ccoord);	
  Array tr_cccoord = cccoord.transpose();
  Array inv_tr_cccoord = tr_cccoord.inverse();
	
  Array cccon_vec(dim);
  cccon_vec.setArray(ccon_vec);
  Array tr_cccon_vec = cccon_vec.transpose();
	
  Array pphase_frac = inv_tr_cccoord * tr_cccon_vec; // this is a column vector
	
  double sum = 0.0;
	
  for (int i=0;i<dim;i++){
    phase_frac.push_back(pphase_frac.elem(i+1,1));
    sum = sum + pphase_frac.elem(i+1,1);
  }
  phase_frac.push_back(1.0-sum);
	
  sum = 0.0;
  for (int i=0;i<dim;i++){
    sum = sum + normal_vec[i]*arr.coordinate[i];
  }
	
  // en_facet = ( -d - ax - by )/c;
  double en_facet = (-offset - sum)/normal_vec[dim];
  delE_from_facet = arr.fenergy - en_facet;
  /*// norm_dist_from_facet = mod(ax1+by1+cz1+d)/sqrt(a^2+b^2+c^2);
    double sum1 = 0.0;
    for(int i=0;i<dim_whole;i++){
    sum1+=normal_vec[i]*arr.coordinate[i];
    }
    sum1+=offset;
    double sum2 =0.0;
    for(int i=0;i<dim_whole;i++){
    sum1+=normal_vec[i]*normal_vec[i];
    }
    norm_dist_from_facet = fabs(sum1)/sqrt(sum2);	*/
	
	
  for(int i=0;i<dim_whole;i++){		
    if((phase_frac[i] < 0.0) || (phase_frac[i] > 1.0)) return(false);
  }	
  return (true);	
	
}  // end of find_endiff 
//************************************************************

//************************************************************
void facet::get_mu(){   // added by jishnu
	
  int no_of_comp = corner.size()-1;
  for (int i=0;i<no_of_comp;i++){
    double value = (-offset - normal_vec[i])/normal_vec[corner.size()];
    mu.push_back(value);
  }
	
} // end of get_mu
//************************************************************
//************************************************************
void hull::sort_conc(){
	
  if(point.size() == 0){
    cout << "no points in hull object \n";
    cout << "quitting sort_conc \n";
    return;
  }
	
  //sort the last column first and go left in the coordinate vector of point[]
	
  for(int c=point[0].coordinate.size()-2; c >= 0; c--){
    for(int i=0; i<point.size(); i++){
      for(int j=i+1; j<point.size(); j++){
	if(point[i].coordinate[c] > point[j].coordinate[c]){
	  arrangement tarrange = point[j];
	  point[j]=point[i];
	  point[i]=tarrange;
	}
      }
    }
  }
	
	
}
//************************************************************
void hull::write_hull(){
	
  string hull_file="hull";	
  ofstream hull;
  hull.open(hull_file.c_str());
  if(!hull){
    cout << "cannot open hull file.\n";
    return;
  }
  hull << "# formation_energy        meaningless          concentrations            meaningless            name \n";
  for(int i=0; i<point.size(); i++){
    hull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << point[i].coordinate[point[i].coordinate.size()-1];
    hull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) << 0;
    for(int j=0;j<(point[i].coordinate.size()-1);j++){
      hull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(20)<<setprecision(9) <<  point[i].coordinate[j];
    }
    hull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0) << 0;
    hull <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0) << point[i].name << "\n";
  }
	
  hull.close();
	
	
  string facet_file="facet";	
  ofstream facett;
  facett.open(facet_file.c_str());
  if(!facett){
    cout << "cannot open facet file.\n";
    return;
  }
	
  for(int i=0;i<face.size();i++){
    // cout << face[i].normal_vec[0] << "   "<< face[i].normal_vec[1]<< "     "<< face[i].offset << "  \n";
    for(int j=0;j<face[i].corner.size();j++){	
      facett <<setiosflags(ios::right)<<setiosflags(ios::fixed)<<setw(15)<<setprecision(0)<< face[i].corner[j].name;
    }
    facett << "\n";
  }
	
  facett.close();
	
	
}
//************************************************************

//************************************************************

void chempot::initialize(concentration conc){
  for(int i=0; i<conc.compon.size(); i++){
    vector<double> tm;
    vector<specie> tcompon;
    for(int j=0; j<conc.compon[i].size(); j++){
      tm.push_back(0);
      tcompon.push_back(conc.compon[i][j]);
    }
    m.push_back(tm);
    compon.push_back(tcompon);
  }
}
//************************************************************

void chempot::initialize(vector<vector< specie > > init_compon){
  for(int i=0; i<init_compon.size(); i++){
    vector<double> tm;
    vector<specie> tcompon;
    for(int j=0; j<init_compon[i].size(); j++){
      tm.push_back(0);
      tcompon.push_back(init_compon[i][j]);
    }
    m.push_back(tm);
    compon.push_back(tcompon);
  }
}



//************************************************************

void chempot::set(facet face){
  int k=0;
  for(int i=0; i<m.size(); i++){
    for(int j=0; j<m[i].size()-1; j++){
      m[i][j]=face.mu[k];
      k++;
    }
    m[i][m[i].size()-1]=0.0;
  }
}


//************************************************************

void chempot::increment(chempot muinc){
  if(m.size() !=muinc.m.size()){
    cout << "Trying to increment chemical potential with wrong dimensioned increment \n";
    return;
  }
  for(int i=0; i<m.size(); i++){
    if(m[i].size() !=muinc.m[i].size()){
      cout << "Trying to increment chemical potential with wrong dimensioned increment \n";
      return;
    }

    for(int j=0; j<m[i].size(); j++){
      m[i][j]=m[i][j]+muinc.m[i][j];
    }
  }

}


//************************************************************

void chempot::print(ostream &stream){
  for(int i=0; i<m.size(); i++){
    for(int j=0; j<m[i].size(); j++){
      stream << m[i][j] << "  ";
    }
  }
}


//************************************************************

void chempot::print_compon(ostream &stream){
  for(int i=0; i<compon.size(); i++){
    for(int j=0; j<compon[i].size(); j++){
      compon[i][j].print(stream);
      stream << "  ";
    }
  }
}




//************************************************************

void trajectory::initialize(concentration conc){
  Rx.clear();
  Ry.clear();
  Rz.clear();
  R2.clear();
  spin.clear();
  for(int i=0; i<conc.compon.size(); i++){
    vector<double> tR;
    vector<int> tspin;
    vector<specie> telements;
    for(int j=0; j<conc.compon[i].size(); j++){
      tR.push_back(0.0);
      tspin.push_back(conc.compon[i][j].spin);
      telements.push_back(conc.compon[i][j]);
    }
    Rx.push_back(tR);
    Ry.push_back(tR);
    Rz.push_back(tR);
    R2.push_back(tR);
    spin.push_back(tspin);
    elements.push_back(telements);
  }
}


//************************************************************

void trajectory::set_zero(){
  for(int i=0; i<Rx.size(); i++){
    for(int j=0; j<Rx[i].size(); j++){
      Rx[i][j]=0.0;
      Ry[i][j]=0.0;
      Rz[i][j]=0.0;
      R2[i][j]=0.0;
    }
  }


}


//************************************************************

void trajectory::increment(trajectory R){

  if(Rx.size() != R.Rx.size()){
    cout << "incompatibility in trajectory incrementer \n";
    return;
  }

  for(int i=0; i<Rx.size(); i++){

    if(Rx[i].size() != R.Rx[i].size()){
      cout << "incompatibility in trajectory incrementer \n";
      return;
    }

    for(int j=0; j<Rx[i].size(); j++){
      Rx[i][j]=Rx[i][j]+R.Rx[i][j];
      Ry[i][j]=Ry[i][j]+R.Ry[i][j];
      Rz[i][j]=Rz[i][j]+R.Rz[i][j];
      R2[i][j]=R2[i][j]+R.R2[i][j];
    }
  }
}


//************************************************************

void trajectory::normalize(double D){
  for(int i=0; i<spin.size(); i++){
    for(int j=0; j<spin[i].size(); j++){
      R2[i][j]=R2[i][j]/D;
    }
  }
}


//************************************************************

void trajectory::normalize(concentration conc){
  for(int i=0; i<spin.size(); i++){
    for(int j=0; j<spin[i].size(); j++){
      if(abs(conc.occup[i][j]) > tol){
	R2[i][j]=R2[i][j]/conc.occup[i][j];
      }
    }
  }
}


//************************************************************

void trajectory::print(ostream &stream){
  for(int i=0; i<R2.size(); i++){
    for(int j=0; j<R2[i].size(); j++){
      stream << R2[i][j] << "  ";
    }
  }
}



//************************************************************

void trajectory::print_elements(ostream &stream){
  for(int i=0; i<elements.size(); i++){
    for(int j=0; j<elements[i].size(); j++){
      elements[i][j].print(stream);
      stream << "  ";
    }
  }
}



//************************************************************

Monte_Carlo::Monte_Carlo(structure in_prim, structure in_struc, multiplet in_basiplet, int idim, int jdim, int kdim){

  prim=in_prim;
  basiplet=in_basiplet;
  di=idim; dj=jdim; dk=kdim;

  prim.get_trans_mat();
  prim.update_lat();
  prim.update_struc();

  //check whether the monte carlo cell dimensions are compatible with those of init_struc
  //if not, use the new ones suggested by compatible

  int ndi,ndj,ndk;
  if(!compatible(in_struc,ndi,ndj,ndk)){
    di=ndi; dj=ndj; dk=ndk;
    cout << "New Monte Carlo cell dimensions have been chosen to make the cell\n";
    cout << "commensurate with the initial configuration.\n";
    cout << "The dimensions now are " << di << " " << dj << " " << dk << "\n";
    cout << "\n";
  }

  nuc=di*dj*dk;

  si=6*di; sj=6*dj; sk=6*dk;
  arrayi = new int [2*di];
  arrayj = new int [2*dj];
  arrayk = new int [2*dk];
  for(int i=0; i<di; ++i){arrayi[i] = i; arrayi[i+di] = i;}
  for(int j=0; j<dj; ++j){arrayj[j] = j; arrayj[j+dj] = j;}
  for(int k=0; k<dk; ++k){arrayk[k] = k; arrayk[k+dk] = k;}

  idum=time(NULL);

  //collect basis sites and determine # of basis sites bd
  collect_basis();
  db=basis.size();



  ind1=di*dj*dk; ind2=dj*dk; ind3=dk;
  nmcL=di*dj*dk*db;
  mcL= new int[nmcL];
  ltoi = new int[nmcL];
  ltoj = new int[nmcL];
  ltok = new int[nmcL];
  ltob = new int[nmcL];

  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
	for(int b=0; b<db; b++){
	  int l=index(i,j,k,b);
	  ltoi[l]=i;
	  ltoj[l]=j;
	  ltok[l]=k;
	  ltob[l]=b;
	}
      }
    }
  }

  //initialize the concentration variables
  conc.collect_components(prim);
  collect_sublat();
  num_atoms=conc;
  num_hops=conc;
  AVconc=conc;
  AVsublat_conc=sublat_conc;
  AVnum_atoms=num_atoms;
  assemble_conc_basis_links();

  //initialize the susceptibility and thermodynamic factor variables
  Susc.initialize(conc);
  AVSusc.initialize(conc);

  cout << "Dimensions of Susc are = " << Susc.f.size() << " and " << Susc.f[0].size() << "\n";

  //initialize correlation vector

  for(int i=0; i<basiplet.orb.size(); i++){
    for(int j=0; j<basiplet.orb[i].size(); j++){
      AVcorr.push_back(0.0);
    }
  }

  cout << "AVcorr initialized with " << AVcorr.size() << " elements.\n";

  //initialize the occupation variables in mcL array with the arrangement in init_struc
  initialize(in_struc);

  generate_ext_monteclust(basis, basiplet, montiplet);
  generate_eci_arrays();
}




//************************************************************
void Monte_Carlo::collect_basis(){
  for(int na=0; na<prim.atom.size(); na++){
    if(prim.atom[na].compon.size() >= 2) basis.push_back(prim.atom[na]);
  }
  for(int nb=0; nb<basis.size(); nb++){
    basis[nb].assign_spin();
  }
}


//************************************************************

void Monte_Carlo::assemble_conc_basis_links(){

  //first make the basis to concentration vector basis_to_conc
  for(int i=0; i<basis.size(); i++){
    basis_to_conc.push_back(-1);
    for(int j=0; j<conc.compon.size(); j++){
      if(compare(basis[i].compon,conc.compon[j]))basis_to_conc[i]=j;
    }
    if(basis_to_conc[i] == -1){
      cout << "incompatibility between the basis and the concentration object\n";
      cout << "quitting assemble_conc_links() \n";
      return;
    }
  }

  //next make the conc_to_basis vector
  for(int i=0; i<conc.compon.size(); i++){
    vector<int> tconc_to_basis;
    for(int j=0; j<basis.size(); j++){
      if(compare(basis[j].compon,conc.compon[i])) tconc_to_basis.push_back(j);
    }
    if(tconc_to_basis.size() == 0){
      cout << "incompatibility between the basis and the concentration object\n";
      cout << "quitting assemble_conc_links() \n";
      return;
    }
    conc_to_basis.push_back(tconc_to_basis);
  }
}





//************************************************************
void Monte_Carlo::collect_sublat(){

  //first generate the orbits of non equivalent points
  vector<orbit> points;

  for(int i=0; i< basis.size(); i++){
    //check whether the basis site already exists in an orbit
    bool found = false;
    ////////////////////////////////////////////////////////////////////////////////
    //cout << "points size: " << points.size() << "\n";
    ////////////////////////////////////////////////////////////////////////////////
    for(int np=0; np<points.size(); np++){
      for(int ne=0; ne<points[np].equiv.size(); ne++){
	////////////////////////////////////////////////////////////////////////////////
	//swoboda
	//cout << "\ncompare basis and points\n";
	//for(int x=0; x<3; x++){
	//    cout << "basis[" << i << "]: " << basis[i].fcoord[x] << "\tp[" << np << "]equiv[" << ne << "]: " <<
	//    points[np].equiv[ne].point[0].fcoord[x] << "\n";
	//}
	//cout << "\n";
	////////////////////////////////////////////////////////////////////////////////
	if(compare(basis[i],points[np].equiv[ne].point[0])) found = true;
      }
    }
    if(!found){
      cluster tclust;
      tclust.point.push_back(basis[i]);
      orbit torb;
      torb.equiv.push_back(tclust);
      get_equiv(torb,prim.factor_group);
      points.push_back(torb);
    }
  }


  cout << " THE NUMBER OF DISTINCT POINT CLUSTERS ARE \n";
  cout << points.size() << "\n";



  for(int i=0; i<points.size(); i++){
    sublat_conc.compon.push_back(points[i].equiv[0].point[0].compon);
    vector< double> toccup;
    for(int j=0; j<points[i].equiv[0].point[0].compon.size(); j++) toccup.push_back(0.0);
    sublat_conc.occup.push_back(toccup);
    sublat_conc.mu.push_back(toccup);
  }

  // fill the basis_to_sublat vector

  for(int i=0; i<basis.size(); i++){
    bool mapped= false;
    for(int j=0; j<points.size(); j++){
      for(int k=0; k< points[j].equiv.size(); k++){
	if(compare(basis[i],points[j].equiv[k].point[0])){
	  mapped = true;
	  basis_to_sublat.push_back(j);
	}
      }
    }
    if(!mapped){
      cout << " unable to map a basis site to a sublattice site \n";
    }
  }

  // fill the sublat_to_basis double vector

  for(int j=0; j<points.size(); j++){
    vector<int> tbasis;
    for(int i=0; i<basis.size(); i++){
      if(basis_to_sublat[i] == j) tbasis.push_back(i);
    }
    sublat_to_basis.push_back(tbasis);
  }

}





//************************************************************

void Monte_Carlo::update_mu(chempot mu){
  //takes the mu's and puts them in the right spots in basis
  for(int i=0; i<basis.size(); i++){
    basis[i].mu.clear();
    for(int j=0; j<basis[i].compon.size(); j++){
      basis[i].mu.push_back(mu.m[basis_to_conc[i]][j]);
    }
    basis[i].assemble_flip();
  }

}




//************************************************************
void Monte_Carlo::invert_index(int &i, int &j, int &k, int &b, int l){
  k = l % dk;
  j = ((l - k) / dk) % dj;
  i = ((l - k - j*dk) / (dj*dk)) % di ;
  b = (l - k - j*dk - i*dj*dk) / (dk*dj*di) ;
}

//************************************************************

void Monte_Carlo::generate_eci_arrays(){

  //arrays that need to be constructed:
  //     - eci
  //     - multiplicity
  //     - number of points
  //     - shift array
  //     - start_end array (contains initial and final indices for the different bases in all above arrays)

  //first work out the dimensions so that arrays with appropriate lengths can be dimensioned
  //     for each basis site, go through the clusters and find those with non-zero eci - count the number of them
  //     for each basis site, go through all sites and exponents of the clusters and count the number of them

  neci=0;
  ns=0;
  nse=4*db;

  // determine neci and ns

  for(int nm=0; nm<montiplet.size(); nm++){
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
        if(abs(montiplet[nm].orb[np][no].eci) >= 0.000001){
          neci++;
          for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
            for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
              ns=ns+(montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1)*4;
            }
          }
        }
      }
    }
  }
  //allocate memory for s eci, mult, nump and startend ( determined by number of basis sites)

  s= new int[ns];
  eci= new double[neci];
  nums= new int[neci];
  nump= new int[neci];
  mult= new int[neci];
  startend= new int[nse];

  //fill up all the arrays by going through them again

  //first assign the empty eci

  if(montiplet.size()>=1){
    if(montiplet[0].orb.size()>=1){
      if(montiplet[0].orb[0].size()>=1){
        eci_empty=montiplet[0].orb[0][0].eci;
      }
    }
  }
  //cout<<"empty eci is assigned \n";

  int i=0;
  int j=0;
  for(int nm=0; nm<montiplet.size(); nm++){
    startend[nm*4+0]=i;
    startend[nm*4+2]=j;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
        if(fabs(montiplet[nm].orb[np][no].eci) >= 0.000001){
          eci[i]=montiplet[nm].orb[np][no].eci;
          mult[i]=montiplet[nm].orb[np][no].equiv.size();
          nums[i]=montiplet[nm].orb[np][no].equiv[0].point.size();
          if(nums[i] == 0){
            cout << "Serious problem: a cluster has zero points\n";
            cout << "There will be problems when calculating the energy \n";
          }
          nump[i]=0;
          for(int n=0; n<montiplet[nm].orb[np][no].equiv[0].point.size(); n++){
            nump[i]=nump[i]+montiplet[nm].orb[np][no].equiv[0].point[n].bit+1;
          }
          i++;

          for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
            for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
              for(int nn=0; nn<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; nn++){
                for(int nnn=0; nnn<4; nnn++){
                  s[j]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[nnn];
                  j++;
                }
              }
            }
          }

        }
      }
    }
    startend[nm*4+1]=i;
    startend[nm*4+3]=j;
  }
}


//************************************************************

void Monte_Carlo::write_point_energy(ofstream &out){

  // for each basis site, get all points that are accessed

  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  vector< vector< vector< mc_index > > > collect;

  for(int nm=0; nm< montiplet.size(); nm++){
    vector< vector< mc_index > > tcollect;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
	if(abs(montiplet[nm].orb[np][no].eci) >= 0.000001){
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
	    for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
	      vector <mc_index> tempcollect;
	      mc_index tpoint;
	      tpoint.basis_flag=montiplet[nm].orb[np][no].equiv[ne].point[n].basis_flag;
	      for(int i=0; i<4; i++){
		tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
		tpoint.num_specie=montiplet[nm].orb[np][no].equiv[ne].point[n].compon.size();
	      }
	      bool already_present=false;
	      for(int i=0; i<tcollect.size(); i++){
		if(compare(tpoint,tcollect[i][0])){
		  already_present=true;
		  //adding check to ensure size of collect[] is correct
		  if(tcollect[i].size() < (montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1)){
		    for(int j=tcollect[i].size(); j<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; j++){
		      mc_index tempoint;
		      tempoint=tpoint;
		      string j_num;
		      string p_num;
		      tempoint.name="p";
		      int_to_string(i,p_num,10);
		      tempoint.name.append(p_num);
		      int_to_string(j,j_num,10);
		      tempoint.name.append(j_num);
		      tcollect[i].push_back(tempoint);
		    }
		  }
		  break;
		}
	      }
	      if(already_present == false){
		tpoint.l=tcollect.size();
		string p_num;
		int_to_string(tpoint.l,p_num,10);
		tpoint.name="p";
		tpoint.name.append(p_num);
		for(int j=0; j<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; j++){
		  mc_index tempoint;
		  tempoint=tpoint;
		  string j_num;
		  int_to_string(j,j_num,10);
		  tempoint.name.append(j_num);
		  tempcollect.push_back(tempoint);
		}
		tcollect.push_back(tempcollect);
	      }
	    }
	  }
	}
      }
    }
    collect.push_back(tcollect);
  }

  ////////////////////////////////////////////////////////////////////////////////

  out << "double Monte_Carlo::pointenergy(int i, int j, int k, int b){\n";

  out << "  double energy = 0.0;\n";
  //  if(montiplet.size() > 0){
  //    if(montiplet[0].orb.size() > 0){
  //      if(montiplet[0].orb[0].size() > 0) out << montiplet[0].orb[0][0].eci << ";\n";
  //      else out << "0.0;\n";
  //    }
  //    else out << "0.0;\n";
  //  }
  //  else out << "0.0;\n";



  out << "  int l; " << "\n";

  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  for(int b=0; b<collect.size(); b++){
    out << "  if(b == " << b << "){\n";

    for(int n=0; n<collect[b].size(); n++){
      for(int m=0; m<collect[b][n].size(); m++){
	if(collect[b][n][m].basis_flag=='1'){
	  //for occ basis
	  if(collect[b][n][m].num_specie%2!=0){
	    int num;
	    num=(collect[b][n][m].num_specie-1)/2;
	    out << "     l=index(i";
	    if(collect[b][n][m].shift[0] == 0) out << ",j";
	    if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	    if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	    if(collect[b][n][m].shift[1] == 0) out << ",k";
	    if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	    if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	    if(collect[b][n][m].shift[2] == 0) out << ",";
	    if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	    if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	    out << collect[b][n][m].shift[3] << ");\n";

	    string n_num,m_num;
	    string n_name,t_name;
	    int_to_string(n,n_num,10);
	    n_name="num";
	    t_name="t";
	    n_name.append(n_num);
	    t_name.append(n_num);
	    int_to_string(m,m_num,10);
	    n_name.append(m_num);
	    t_name.append(m_num);

	    //	    out << "     int " << t_name << "=" << m+1 << ";\n";
	    //	    out << "     int " << n_name << "=" << num << ";\n";
	    //                      out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //			n_name << "+" << t_name << "-1)+0.9));\n";
	    out << "     double " << collect[b][n][m].name << "= 0.5*mcL[l]+0.5;\n"; 
	  }
	  else{
	    int num;
	    num=collect[b][n][m].num_specie/2;
	    out << "     l=index(i";
	    if(collect[b][n][m].shift[0] == 0) out << ",j";
	    if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	    if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	    if(collect[b][n][m].shift[1] == 0) out << ",k";
	    if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	    if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	    if(collect[b][n][m].shift[2] == 0) out << ",";
	    if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	    if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	    out << collect[b][n][m].shift[3] << ");\n";

	    string m_num,n_num;
	    string n_name,t_name;
	    int_to_string(n,n_num,10);
	    n_name="num";
	    t_name="t";
	    n_name.append(n_num);
	    t_name.append(n_num);
	    int_to_string(m,m_num,10);
	    n_name.append(m_num);
	    t_name.append(m_num);

	    //	    out << "     int " << t_name << "=" << m+1 << ";\n";
	    //	    out << "     int " << n_name << "=" << num << ";\n";
	    //	    if(m+1<=num) out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //			   n_name << "+" << t_name << "-1)+0.9));\n";
	    //	    else out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //		   n_name << "+" << t_name << ")+.9));\n";
	    out << "     double " << collect[b][n][m].name << "=0.5*mcL[l]+0.5;\n"; 
	  }
	}

	else{
	  //for spin basis
	  out << "     l=index(i";
	  if(collect[b][n][m].shift[0] == 0) out << ",j";
	  if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	  if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	  if(collect[b][n][m].shift[1] == 0) out << ",k";
	  if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	  if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	  if(collect[b][n][m].shift[2] == 0) out << ",";
	  if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	  if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	  out << collect[b][n][m].shift[3] << ");\n";
	  out << "     int " << collect[b][n][m].name << "=mcL[l]";
	  if(m==0) out << ";\n";
	  if(m>0){
	    for(int mm=0; mm<m; mm++){
	      out << "*mcL[l]";
	      if(mm==m-1) out << ";\n";
	    }
	  }
	}
      }
    }

    out << "\n";

    out << "     energy = energy";

    int nm=b;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
	if(fabs(montiplet[nm].orb[np][no].eci) >= 0.000001){
	  if(montiplet[nm].orb[np][no].eci < 0.0) out << montiplet[nm].orb[np][no].eci << "*(";
	  if(montiplet[nm].orb[np][no].eci > 0.0) out << "+" << montiplet[nm].orb[np][no].eci << "*(";
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
	    if(ne > 0) out << "+";
	    for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
	      mc_index tpoint;
	      for(int i=0; i<4; i++){
		tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
	      }
	      for(int i=0; i<collect[b].size(); i++){
		if(compare(tpoint,collect[b][i][0])){
		  int j;
		  j=montiplet[nm].orb[np][no].equiv[ne].point[n].bit;
		  //cout << "b: " << b << "\ti: " << i << "\tj: " << j << "\tname: " << collect[b][i][j].name << "\n";
		  out << collect[b][i][j].name;
		  if(n == montiplet[nm].orb[np][no].equiv[ne].point.size()-1){
		    out << "";
		  }
		  else out << "*";
		  break;
		}
	      }
	    }
	  }
	  out << ")";
	}
      }
    }
    out << ";\n";
    out << "     return energy;\n";
    out << "  }\n";
    out << "\n";
    out << "\n";
  }

  out << "}\n";

}





//************************************************************

void Monte_Carlo::write_normalized_point_energy(ofstream &out){

  // for each basis site, get all points that are accessed

  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  vector< vector< vector< mc_index > > > collect;

  for(int nm=0; nm< montiplet.size(); nm++){
    vector< vector< mc_index > > tcollect;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
	if(abs(montiplet[nm].orb[np][no].eci) >= 0.000001){
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
	    for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
	      vector <mc_index> tempcollect;
	      mc_index tpoint;
	      tpoint.basis_flag=montiplet[nm].orb[np][no].equiv[ne].point[n].basis_flag;
	      for(int i=0; i<4; i++){
		tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
		tpoint.num_specie=montiplet[nm].orb[np][no].equiv[ne].point[n].compon.size();
	      }
	      bool already_present=false;
	      for(int i=0; i<tcollect.size(); i++){
		if(compare(tpoint,tcollect[i][0])){
		  already_present=true;
		  //adding check to ensure size of collect[] is correct
		  if(tcollect[i].size() < (montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1)){
		    for(int j=tcollect[i].size(); j<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; j++){
		      mc_index tempoint;
		      tempoint=tpoint;
		      string j_num;
		      string p_num;
		      tempoint.name="p";
		      int_to_string(i,p_num,10);
		      tempoint.name.append(p_num);
		      int_to_string(j,j_num,10);
		      tempoint.name.append(j_num);
		      tcollect[i].push_back(tempoint);
		    }
		  }
		  break;
		}
	      }
	      if(already_present == false){
		tpoint.l=tcollect.size();
		string p_num;
		int_to_string(tpoint.l,p_num,10);
		tpoint.name="p";
		tpoint.name.append(p_num);
		for(int j=0; j<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; j++){
		  mc_index tempoint;
		  tempoint=tpoint;
		  string j_num;
		  int_to_string(j,j_num,10);
		  tempoint.name.append(j_num);
		  tempcollect.push_back(tempoint);
		}
		tcollect.push_back(tempcollect);
	      }
	    }
	  }
	}
      }
    }
    collect.push_back(tcollect);
  }

  ////////////////////////////////////////////////////////////////////////////////

  out << "double Monte_Carlo::normalized_pointenergy(int i, int j, int k, int b){\n";

  out << "  double energy = 0.0;\n";
  //  if(montiplet.size() > 0){
  //    if(montiplet[0].orb.size() > 0){
  //      if(montiplet[0].orb[0].size() > 0) out << montiplet[0].orb[0][0].eci << ";\n";
  //      else out << "0.0;\n";
  //    }
  //    else out << "0.0;\n";
  //  }
  //  else out << "0.0;\n";



  out << "  int l; " << "\n";

  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  for(int b=0; b<collect.size(); b++){
    out << "  if(b == " << b << "){\n";

    for(int n=0; n<collect[b].size(); n++){
      for(int m=0; m<collect[b][n].size(); m++){
	if(collect[b][n][m].basis_flag=='1'){
	  //for occ basis
	  if(collect[b][n][m].num_specie%2!=0){
	    int num;
	    num=(collect[b][n][m].num_specie-1)/2;
	    out << "     l=index(i";
	    if(collect[b][n][m].shift[0] == 0) out << ",j";
	    if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	    if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	    if(collect[b][n][m].shift[1] == 0) out << ",k";
	    if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	    if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	    if(collect[b][n][m].shift[2] == 0) out << ",";
	    if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	    if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	    out << collect[b][n][m].shift[3] << ");\n";

	    string n_num,m_num;
	    string n_name,t_name;
	    int_to_string(n,n_num,10);
	    n_name="num";
	    t_name="t";
	    n_name.append(n_num);
	    t_name.append(n_num);
	    int_to_string(m,m_num,10);
	    n_name.append(m_num);
	    t_name.append(m_num);

	    //	    out << "     int " << t_name << "=" << m+1 << ";\n";
	    //	    out << "     int " << n_name << "=" << num << ";\n";
	    //                      out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //			n_name << "+" << t_name << "-1)+0.9));\n";
	    out << "     double " << collect[b][n][m].name << "=0.5*mcL[l]+0.5;\n";
	  }
	  else{
	    int num;
	    num=collect[b][n][m].num_specie/2;
	    out << "     l=index(i";
	    if(collect[b][n][m].shift[0] == 0) out << ",j";
	    if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	    if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	    if(collect[b][n][m].shift[1] == 0) out << ",k";
	    if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	    if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	    if(collect[b][n][m].shift[2] == 0) out << ",";
	    if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	    if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	    out << collect[b][n][m].shift[3] << ");\n";

	    string m_num,n_num;
	    string n_name,t_name;
	    int_to_string(n,n_num,10);
	    n_name="num";
	    t_name="t";
	    n_name.append(n_num);
	    t_name.append(n_num);
	    int_to_string(m,m_num,10);
	    n_name.append(m_num);
	    t_name.append(m_num);

	    //	    out << "     int " << t_name << "=" << m+1 << ";\n";
	    //	    out << "     int " << n_name << "=" << num << ";\n";
	    //	    if(m+1<=num) out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //			   n_name << "+" << t_name << "-1)+0.9));\n";
	    //	    else out << "     int " << collect[b][n][m].name << "=floor(1/(abs(mcL[l]-" <<
	    //		   n_name << "+" << t_name << ")+.9));\n";
	    out << "     double " << collect[b][n][m].name << "=0.5*mcL[l]+0.5;\n";
	  }
	}

	else{
	  //for spin basis
	  out << "     l=index(i";
	  if(collect[b][n][m].shift[0] == 0) out << ",j";
	  if(collect[b][n][m].shift[0] < 0) out << collect[b][n][m].shift[0] << ",j";
	  if(collect[b][n][m].shift[0] > 0) out << "+" << collect[b][n][m].shift[0] << ",j";

	  if(collect[b][n][m].shift[1] == 0) out << ",k";
	  if(collect[b][n][m].shift[1] < 0) out << collect[b][n][m].shift[1] << ",k";
	  if(collect[b][n][m].shift[1] > 0) out << "+" << collect[b][n][m].shift[1] << ",k";

	  if(collect[b][n][m].shift[2] == 0) out << ",";
	  if(collect[b][n][m].shift[2] < 0) out << collect[b][n][m].shift[2] << ",";
	  if(collect[b][n][m].shift[2] > 0) out << "+" << collect[b][n][m].shift[2] << ",";

	  out << collect[b][n][m].shift[3] << ");\n";
	  out << "     int " << collect[b][n][m].name << "=mcL[l]";
	  if(m==0) out << ";\n";
	  if(m>0){
	    for(int mm=0; mm<m; mm++){
	      out << "*mcL[l]";
	      if(mm==m-1) out << ";\n";
	    }
	  }
	}
      }
    }

    out << "\n";

    out << "     energy = energy";

    int nm=b;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
	if(fabs(montiplet[nm].orb[np][no].eci) >= 0.000001){
	  double npoints=montiplet[nm].orb[np][no].equiv[0].point.size();
	  if(montiplet[nm].orb[np][no].eci < 0.0) out << montiplet[nm].orb[np][no].eci/npoints << "*(";
	  if(montiplet[nm].orb[np][no].eci > 0.0) out << "+" << montiplet[nm].orb[np][no].eci/npoints << "*(";
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
	    if(ne > 0) out << "+";
	    for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
	      mc_index tpoint;
	      for(int i=0; i<4; i++){
		tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
	      }
	      for(int i=0; i<collect[b].size(); i++){
		if(compare(tpoint,collect[b][i][0])){
		  int j;
		  j=montiplet[nm].orb[np][no].equiv[ne].point[n].bit;
		  //cout << "b: " << b << "\ti: " << i << "\tj: " << j << "\tname: " << collect[b][i][j].name << "\n";
		  out << collect[b][i][j].name;
		  if(n == montiplet[nm].orb[np][no].equiv[ne].point.size()-1){
		    out << "";
		  }
		  else out << "*";
		  break;
		}
	      }
	    }
	  }
	  out << ")";
	}
      }
    }
    out << ";\n";
    out << "     return energy;\n";
    out << "  }\n";
    out << "\n";
    out << "\n";
  }

  out << "}\n";

}


//************************************************************ 

// Routine adapted by John from Monte_Carlo::write_point_energy

void Monte_Carlo::write_point_corr(ofstream &out){


  // for each basis site, get all points that are accessed (points in all clusters that contain the basis site, and also have non-zero eci)

  out << "\n \n//************************************************************ \n \n";

  vector< vector< mc_index > > collect;

  for(int nm=0; nm<montiplet.size(); nm++){
    vector<mc_index> tcollect;
    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
	if(abs(montiplet[nm].orb[np][no].eci) >= 0.000001){
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
            for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){
	      mc_index tpoint;
	      for(int i=0; i<4; i++) tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
	      
	      // check whether tpoint is already among the list in tcollect
	      
	      bool already_present=false;
	      for(int i=0; i<tcollect.size(); i++){
		if(compare(tpoint,tcollect[i])){
		  already_present=true;
		  break;
		}
	      }
	      if(already_present == false){
		tpoint.l=tcollect.size();
		string p_num;
		int_to_string(tpoint.l,p_num,10);
		tpoint.name="p";
		tpoint.name.append(p_num);
		tcollect.push_back(tpoint);
	      }
	    }
          }
	}
      }
    }
    collect.push_back(tcollect);
  }


  // Begin writing Mote_Carlo::pointcorr

  out << "void Monte_Carlo::pointcorr(int i, int j, int k, int b){\n";

  
  out << "  int l; " << "\n";
  
  //Write the variable assignments -> spins of relevant sites are stored in doubles with names of form p0, p1, etc.
  
  for(int b=0; b< collect.size(); b++){
    out << "  if(b == " << b << "){\n";
    for(int n=0; n<collect[b].size(); n++){
      out << "     l=index(i";
      if(collect[b][n].shift[0] == 0) out << ",j";
      if(collect[b][n].shift[0] < 0) out << collect[b][n].shift[0] << ",j";
      if(collect[b][n].shift[0] > 0) out << "+" << collect[b][n].shift[0] << ",j";

      if(collect[b][n].shift[1] == 0) out << ",k";
      if(collect[b][n].shift[1] < 0) out << collect[b][n].shift[1] << ",k";
      if(collect[b][n].shift[1] > 0) out << "+" << collect[b][n].shift[1] << ",k";
      
      if(collect[b][n].shift[2] == 0) out << ",";
      if(collect[b][n].shift[2] < 0) out << collect[b][n].shift[2] << ",";
      if(collect[b][n].shift[2] > 0) out << "+" << collect[b][n].shift[2] << ",";

      out << collect[b][n].shift[3] << "); \n";
      out << "     double " << collect[b][n].name << "=mcL[l]; \n";
    }


    //write out the correlation formulas in terms of the pxx

    out << "\n";


    int nm=b;

    //First add the empty cluster
    out << "     AVcorr[0]+=1.0/" << collect.size() << ";\n";

    for(int np=1; np<montiplet[nm].orb.size(); np++){
      for(int no=0; no<montiplet[nm].orb[np].size(); no++){
        if(fabs(montiplet[nm].orb[np][no].eci) >= 0.000001){

	  int mult_num=0;
	  bool break_flag=false;
	  for(int i=0; i<basiplet.orb.size(); i++){
	    for(int j=0; j<basiplet.orb[i].size(); j++){
	      if(montiplet[nm].index[np][no]==mult_num){
		mult_num=basiplet.orb[i][j].equiv.size();
		break_flag=true;
	      }
	      if(break_flag) break;
	      mult_num++;
	    }
	    if(break_flag) break;
	  }
		

	  out << "     AVcorr[" << montiplet[nm].index[np][no] << "]+=(";
	  
	  for(int ne=0; ne<montiplet[nm].orb[np][no].equiv.size(); ne++){
	    if(ne > 0) out << "+";
	    for(int n=0; n<montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){

	      //find the name of the point
	      mc_index tpoint;
	      for(int i=0; i<4; i++) tpoint.shift[i]=montiplet[nm].orb[np][no].equiv[ne].point[n].shift[i];
	      for(int i=0; i<collect[b].size(); i++){
		if(compare(tpoint,collect[b][i])){
		  for(int j=0; j<montiplet[nm].orb[np][no].equiv[ne].point[n].bit+1; j++){
		    out << collect[b][i].name; 
		    if(n == montiplet[nm].orb[np][no].equiv[ne].point.size()-1 && j == montiplet[nm].orb[np][no].equiv[ne].point[n].bit) out << "";
		    else out << "*";																	   
		  }
		  break;
		}
	      }
	      
	    }	    
	  }
	  out << ")";

	  // Divide by multiplicity and/or size of each cluster
	  if(mult_num>1)
	    out << "/" << mult_num*np;
	  else if(np>1)
	    out << "/" << np;
	  out << ";\n";

	}
      }
    }
    out << "     return;\n";
    out << "  }\n";
    out << "\n";
    out << "\n";
  }
  
  out << "}\n";
  
}


//************************************************************

void Monte_Carlo::write_monte_h(string class_file){
  ifstream in;
  in.open(class_file.c_str());
	
  if(!in){
    cout << "cannot open the " << class_file << " file \n";
    cout << "no monte.h created \n";
    return;
  }
	
  ofstream out;
  out.open("monte.h");
	
  if(!out){
    cout << "unable to create/open monte.h\n";
    cout << "no monte.h created \n";
    return;
  }
	
  bool point_energy_written=false;
  string line;
  while(getline(in,line) && !point_energy_written){
    string check=line.substr(0,31);
    if(check == "double Monte_Carlo::pointenergy"){
      in.close();
      write_point_energy(out);
      out << "\n";
      out << "\n";
      write_normalized_point_energy(out);
      out << "\n";
      out << "\n";
      write_point_corr(out);
      //out << "\n";   // added by jishnu
      //out << "\n";   // added by jishnu
      //write_environment_bool_table(out);   // added by jishnu
      //out << "\n";   // added by jishnu
      //out << "\n";   // added by jishnu
      //write_evaluate_bool(out);   // added by jishnu
      point_energy_written=true;
      out.close();
    }
    else{
      out << line;
      out << "\n";
    }
  }
	
  if(!point_energy_written){
    in.close();
    write_point_energy(out);
    out << "\n";
    out << "\n";
    write_normalized_point_energy(out);
    out << "\n";
    out << "\n";
    write_point_corr(out);
    //out << "\n";   // added by jishnu
    //out << "\n";   // added by jishnu
    //write_environment_bool_table(out);   // added by jishnu
    //out << "\n";   // added by jishnu
    //out << "\n";   // added by jishnu
    //write_evaluate_bool(out);   // added by jishnu
    point_energy_written=true;
    out.close();
  }
	
  in.close();
  return;
	
	
}


//************************************************************
//************************************************************

void Monte_Carlo::write_monte_xyz(ostream &stream){
  //first update the monte carlo structure to have the correct occupancies
  //then print the monte carlo structure
  for(int i=0; i<Monte_Carlo_cell.atom.size(); i++){
    int l=Monte_Carlo_cell.atom[i].bit;
    for(int j=0; j<Monte_Carlo_cell.atom[i].compon.size(); j++){
      if(mcL[l] == Monte_Carlo_cell.atom[i].compon[j].spin) Monte_Carlo_cell.atom[i].occ=Monte_Carlo_cell.atom[i].compon[j];
    }
  }

  Monte_Carlo_cell.write_struc_xyz(stream, conc);

}

//************************************************************
//Added by Aziz
//************************************************************

void Monte_Carlo::write_monte_poscar(ostream &stream){
  //first update the monte carlo structure to have the correct occupancies
  //then print the monte carlo structure
  for(int i=0; i<Monte_Carlo_cell.atom.size(); i++){
    int l=Monte_Carlo_cell.atom[i].bit;
    for(int j=0; j<Monte_Carlo_cell.atom[i].compon.size(); j++){
      if(mcL[l] == Monte_Carlo_cell.atom[i].compon[j].spin) Monte_Carlo_cell.atom[i].occ=Monte_Carlo_cell.atom[i].compon[j];
    }
  }

  Monte_Carlo_cell.write_struc_poscar(stream);

}





//************************************************************

bool Monte_Carlo::compatible(structure struc, int &ndi, int &ndj, int &ndk){

  // checks whether the given MC-dimensions (di,dj,dk) are compatible with the struc unit cell
  // if not, new ones (ndi,ndj,ndk which are the smallest ones just larger than di,dj,dk) are
  // suggested that are compatible


  double mclat[3][3],struclat[3][3],inv_struclat[3][3],strucTOmc[3][3];

  //since struc could correspond to a relaxed structure, its lat[][] may not be integer
  //multiples of the prim.lat[][]
  //therefore, we first find the closest integer multiples of prim.lat[][]
  //and then calculate the struclat[][] as these integer multiples of prim.lat[][]

  struc.generate_slat(prim);   // finds the closest integer multiples of struc.lat in terms of prim.lat and puts it into struc.slat

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      struclat[i][j]=0.0;
      for(int k=0; k<3; k++){
        struclat[i][j]=struclat[i][j]+struc.slat[i][k]*prim.lat[k][j];
      }
    }
  }

  inverse(struclat,inv_struclat);

  ndi=di; ndj=dj; ndk=dk;

  int int_rows;
  do{

    for(int j=0; j<3; j++){
      mclat[0][j]=ndi*prim.lat[0][j];
      mclat[1][j]=ndj*prim.lat[1][j];
      mclat[2][j]=ndk*prim.lat[2][j];
    }

    //Determine the matrix that relates the Monte Carlo cell (mclat) to the struc cell (struclat)

    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        strucTOmc[i][j]=0.0;
	for(int k=0; k<3; k++){
          strucTOmc[i][j]=strucTOmc[i][j]+mclat[i][k]*inv_struclat[k][j];
        }
      }
    }

    //Test whether the rows of strucTOmc[][] are all integers (within numerical noise)
    //if not, increment the mc dimensions corresponding to the non-integer rows

    int_rows=0;
    if(!is_integer(strucTOmc[0])) ndi++;
    else int_rows++;
    if(!is_integer(strucTOmc[1])) ndj++;
    else int_rows++;
    if(!is_integer(strucTOmc[2])) ndk++;
    else int_rows++;

  }while(int_rows !=3);

  if(ndi == di && ndj == dj && ndk == dk) return true;
  else return false;

}

//************************************************************

void Monte_Carlo::initialize(structure init_struc){

  init_struc.generate_slat(prim);
  init_struc.map_on_expanded_prim_basis(prim);


  //create the Monte_Carlo_cell structure

  for(int i=0; i<200; i++) Monte_Carlo_cell.title[i]=prim.title[i];
  Monte_Carlo_cell.scale=prim.scale;
  for(int j=0; j<3; j++){
    Monte_Carlo_cell.lat[0][j]=di*prim.lat[0][j];
    Monte_Carlo_cell.lat[1][j]=dj*prim.lat[1][j];
    Monte_Carlo_cell.lat[2][j]=dk*prim.lat[2][j];
  }

  //determine the relation between the Monte_Carlo_cell.lat[][] and the init_struc.lat[][]
  //this goes into Monte_Carlo_cell.slat[][]
  //then expand the init_struc and fill up the Monte_Carlo_cell with it
  //It is assumed that the Monte_Carlo_cell and init_struc are compatible with each other

  Monte_Carlo_cell.generate_slat(init_struc);
  Monte_Carlo_cell.expand_prim_basis(init_struc);


  //for each site within the Monte_Carlo_cell (with more than two components) determine the shift indices (indicates indices of unit cell and basis)
  //to do that we want the Monte_Carlo_cell coordinates in the primitive cell coordinate system

  Monte_Carlo_cell.generate_slat(prim);
  Monte_Carlo_cell.get_trans_mat();

  ////////////////////////////////////////////////////////////////////////////////
  //swoboda
  //for(int na=0; na<2000; na++){ //Monte_Carlo_cell.atom.size(); na++){
  //    cout << "atom[" << na << "]: " << Monte_Carlo_cell.atom[na].bit << "\t";
  //}
  //cout << "\nsize of Monte Carlo cell: " << Monte_Carlo_cell.atom.size() << "\n";
  ////////////////////////////////////////////////////////////////////////////////

  for(int na=0; na<Monte_Carlo_cell.atom.size(); na++){
    if(Monte_Carlo_cell.atom[na].compon.size() >= 2){
      atompos hatom=Monte_Carlo_cell.atom[na];
      conv_AtoB(Monte_Carlo_cell.StoP,Monte_Carlo_cell.atom[na].fcoord,hatom.fcoord);
      get_shift(hatom,basis);

      //assign the atom[na].bit variable the index within the mcL array for this site

      int i=hatom.shift[0];
      int j=hatom.shift[1];
      int k=hatom.shift[2];
      int b=hatom.shift[3];

      int l=index(i,j,k,b);
      Monte_Carlo_cell.atom[na].bit=l;
      
      mcL[l]=Monte_Carlo_cell.atom[na].occ.spin;
    }
  }
}



//************************************************************

void Monte_Carlo::initialize(concentration in_conc){

  //visit each site of the cell
  //randomly assign a spin with the probability in conc to the Monte Carlo cell

  for(int i=0; i < di; i++){
    for(int j=0; j < dj; j++){
      for(int k=0; k < dk; k++){
	for(int b=0; b<db; b++){
	  //check whether the basis site is regular (as opposed to activated)
	  //SPECIFIC FOR LITIS2
	  if(basis[b].bit == 0){
	    //if(basis[b].basis_flag == '0'){
	    int l=index(i,j,k,b);
	    //get the basis to concentration link to determine which spins can occupy that site
	    int c=basis_to_conc[b];

	    double p=ran0(idum);
	    double sum=0.0;
	    for(int d=0; d<in_conc.compon[c].size(); d++){
	      sum=sum+in_conc.occup[c][d];
	      if(p <= sum){
		mcL[l]=in_conc.compon[c][d].spin;
		break;
	      }
	    }
	  }
	}
      }
    }
  }
  //calculate the new concentration
  calc_concentration();
  calc_sublat_concentration();
  calc_num_atoms();
}

//************************************************************

void Monte_Carlo::initialize_1vac(concentration in_conc){
  ////////////////////////////////////////////////////////////////////////////////
  //added by Ben Swoboda
  //initialize changed so that only 1 site in structure becomes a vacancy.

  int count=0;

  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
	for(int b=0; b<db; b++){
	  //check if basis site is a lattice site
	  if(basis[b].bit==0){
	    int l=index(i,j,k,b);
	    int c=basis_to_conc[b];
	    if(count==0){
	      mcL[l]=in_conc.compon[0][1].spin;
	      count++;
	    }
	    else {
	      mcL[l]=in_conc.compon[0][0].spin;
	    }
	  }
	}
      }
    }
  }

  //calculate the new concentration
  calc_concentration();
  calc_sublat_concentration();
  calc_num_atoms();

}

//************************************************************
//Added by Aziz : Beginning
//************************************************************

void Monte_Carlo::initialize_1_specie(double in_conc){

  //This routine fills the bulk to the specified concentration in_conc  
  //This routine assumes that there is only 1 type of specie in the system
  //DO NOT USE FOR MULTI-SPECIES SYSTEMS 
  
  //visit each site of the cell
  //randomly assign a spin with the probability in conc to the Monte Carlo cell


  for(int i=0; i < di; i++){
    for(int j=0; j < dj; j++){
      for(int k=0; k < dk; k++){
	    for(int b=0; b<db; b++){
	      //check whether the basis site is regular (as opposed to activated)
	      //SPECIFIC FOR LITIS2
	      if(basis[b].bit == 0){
	      //if(basis[b].basis_flag == '0'){
	        int l=index(i,j,k,b);


	        double p=ran0(idum);

            if(p <= in_conc){
		      mcL[l]=1;
	        }
		    else { //Empty site that is not in vacuum
		      mcL[l]=-1; 
	        }
	      }
        }
       }
    }
  }
  //calculate the new concentration
  calc_concentration(); 
  calc_num_atoms();
}



//************************************************************
//Added by Aziz : End
//************************************************************

//************************************************************

////////////////////////////////////////////////////////////////////////////////
//************************************************************

void Monte_Carlo::calc_energy(double &energy){
  energy=0.0;
  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
	energy=energy+eci_empty;
        for(int b=0; b<db; b++){
          double temp=normalized_pointenergy(i,j,k,b);
          energy=energy+temp;
        }
      }
    }
  }
}



//************************************************************

void Monte_Carlo::calc_num_atoms(){
  //determines how many of each component there are in the monte carlo cell

  for(int i=0; i<num_atoms.compon.size(); i++){
    for(int j=0; j<num_atoms.compon[i].size(); j++){
      num_atoms.occup[i][j]=0.0;
    }
  }
  //go through the Monte Carlo lattice
  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
        for(int b=0; b<db; b++){
          //determine the concentration unit
          int c=basis_to_conc[b];
          int l=index(i,j,k,b);
          for(int m=0; m<num_atoms.compon[c].size(); m++){
            if(mcL[l] == num_atoms.compon[c][m].spin) num_atoms.occup[c][m]++;
          }
        }
      }
    }
  }

}

//************************************************************

void Monte_Carlo::calc_concentration(){

  for(int i=0; i<conc.compon.size(); i++){
    for(int j=0; j<conc.compon[i].size(); j++){
      conc.occup[i][j]=0.0;
    }
  }
  //go through the Monte Carlo lattice
  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
        for(int b=0; b<db; b++){
	  //only calculate the concentration over the regular sites
	  //SPECIFIC FOR LITIS2
	  if(basis[b].bit == 0){
	    // if(basis[b].basis_flag == '0'){
	    //determine the concentration unit
	    int c=basis_to_conc[b];
	    int l=index(i,j,k,b);
	    for(int m=0; m<conc.compon[c].size(); m++){
	      if(mcL[l] == conc.compon[c][m].spin) conc.occup[c][m]++;
	    }
	  }
        }
      }
    }
  }


  for(int i=0; i<conc.compon.size(); i++){
    //number of sublattices with the i'th concentration unit
    int n=0;
    for(int ii=0; ii<conc_to_basis[i].size(); ii++)
      //SPECIFIC FOR LITIS2
      // if(basis[conc_to_basis[i][ii]].basis_flag == '0') n++;
      if(basis[conc_to_basis[i][ii]].bit == 0) n++;

    //number of sites in the crystal with the i'th concentration unit
    n=n*di*dj*dk;
    if(n != 0){
      for(int j=0; j<conc.compon[i].size(); j++){
	conc.occup[i][j]=conc.occup[i][j]/n;
      }
    }
  }
}


//************************************************************

void Monte_Carlo::calc_sublat_concentration(){

  for(int i=0; i<sublat_conc.compon.size(); i++){
    for(int j=0; j<sublat_conc.compon[i].size(); j++){
      sublat_conc.occup[i][j]=0.0;
    }
  }
  //go through the Monte Carlo lattice
  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
        for(int b=0; b<db; b++){
	  //only calculate the concentration over the regular sites
	  if(basis[b].bit == 0){
	    //determine the concentration unit
	    int c=basis_to_sublat[b];
	    int l=index(i,j,k,b);
	    for(int m=0; m<sublat_conc.compon[c].size(); m++){
	      if(mcL[l] == sublat_conc.compon[c][m].spin) sublat_conc.occup[c][m]++;
	    }
	  }
        }
      }
    }
  }


  for(int i=0; i<sublat_conc.compon.size(); i++){
    //number of sublattices with the i'th concentration unit
    int n=0;
    for(int ii=0; ii<sublat_to_basis[i].size(); ii++)
      if(basis[sublat_to_basis[i][ii]].bit == 0)n++;

    //number of sites in the crystal with the i'th concentration unit
    n=n*di*dj*dk;
    if(n != 0){
      for(int j=0; j<sublat_conc.compon[i].size(); j++){
	sublat_conc.occup[i][j]=sublat_conc.occup[i][j]/n;
      }
    }
  }
}





//************************************************************

void Monte_Carlo::update_num_hops(int l, int ll, int b, int bb){
  int c=basis_to_conc[b];
  int cc=basis_to_conc[bb];

  for(int i=0; i<num_hops.compon[c].size(); i++){
    if(num_hops.compon[c][i].spin == mcL[l]) num_hops.occup[c][i]++;
  }

  for(int i=0; i<num_hops.compon[cc].size(); i++){
    if(num_hops.compon[cc][i].spin == mcL[ll]) num_hops.occup[cc][i]++;
  }
  return;

}




//************************************************************

double Monte_Carlo::calc_grand_canonical_energy(chempot mu){
  double energy;
  calc_energy(energy);
  calc_num_atoms();

  for(int i=0; i<num_atoms.compon.size(); i++){
    for(int j=0; j<num_atoms.compon[i].size(); j++){
      energy=energy-mu.m[i][j]*num_atoms.occup[i][j];
    }
  }
  return energy;

}



//************************************************************

void Monte_Carlo::grand_canonical(double beta, chempot mu, int n_pass, int n_equil_pass){
  int i,j,k,b,l;
  fluctuation uncorr_susc;
  uncorr_susc.initialize(conc);

  if(n_pass <=n_equil_pass){
    cout << "Npass must be larger than Nequil\n";
    cout << "Quitting grand_canonical()\n";
    exit(1);
  }

  //initialize all the thermodynamic averages to zero
  AVenergy=0.0;
  heatcap=0.0;
  AVconc.set_zero();
  AVnum_atoms.set_zero();
  AVsublat_conc.set_zero();
  AVSusc.set_zero();
  flipfreq=0.0;

  if(corr_flag){
    for(int n=0; n<AVcorr.size(); n++){
      AVcorr[n]=0.0;
    }
  }

  //copy mu into the basis objects and update the flip and dmu
  update_mu(mu);

  double grand_energy=calc_grand_canonical_energy(mu);

  cout << "grand canonical energy = " << grand_energy/nuc << "\n";

  for(int n=0; n<n_pass; n++){
    //pick nmcL lattice sites at random
    for(int nn=0; nn<nmcL; nn++){
      i=int(di*ran0(idum));
      j=int(dj*ran0(idum));
      k=int(dk*ran0(idum));
      b=int(db*ran0(idum));
      l=index(i,j,k,b);


      //determine index of current occupant at site l
      int co=basis[b].iflip(mcL[l]);

      //pick a flip event
      int f=int(ran0(idum)*basis[b].flip[co].size());

      int tspin=mcL[l];

      //get point energy before the spin flip
      double en_before=pointenergy(i,j,k,b)-basis[b].mu[co];

      //get point energy after the spin flip
      mcL[l]=basis[b].flip[co][f];
      int no=basis[b].iflip(mcL[l]);
      double en_after=pointenergy(i,j,k,b)-basis[b].mu[no];

      double delta_energy;
      delta_energy=en_after-en_before;

      if(delta_energy < 0.0 || exp(-delta_energy*beta) >=ran0(idum)){
        flipfreq++;
        grand_energy=grand_energy+delta_energy;
      }
      else{
        mcL[l]=tspin;
      }
    }

    if(n >= n_equil_pass){
      if(corr_flag){
	//Update average correlations
	for(i=0; i<di; i++){
	  for(j=0; j<dj; j++){
	    for(k=0; k<dk; k++){
	      for(b=0; b<db; b++){
		pointcorr(i, j, k, b);
	      }
	    }
	  }
	}
      }
      AVenergy=AVenergy+grand_energy;
      heatcap=heatcap+(grand_energy*grand_energy);
      calc_num_atoms();
      calc_concentration();
      calc_sublat_concentration();
      AVconc.increment(conc);
      AVsublat_conc.increment(sublat_conc);
      AVnum_atoms.increment(num_atoms);
      Susc.evaluate(num_atoms);
      AVSusc.increment(Susc);
    }

  }

  AVenergy=AVenergy/(n_pass-n_equil_pass);
  if(corr_flag){
    for(int n=0; n<AVcorr.size(); n++){
      AVcorr[n]=AVcorr[n]/(di*dj*dk*(n_pass-n_equil_pass));
    }
  }
  heatcap=heatcap/(n_pass-n_equil_pass);
  heatcap=(heatcap-(AVenergy*AVenergy))*(beta*beta)*kb;
  AVconc.normalize(n_pass-n_equil_pass);
  AVsublat_conc.normalize(n_pass-n_equil_pass);
  AVnum_atoms.normalize(n_pass-n_equil_pass);
  uncorr_susc.evaluate(AVnum_atoms);
  AVSusc.normalize(n_pass-n_equil_pass);
  AVSusc.decrement(uncorr_susc);
  AVSusc.normalize(1.0/beta);
  flipfreq=flipfreq/(n_pass*nmcL);

}

//************************************************************

void Monte_Carlo::canonical_1_species(double beta,int Temp, int n_pass, int n_equil_pass, int n_pass_output){


  //This subroutine performs canonical Monte Carlo simulations in the specific case of an A/B system


  int i_1,j_1,k_1,b_1,l_1; //Position indices for site 1
  int i_2,j_2,k_2,b_2,l_2; //Position indices for site 2
  string outstruc_file_title, outenergy_file_title;
  ofstream outstruc_file,outenergy_file;
  string step;
  double energy, energy_intermediate_configuration;

  fluctuation uncorr_susc;
  uncorr_susc.initialize(conc);

  if(n_pass <=n_equil_pass){
    cout << "Npass must be larger than Nequil\n";
    cout << "Quitting canonical_1_species \n";
    exit(1);
  }

  //initialize all the thermodynamic averages to zero
  AVenergy=0.0;
  heatcap=0.0;
  AVconc.set_zero();
  AVnum_atoms.set_zero();
  AVsublat_conc.set_zero();
  AVSusc.set_zero();
  flipfreq=0.0;
  

  if(corr_flag){
    for(int n=0; n<AVcorr.size(); n++){
      AVcorr[n]=0.0;
    }
  }


  //Calculating the total energy 

  calc_energy(energy); 

  //Outputting the value to the screen
  cout << "The energy per unit cell is = " << energy/nuc << "\n";
 
  
  bool right_combination=false; //This variable tells us if a valid combination of sites has been found (so that the concentration is maintained)
  


  for(int n=0; n<n_pass; n++){
    //n_pass passes are made. For each pass, pick nmcL canonical pairs of sites (at random)
    
    if (n >n_equil_pass && n%n_pass_output==0 ) {
        
        outstruc_file_title="structure.";
	    step="";
	    int_to_string(Temp,step,10);
	    outstruc_file_title.append(step);
	    outstruc_file_title.append(".");
	    step="";
	    int_to_string(n,step,10);
	    outstruc_file_title.append(step);
        outstruc_file_title.append(".xyz");
	    outstruc_file.open(outstruc_file_title.c_str());
	    write_monte_xyz(outstruc_file);
	    outstruc_file.close();
	    
        outstruc_file_title="structure.";
	    step="";
	    int_to_string(Temp,step,10);
	    outstruc_file_title.append(step);
	    outstruc_file_title.append(".");
	    step="";
	    int_to_string(n,step,10);
	    outstruc_file_title.append(step);
        outstruc_file_title.append(".POSCAR");
	    outstruc_file.open(outstruc_file_title.c_str());
	    write_monte_poscar(outstruc_file);
	    outstruc_file.close();
	    
        outenergy_file_title="structure.";
	    step="";
	    int_to_string(Temp,step,10);
	    outenergy_file_title.append(step);
	    outenergy_file_title.append(".");
	    step="";
	    int_to_string(n,step,10);
	    outenergy_file_title.append(step);
        outenergy_file_title.append(".energy");
	    outenergy_file.open(outenergy_file_title.c_str());
	    calc_energy(energy_intermediate_configuration);
	    energy_intermediate_configuration=energy_intermediate_configuration/(db*nuc);
	    outenergy_file << energy_intermediate_configuration << " " ;
	    outenergy_file.close();
      }
	
	
    for(int nn=0; nn<nmcL; nn++){

      //Loop until a valid pair of sites is selected
      //The boolean variable "right_combination" tells us if a valid combination has been identified
      

      right_combination = false;

      while (right_combination==false) {
        i_1=int(di*ran0(idum));
        j_1=int(dj*ran0(idum));
        k_1=int(dk*ran0(idum));
        b_1=int(db*ran0(idum));  
        l_1=index(i_1,j_1,k_1,b_1);

        i_2=int(di*ran0(idum));
        j_2=int(dj*ran0(idum));
        k_2=int(dk*ran0(idum));
        b_2=int(db*ran0(idum)); 
        l_2=index(i_2,j_2,k_2,b_2);

        //Testing if a valid combination (electron/hole or lithium-ion/vacancy pair) has been identified
	    if (  ( (mcL[l_2]==1 && mcL[l_1]==-1) || (mcL[l_2]==-1 && mcL[l_1]==1) ) ){
	      right_combination=true;	 
        }
		
      }
	  

	  //Storing the initial spin states 
      int tspin_2=mcL[l_2];  
	  int tspin_1=mcL[l_1];

	  //To calculate the difference in energy due to the 2 flips, we  first perform
	  // the flip on site 2 and then the flip on site 1. We calculate the pointenergy
	  //of each site, before and after the flip.

	  //Site 2 flip
      //get point energy before the flip of site 2
       double en_before_flip_2=pointenergy(i_2,j_2,k_2,b_2);

       //get point energy after the flip of site 2
       mcL[l_2]=-mcL[l_2];
	   double en_after_flip_2=pointenergy(i_2,j_2,k_2,b_2);

      double delta_energy_flip_2;
      delta_energy_flip_2=en_after_flip_2-en_before_flip_2;


	  //Site 1 flip
      //get point energy before the flip of site 1
	  double en_before_flip_1=pointenergy(i_1,j_1,k_1,b_1);

      //Get point energy after the flip of site 1
      mcL[l_1]=-mcL[l_1];
	  double en_after_flip_1=pointenergy(i_1,j_1,k_1,b_1);

      double delta_energy_flip_1;
      delta_energy_flip_1=en_after_flip_1-en_before_flip_1;

	  //Calculating the total change in energy
	  double delta_energy=delta_energy_flip_2+delta_energy_flip_1;

      if(delta_energy < 0.0 || exp(-delta_energy*beta) >=ran0(idum)){
        flipfreq++;
        energy=energy+delta_energy;
      }
      else{
        mcL[l_2]=tspin_2;
	    mcL[l_1]=tspin_1;
      }
    }

	//Calculating the correlations
    //(Note : We have to check if these correlation still apply in the canonical space)

    if(n >= n_equil_pass){
      if(corr_flag){
	//Update average correlations
	for(int i=0; i<di; i++){
	  for(int j=0; j<dj; j++){
	    for(int k=0; k<dk; k++){
	      for(int b=0; b<db; b++){
		pointcorr(i, j, k, b);
	      }
	    }
	  }
	}
     }
      AVenergy=AVenergy+energy;
      heatcap=heatcap+(energy*energy);
      calc_num_atoms();
      calc_concentration();
      calc_sublat_concentration();
      AVconc.increment(conc);
      AVsublat_conc.increment(sublat_conc);
      AVnum_atoms.increment(num_atoms);
      Susc.evaluate(num_atoms);
      AVSusc.increment(Susc);
    }

  }

  AVenergy=AVenergy/(n_pass-n_equil_pass);
  if(corr_flag){
    for(int n=0; n<AVcorr.size(); n++){
      AVcorr[n]=AVcorr[n]/(di*dj*dk*(n_pass-n_equil_pass));
    }
  }
  heatcap=heatcap/(n_pass-n_equil_pass);
  heatcap=(heatcap-(AVenergy*AVenergy))*(beta*beta)*kb;
  AVconc.normalize(n_pass-n_equil_pass);
  AVsublat_conc.normalize(n_pass-n_equil_pass);
  AVnum_atoms.normalize(n_pass-n_equil_pass);
  uncorr_susc.evaluate(AVnum_atoms);
  AVSusc.normalize(n_pass-n_equil_pass);
  AVSusc.decrement(uncorr_susc);
  AVSusc.normalize(1.0/beta);
  flipfreq=flipfreq/(n_pass*nmcL);
  

}


//************************************************************


// Code Edited - function added by John Thomas                                                                                                                                                                                                                                     
double Monte_Carlo::lte(double beta, chempot mu){
	
  int l;
	
  update_mu(mu);
  double phi_expanded=calc_grand_canonical_energy(mu);
	
  for(int b=0; b<db; b++){
    for(int i=0; i<di; i++){
      for(int j=0; j<dj; j++){
	for(int k=0; k<dk; k++){
					
	  l=index(i,j,k,b);
					
	  int co=basis[b].iflip(mcL[l]);
					
	  for(int f=0; f<basis[b].flip[co].size(); f++){
						
	    int tspin=mcL[l];
	    double en_before=pointenergy(i,j,k,b)-basis[b].mu[co];
	    //cout << "en_before = "<< en_before << "\n";
						
	    mcL[l]=basis[b].flip[co][f];
	    int no=basis[b].iflip(mcL[l]);
	    double en_after=pointenergy(i,j,k,b)-basis[b].mu[no];
	    //cout << "en_after = "<< en_after << "\n";
						
	    double delta_energy;
	    delta_energy=en_after-en_before;
	    mcL[l]=tspin;
						
	    if(delta_energy<0){
	      cout << "Configuration is not a ground state at current chemical potential.  Please re-initialize with correct ground state. \n";
	      exit(1); 
	    }
	    phi_expanded=phi_expanded-exp(-delta_energy*beta)/beta;
	  }
	}
      }
    }
  }
	
  phi_expanded=phi_expanded/nuc;
  cout << "Finished Low Temperature Expansion.  Free energy at initial temperature is " << phi_expanded << "\n";
  return phi_expanded;
}
//\end edited code                                                                                                                                                                                                                                                                 

//************************************************************

void Monte_Carlo::initialize_kmc(){
	
  //for each basis site, find the corresponding multiplet of montiplet
  //push back the nearest neighbor clusters into endpoints
	
  char name[2];
  name[0]='V';
  name[1]='a';
	
	
	
  //read in the hop clusters
  //identify final states and activated states
  //generate basiplet and then montiplet for each site
  //collect other information - shift vectors, hop vectors etc.
	
  string hop_file="hops";
	
  ifstream in;
  in.open(hop_file.c_str());
  if(!in){
    cout << hop_file << " cannot be read \n";
    cout << "not initializing kinetic Monte Carlo \n";
    return;
  }
	
  int nhops;
  vector<orbit> hoporb;
  char buff[200];
  in >> nhops;
  in.getline(buff,199);
  for(int i=0; i< nhops; i++){
    cluster tclus;
    int np;
    in >> np >> np;
    in.getline(buff,199);
    if(np != 2 && np != 3){
      cout << "hop cluster in " << hop_file << " is not legitimate \n";
      cout << "use only 2 or 3 point clusters \n";
      cout << "not initializing kinetic Monte Carlo \n";
      return;
    }
		
    tclus.readf(in,np);
    //for each point indicate whether they are part of the regular lattice or are an activated
    //state
    //the bit in atompos is used to indicate whether the point is regular or activated
    // regular: bit = 0 ; activated: bit = 1.
    //if there are only two sites in the hop, then all sites are regular
    //if there are three sites in the hop, then the middle one is the activated
		
    if(tclus.point.size() == 2){
      for(int j=0; j<2; j++){
	tclus.point[j].bit=0;
      }
    }
		
    if(tclus.point.size() == 3){
      for(int j=0; j<3; j++){
	tclus.point[j].bit=0;
      }
      tclus.point[1].bit=1;
    }
		
    orbit torb;
    torb.equiv.push_back(tclus);
    hoporb.push_back(torb);
		
  }
	
  in.close();
	
  //for each point in the hop clusters, we need to compare to basis[] and assign those points the allowed components etc
	
	
  for(int i=0; i< hoporb.size(); i++){
    for(int j=0; j< hoporb[i].equiv.size(); j++){
      for(int k=0; k<hoporb[i].equiv[j].point.size(); k++){
	bool mapped=false;
	for(int l=0; l< basis.size(); l++){
	  int trans[3];
	  if(compare(basis[l].fcoord,hoporb[i].equiv[j].point[k].fcoord,trans)){
	    hoporb[i].equiv[j].point[k].compon.clear();
	    for(int m=0; m < basis[l].compon.size(); m++){
	      hoporb[i].equiv[j].point[k].compon.push_back(basis[l].compon[m]);
	    }
	    mapped=true;
	  }
	}
	if(!mapped){
	  cout << " a point from the hop cluster was not mapped onto \n";
	  cout << "a basis site of your monte carlo system \n";
	}
      }
    }
  }
	
	
	
	
  //for each hop cluster we need to get the orbit using the prim symmetry operations
	
  for(int i=0; i< hoporb.size(); i++){
    get_equiv(hoporb[i],prim.factor_group);
  }
	
  //next get a montiplet type structure (hop clusters radiating out of each basis site using
  //get_ext_montiplet()
  //first we need to copy the vector of orbits into a multiplet
	
  multiplet hoptiplet;
	
  //fill up the empty and point slots in the hoptiplet first
  for(int i=0; i<2; i++){
    vector<orbit>torbvec;
    hoptiplet.orb.push_back(torbvec);
  }
	
  //check whether any of the hop clusters have two sites
  for(int i=0; i < hoporb.size(); i++){
    if(hoporb[i].equiv[0].point.size() == 2){
      if(hoptiplet.orb.size()<3){
	vector<orbit> torbvec;
	torbvec.push_back(hoporb[i]);
	hoptiplet.orb.push_back(torbvec);
      }
      else{
	hoptiplet.orb[2].push_back(hoporb[i]);
      }
    }
  }
  //if no hop clusters with two sites are present pushback an empty
  if(hoptiplet.orb.size()<3){
    vector<orbit>torbvec;
    hoptiplet.orb.push_back(torbvec);
  }
	
  //check whether any of the hop clusters have three sites
  for(int i=0; i < hoporb.size(); i++){
    if(hoporb[i].equiv[0].point.size() ==3){
      if(hoptiplet.orb.size()<4){
	vector<orbit> torbvec;
	torbvec.push_back(hoporb[i]);
	hoptiplet.orb.push_back(torbvec);
      }
      else{
	hoptiplet.orb[3].push_back(hoporb[i]);
      }
    }
  }
	
  //for each basis site get a montiplet
	
  vector<multiplet> montihoptiplet;
	
  generate_ext_monteclust(basis,hoptiplet,montihoptiplet);
	
  //now construct the hop class
	
  jumps.clear();
	
  if(basis.size() != montihoptiplet.size()){
    cout << "mismatch between the size of basis and montihoptiplet \n";
    cout << "not initializing kinetic Monte Carlo \n";
    return;
  }
	
  for(int i=0; i<basis.size(); i++){
    //make an empty vector of hops
    vector<hop> tjumpvec;
    jumps.push_back(tjumpvec);
		
    //determine whether this basis site is an endpoint of any hop - as opposed to only serving as an activated state
    bool endpoint=false;
    for(int np=2; np<montihoptiplet[i].orb.size(); np++){
      for(int n=0; n<montihoptiplet[i].orb[np].size(); n++){
	for(int l=0; l<montihoptiplet[i].orb[np][n].equiv[0].point.size(); l++){
	  if(compare(basis[i].fcoord,montihoptiplet[i].orb[np][n].equiv[0].point[l].fcoord) &&
	     compare(basis[i].compon,montihoptiplet[i].orb[np][n].equiv[0].point[l].compon)){
	    if(montihoptiplet[i].orb[np][n].equiv[0].point[l].bit == 0 ) endpoint = true;
						
	  }
	}
      }
    }
		
    if(endpoint){
      for(int np=2; np<montihoptiplet[i].orb.size(); np++){
	for(int n=0; n<montihoptiplet[i].orb[np].size(); n++){
					
	  //check whether this particular cluster contains the basis as an endpoint and not an activated state
	  bool endpoint2 = false;
	  for(int l=0; l<montihoptiplet[i].orb[np][n].equiv[0].point.size(); l++){
	    if(compare(basis[i].fcoord,montihoptiplet[i].orb[np][n].equiv[0].point[l].fcoord) &&
	       compare(basis[i].compon,montihoptiplet[i].orb[np][n].equiv[0].point[l].compon)){
	      if(montihoptiplet[i].orb[np][n].equiv[0].point[l].bit == 0 ) endpoint2 = true;
	    }
	  }
					
	  if(endpoint2){
	    hop tjump;
	    tjump.b=i;          // assign the basis index for this jump object
	    tjump.initial=basis[i];
	    tjump.vac_spin_init=basis[i].get_spin(name);
						
	    for(int ne=0; ne<montihoptiplet[i].orb[np][n].equiv.size(); ne++){
							
	      tjump.endpoints.push_back(montihoptiplet[i].orb[np][n].equiv[ne]);
	      tjump.endpoints[ne].get_cart(prim.FtoC);
	      for(int k=0; k<tjump.endpoints[ne].point.size(); k++){
		//for each of the neighbors, get the cartesian coordinates
		if(!compare(tjump.endpoints[ne].point[k],basis[i])){
		  if(tjump.endpoints[ne].point[k].bit == 0){
		    mc_index tfinal;
		    int tvac_spin;
		    vec tjump_vec;
		    double tleng;
		    for(int l=0; l<3; l++){
		      tjump_vec.fcoord[l]=tjump.endpoints[ne].point[k].fcoord[l]-basis[i].fcoord[l];
		      tjump_vec.ccoord[l]=tjump.endpoints[ne].point[k].ccoord[l]-basis[i].ccoord[l];
		      tjump_vec.ccoord[l]=tjump_vec.ccoord[l]*(1.0e-8);   //convert to cm
		      tjump_vec.frac_on=true;
		      tjump_vec.cart_on=true;
		    }
		    tleng=tjump_vec.calc_dist();
		    for(int l=0; l<4; l++){
		      tfinal.shift[l]=tjump.endpoints[ne].point[k].shift[l];
		    }
		    tvac_spin=tjump.endpoints[ne].point[k].get_spin(name);
		    tjump.jump_vec.push_back(tjump_vec);
		    tjump.jump_leng.push_back(tleng);
		    tjump.final.push_back(tfinal);
		    tjump.vac_spin.push_back(tvac_spin);
		  }
		  if(tjump.endpoints[ne].point[k].bit == 1){
		    mc_index tactivated;
		    for(int l=0; l<4; l++){
		      tactivated.shift[l]=tjump.endpoints[ne].point[k].shift[l];
		    }
		    tjump.activated.push_back(tactivated);
		  }
		}
	      }
	    }
	    bool clear=true;
	    tjump.get_reach(montiplet,clear,basis);
	    jumps[i].push_back(tjump);
	  }
	}
      }
    }
  }
	
  cout << "about to enter extend_reach() \n";
	
  extend_reach();
	
  cout << "just passed extend_reach() \n";
	
  for(int i=0; i < jumps.size(); i++){
    cout << "HOP INFO FOR BASIS SITE  i = " << i << "\n";
    for(int j=0; j < jumps[i].size(); j++){
      cout << " HOP j = " << j << "\n";
      jumps[i][j].print_hop_info(cout);
    }
  }
	
  //determine the length of the prob and cprob arrays
  np=0;
  for(int i=0; i<jumps.size(); i++){
    for(int j=0; j<jumps[i].size(); j++){
      np=np+jumps[i][j].endpoints.size();
    }
  }
  np=np*nuc;
	
  prob = new double[np];
  cprob = new double[np];
  ltosp = new int[nmcL];
  ptol = new int[np];
  ptoj = new int[np];  // probability to jump
  ptoh = new int[np];   //probability to hop
	
		
  if(basis.size() != jumps.size()){
    cout << "number of basis sites is not equal to number of jump vectors\n";
    cout << "some kind of error occurred \n";
    cout << "NOT CONTINUING WITH initialize_kmc()\n";
    return;
  }
	
  //figure out which basis sites serve as a regular site at least once (store that info as bit=0, otherwise bit=1)
  for(int i=0; i<basis.size(); i++){
    if(jumps[i].size() == 0) basis[i].bit =1;
    else basis[i].bit =0;
  }
	
	
  Rx = new double[nmcL];
  Ry = new double[nmcL];
  Rz = new double[nmcL];
	
  int p=0;
  for(int l=0; l<nmcL; l++){
    if(basis[ltob[l]].bit == 0){
      ltosp[l]=p;
      for(int j=0; j<jumps[ltob[l]].size(); j++){
	for(int h=0; h<jumps[ltob[l]][j].endpoints.size(); h++){
	  ptol[p]=l;
	  ptoj[p]=j;
	  ptoh[p]=h;
	  p++;
	}
      }
    }
  }
	
  cout << "WE HAVE FILLED UP ALL THE HOP ARRAYS \n";
  cout << " p = " << p << "\n";
  cout << " np = " << np << "\n";
	
	
  hop_leng=0.0;
  for(int i=0; i<jumps.size(); i++){
    for(int j=0; j<jumps[i].size(); j++){
      for(int k=0; k<jumps[i][j].jump_leng.size(); k++){
	if(jumps[i][j].jump_leng[k] > hop_leng) hop_leng = jumps[i][j].jump_leng[k];
      }
    }
  }
	
  cout << "The maximum hop distance is " << hop_leng << "\n";
	
  if(hop_leng == 0.0){
    cout << "the maximum hop length is zero \n";
    cout << "you will have problems in your kmc simulation \n";
    cout << "check your hop file \n";
  }


  R.initialize(conc);
  kinL.initialize(conc);
  Dtrace.initialize(conc);
  corrfac.initialize(conc);
  AVkinL.initialize(conc);
  AVDtrace.initialize(conc);
  AVcorrfac.initialize(conc);

}

//************************************************************
// goes through the jump structure and for every hop reach site, adds additional reach sites that include sites that 
// each site in the existing reach can hop to

void Monte_Carlo::extend_reach(){

  for(int i=0; i<jumps.size(); i++){
    for(int j=0; j<jumps[i].size(); j++){
      for(int h=0; h < jumps[i][j].reach.size(); h++){
	int rs=jumps[i][j].reach[h].size();
	for(int r=0; r<rs; r++){
	  int b=jumps[i][j].reach[h][r].shift[3];
	  //for every hop from basis site b make the shift 
	  for(int ht=0; ht < jumps[b].size(); ht++){
	    for(int hh=0; hh<jumps[b][ht].endpoints.size(); hh++){
	      mc_index treach;
	      treach.shift[0]=jumps[i][j].reach[h][r].shift[0]+jumps[b][ht].final[hh].shift[0];
	      treach.shift[1]=jumps[i][j].reach[h][r].shift[1]+jumps[b][ht].final[hh].shift[1];
	      treach.shift[2]=jumps[i][j].reach[h][r].shift[2]+jumps[b][ht].final[hh].shift[2];
	      treach.shift[3]=jumps[b][ht].final[hh].shift[3];

	      //check whether this treach point already exists in reach[h]
	      bool isnew=true;
	      for(int n=0; n<jumps[i][j].reach[h].size(); n++){
		if(compare(jumps[i][j].reach[h][n],treach)){
		  isnew=false;
		  break;
		}
	      }
	      if(isnew) jumps[i][j].reach[h].push_back(treach);


	    }
	  }
	}
      }
    }
  }


}


//************************************************************

void Monte_Carlo::get_hop_prob(int i, int j, int k, int b, double beta){
  double nu=1.0e13;
	
  if(jumps[b].size() == 0) return;
	
  int l=index(i,j,k,b);
	
  if(mcL[l] == -1){
    //  if(mcL[l] == jumps[b][0].vac_spin_init){
		
    int cumh=0;
    for(int ht=0; ht < jumps[b].size(); ht++){
      for(int h=0; h<jumps[b][ht].endpoints.size(); h++){
	int ii=i+jumps[b][ht].final[h].shift[0];
	int jj=j+jumps[b][ht].final[h].shift[1];
	int kk=k+jumps[b][ht].final[h].shift[2];
	int bb=jumps[b][ht].final[h].shift[3];
				
	int ll=index(ii,jj,kk,bb);
	int pp=ltosp[l]+cumh;
	cumh++;
				
	//check whether this endpoint is occupied
				
	if(mcL[ll] == 1){
	  //	if(mcL[ll] != jumps[b][ht].vac_spin[h]){
	  double barrier=calc_barrier(i,j,k,b,ii,jj,kk,bb,l,ll,ht,h);
	  if(barrier > 0) prob[pp]=nu*exp(-beta*barrier);
	  else prob[pp]=nu;
	}
	else{ //the endpoint h is not occupied
	  prob[pp]=0.0;
	}
      }  
    }
  }
  else{  // the site is not vacant
    int cumh=0;
    for(int ht=0; ht < jumps[b].size(); ht++){
      for(int h=0; h<jumps[b][ht].endpoints.size(); h++){
	int pp=ltosp[l]+cumh;
	cumh++;
	prob[pp]=0.0;
      }
    }
  }
}


//************************************************************

double Monte_Carlo::calc_barrier(int i, int j, int k, int b, int ii, int jj, int kk, int bb, int l, int ll, int ht, int h){
  double barrier,penA1,penA2,penB1,penB2;
  double kra=0.28;
	
  if(jumps[b][ht].activated.size() > 0){   // intermediate activated state
		
    int iii=i+jumps[b][ht].activated[h].shift[0];
    int jjj=j+jumps[b][ht].activated[h].shift[1];
    int kkk=k+jumps[b][ht].activated[h].shift[2];
    int bbb=jumps[b][ht].activated[h].shift[3];
		
    int lll=index(iii,jjj,kkk,bbb);
		
		
    penA1=pointenergy(ii,jj,kk,bb);		
    penB1=pointenergy(iii,jjj,kkk,bbb);
		
    mcL[ll]=-mcL[ll];	
    mcL[lll]=-mcL[lll];
    penA2=pointenergy(ii,jj,kk,bb);
    penB2=pointenergy(iii,jjj,kkk,bbb);
		
    mcL[ll]=-mcL[ll];
    mcL[lll]=-mcL[lll];
		
    barrier=(penA2+penB2-penA1-penB1);
    return barrier;
  }
  else{   // no intermediate activated state
    penB1=pointenergy(i,j,k,b);
    penB2=pointenergy(ii,jj,kk,bb);
		
    mcL[l]=-mcL[l];
    mcL[ll]=-mcL[ll];
		
    penA1=pointenergy(i,j,k,b);
    penA2=pointenergy(ii,jj,kk,bb);
		
    mcL[l]=-mcL[l];
    mcL[ll]=-mcL[ll];
		
		
    barrier=0.5*(penA1+penA2-penB1-penB2)+kra;
    return barrier;
		
  }
}





//************************************************************

void Monte_Carlo::initialize_prob(double beta){

  for(int i=0; i<di; i++){
    for(int j=0; j<dj; j++){
      for(int k=0; k<dk; k++){
	for(int b=0; b<db; b++){
	  if(basis[b].bit == 0){
	    get_hop_prob(i,j,k,b,beta);
	  }
	}
      }
    }
  }
}


//************************************************************

int Monte_Carlo::pick_hop(){
  //get cummulative hop array
  //pick the interval in which ran0(idum) falls
  //the event that occurs is the index bounding the interval from above

  cprob[0]=prob[0];
  for(int p=1; p<np; p++) cprob[p]=cprob[p-1]+prob[p];

  tot_prob=cprob[np-1];

  double ran=ran0(idum)*tot_prob;

  int il=np-1;
  int mm=il;
  int ir=-1;

  while(il-ir > 1){
    double interv=il-ir;
    int mid=int(ceil(interv/2.0))+ir;

    if(cprob[mid] > ran) il=mid;
    else ir=mid;
  }

  return il;

}


//************************************************************

void Monte_Carlo::update_prob(int i, int j, int k, int b, int ht, int h, double beta){
  for(int n=0; n<jumps[b][ht].reach[h].size(); n++){
    int ii=i+jumps[b][ht].reach[h][n].shift[0];
    int jj=j+jumps[b][ht].reach[h][n].shift[1];
    int kk=k+jumps[b][ht].reach[h][n].shift[2];
    int bb=jumps[b][ht].reach[h][n].shift[3];
    get_hop_prob(ii,jj,kk,bb,beta);
  }
}




//************************************************************

void Monte_Carlo::kinetic(double beta, double n_pass, double n_equil_pass){

  //do some initializations of R-vectors and fluctuation variables
  for(int l=0; l<nmcL; l++){
    Rx[l]=0.0;
    Ry[l]=0.0;
    Rz[l]=0.0;
  }

  double tRx,tRy,tRz;
  double kmc_time=0.0;

  R.set_zero();
  kinL.set_zero();
  Dtrace.set_zero();
  corrfac.set_zero();
  num_hops.set_zero();
  AVkinL.set_zero();
  AVDtrace.set_zero();
  AVcorrfac.set_zero();

  calc_num_atoms();

  cout << "The number of atoms are \n";
  num_atoms.print_concentration(cout);
  cout << "\n";

  calc_concentration();


  initialize_prob(beta);
  

  for(int n=0; n<n_pass; n++){
    for(int nn=0; nn<nmcL; nn++){

      int p=pick_hop();

      //determine site of the hop and the hop for that site
      int l=ptol[p];
      int h=ptoh[p];
      int ht=ptoj[p];

      int i=ltoi[l];
      int j=ltoj[l];
      int k=ltok[l];
      int b=ltob[l];


      //determine the end point of the hop

      int ii=i+jumps[b][ht].final[h].shift[0];
      int jj=j+jumps[b][ht].final[h].shift[1];
      int kk=k+jumps[b][ht].final[h].shift[2];
      int bb=jumps[b][ht].final[h].shift[3];


      int ll=index(ii,jj,kk,bb); 


      int temp=mcL[l];
      mcL[l]=mcL[ll];
      mcL[ll]=temp;

      tRx=Rx[l];
      tRy=Ry[l];
      tRz=Rz[l];

      Rx[l]=Rx[ll]-jumps[b][ht].jump_vec[h].ccoord[0];
      Ry[l]=Ry[ll]-jumps[b][ht].jump_vec[h].ccoord[1];
      Rz[l]=Rz[ll]-jumps[b][ht].jump_vec[h].ccoord[2];

      Rx[ll]=tRx+jumps[b][ht].jump_vec[h].ccoord[0];
      Ry[ll]=tRy+jumps[b][ht].jump_vec[h].ccoord[1];
      Rz[ll]=tRz+jumps[b][ht].jump_vec[h].ccoord[2];

      kmc_time=kmc_time-log(ran0(idum))/tot_prob;

      update_num_hops(l,ll,b,bb);

      update_prob(i,j,k,b,ht,h,beta);
    }


    if(n > n_equil_pass){

      collect_R();

      kinL.evaluate(R);
      Dtrace=R;
      corrfac=R;
      Dtrace.normalize(num_atoms);

      double norm1=kmc_time*6.0;
      double norm2=hop_leng*hop_leng;

      kinL.normalize(norm1);
      double size=nuc;
      kinL.normalize(size,num_atoms);
      Dtrace.normalize(norm1);
      corrfac.normalize(norm2);
      corrfac.normalize(num_hops);


      AVkinL.increment(kinL);
      AVDtrace.increment(Dtrace);
      AVcorrfac.increment(corrfac);

    }

  }


  double norm3=n_pass-n_equil_pass;

  AVkinL.normalize(norm3);
  AVDtrace.normalize(norm3);
  AVcorrfac.normalize(norm3);


}



//************************************************************

void Monte_Carlo::collect_R(){



  for(int i=0; i<R.spin.size(); i++){

    for(int j=0; j<R.spin[i].size(); j++){
      R.Rx[i][j]=0.0;
      R.Ry[i][j]=0.0;
      R.Rz[i][j]=0.0;
      R.R2[i][j]=0.0;
      for(int l=0; l<nmcL; l++){
	//if(basis[ltob[l]].basis_flag == '0'){
	if(basis[ltob[l]].bit == 0){
	  if(mcL[l] == R.spin[i][j]){
	    R.Rx[i][j]=R.Rx[i][j]+Rx[l];
	    R.Ry[i][j]=R.Ry[i][j]+Ry[l];
	    R.Rz[i][j]=R.Rz[i][j]+Rz[l];
	    R.R2[i][j]=R.R2[i][j]+Rx[l]*Rx[l]+Ry[l]*Ry[l]+Rz[l]*Rz[l];
	  }
	}
      }
    }
  }
}


//************************************************************

void hop::get_reach(vector<multiplet> montiplet, bool clear, vector<atompos> basis){
	
  if(clear) reach.clear();
  for(int i=0; i < endpoints.size(); i++){
    vector < mc_index > treach;
    for(int j=0; j < endpoints[i].point.size(); j++){
      for(int k=0; k< montiplet.size(); k++){
	if(montiplet[k].orb[1].size()>1 || montiplet[k].orb[1][0].equiv.size() > 1 || montiplet[k].orb[1][0].equiv[0].point.size() > 1){
	  cout << "we have an unanticipated montiplet structure\n";
	  cout << "leaving get_reach and reach vector is not constructed \n";
	  cout << "expect errors\n";
	  return;
	}
	int trans[3];
	if(compare(endpoints[i].point[j],montiplet[k].orb[1][0].equiv[0].point[0],trans)){
	  for(int ii=0; ii<montiplet[k].orb.size(); ii++){
	    for(int jj=0; jj<montiplet[k].orb[ii].size(); jj++){
	      if(abs(montiplet[k].orb[ii][jj].eci) > 0.000000001){
		for(int kk=0; kk<montiplet[k].orb[ii][jj].equiv.size(); kk++){
		  for(int ll=0; ll<montiplet[k].orb[ii][jj].equiv[kk].point.size(); ll++){
		    mc_index tsite;
		    //first check whether the site is regular or activated before pushing back
		    if(basis[montiplet[k].orb[ii][jj].equiv[kk].point[ll].shift[3]].bit == 0){
		      for(int n=0; n<3; n++) tsite.shift[n]=montiplet[k].orb[ii][jj].equiv[kk].point[ll].shift[n]+trans[n];
		      tsite.shift[3]=montiplet[k].orb[ii][jj].equiv[kk].point[ll].shift[3];
		      if(new_mc_index(treach,tsite)) treach.push_back(tsite);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    reach.push_back(treach);
  }
}

//************************************************************

void hop::print_hop_info(ostream &stream){
  stream << "hop for basis site " << b << "\n";
  stream << "The vacancy spin of the initial site is " << vac_spin_init << "\n";
  stream << "The initial site is: \n";
  initial.print(stream);
  stream << "\n";
  stream << "NUMBER of hops = " << endpoints.size() << "\n";
  stream << "The clusters corresponding to each hop are: \n";
  for(int i=0; i<endpoints.size(); i++){
    stream << "Cluster " << i << "\n";
    endpoints[i].print(stream);
  }
  stream << "\n";
  stream << "\n";
  for(int i=0; i<jump_vec.size(); i++){
    stream << "jump vector for hop " << i << "\n";
    jump_vec[i].print_cart(stream);
  }
  stream << "\n";
  stream << "\n";

  stream << "Shifts of the final states of the hops \n";
  for(int i=0; i<final.size(); i++){
    stream << "final state " << i << "  ";
    final[i].print(stream);
    stream << "\n";
  }
  stream << "\n";
  stream << "\n";

  stream << "Reach of the hop \n";
  for(int i=0; i< reach.size(); i++){
    stream << "Reach for hop " << i << "\n";
    for(int j=0; j< reach[i].size(); j++){
      reach[i][j].print(stream);
      stream << "\n";
    }
    stream << "\n";
  }
  stream << "\n";
  stream << "\n";
  
}


//************************************************************

void mc_index::print(ostream &stream){
  for(int i=0; i<4; i++) stream << shift[i] << "   ";

}




//************************************************************

void fluctuation::initialize(concentration conc){
  //dimensions the fluctuation object to be compatible with the concentration object conc
  f.clear();
  compon.clear();
  for(int i=0; i<conc.compon.size(); i++){
    vector<double> tf;
    vector< specie > tcompon;
    for(int j=0; j<conc.compon[i].size(); j++){
      tcompon.push_back(conc.compon[i][j]);
      for(int k=i; k<conc.compon.size(); k++){
        for(int l=0; l<conc.compon[k].size(); l++){
	  if(k!=i || l>=j)
            tf.push_back(0.0);
	}
      }
    }
    f.push_back(tf);
    compon.push_back(tcompon);
  }

}



//************************************************************

void fluctuation::set_zero(){
  for(int i=0; i<f.size(); i++){
    for(int j=0; j<f[i].size(); j++){
      f[i][j]=0.0;
    }
  }
}


//************************************************************

void fluctuation::evaluate(concentration conc){

  if(conc.compon.size() != f.size()){
    cout << "concentration and fluctuation variables are not compatible\n";
    cout << "no update of fluctuation\n";
    return;
  }
  for(int i=0; i<conc.compon.size(); i++){
    int m=0;
    for(int j=0; j<conc.compon[i].size(); j++){
      for(int k=i; k<conc.compon.size(); k++){
        for(int l=0; l<conc.compon[k].size(); l++){
	  if(k!=i || l>=j){	
	    f[i][m]=conc.occup[i][j]*conc.occup[k][l];
	    m++;
	  }
	}
      }
    }
  }

}


//************************************************************

void fluctuation::evaluate(trajectory R){

  if(R.Rx.size() != f.size()){
    cout << "trajectory and fluctuation variables are not compatible\n";
    cout << "no update of fluctuation\n";
    return;
  }
  for(int i=0; i<R.Rx.size(); i++){
    int m=0;
    for(int j=0; j<R.Rx[i].size(); j++){
      for(int k=i; k<R.Rx.size(); k++){
        for(int l=0; l<R.Rx[k].size(); l++){
	  if(k!=i || l>=j){
	    f[i][m]=R.Rx[i][j]*R.Rx[k][l]+R.Ry[i][j]*R.Ry[k][l]+R.Rz[i][j]*R.Rz[k][l];
	    m++;
	  }
	}
      }
    }
  }

}



//************************************************************

void fluctuation::increment(fluctuation FF){

  if(f.size() != FF.f.size()){
    cout << "fluctuation variables are not compatible\n";
    cout << "cannot update \n";
    return;
  }

  for(int i=0; i<f.size(); i++){

    if(f[i].size() != FF.f[i].size()){
      cout << "fluctuation variables are not compatible\n";
      cout << "cannot update \n";
      return;
    }

    for(int j=0; j<f[i].size(); j++){
      f[i][j]=f[i][j]+FF.f[i][j];
    }
  }
}


//************************************************************

void fluctuation::decrement(fluctuation FF){

  if(f.size() != FF.f.size()){
    cout << "fluctuation variables are not compatible\n";
    cout << "cannot update \n";
    return;
  }

  for(int i=0; i<f.size(); i++){

    if(f[i].size() != FF.f[i].size()){
      cout << "fluctuation variables are not compatible\n";
      cout << "cannot update \n";
      return;
    }

    for(int j=0; j<f[i].size(); j++){
      f[i][j]=f[i][j]-FF.f[i][j];
    }
  }
}




//************************************************************

void fluctuation::normalize(double n){
  for(int i=0; i<f.size(); i++){
    for(int j=0; j<f[i].size(); j++){
      f[i][j]=f[i][j]/n;
    }
  }
}


//************************************************************

void fluctuation::normalize(double n, concentration conc){
  //for single component diffusion only normalizes the diagonal elements and sets the off-diagonal terms to zero
  if(compon.size() == 1 && compon[0].size() == 2){
    int m=0;
    for(int j=0; j< compon[0].size(); j++){
      for(int k=j; k< compon[0].size(); k++){
	if(j == k){
	  if(abs(conc.occup[0][j]) > tol) f[0][m]=f[0][m]/conc.occup[0][j];
	}
	else{
	  f[0][m]=0.0;
	}
	m++;
      }
    }
  }
  else{
    normalize(n);
  }

}



//************************************************************

void fluctuation::print(ostream &stream){
  if(compon.size() == 1 && compon[0].size() == 2){
    //only print the diagonal elements
    int m=0;
    for(int j=0; j<compon[0].size(); j++){
      for(int k=j; k<compon[0].size(); k++){
	if(j == k) stream << f[0][m] << "  ";
      }
    }
    m++;
  }
  
  else{
    for(int i=0; i<f.size(); i++){
      for(int j=0; j<f[i].size(); j++){
	stream << f[i][j] << "  ";
      }
    }
  }
  
}




//************************************************************

void fluctuation::print_elements(ostream &stream){

  if(compon.size() == 1 && compon[0].size() == 2){
    for(int j=0; j<compon[0].size(); j++){
      compon[0][j].print(stream);
      stream << "_";
      compon[0][j].print(stream);
      stream << "  ";
    }
  }
  else{
    for(int i=0; i<compon.size(); i++){
      for(int j=0; j<compon[i].size(); j++){
	for(int k=i; k<compon.size(); k++){
	  for(int l=0; l<compon[k].size(); l++){
	    if(k!=i || l>=j){
	      compon[i][j].print(stream);
	      stream << "_";
	      compon[k][l].print(stream);
	      stream << "  ";
	    }
	  }
	}
      }
    }
  }
}

//************************************************************

////////////////////////////////////////////////////////////////////////////////
//************************************************************
//added by Ben Swoboda
void get_clust_func(atompos atom1, atompos atom2, double &clust_func){
  int i,index;
  atom1.basis_vec.clear();

  index=atom2.bit;

  //determine which basis is to be used and store values in basis vector
  if(atom1.basis_flag == '1'){
    for(i=0; i<atom1.p_vec.size(); i++){
      atom1.basis_vec.push_back(atom1.p_vec[i]);
    }
    //        index++;
  }
  else{
    for(i=0; i<atom1.spin_vec.size(); i++){
      atom1.basis_vec.push_back(atom1.spin_vec[i]);
    }
  }

  clust_func=clust_func*atom1.basis_vec[index];

  basis_type=atom1.basis_flag;

  return;

}
//************************************************************
////////////////////////////////////////////////////////////////////////////////


//************************************************************

void calc_correlations(structure struc, multiplet super_basiplet, arrangement &conf){
  int nm,no,nc,np,na,i,j,k;
  double tcorr,clust_func;

  //push back the correlation for the empty cluster
  tcorr=1.0;
  conf.correlations.push_back(tcorr);

  for(nm=1; nm<super_basiplet.orb.size(); nm++){
    for(no=0; no<super_basiplet.orb[nm].size(); no++){
      tcorr=0;
      if(super_basiplet.orb[nm][no].equiv.size()==0) {
	cout << "something screwed up, no cluster in your orbit \n";
	exit(1);
      }
      for(nc=0; nc<super_basiplet.orb[nm][no].equiv.size(); nc++){
	clust_func=1;
	for(np=0; np<super_basiplet.orb[nm][no].equiv[nc].point.size(); np++){
	  for(na=0; na < struc.atom.size(); na++){
	    if(compare(struc.atom[na].fcoord,super_basiplet.orb[nm][no].equiv[nc].point[np].fcoord)){
	      //	      for(i=0; i<=super_basiplet.orb[nm][no].equiv[nc].point[np].bit; i++)
	      //		clust_func=clust_func*struc.atom[na].occ.spin;
	      //	      break;
	      ////////////////////////////////////////////////////////////////////////////////
	      //added by Ben Swoboda
	      //call get_clust_func => returns cluster function using desired basis
	      get_clust_func(struc.atom[na], super_basiplet.orb[nm][no].equiv[nc].point[np], clust_func);
	      //cout << "spin: " << struc.atom[na].occ.spin << "\tbit: " << super_basiplet.orb[nm][no].equiv[nc].point[np].bit <<
	      //"\tspecie: " << struc.atom[na].occ.name << "\n";
	      break;
	      ////////////////////////////////////////////////////////////////////////////////
	    }
	  }
	  if(na == struc.atom.size()){
	    cout << "have not mapped a cluster point on the crystal \n";
	    cout << "inside of calc_correlations \n";
	  }
	}
	tcorr=tcorr+clust_func;
      }
      tcorr=tcorr/super_basiplet.orb[nm][no].equiv.size();
      conf.correlations.push_back(tcorr);
    }
  }

  return;

}


//************************************************************
void get_super_basis_vec(structure &superstruc, vector < vector < vector < int > > > &super_basis_vec){
  int ns,i,j;

  /*Populate super_basis_vec, which contains the spin basis (including exponentiation or bit 
    differentiation for ternary and higher order systems).  The super basis vector is ordered as follows:
    atom -> component specie -> basis components.  Currently supports spin and occupation bases, but
    will have to be edited to accomodate new bases.  Can be made more general if the data structures
    for storing occupation variable bases are generalized. */

  super_basis_vec.clear();

  //cout << "Printing super_basis_vec:  \n";
  for(i=0; i<superstruc.atom.size(); i++){
    vector < vector < int > > tcompon_vec;
    for(j=0; j<superstruc.atom[i].compon.size(); j++){
      superstruc.atom[i].occ=superstruc.atom[i].compon[j];
      vector<int> tbasis_vec;
      if(superstruc.atom[i].basis_flag=='1' && superstruc.atom[i].compon.size()>1){
	//cout << "\n Atom " << i << ": ";
	//cout << "Occupation basis ";
	get_basis_vectors(superstruc.atom[i]);
	for(int k=0; k<superstruc.atom[i].p_vec.size(); k++){
	  tbasis_vec.push_back(superstruc.atom[i].p_vec[k]);
	  //cout << "  " << tbasis_vec.back();
	}
      }
      else if(superstruc.atom[i].compon.size()>1){
	//cout << "\n Atom " << i << ": ";
	//cout << "Spin basis ";
	get_basis_vectors(superstruc.atom[i]);
	for(int k=0; k<superstruc.atom[i].spin_vec.size(); k++){
	  tbasis_vec.push_back(superstruc.atom[i].spin_vec[k]);
	  //cout << "  " << tbasis_vec.back() ;
	}
      }
      tcompon_vec.push_back(tbasis_vec);
    }
    super_basis_vec.push_back(tcompon_vec);
    superstruc.atom[i].occ=superstruc.atom[i].compon[0];
  }
  return;
}



//************************************************************

void get_corr_vector(structure &struc, multiplet &super_basiplet, vector< vector< vector< vector< int > > > > &corr_vec){
  int nm,no,nc,np,na,i,j,k;
  /*Go through super_basiplet and determine the basis functions for each correlation.  Store the site indeces and bit orderings
    for each local basis function of each correlation.  corr_vec is ordered as follows:
    basis function (as numbered in BCLUST) -> local equivalent basis functions -> atomic sites -> site number and bit/exponent value
  */

  corr_vec.clear();
  for(nm=1; nm<super_basiplet.orb.size(); nm++){
    for(no=0; no<super_basiplet.orb[nm].size(); no++){

      if(super_basiplet.orb[nm][no].equiv.size()==0) {
	cout << "something screwed up, no cluster in your orbit \n";
	exit(1);
      }
      vector < vector < vector < int > > > tcorr_vec;
      for(nc=0; nc<super_basiplet.orb[nm][no].equiv.size(); nc++){
	vector< vector< int > > func_vec;
	for(np=0; np<super_basiplet.orb[nm][no].equiv[nc].point.size(); np++){
	  for(na=0; na < struc.atom.size(); na++){
	    if(compare(struc.atom[na].fcoord,super_basiplet.orb[nm][no].equiv[nc].point[np].fcoord)){
	      vector< int > bit_vec;
	      bit_vec.push_back(na); //push back site index
	      bit_vec.push_back(super_basiplet.orb[nm][no].equiv[nc].point[np].bit); //push back bit value
	      func_vec.push_back(bit_vec);
	      break;
	    }
	  }
	  if(na == struc.atom.size()){
	    cout << "have not mapped a cluster point on the crystal \n";
	    cout << "inside of calc_correlations \n";
	  }
	}
	tcorr_vec.push_back(func_vec);
      }
      corr_vec.push_back(tcorr_vec);
    }
  }
  return;

}



//************************************************************

bool new_conf(arrangement &conf,superstructure &superstruc){

  for(int i=0; i<superstruc.conf.size(); i++){
    int j=0;
    while(j < superstruc.conf[i].correlations.size() &&
	  abs(conf.correlations[j]-superstruc.conf[i].correlations[j]) < tol){
      j++;
    }
    if(j == superstruc.conf[i].correlations.size()) return false;
  }
  return true;
}



//************************************************************

bool new_conf(arrangement &conf,vector<superstructure> &superstruc){

  for(int i=0; i<superstruc.size(); i++){
    if(!new_conf(conf,superstruc[i])) return false;
  }
  return true;
}




//************************************************************

void generate_ext_clust(structure struc, int min_num_compon, int max_num_points,
			vector<double> max_radius, multiplet &clustiplet){

  int i,j,k,m,n,np,nc;
  int dim[3];
  int num_basis=0;
  vector<atompos> basis;
  vector<atompos> gridstruc;



  //first get the basis sites

  {           // make the basis from which the sites for the local clusters are to be picked
    for(i=0; i<struc.atom.size(); i++)
      if(struc.atom[i].compon.size() >= min_num_compon){
	basis.push_back(struc.atom[i]);
	num_basis++;
      }
  }           // end of the basis generation



  //make the empty cluster
  //then starting from the empty cluster, start building multipoint clusters


  //make the empty cluster

  {          // beginning of the point cluster generation block
    vector<orbit> torbvec;
    cluster tclust;
    tclust.max_leng=0;
    tclust.min_leng=0;
    orbit torb;
    torb.equiv.push_back(tclust);
    torbvec.push_back(torb);
    clustiplet.orb.push_back(torbvec);
  }          // end of the empty cluster generation block


  //make a sphere with max_radius[2] and collect all crystal sites within that

  lat_dimension(struc.lat,max_radius[2],dim);



  for(i=-dim[0]; i<=dim[0]; i++){
    for(j=-dim[1]; j<=dim[1]; j++){
      for(k=-dim[2]; k<=dim[2]; k++){
        vec tlatt;

        tlatt.ccoord[0]=i*struc.lat[0][0]+j*struc.lat[1][0]+k*struc.lat[2][0];
        tlatt.ccoord[1]=i*struc.lat[0][1]+j*struc.lat[1][1]+k*struc.lat[2][1];
        tlatt.ccoord[2]=i*struc.lat[0][2]+j*struc.lat[1][2]+k*struc.lat[2][2];
        tlatt.cart_on=true;
        conv_AtoB(struc.CtoF,tlatt.ccoord,tlatt.fcoord);
        tlatt.frac_on=true;

        for(m=0; m<num_basis; m++){
          atompos tatom;
	  tatom=basis[m];
          for(int ii=0; ii<3; ii++){
            tatom.ccoord[ii]=basis[m].ccoord[ii]+tlatt.ccoord[ii];
            tatom.fcoord[ii]=basis[m].fcoord[ii]+tlatt.fcoord[ii];
          }

          //get distance to closest basis site in the unit cell at the origin

          double min_dist=1e20;
          for(n=0; n<num_basis; n++){
	    double temp[3];
	    double dist=0.0;
	    for(int ii=0; ii<3; ii++){
	      temp[ii]=tatom.ccoord[ii]-basis[n].ccoord[ii];
	      dist=dist+temp[ii]*temp[ii];
	    }
	    dist=sqrt(dist);
	    if(dist < min_dist)min_dist=dist;
          }
          if(min_dist < max_radius[2]) {
            gridstruc.push_back(tatom);
          }
        }
      }
    }
  }



  //for each cluster of the previous size, add points from gridstruc
  //   - see if the new cluster satisfies the size requirements
  //   - see if it is new
  //   - generate all its equivalents

  for(np=1; np<=max_num_points; np++){
    vector<orbit> torbvec;
    for(nc=0; nc<clustiplet.orb[np-1].size(); nc++){
      for(n=0; n<gridstruc.size(); n++){
        cluster tclust;
        atompos tatom;

        tatom=gridstruc[n];

        if(clustiplet.orb[np-1][nc].equiv.size() == 0){
          cout << "something screwed up \n";
          exit(1);
        }


        tclust=clustiplet.orb[np-1][nc].equiv[0];

        tclust.point.push_back(tatom);

        tclust.get_dimensions();

	if(tclust.point.size() == 1 && new_clust(tclust,torbvec)){
          orbit torb;
          torb.equiv.push_back(tclust);

          get_equiv(torb,struc.factor_group);

          torbvec.push_back(torb);
        }
	else{
	  if(tclust.max_leng < max_radius[np] && tclust.min_leng > tol && new_clust(tclust,torbvec)){
	    orbit torb;
	    torb.equiv.push_back(tclust);
	    get_equiv(torb,struc.factor_group);
	    torbvec.push_back(torb);
	  }
	}
      }
    }
    clustiplet.orb.push_back(torbvec);
    clustiplet.sort(np);
  }


}




//************************************************************


//************************************************************
void read_junk(istream &stream){  // added by jishnu //to skip the reading till the new line
  char j;
  do{
    stream.get(j);		
  }while(!(stream.peek()=='\n'));
}
//************************************************************

void generate_loc_clust(structure struc, int min_num_compon, int max_num_points,
			vector<double> max_radius, multiplet &loc_clustiplet, cluster clust){

  int i,j,k,m,n,np,nc;
  int dim[3];
  int num_basis=0;
  vector<atompos> basis;
  vector<atompos> gridstruc;

  //first get the basis sites

  {           // make the basis from which the sites for the local clusters are to be picked
    for(i=0; i<struc.atom.size(); i++)
      if(struc.atom[i].compon.size() >= min_num_compon){
	basis.push_back(struc.atom[i]);
	num_basis++;
      }
  }           // end of the basis generation


  //now generate the local clusters emanating from input cluster

  //the first multiplet is simply the cluster itself

  {          // beginning of the point cluster generation block
    vector<orbit> torbvec;
    orbit torb;
    torb.equiv.push_back(clust);
    torbvec.push_back(torb);

    loc_clustiplet.orb.push_back(torbvec);
  }          // end of the empty cluster generation block




  //make a sphere with max_radius[2]+clust.max_leng and collect all crystal basis sites within that

  lat_dimension(struc.lat,max_radius[2]+clust.max_leng,dim);

  for(i=-dim[0]; i<=dim[0]; i++){
    for(j=-dim[1]; j<=dim[1]; j++){
      for(k=-dim[2]; k<=dim[2]; k++){
        vec tlatt;

        tlatt.ccoord[0]=i*struc.lat[0][0]+j*struc.lat[1][0]+k*struc.lat[2][0];
        tlatt.ccoord[1]=i*struc.lat[0][1]+j*struc.lat[1][1]+k*struc.lat[2][1];
        tlatt.ccoord[2]=i*struc.lat[0][2]+j*struc.lat[1][2]+k*struc.lat[2][2];
        tlatt.cart_on=true;
        conv_AtoB(struc.CtoF,tlatt.ccoord,tlatt.fcoord);
        tlatt.frac_on=true;

        for(m=0; m<num_basis; m++){
          atompos tatom;
	  tatom=basis[m];
          for(int ii=0; ii<3; ii++){
            tatom.ccoord[ii]=basis[m].ccoord[ii]+tlatt.ccoord[ii];
            tatom.fcoord[ii]=basis[m].fcoord[ii]+tlatt.fcoord[ii];
          }

          //get distance to the site of interest atom_num

	  double min_dist=1.0e20;
	  for(n=0; n<clust.point.size(); n++){
	    double temp[3];
	    double dist=0.0;
	    for(int ii=0; ii<3; ii++){
	      temp[ii]=tatom.ccoord[ii]-clust.point[n].ccoord[ii];
	      dist=dist+temp[ii]*temp[ii];
	    }
	    dist=sqrt(dist);
	    if(dist < min_dist)min_dist=dist;
	  }
	  if(min_dist < max_radius[2]) {
	    gridstruc.push_back(tatom);
	  }
	}
      }
    }
  }

  //for each cluster of the previous size, add points from gridstruc
  //   - see if the new cluster satisfies the size requirements
  //   - see if it is new
  //   - generate all its equivalents

  for(np=1; np<=max_num_points-1; np++){
    vector<orbit> torbvec;

    for(nc=0; nc<loc_clustiplet.orb[np-1].size(); nc++){
      for(n=0; n<gridstruc.size(); n++){
        cluster tclust;
        atompos tatom;

        tatom=gridstruc[n];

        if(loc_clustiplet.orb[np-1][nc].equiv.size() == 0){
          cout << "something screwed up \n";
          exit(1);
        }

        tclust=loc_clustiplet.orb[np-1][nc].equiv[0];
        tclust.point.push_back(tatom);

        tclust.get_dimensions();

        if(tclust.max_leng < max_radius[np+1] && tclust.min_leng > tol && new_loc_clust(tclust,torbvec)){
          orbit torb;
          torb.equiv.push_back(tclust);
          get_loc_equiv(torb,clust.clust_group);
          torbvec.push_back(torb);
        }
      }
    }
    loc_clustiplet.orb.push_back(torbvec);
    loc_clustiplet.sort(np);
  }

  cout << "\n";
  cout << "LOCAL CLUSTER \n";
  cout << "\n";

  for(np=1; np<=max_num_points-1; np++)
    loc_clustiplet.print(cout);

}


//************************************************************

void calc_clust_symmetry(structure struc, cluster &clust){
  int pg,na,i,j,k,n,m,num_suc_maps;
  atompos hatom;
  double shift[3],temp[3];
  sym_op tclust_group;
  cluster tclust;


  //all symmetry operations are done within the fractional coordinate system
  //since translations back into the unit cell are straightforward


  //apply a point group operation to the cluster
  //then see if a translation (shift) maps the transformed cluster onto the original cluster
  //if so test whether this point group operation and translation maps the crystal onto itself


  for(pg=0; pg<struc.point_group.size(); pg++){

    tclust=clust.apply_sym(struc.point_group[pg]);

    for(i=0; i<clust.point.size(); i++){
      if(compare(clust.point[0].compon, tclust.point[i].compon)){

	for(j=0; j<3; j++) shift[j]=clust.point[0].fcoord[j]-tclust.point[i].fcoord[j];

	//check whether the cluster maps onto itself

	num_suc_maps=0;
	for(n=0; n<clust.point.size(); n++){
	  for(m=0; m<clust.point.size(); m++){
	    if(compare(clust.point[n].compon, clust.point[m].compon)){
	      for(j=0; j<3; j++) temp[j]=clust.point[n].fcoord[j]-tclust.point[m].fcoord[j]-shift[j];

	      k=0;
	      for(j=0; j<3; j++)
		if(abs(temp[j]) < 0.00005 ) k++;
	      if(k==3)num_suc_maps++;
	    }
	  }
	}

	if(num_suc_maps == clust.point.size()){
	  //the cluster after transformation and translation maps onto itself
	  //now check whether the rest of the crystal maps onto itself

	  //apply the point group to the crystal

	  vector<atompos> tatom;
	  for(na=0; na<struc.atom.size(); na++){
	    hatom=struc.atom[na].apply_sym(struc.point_group[pg]);
	    tatom.push_back(hatom);
	  }
	  //check whether after translating with shift it maps onto itself
	  num_suc_maps=0;
	  for(n=0; n<struc.atom.size(); n++){
	    for(m=0; m<struc.atom.size(); m++){
	      if(compare(struc.atom[n].compon,struc.atom[m].compon)){
		for(j=0; j<3; j++) temp[j]=struc.atom[n].fcoord[j]-tatom[m].fcoord[j]-shift[j];
		within(temp);

		k=0;
		for(j=0; j<3; j++)
		  if(abs(temp[j]) < 0.00005 ) k++;
		if(k==3)num_suc_maps++;
	      }
	    }
	  }

	  if(num_suc_maps == struc.atom.size()){

	    //check whether the symmetry operation already exists in the factorgroup array

	    int ll=0;
	    for(int ii=0; ii<clust.clust_group.size(); ii++)
	      if(compare(struc.point_group[pg].fsym_mat,clust.clust_group[ii].fsym_mat)
		 && compare(shift,clust.clust_group[ii].ftau) )break;
	      else ll++;

	    // if the symmetry operation is new, add it to the clust_group vector
	    // and update all info about the sym_op object

	    if(clust.clust_group.size() == 0 || ll == clust.clust_group.size()){
	      tclust_group.frac_on=false;
	      tclust_group.cart_on=false;
	      for(int jj=0; jj<3; jj++){
		tclust_group.ftau[jj]=shift[jj];
		for(int kk=0; kk<3; kk++){
		  tclust_group.fsym_mat[jj][kk]=struc.point_group[pg].fsym_mat[jj][kk];
		  tclust_group.lat[jj][kk]=struc.lat[jj][kk];
		}
	      }
	      tclust_group.frac_on=true;
	      tclust_group.update();
	      clust.clust_group.push_back(tclust_group);
	    }
	  }
	  tatom.clear();
	}
      }
    }
  }
  return;
}

//************************************************************

void generate_ext_basis_environ(structure struc, multiplet clustiplet, multiplet &basiplet){  // jishnu
  int np,no,i,j,k;
	
  //go through the clustiplet
  //for each orbit within the clustiplet, take the first of the equivalent clusters,
  //  enumerate each exponent sequence for that cluster
  //     first check if the exponent sequence on the cluster is new compared with already considered sequences
  //     for each new exponent sequence, generate all equivalent basis clusters by doing:
  //                 -for each factor_group
  //                    apply factor_group to the cluster
  //                    determine the cluster group symmetry
  //                    for each clust_group
  //                       apply clust_group to the cluster
  //                       determine if a new cluster
	
	
  if(clustiplet.orb.size() > 0){
    basiplet.orb.push_back(clustiplet.orb[0]);
  }
	
  for(np=1; np<clustiplet.orb.size(); np++){
    vector<orbit> torbvec;
    for(no=0; no<clustiplet.orb[np].size(); no++){
      cluster tclust;
      tclust=clustiplet.orb[np][no].equiv[0];
			
      //enumerate each exponent sequence for this cluster
      int last=0;
			
      while(last == 0){
	tclust.point[0].bit++;
	for(i=0; i<(tclust.point.size()-1); i++){
	  if(tclust.point[i].bit !=0 && tclust.point[i].bit%tclust.point[i].compon.size() == 0){   // changed by jishnu // got rid of '-1' to get all possible bits (considering Vacancies as one component)
	    tclust.point[i+1].bit++;
	    tclust.point[i].bit=0;
	  }
	}
	if(tclust.point[tclust.point.size()-1].bit !=0 &&
	   tclust.point[tclust.point.size()-1].bit%tclust.point[tclust.point.size()-1].compon.size() == 0){ // changed by jishnu // got rid of '-1' to get all possible bits (considering Vacancies as one component)
	  last=last+1;
	  tclust.point[tclust.point.size()-1].bit=0;
	}
				
				
	// check if this cluster already exists
	if(new_clust(tclust, torbvec)){
					
	  //if not apply factor group
					
	  orbit torb;
	  for(int fg=0; fg<struc.factor_group.size(); fg++){
	    cluster tclust1;
	    tclust1=tclust.apply_sym(struc.factor_group[fg]);
	    within(tclust1);
	    tclust1.get_cart(struc.FtoC);
						
	    //determine cluster symmetry
						
	    calc_clust_symmetry(struc,tclust1);
						
	    //apply cluster symmetry and check if already part of current orbit
						
	    for(int cg=0; cg<tclust1.clust_group.size(); cg++){
	      cluster tclust2;
	      tclust2=tclust1.apply_sym(tclust1.clust_group[cg]);
	      tclust2.get_cart(struc.FtoC);
							
	      if(new_clust(tclust2,torb)){
		torb.equiv.push_back(tclust2);
	      }
	    }
	  }
	  if(torb.equiv.size() !=0){
	    torbvec.push_back(torb);
	  }
	}
      }
    }
    basiplet.orb.push_back(torbvec);
  }
  return;
}



////////////////////////////////////////////////////////////////////////////////


//************************************************************

void generate_ext_basis(structure struc, multiplet clustiplet, multiplet &basiplet){
  int np,no,i,j,k;

  //go through the clustiplet
  //for each orbit within the clustiplet, take the first of the equivalent clusters,
  //  enumerate each exponent sequence for that cluster
  //     first check if the exponent sequence on the cluster is new compared with already considered sequences
  //     for each new exponent sequence, generate all equivalent basis clusters by doing:
  //                 -for each factor_group
  //                    apply factor_group to the cluster
  //                    determine the cluster group symmetry
  //                    for each clust_group
  //                       apply clust_group to the cluster
  //                       determine if a new cluster


  if(clustiplet.orb.size() > 0){
    basiplet.orb.push_back(clustiplet.orb[0]);
  }

  for(np=1; np<clustiplet.orb.size(); np++){
    vector<orbit> torbvec;
    for(no=0; no<clustiplet.orb[np].size(); no++){
      cluster tclust;
      tclust=clustiplet.orb[np][no].equiv[0];

      //enumerate each exponent sequence for this cluster
      int last=0;

      while(last == 0){
	tclust.point[0].bit++;
	for(i=0; i<(tclust.point.size()-1); i++){
	  if(tclust.point[i].bit !=0 && tclust.point[i].bit%(tclust.point[i].compon.size()-1) == 0){
	    tclust.point[i+1].bit++;
	    tclust.point[i].bit=0;
	  }
	}
	if(tclust.point[tclust.point.size()-1].bit !=0 &&
	   tclust.point[tclust.point.size()-1].bit%(tclust.point[tclust.point.size()-1].compon.size()-1) == 0){
	  last=last+1;
	  tclust.point[tclust.point.size()-1].bit=0;
	}


	// check if this cluster already exists
	if(new_clust(tclust, torbvec)){

	  //if not apply factor group

	  orbit torb;
	  for(int fg=0; fg<struc.factor_group.size(); fg++){
	    cluster tclust1;
	    tclust1=tclust.apply_sym(struc.factor_group[fg]);
	    within(tclust1);
	    tclust1.get_cart(struc.FtoC);

	    //determine cluster symmetry

	    calc_clust_symmetry(struc,tclust1);

	    //apply cluster symmetry and check if already part of current orbit

	    for(int cg=0; cg<tclust1.clust_group.size(); cg++){
	      cluster tclust2;
	      tclust2=tclust1.apply_sym(tclust1.clust_group[cg]);
	      tclust2.get_cart(struc.FtoC);

	      if(new_clust(tclust2,torb)){
		torb.equiv.push_back(tclust2);
	      }
	    }
	  }
	  if(torb.equiv.size() !=0){
	    torbvec.push_back(torb);
	  }
	}
      }
    }
    basiplet.orb.push_back(torbvec);
  }
  return;
}



////////////////////////////////////////////////////////////////////////////////
//added by anton - filters a multiplet for clusters containing just one activated site (with occupation basis = 1)

//************************************************************
void filter_activated_clust(multiplet clustiplet, multiplet &activatedclustiplet){


  //clear activatedclustiplet
  activatedclustiplet.orb.clear();
  activatedclustiplet.size.clear();
  activatedclustiplet.order.clear();
  activatedclustiplet.index.clear();
  activatedclustiplet.subcluster.clear();

  //copy the empty cluster into activatedclustiplet

  if(clustiplet.orb.size() >= 1){
    activatedclustiplet.orb.push_back(clustiplet.orb[0]);
  }
  

  for(int i=1; i<clustiplet.orb.size(); i++){
    vector<orbit> torb;
    for(int j=0; j<clustiplet.orb[i].size(); j++){
      int num_activated=0;
      for(int k=0; k<clustiplet.orb[i][j].equiv[0].point.size(); k++){
	if(clustiplet.orb[i][j].equiv[0].point[k].basis_flag == '1') {
	  num_activated++;
	}
      }
      if(num_activated == 1){
	torb.push_back(clustiplet.orb[i][j]);
      }
    }
    activatedclustiplet.orb.push_back(torb);
  }
  
  
}

////////////////////////////////////////////////////////////////////////////////

//added by anton
//************************************************************

void merge_multiplets(multiplet clustiplet1, multiplet clustiplet2, multiplet &clustiplet3){

  //merge clustiplet1 with clustiplet2 and put it in clustiplet3


  for(int np=0; np<clustiplet1.orb.size(); np++){
    clustiplet3.orb.push_back(clustiplet1.orb[np]);
  }


  for(int np=1; np<clustiplet2.orb.size(); np++){
    if(np> clustiplet3.orb.size()) clustiplet3.orb.push_back(clustiplet2.orb[np]);
    else{
      //add only orbits from clustiplet2.orb[np] that are new
      for(int i=0; i<clustiplet2.orb[np].size(); i++){
	bool isnew = true;
	for(int j=0; j<clustiplet3.orb[np].size(); j++){
	  if(compare(clustiplet2.orb[np][i],clustiplet3.orb[np][j])){
	    isnew = false;
	    break;
	  }
	}
	if(isnew) clustiplet3.orb[np].push_back(clustiplet2.orb[np][i]);
      }
    }
  }

  for(int np=1; np < clustiplet3.orb.size(); np++){
    clustiplet3.sort(np);
  }

}




//************************************************************

void generate_ext_monteclust(vector<atompos> basis, multiplet basiplet, vector<multiplet> &montiplet){


  //takes the basiplet and for each basis site in prim and enumerates all clusters that go through that basis site
  //then it determines the shift indices for each point of these clusters

  //for each basis site (make sure the basis site is within the primitive unit cell)
  //go through the basiplet and for each cluster within an orbit
  //        translate the cluster such that each point has been within the primitive unit cell
  //        if the resulting cluster has a point that coincides with the basis site,
  //        add it to the orbit for that cluster type in montiplet
  //        also determine the shift indices for this cluster

  for(int na=0; na<basis.size(); na++){   //select an atompos from prim, LOOP in atom

    multiplet tmontiplet;
    int clust_count=0;           //Edited by John -> Track mapping of clusters to basis sites   
    //make sure the basis atom is within the primitive unit cell
    within(basis[na]);

    //first make an empty cluster and add it to tmontiplet
    {
      cluster tclust;
      orbit torb;  // has a non-equivalent cluster, a colection of equivalent clusters, and ECI value
      torb.eci=basiplet.orb[0][0].eci;  // ECI of empty cluster
      torb.equiv.push_back(tclust);     // push empty cluster
      vector<orbit>torbvec;
      torbvec.push_back(torb);
      tmontiplet.orb.push_back(torbvec);

      vector<int> tind_vec;               //Edited code                                                                                                                                                                         
      tind_vec.push_back(clust_count);    // Edited code                                                                                                                                                                        
      tmontiplet.index.push_back(tind_vec);  //Edited code                                                                                                                                                                      
      clust_count++; // Edited code                        
    }

    //go through each cluster of basiplet
    for(int np=1; np<basiplet.orb.size(); np++){    // select one row of orb table, from row 1, LOOP in orb
      vector<orbit> torbvec;
      vector<int> tind_vec;  // Edited code 

      for(int no=0; no<basiplet.orb[np].size(); no++){  // select each clumn in a given row of orb table, LOOP orbit
	orbit torb;
	torb.eci=basiplet.orb[np][no].eci;
	bool found=false;
	for(int neq=0; neq<basiplet.orb[np][no].equiv.size(); neq++){  // LOOP equiv.size()
	  cluster tclust=basiplet.orb[np][no].equiv[neq];   // select an equivalent cluster

	  //for each point of the cluster translate the cluster so that point lies
	  //within the primitive unit cell

	  for(int n=0; n<tclust.point.size(); n++){ //LOOP n<np, np is the row index of orb table, also is the size of cluster, 1 means point
	    within(tclust,n);      //translate nth point of a cluster into prim cell, other points translate accordingly
	    //check whether the basis site basis[na] belongs to this cluster
	    if(compare(tclust.point[n].fcoord, basis[na].fcoord)){  //if the selected point belongs to a given basis point, then push this cluster
	      //add the cluster to the orbit
	      torb.equiv.push_back(tclust);  // if the given basis site basis[na] belongs this cluster, then push this cluster
	      found=true;
	    }
	  }  // LOOP n<np


	}    // LOOP equiv.size()
	if(found){
	  torbvec.push_back(torb);
	  tind_vec.push_back(clust_count); //Edited Code -- map index of cluster in basiplet to index of cluster in montiplet
	}
	clust_count++;
      }  // LOOP orbit
      tmontiplet.orb.push_back(torbvec);
      tmontiplet.index.push_back(tind_vec);  // Edited code   
    }   //LOOP in orb
    montiplet.push_back(tmontiplet);
  } // LOOP atom
  //after the above part, each basis atompos has a orb table.



  //work out the shift tables in each cluster object
  //these tell us which basis the point of the cluster belongs to and the coordinates of the unit cell

  for(int nm=0; nm < montiplet.size(); nm++){  // LOOP start
    for(int np=0; np < montiplet[nm].orb.size(); np++){
      for(int no=0; no < montiplet[nm].orb[np].size(); no++){
        for(int ne=0; ne < montiplet[nm].orb[np][no].equiv.size(); ne++){
          for(int n=0; n < montiplet[nm].orb[np][no].equiv[ne].point.size(); n++){

            get_shift(montiplet[nm].orb[np][no].equiv[ne].point[n], basis);

          }
        }
      }
    }
  }  // LOOP end
}






//************************************************************
//Useful functions
//************************************************************

void double_to_string(double n, string &a, int dec_places){
  //Only works for base 10 numbers
  double nn;
  int i;
  if(n<0)
    a.push_back('-');
  n=abs(n);
  nn=floor(n);
  i=int(nn);
  int_to_string(i, a, 10);
  if(dec_places>0)
    a.push_back('.');
  while(dec_places>0){
    n=10*(n-nn);
    nn=floor(n);
    i=int(nn);
    int_to_string(i, a, 10);
    dec_places--;	
  }
  return;
}

//************************************************************


void int_to_string(int i, string &a, int base){
  int ii=i;
  string aa;

  if(ii==0){
    a.push_back(ii+48);
    return;
  }

  if(ii<0)a.push_back('-');
  ii=abs(ii);
  
  int remain=ii%base;

  while(ii > 0){
    aa.push_back(ii%base+48);
    ii=(ii-remain)/base;
    remain=ii%base;
  }
  for(ii=aa.size()-1; ii>=0; ii--){
    a.push_back(aa[ii]);
  }
  return;
}


//************************************************************

double determinant(double mat[3][3]){
  return mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])-
    mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])+
    mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);
}

//************************************************************

void inverse(double mat[3][3], double invmat[3][3]){
  double det=determinant(mat);
  invmat[0][0]=(+mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])/det;
  invmat[0][1]=(-mat[0][1]*mat[2][2]+mat[0][2]*mat[2][1])/det;
  invmat[0][2]=(+mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1])/det;
  invmat[1][0]=(-mat[1][0]*mat[2][2]+mat[1][2]*mat[2][0])/det;
  invmat[1][1]=(+mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0])/det;
  invmat[1][2]=(-mat[0][0]*mat[1][2]+mat[0][2]*mat[1][0])/det;
  invmat[2][0]=(+mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])/det;
  invmat[2][1]=(-mat[0][0]*mat[2][1]+mat[0][1]*mat[2][0])/det;
  invmat[2][2]=(+mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])/det;
  return;
}

//************************************************************
//mat3=mat1*mat2
void matrix_mult(double mat1[3][3], double mat2[3][3], double mat3[3][3]){
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      mat3[i][j]=0.0;
      for(int k=0; k<3; k++){
	mat3[i][j]=mat3[i][j]+mat1[i][k]*mat2[k][j];
      }
    }
  }
}


//************************************************************
//Added by John
//Given vector vec1, find a perpendicular vector vec2
void get_perp(double vec1[3], double vec2[3]){
  
  for(int i=0; i<3; i++){
    vec2[i]=0.0;
  }
  for(int i=0; i<3; i++){
    if(abs(vec1[i])<tol){ 
      vec2[i]=1.0;
      return;
    }
  }

  vec2[0]=vec2[1]=1.0;
  vec2[2]=-(vec1[0]+vec1[1])/vec1[2];
  normalize(vec2,1.0);
  return;
}
//\End Addition

//************************************************************
//Added by John
//Given vectors vec1 and vec2, find perpendicular vector vec3
void get_perp(double vec1[3], double vec2[3], double vec3[3]){
  for(int i=0; i<3; i++)
    vec3[i]=vec1[(i+1)%3]*vec2[(i+2)%3]-vec1[(i+2)%3]*vec2[(i+1)%3];
  if(normalize(vec3, 1.0))
    return;
  else get_perp(vec1, vec3);
}
//\End Addition

//************************************************************
//Added by John
bool normalize(double vec1[3], double length){
  double tmag=0.0;
  for(int i=0; i<3; i++){
    tmag+=vec1[i]*vec1[i];
  }
  tmag=sqrt(tmag);
  if(tmag>tol){
    for(int i=0; i<3; i++){
      if(abs(vec1[i])>tol)
	vec1[i]*=length/tmag;
      else vec1[i]=0.0;
    }
    return true;
  }
  else return false;
}
//\End addition
//************************************************************

void coord_trans_mat(double lat[3][3], double FtoC[3][3], double CtoF[3][3]){
  int i,j;

  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      FtoC[i][j]=lat[j][i];

  inverse(FtoC,CtoF);
}


//************************************************************

bool compare(double mat1[3][3], double mat2[3][3]){
  int i,j,k;

  k=0;
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      if(abs(mat1[i][j]-mat2[i][j]) < tol) k++;

  if(k == 9) return true;
  else return false;
}


//************************************************************

bool compare(double vec1[3], double vec2[3]){
  int i,k;

  k=0;
  for(i=0; i<3; i++)
    if(abs(vec1[i]-vec2[i]) < tol) k++;

  if(k == 3) return true;
  else return false;
}


//************************************************************

bool compare(double vec1[3], double vec2[3], int trans[3]){

  double ftrans[3];
  for(int i=0; i<3; i++){
    ftrans[i]=vec1[i]-vec2[i];
    trans[i]=0;
    //check whether all elements of ftrans[3] are within a tolerance of an integer
    while(abs(ftrans[i]) > tol && ftrans[i] > 0){
      ftrans[i]=ftrans[i]-1.0;
      trans[i]=trans[i]+1;
    }
    while(abs(ftrans[i]) > tol && ftrans[i] < 0){
      ftrans[i]=ftrans[i]+1.0;
      trans[i]=trans[i]-1;
    }
    if(abs(ftrans[i] > tol)) return false;
  }
  return true;
}


//************************************************************

bool compare(char name1[2], char name2[2]){
  for(int i=0; i<2; i++)
    if(name1[i] != name2[i]) return false;
  return true;
}


//************************************************************

bool compare(specie compon1, specie compon2){   // changed by jishnu
  int i,k;

  //k=0;
  //for(i=0; i<2; i++)
  //  if(compon1.name[i] == compon2.name[i]) k++;
  //if(k == 2) return true;
  if(compon1.name.compare(compon2.name) == 0) return true;
  else return false;
}


//************************************************************

bool compare(vector<specie> compon1, vector<specie> compon2){
  int i,j,k;
  int num_suc_maps,ll;

  if(compon1.size() != compon2.size()) return false;

  num_suc_maps=0;
  for(i=0; i<compon1.size(); i++){
    for(j=0; j<compon2.size(); j++){
      if(compare(compon1[i],compon2[j])) num_suc_maps++;
    }
  }
  if(num_suc_maps == compon1.size()) return true;
  else return false;
}


//************************************************************

bool compare(atompos &atom1, atompos &atom2){
  int i,j,k,l;

  k=0;
  l=0;
  for(i=0; i<3; i++)
    if(abs(atom1.fcoord[i]-atom2.fcoord[i]) < tol) k++;
    else return false;

  if(k == 3 && compare(atom1.compon, atom2.compon) &&
     atom1.bit == atom2.bit) return true;
  else return false;

}

//************************************************************

bool compare_just_coordinates(atompos &atom1, atompos &atom2){
  int i,j,k,l;
	
  k=0;
  l=0;
  for(i=0; i<3; i++)
    if(abs(atom1.fcoord[i]-atom2.fcoord[i]) < tol) k++;
    else return false;
	
  //if(k == 3 && compare(atom1.compon, atom2.compon) ) return true;
  if(k == 3) return true;
  else return false;
	
}



//************************************************************

bool compare(atompos atom1, atompos atom2, int trans[3]){

  double ftrans[3];
  for(int i=0; i<3; i++){
    ftrans[i]=atom1.fcoord[i]-atom2.fcoord[i];
    trans[i]=0;
    //check whether all elements of ftrans[3] are within a tolerance of an integer
    while(abs(ftrans[i]) > tol && ftrans[i] > 0){
      ftrans[i]=ftrans[i]-1.0;
      trans[i]=trans[i]+1;
    }
    while(abs(ftrans[i]) > tol && ftrans[i] < 0){
      ftrans[i]=ftrans[i]+1.0;
      trans[i]=trans[i]-1;
    }
    if(abs(ftrans[i] > tol)) return false;
  }
  return true;
}
//************************************************************

bool compare_just_coordinates(cluster &clust1, cluster &clust2){

  ////////////////////////////////////////////////////////////////////////////////
  //added by anton
  if(clust1.point.size() != clust2.point.size()) return false;

  int k=0;
  for(int np1=0; np1<clust1.point.size(); np1++){
    for(int np2=0; np2<clust2.point.size(); np2++){
      if(compare_just_coordinates(clust1.point[np1],clust2.point[np2])) k++;
    }
  }

  if(k == clust1.point.size()) return true;
  else return false;
}



//************************************************************

//************************************************************

bool compare(cluster &clust1, cluster &clust2){

  ////////////////////////////////////////////////////////////////////////////////
  //added by anton
  if(clust1.point.size() != clust2.point.size()) return false;

  int k=0;
  for(int np1=0; np1<clust1.point.size(); np1++){
    for(int np2=0; np2<clust2.point.size(); np2++){
      if(compare(clust1.point[np1],clust2.point[np2])) k++;
    }
  }

  if(k == clust1.point.size()) return true;
  else return false;
}



//************************************************************
////////////////////////////////////////////////////////////////////////////////
//added by anton
bool compare(orbit orb1, orbit orb2){

  if(orb1.equiv.size() != orb2.equiv.size()) return false;

  int k=0; 
  for(int ne1=0; ne1<orb1.equiv.size(); ne1++){
    for(int ne2=0; ne2<orb2.equiv.size(); ne2++){
      if(compare(orb1.equiv[ne1],orb2.equiv[ne2])) k++; 
    }
  }

  if(k == orb1.equiv.size()) return true;
  else return false;
}




////////////////////////////////////////////////////////////////////////////////




//************************************************************

bool compare(vector<double> vec1, vector<double> vec2){
  if(vec1.size() != vec2.size()) return false;
  for(int i=0; i<vec1.size(); i++)
    if(abs(vec1[i]-vec2[i]) > tol) return false;

  return true;

}


//************************************************************

bool compare(concentration conc1, concentration conc2){
  if(conc1.compon.size() != conc2.compon.size()) return false;
  for(int i=0; i<conc1.compon.size(); i++){
    if(!compare(conc1.compon[i],conc2.compon[i])) return false;
    if(!compare(conc1.occup[i],conc2.occup[i]))return false;
  }

  return true;

}


//************************************************************

bool compare(mc_index m1, mc_index m2){
  for(int i=0; i<4; i++)
    if(m1.shift[i] != m2.shift[i]) return false;

  return true;

}





//************************************************************

bool new_mc_index(vector<mc_index> v1, mc_index m2){
  for(int i=0; i<v1.size(); i++)
    if(compare(v1[i],m2)) return false;

  return true;
}


//************************************************************

bool is_integer(double vec[3]){
  int j,k;
  k=0;
  for(j=0; j<3; j++)
    if(abs(vec[j]-ceil(vec[j])) < tol || abs(vec[j]-floor(vec[j])) < tol) k++;

  if(k == 3) return true;
  else return false;
}


//************************************************************

bool is_integer(double mat[3][3]){
  int i,j,k;
  k=0;
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      if(abs(mat[i][j]-ceil(mat[i][j])) < tol || abs(mat[i][j]-floor(mat[i][j])) < tol) k++;

  if(k == 9) return true;
  else return false;
}




//************************************************************

void within(double fcoord[3]){
  int i;
  for(i=0; i<3; i++){
    while(fcoord[i] < 0.0)fcoord[i]=fcoord[i]+1.0;
    while(fcoord[i] >0.99999)fcoord[i]=fcoord[i]-1.0;
  }
  return;
}


//************************************************************

void within(atompos &atom){
  int i;
  for(i=0; i<3; i++){
    while(atom.fcoord[i] < 0.0)atom.fcoord[i]=atom.fcoord[i]+1.0;
    while(atom.fcoord[i] >0.99999)atom.fcoord[i]=atom.fcoord[i]-1.0;
  }
  return;
}

//************************************************************
//added by Ben Swoboda
// used to translate all PRIM coordinates to unit cell

void within(structure &struc){
  int i;
  for(i=0; i<struc.atom.size(); i++){
    within(struc.atom[i]);
  }
  return;
}

//************************************************************
// translates a cluster so that its first point is within the unit cell

void within(cluster &clust){
  int i,np;
  for(i=0; i<3; i++){
    while(clust.point[0].fcoord[i] < 0.0){
      clust.point[0].fcoord[i]=clust.point[0].fcoord[i]+1.0;
      for(np=1; np<clust.point.size(); np++)
	clust.point[np].fcoord[i]=clust.point[np].fcoord[i]+1.0;
    }
    while(clust.point[0].fcoord[i] >0.99999){
      clust.point[0].fcoord[i]=clust.point[0].fcoord[i]-1.0;
      for(np=1; np<clust.point.size(); np++)
	clust.point[np].fcoord[i]=clust.point[np].fcoord[i]-1.0;
    }
  }
  return;
}


//************************************************************
// translates a cluster so that its nth  point is within the unit cell

void within(cluster &clust, int n){
  int i,np;
  for(i=0; i<3; i++){
    while(clust.point[n].fcoord[i] < 0.0){
      for(np=0; np<clust.point.size(); np++)
	clust.point[np].fcoord[i]=clust.point[np].fcoord[i]+1.0;
    }
    while(clust.point[n].fcoord[i] > 0.99999){
      for(np=0; np<clust.point.size(); np++)
	clust.point[np].fcoord[i]=clust.point[np].fcoord[i]-1.0;
    }
  }
  return;
}




//************************************************************
/*
  This function, given a matrix with as rows the cartesian coordinates
  of the unit cell vectors defining a lattice, returns the lengths of the
  unit cell vectors (in latparam) and the angles between the vectors (latparam)
*/
//************************************************************

void latticeparam(double lat[3][3], double latparam[3], double latangle[3])
{
  int i,j;
  double temp;

  //calculate the a,b,c lattice parameters = length of the unit cell vectors

  for(i=0; i<3; i++){
    latparam[i]=0.0;
    for(j=0; j<3; j++)latparam[i]=latparam[i]+lat[i][j]*lat[i][j];
    latparam[i]=sqrt(latparam[i]);
  }

  //calculate the angles between the unit cell vectors

  for(i=0; i<3; i++){
    latangle[i]=0.0;
    for(j=0; j<3; j++) latangle[i]=latangle[i]+lat[(i+1)%3][j]*lat[(i+2)%3][j];
    temp=latangle[i]/(latparam[(i+1)%3]*latparam[(i+2)%3]);

    //make sure numerical errors don't place the arguments outside of the [-1,1] interval

    if((temp-1.0) > 0.0)temp=1.0;
    if((temp+1.0) < 0.0)temp=-1.0;
    latangle[i]=(180.0/3.141592654)*acos(temp);
  }

  return;
}


//************************************************************
/*
  This function, given a matrix with as rows the cartesian coordinates
  of the unit cell vectors defining a lattice, returns the lengths of the
  unit cell vectors (in latparam) and the angles between the vectors (latparam)
  It also determines which vector is largest, smallest and in between in length.
*/
//************************************************************

void latticeparam(double lat[3][3], double latparam[3], double latangle[3], int permut[3])
{
  int i,j;
  double temp;

  //calculate the a,b,c lattice parameters = length of the unit cell vectors

  for(i=0; i<3; i++){
    latparam[i]=0.0;
    for(j=0; j<3; j++)latparam[i]=latparam[i]+lat[i][j]*lat[i][j];
    latparam[i]=sqrt(latparam[i]);
  }

  //calculate the angles between the unit cell vectors

  for(i=0; i<3; i++){
    latangle[i]=0.0;
    for(j=0; j<3; j++) latangle[i]=latangle[i]+lat[(i+1)%3][j]*lat[(i+2)%3][j];
    temp=latangle[i]/(latparam[(i+1)%3]*latparam[(i+2)%3]);

    //make sure numerical errors don't place the arguments outside of the [-1,1] interval

    if((temp-1.0) > 0.0)temp=1.0;
    if((temp+1.0) < 0.0)temp=-1.0;
    latangle[i]=(180.0/3.141592654)*acos(temp);
  }

  int imin,imax,imid;
  double min,max;

  max=min=latparam[0];
  imax=imin=0;

  for(i=0; i<3; i++){
    if(max <= latparam[i]){
      max=latparam[i];
      imax=i;
    }
    if(min > latparam[i]){
      min=latparam[i];
      imin=i;
    }
  }

  for(i=0; i<3; i++) if(i != imin && i !=imax)imid=i;

  //if all lattice parameters are equal length, numerical noise may cause imin=imax

  if(imin == imax)
    for(i=0; i<3; i++) if(i != imin && i !=imid)imax=i;

  permut[0]=imin;
  permut[1]=imid;
  permut[2]=imax;

  return;

}



//************************************************************
/*
  This function, given a lattice with vectors lat[3][3], finds the
  dimensions along the unit cell vectors such that a sphere of given radius
  fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
  centered at the origin.

  The algorithm works by getting the normal (e.g. n1) to each pair of lattice
  vectors (e.g. a2, a3), scaling this normal to have length radius and
  then projecting this normal parallel to the a2,a3 plane onto the
  remaining lattice vector a1. This will tell us the number of a1 vectors
  needed to make a grid to encompass the sphere.
*/
//************************************************************

void lat_dimension(double lat[3][3], double radius, int dim[3]){
  int i,j,k;
  double inv_lat[3][3],normals[3][3],length[3];
  double frac_normals[3][3];

  //get the normals to pairs of lattice vectors of length radius

  for(i=0; i<3; i++){
    for(j=0; j<3; j++)normals[i][j]=lat[(i+1)%3][(j+1)%3]*lat[(i+2)%3][(j+2)%3]-
			lat[(i+1)%3][(j+2)%3]*lat[(i+2)%3][(j+1)%3];

    length[i]=0;
    for(j=0; j<3; j++)
      length[i]=length[i]+normals[i][j]*normals[i][j];
    length[i]=sqrt(length[i]);

    for(j=0; j<3; j++)normals[i][j]=radius*normals[i][j]/length[i];

  }

  //get the normals in the coordinates system of the lattice vectors

  inverse(lat,inv_lat);


  for(i=0; i<3; i++){
    for(j=0; j<3; j++){
      frac_normals[i][j]=0;
      for(k=0; k<3; k++)
	frac_normals[i][j]=frac_normals[i][j]+inv_lat[k][j]*normals[i][k];
    }
  }


  //the diagonals of frac_normal contain the dimensions of the lattice grid that
  //encompasses a sphere of radius = radius

  for(i=0; i<3; i++) dim[i]=(int)ceil(abs(frac_normals[i][i]));

  return;
}



//************************************************************

void conv_AtoB(double AtoB[3][3], double Acoord[3], double Bcoord[3]){
  int i,j;

  for(i=0; i<3; i++){
    Bcoord[i]=0.0;
    for(j=0; j<3 ; j++) Bcoord[i]=Bcoord[i]+AtoB[i][j]*Acoord[j];
  }
}



//************************************************************

double distance(atompos atom1,atompos atom2){
  double dist=0.0;
  for(int i=0; i<3; i++){
    dist=dist+(atom1.ccoord[i]-atom2.ccoord[i])*(atom1.ccoord[i]-atom2.ccoord[i]);
  }
  dist=sqrt(dist);
  return dist;
}


//************************************************************
//this routine starts from the first cluster in the orbit and
//generates all equivalent clusters by applying the factor_group
//symmetry operations

void get_equiv(orbit &orb, vector<sym_op> &op){
  int fg;
  cluster tclust1;


  if(orb.equiv.size() == 0){
    cout << "No cluster present \n";
    exit(1);
  }

  tclust1=orb.equiv[0];

  orb.equiv.clear();

  for(fg=0; fg < op.size(); fg++){
    cluster tclust2;
    tclust2=tclust1.apply_sym(op[fg]);
    within(tclust2);
    tclust2.get_cart(op[fg].FtoC);

    if(new_clust(tclust2,orb)){
      orb.equiv.push_back(tclust2);
    }

  }
}


//************************************************************
//checks to see whether a cluster clust already belongs to an orbit of clusters

bool new_clust(cluster clust, orbit &orb){
  int np,ne;


  if(orb.equiv.size() == 0) return true;

  if(clust.point.size() != orb.equiv[0].point.size()){
    cout << " comparing clusters with different number of points \n";
    cout << " the new_clust function was not meant for that \n";
    exit(1);
  };

  for(np=0; np<clust.point.size(); np++){
    cluster tclust;
    tclust=clust;
    within(tclust,np);
    for(ne=0; ne<orb.equiv.size(); ne++){
      if(compare(tclust,orb.equiv[ne])) return false;
    }
  }
  return true;
}



//************************************************************

bool new_clust(cluster clust, vector<orbit> &orbvec){
  int nc,non_match;

  if(orbvec.size() == 0) return true;

  if(clust.point.size() != orbvec[0].equiv[0].point.size()){
    cout << " comparing clusters with different number of points \n";
    cout << " the new_clust function was not meant for that \n";
    exit(1);
  };


  non_match=0;
  for(nc=0; nc<orbvec.size(); nc++){
    if(abs(orbvec[nc].equiv[0].max_leng - clust.max_leng) < 0.0001 &&
       abs(orbvec[nc].equiv[0].min_leng - clust.min_leng) < 0.0001){
      if(new_clust(clust,orbvec[nc]))non_match++;
    }
    else non_match++;
  }

  if(non_match == orbvec.size())return true;
  else return false;

}



//************************************************************
//this routine starts from the first cluster in the orbit and
//generates all equivalent clusters by applying the site_point_group
//symmetry operations

void get_loc_equiv(orbit &orb, vector<sym_op> &op){
  int g;
  cluster tclust1;


  if(orb.equiv.size() == 0){
    cout << "No cluster present \n";
    exit(1);
  }

  tclust1=orb.equiv[0];


  orb.equiv.clear();

  for(g=0; g < op.size(); g++){
    cluster tclust2;
    tclust2=tclust1.apply_sym(op[g]);
    tclust2.get_cart(op[g].FtoC);

    if(new_loc_clust(tclust2,orb)){
      orb.equiv.push_back(tclust2);
    }
  }
}


//************************************************************
//checks to see whether a cluster clust already belongs to an orbit of clusters

bool new_loc_clust(cluster clust, orbit orb){
  int np,ne;


  if(orb.equiv.size() == 0) return true;

  if(clust.point.size() != orb.equiv[0].point.size()){
    cout << " comparing clusters with different number of points \n";
    cout << " the new_clust function was not meant for that \n";
    exit(1);
  };

  for(ne=0; ne<orb.equiv.size(); ne++){
    if(compare(clust,orb.equiv[ne])) return false;
  }
  return true;
}



//************************************************************

bool new_loc_clust(cluster clust, vector<orbit> orbvec){
  int nc,non_match;


  if(orbvec.size() == 0) return true;

  if(clust.point.size() != orbvec[0].equiv[0].point.size()){
    cout << " comparing clusters with different number of points \n";
    cout << " the new_clust function was not meant for that \n";
    exit(1);
  };


  non_match=0;
  for(nc=0; nc<orbvec.size(); nc++){
    if(abs(orbvec[nc].equiv[0].max_leng-clust.max_leng) < 0.0001 &&
       abs(orbvec[nc].equiv[0].min_leng-clust.min_leng) < 0.0001){
      if(new_loc_clust(clust,orbvec[nc]))non_match++;
    }
    else non_match++;
  }

  if(non_match == orbvec.size())return true;
  else return false;

}


//************************************************************

void structure::read_species(){
	
  double temp;
  string tstring;
  vector <string> names;
  vector <double> masses;
  vector <double> magmoms;
  vector <double> Us;
  vector <double> Js;
	
  ifstream in;
  if(scandirectory(".","species")) in.open("species");
  else if(scandirectory(".","SPECIES")) in.open("SPECIES");
  else{
    cout << "No SPECIES file in the current directory \n";
    return;
  }
  if(!in){
    cout << "cannot open species\n";
    return;
  }
	
  while(tstring != "mass") {
    in >> tstring;
    if (tstring != "mass") names.push_back(tstring);
  } 
	
  for(int i = 0; i < names.size(); i++) {        
    in >>  temp;
    masses.push_back(temp); 
  }   
  in >> tstring;    
  for(int i = 0; i < names.size(); i++) {        
    in >>  temp;
    magmoms.push_back(temp); 
  }   
  in >> tstring;    
  for(int i = 0; i < names.size(); i++) {        
    in >>  temp;
    Us.push_back(temp); 
  }   
  in >> tstring;    
  for(int i = 0; i < names.size(); i++) {        
    in >>  temp;
    Js.push_back(temp); 
  }   
	
	
  for(int i=0;i<atom.size();i++) {
    for(int j=0;j<atom[i].compon.size();j++) {
      if(!(atom[i].compon[j].name.compare("Va") == 0)) {
	for(int k=0;k<names.size();k++) {
	  if(atom[i].compon[j].name.compare(names[k]) == 0) {
	    atom[i].compon[j].mass = masses[k];
	    atom[i].compon[j].magmom = magmoms[k];
	    atom[i].compon[j].U = Us[k];
	    atom[i].compon[j].J = Js[k];
	  }
	}
      }
    }
  }
	
	
	
  in.close();
	
  return;
}


//************************************************************

//************************************************************

void read_cspecs(vector<double> &max_radius){
  char buff[200];
  int dummy;
  double radius;
	
  ifstream in;
  if(scandirectory(".","cspecs")) in.open("cspecs");
  else if(scandirectory(".","CSPECS")) in.open("CSPECS");
  else cout << "No CSPECS file in the current directory \n";
  if(!in){
    cout << "cannot open cspecs\n";
    return;
  }
	
  radius=0.0;
  max_radius.push_back(radius);
  max_radius.push_back(radius);
	
  in.getline(buff,199);
  in.getline(buff,199);
  while(in >> dummy >> radius)
    max_radius.push_back(radius);
	
  in.close();
	
  return;
}


//************************************************************

void write_clust(multiplet clustiplet, string out_file){
  int num_clust=0;
  for(int i=1; i<clustiplet.orb.size(); i++){
    num_clust=num_clust+clustiplet.orb[i].size();
  }

  ofstream out;
  out.open(out_file.c_str());
  if(!out){
    cout << "cannot open " << out_file << "\n";
    return;
  }

  out << num_clust << "\n";
  num_clust=0;
  for(int i=1; i<clustiplet.orb.size(); i++){
    for(int j=0; j<clustiplet.orb[i].size(); j++){
      num_clust++;
      out << num_clust << "  " << clustiplet.orb[i][j].equiv[0].point.size() << "  "
	  << clustiplet.orb[i][j].equiv.size() << "  0  max length "
	  << clustiplet.orb[i][j].equiv[0].max_leng << "\n";
      clustiplet.orb[i][j].equiv[0].print(out);
    }
  }
  out.close();
}


//************************************************************

void write_fclust(multiplet clustiplet, string out_file){
  int num_clust=0;
  for(int i=1; i<clustiplet.orb.size(); i++){
    num_clust=num_clust+clustiplet.orb[i].size();
  }

  ofstream out;
  out.open(out_file.c_str());
  if(!out){
    cout << "cannot open " << out_file <<"\n";
    return;
  }

  out << num_clust << "\n";
  num_clust=0;
  for(int i=1; i<clustiplet.orb.size(); i++){
    for(int j=0; j<clustiplet.orb[i].size(); j++){
      num_clust++;
      out << " Orbit number " << num_clust << "\n";
      out << num_clust << "  " << clustiplet.orb[i][j].equiv[0].point.size() << "  "
	  << clustiplet.orb[i][j].equiv.size() << "  0  max length "
	  << clustiplet.orb[i][j].equiv[0].max_leng << "\n";
      clustiplet.orb[i][j].print(out);
    }
  }
  out.close();
}


//************************************************************

void write_scel(vector<structure> suplat){

  ofstream out;
  out.open("SCEL");
  if(!out){
    cout << "cannot open SCEL \n";
    return;
  }

  out << suplat.size() << "\n";
  for(int i=0; i<suplat.size(); i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	out.precision(5);out.width(5);
	out << suplat[i].slat[j][k] << " ";
      }
      out << "  ";
    }
    out << " volume = " << determinant(suplat[i].slat) << "\n";
  }
  out.close();
}


//************************************************************

void read_scel(vector<structure> &suplat, structure prim){
  suplat.clear();
  int num_scel;
  char buff[200];
  ifstream in;
  in.open("SCEL");
  if(!in){
    cout << "cannot open SCEL \n";
    return;
  }
  in >> num_scel;
  for(int n=0; n<num_scel; n++){
    structure tsuplat;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	in >> tsuplat.slat[i][j];
      }
    }
    in.getline(buff,199);
    tsuplat.scale=prim.scale;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	tsuplat.lat[i][j]=0.0;
	for(int k=0; k<3; k++){
	  tsuplat.lat[i][j]=tsuplat.lat[i][j]+tsuplat.slat[i][k]*prim.lat[k][j];
	}
      }
    }
    suplat.push_back(tsuplat);
  }
  in.close();
}



//************************************************************

bool update_bit(vector<int> max_bit, vector<int> &bit, int &last){

  //given a bit vector, this routine updates the bit by 1
  //if the maximum bit has been reached it returns false
  //last needs to be initialized as zero and bit needs to also be initialized as all zeros

  int bs=bit.size();

  if(last == 0){
    bit[0]++;
    for(int i=0; i<bs-1; i++){
      if(bit[i] !=0 && bit[i]%max_bit[i] == 0){
	bit[i+1]++;
	bit[i]=0;
      }
    }
    if(bit[bs-1] !=0 && bit[bs-1]%max_bit[bs-1] == 0){
      last=last+1;
      bit[bs-1]=0;
    }
    return true;
  }
  else{
    return false;
  }
}


//************************************************************

double ran0(int &idum){
  int IA=16807;
  int IM=2147483647;
  int IQ=127773;
  int IR=2836;
  int MASK=123459876;
  double AM=1.0/IM;

  //minimal random number generator of Park and Miller
  //returns uniform random deviate between 0.0 and 1.0
  //set or rest idum to any  integer value (except the
  //unlikely value MASK) to initialize the sequence: idum must
  //not be altered between calls for successive deviates
  //in a sequence

  int k;
  idum=idum^MASK;     //XOR the two integers
  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if(idum < 0) idum=idum+IM;
  double ran=AM*idum;
  idum=idum^MASK;
  return ran;

}


//************************************************************

//creates the shift vectors

void get_shift(atompos &atom, vector<atompos> basis){

  //first bring atom within the primitive unit cell
  //and document the translations needed for that

  atompos tatom=atom;

  for(int i=0; i<3; i++){
    atom.shift[i]=0;
    while(tatom.fcoord[i] < 0.0){
      atom.shift[i]=atom.shift[i]-1;
      tatom.fcoord[i]=tatom.fcoord[i]+1.0;
    }
    while(tatom.fcoord[i] > 0.99999){
      atom.shift[i]=atom.shift[i]+1;
      tatom.fcoord[i]=tatom.fcoord[i]-1.0;
    }
  }

  //then compare with all basis points and determine which one matches

  int nb=0;
  for(int na=0; na<basis.size(); na++){
    if(basis[na].compon.size() >= 2){
      nb++;
      if(compare(basis[na].fcoord,tatom.fcoord)){
        atom.shift[3]=nb-1;
        break;
      }
    }
  }

}




//************************************************************

// scans the directory 'dirname' to see whether 'filename' resides there
// by Qingchuan Xu

bool scandirectory(string dirname, string filename)
{
  bool exist=false;

  char ch;
  int  n;
  double e0;
  DIR *dir;
  struct dirent *entry;
  char path[PATH_MAX];
  dir = opendir(dirname.c_str());
  if (dir == NULL)
    {
      perror("Error opendir()");
      exit (0);
    }

  while ((entry = readdir(dir)) != NULL && !exist)
    {
      if (entry->d_type == DT_DIR)
        {
          if (strcmp(entry->d_name, ".")&& strcmp(entry->d_name, ".."))
            {
	      if(!fnmatch(filename.c_str(), entry->d_name, 0)) exist=true;
            }
        }
      else if (entry->d_type == DT_REG)
        {
          if (!fnmatch(filename.c_str(), entry->d_name, 0)) exist=true;
        }
    }
  closedir(dir);
  return exist;

}


//************************************************************

// reads the OSZICAR in dirname and extracts the final energy
// by Qingchuan Xu

bool read_oszicar(string dirname, double& e0)
{
  static bool exist=false;
  char ch;
  int  n;
  ifstream readfrom;
  DIR *dir;
  struct dirent *entry;
  char path[PATH_MAX];
  char path1[PATH_MAX];
  dir = opendir(dirname.c_str());
  if (dir == NULL)
    {
      perror("Error opendir()");
      return 0;
    }

  //stop_flag is used to override original recursive behavior (read OSZICAR of lowest level directory).
  //remove all reference to restore previous behavior
  bool stop_flag=true;  
  while ((entry = readdir(dir)) != NULL && stop_flag)
    {
      if (entry->d_type == DT_DIR && !stop_flag)
        {
          if (strcmp(entry->d_name, ".")&& strcmp(entry->d_name, ".."))
            {
              snprintf(path, (size_t) PATH_MAX, "%s/%s", dirname.c_str(),entry->d_name);
              read_oszicar(path, e0);
            }
        }
      else if (entry->d_type == DT_REG || stop_flag)
        {
          if (!fnmatch("OSZICAR", entry->d_name, 0))
            {
	      stop_flag=false;
              exist = true;
              snprintf(path1, (size_t) PATH_MAX, "%s/%s", dirname.c_str(), entry->d_name);
              readfrom.open(path1);
	      do
		{
		  readfrom.get(ch);
		  if (ch=='F')
		    {n=0;
		      do{readfrom.get(ch);
			if(ch=='=')
			  n=n+1;
		      }
		      while (n<2);
		      readfrom>>e0;
		    }
		}
	      while (!readfrom.eof());
	      readfrom.close();

            }
        }
    }
  closedir(dir);
  return exist;
}
// *******************************************************************
// *******************************************************************
bool read_oszicar(string dirname, double& e0, int& count)
{  
  static bool exist=false;
  char ch;
  int  n;
  ifstream readfrom;
  DIR *dir;
  struct dirent *entry;
  char path[PATH_MAX];
  char path1[PATH_MAX];
  dir = opendir(dirname.c_str());
  if (dir == NULL)
    {
      perror("Error opendir()");
      return 0;
    }
  bool stop_flag=true;  //This is just being used to override original recursive behavior.  remove all reference to restore previous behavior
  while ((entry = readdir(dir)) != NULL && stop_flag)
    {
      if (entry->d_type == DT_DIR && !stop_flag)
        {
	  if (strcmp(entry->d_name, ".")&& strcmp(entry->d_name, ".."))
            {
	      snprintf(path, (size_t) PATH_MAX, "%s/%s", dirname.c_str(),entry->d_name);
	      read_oszicar(path, e0, count);
            }
        }
      else if (entry->d_type == DT_REG || stop_flag)
        {
	  if (!fnmatch("OSZICAR", entry->d_name, 0))
            {
	      stop_flag=false;
	      exist = true;
	      snprintf(path1, (size_t) PATH_MAX, "%s/%s", dirname.c_str(), entry->d_name);
	      readfrom.open(path1);
	      count=0;	  
	      do
		{
		  readfrom.get(ch);
		  if (ch=='F')		  
		    {count++;
		      n=0;
		      do{readfrom.get(ch);
			if(ch=='=')
			  n=n+1;
		      }
		      while (n<2);
		      readfrom>>e0;
		    }
		}
	      while (!readfrom.eof());
	      readfrom.close();
				
            }
        }
    }
  closedir(dir);
  return exist;
}
// *******************************************************************

bool read_mc_input(string cond_file, int &n_pass, int &n_equil_pass, int &nx, int &ny, int &nz, double &Tinit, double &Tmin, double &Tmax, double &Tinc, chempot &muinit, chempot &mu_min, chempot &mu_max, vector<chempot> &muinc, int &xyz_step, int &corr_flag, int &temp_chem){

  ifstream in;
  in.open(cond_file.c_str());
  if(!in){
    cout << "cannot open " << cond_file << "\n";
    return false;
  }
  
  char buff[200];
  in >> n_pass;
  in.getline(buff,199);
  in >> n_equil_pass;
  in.getline(buff,199);
  in >> nx;
  in >> ny;
  in >> nz;
  in.getline(buff,199);
  in >> Tinit;
  in.getline(buff,199);
  in >> Tmin;
  in.getline(buff,199);
  in >> Tmax;
  in.getline(buff,199);
  in >> Tinc;
  in.getline(buff,199);

 
  in.getline(buff,199);
  int buffcount=0;
  string tspec;
  double tdouble;
  while(buff[buffcount]!='!'){

    
    while((buff[buffcount]>='A'&&buff[buffcount]<='Z')||(buff[buffcount]>='a'&&buff[buffcount]<='z')){
      tspec.push_back(buff[buffcount]);
      buffcount++;
    }
    if(buff[buffcount]>='-'&&buff[buffcount]<='9'){
      string tstring;
      while(buff[buffcount]>='-'&&buff[buffcount]<='9'){
	tstring.push_back(buff[buffcount]);
	buffcount++;
      }
      tdouble=strtod(tstring.c_str(), NULL);

      for(int ii=0; ii<muinit.compon.size(); ii++){
	for(int jj=0; jj<muinit.compon[ii].size(); jj++){
	  if(!muinit.compon[ii][jj].name.compare(tspec)){
	    //	    muinit.m[ii].erase(jj);
	    muinit.m[ii][jj]=tdouble;
	  }
	}
      }
      tspec.clear();
    }
    else buffcount++;
  }
  buffcount=0;

  in.getline(buff,199);
  while(buff[buffcount]!='!'){
    
    while((buff[buffcount]>='A'&&buff[buffcount]<='Z')||(buff[buffcount]>='a'&&buff[buffcount]<='z')){
      tspec.push_back(buff[buffcount]);
      buffcount++;
    }
    if(buff[buffcount]>='-'&&buff[buffcount]<='9'){
      string tstring;
      while(buff[buffcount]>='-'&&buff[buffcount]<='9'){
	tstring.push_back(buff[buffcount]);
	buffcount++;
      }
      tdouble=strtod(tstring.c_str(), NULL);
      for(int ii=0; ii<mu_min.compon.size(); ii++){
	for(int jj=0; jj<mu_min.compon[ii].size(); jj++){
	  if(!mu_min.compon[ii][jj].name.compare(tspec)){
	    //	    mu_min.m[ii].erase(jj);
	    mu_min.m[ii][jj]=tdouble;
	  }
	}
      }
      tspec.clear();
    }
    else buffcount++;
  } 
  buffcount=0;

  in.getline(buff,199);
  
  while(buff[buffcount]!='!'){
    
    while((buff[buffcount]>='A'&&buff[buffcount]<='Z')||(buff[buffcount]>='a'&&buff[buffcount]<='z')){
      tspec.push_back(buff[buffcount]);
      buffcount++;
    }
    if(buff[buffcount]>='-'&&buff[buffcount]<='9'){
      string tstring;
      while(buff[buffcount]>='-'&&buff[buffcount]<='9'){
	tstring.push_back(buff[buffcount]);
	buffcount++;
      }
      tdouble=strtod(tstring.c_str(), NULL);
      for(int ii=0; ii<mu_max.compon.size(); ii++){
	for(int jj=0; jj<mu_max.compon[ii].size(); jj++){
	  if(!mu_max.compon[ii][jj].name.compare(tspec)){
	    //	    mu_max.m[ii].erase(jj);
	    mu_max.m[ii][jj]=tdouble;
	  }
	}
      }
      tspec.clear();
    }
    else buffcount++;
  }

  buffcount=0;

  
  in.getline(buff,199);

  while(buff[buffcount]!='!'){
    
    while((buff[buffcount]>='A'&&buff[buffcount]<='Z')||(buff[buffcount]>='a'&&buff[buffcount]<='z')){
      tspec.push_back(buff[buffcount]);
      buffcount++;
    }
    if(buff[buffcount]>='-'&&buff[buffcount]<='9'){
      string tstring;
      while(buff[buffcount]>='-'&&buff[buffcount]<='9'){
	tstring.push_back(buff[buffcount]);
	buffcount++;
      }
      tdouble=strtod(tstring.c_str(), NULL);
      chempot tmu;
      tmu.initialize(muinit.compon);
      for(int ii=0; ii<tmu.compon.size(); ii++){
	for(int jj=0; jj<tmu.compon[ii].size(); jj++){
	  if(!tmu.compon[ii][jj].name.compare(tspec)){
	    //	    muinc.m[ii].erase(jj);
	    tmu.m[ii][jj]=tdouble;
	    muinc.push_back(tmu);
	  }
	}
      }
      tspec.clear();
    }
    else buffcount++;
  }
  buffcount=0;
  in >> xyz_step;
  in.getline(buff,199);  
  in >> corr_flag;
  in.getline(buff,199);  
  in >> temp_chem;
  in.getline(buff,199);  
  in.close();
  return true;
}


double Monte_Carlo::pointenergy(int i, int j, int k, int b){
  double energy = 0.0;
  int l; 
  if(b == 0){
     l=index(i,j,k,0);
     int p00=mcL[l];
     l=index(i,j,k,2);
     int p10=mcL[l];
     l=index(i+1,j-1,k,2);
     int p20=mcL[l];
     l=index(i+1,j,k,2);
     int p30=mcL[l];
     l=index(i,j-1,k,2);
     int p40=mcL[l];
     l=index(i+1,j,k,3);
     int p50=mcL[l];
     l=index(i+1,j,k-1,3);
     int p60=mcL[l];
     l=index(i,j,k-1,3);
     int p70=mcL[l];
     l=index(i,j,k,3);
     int p80=mcL[l];
     l=index(i-1,j,k,2);
     int p90=mcL[l];
     l=index(i+2,j-1,k,2);
     int p100=mcL[l];
     l=index(i+2,j,k,2);
     int p110=mcL[l];
     l=index(i-1,j-1,k,2);
     int p120=mcL[l];
     l=index(i+2,j,k,3);
     int p130=mcL[l];
     l=index(i+2,j,k-1,3);
     int p140=mcL[l];
     l=index(i-1,j,k-1,3);
     int p150=mcL[l];
     l=index(i-1,j,k,3);
     int p160=mcL[l];
     l=index(i,j-1,k,1);
     int p170=mcL[l];
     l=index(i,j-1,k-1,1);
     int p180=mcL[l];
     l=index(i,j,k+1,0);
     int p190=mcL[l];
     l=index(i,j,k-1,0);
     int p200=mcL[l];
     l=index(i,j,k-1,2);
     int p210=mcL[l];
     l=index(i+1,j-1,k+1,2);
     int p220=mcL[l];
     l=index(i,j,k+1,2);
     int p230=mcL[l];
     l=index(i+1,j-1,k-1,2);
     int p240=mcL[l];
     l=index(i+1,j,k-1,2);
     int p250=mcL[l];
     l=index(i,j-1,k+1,2);
     int p260=mcL[l];
     l=index(i+1,j,k+1,2);
     int p270=mcL[l];
     l=index(i,j-1,k-1,2);
     int p280=mcL[l];
     l=index(i+3,j,k+1,0);
     int p290=mcL[l];
     l=index(i-3,j,k-1,0);
     int p300=mcL[l];
     l=index(i+3,j,k-1,0);
     int p310=mcL[l];
     l=index(i-3,j,k+1,0);
     int p320=mcL[l];
     l=index(i-1,j,k,0);
     int p330=mcL[l];
     l=index(i+1,j,k,0);
     int p340=mcL[l];
     l=index(i-1,j,k,1);
     int p350=mcL[l];
     l=index(i+1,j,k,1);
     int p360=mcL[l];
     l=index(i-1,j,k-1,1);
     int p370=mcL[l];
     l=index(i+1,j,k-1,1);
     int p380=mcL[l];
     l=index(i+2,j,k,0);
     int p390=mcL[l];
     l=index(i-2,j,k,0);
     int p400=mcL[l];
     l=index(i,j,k-1,1);
     int p410=mcL[l];
     l=index(i,j,k,1);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0300303*(p00*p10+p20*p00+p00*p30+p40*p00)-0.00468447*(p00*p50+p00*p60+p70*p00+p80*p00)-0.00568721*(p00*p90+p100*p00+p00*p110+p120*p00)-0.00344484*(p00*p130+p00*p140+p150*p00+p160*p00)-0.0076456*(p00*p170+p00*p180)+0.00666985*(p00*p190+p200*p00)-0.00379933*(p00*p210+p220*p00+p00*p230+p240*p00+p00*p250+p260*p00+p00*p270+p280*p00)+0.000773675*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0417653*(p00*p330*p10+p340*p00*p30+p20*p40*p00)+0.00188433*(p00*p350*p360+p00*p370*p380+p380*p00*p390+p370*p400*p00+p360*p00*p390+p350*p400*p00)-0.00397897*(p00*p340*p410*p380+p330*p00*p370*p410+p00*p340*p420*p360+p330*p00*p350*p420);
     return energy;
  }


  if(b == 1){
     l=index(i,j,k,1);
     int p00=mcL[l];
     l=index(i,j,k,3);
     int p10=mcL[l];
     l=index(i+1,j+1,k,3);
     int p20=mcL[l];
     l=index(i+1,j,k,3);
     int p30=mcL[l];
     l=index(i,j+1,k,3);
     int p40=mcL[l];
     l=index(i,j,k,2);
     int p50=mcL[l];
     l=index(i,j,k+1,2);
     int p60=mcL[l];
     l=index(i+1,j,k+1,2);
     int p70=mcL[l];
     l=index(i+1,j,k,2);
     int p80=mcL[l];
     l=index(i-1,j,k,3);
     int p90=mcL[l];
     l=index(i+2,j+1,k,3);
     int p100=mcL[l];
     l=index(i+2,j,k,3);
     int p110=mcL[l];
     l=index(i-1,j+1,k,3);
     int p120=mcL[l];
     l=index(i-1,j,k,2);
     int p130=mcL[l];
     l=index(i-1,j,k+1,2);
     int p140=mcL[l];
     l=index(i+2,j,k+1,2);
     int p150=mcL[l];
     l=index(i+2,j,k,2);
     int p160=mcL[l];
     l=index(i,j+1,k,0);
     int p170=mcL[l];
     l=index(i,j+1,k+1,0);
     int p180=mcL[l];
     l=index(i,j,k+1,1);
     int p190=mcL[l];
     l=index(i,j,k-1,1);
     int p200=mcL[l];
     l=index(i,j,k-1,3);
     int p210=mcL[l];
     l=index(i+1,j+1,k+1,3);
     int p220=mcL[l];
     l=index(i,j,k+1,3);
     int p230=mcL[l];
     l=index(i+1,j+1,k-1,3);
     int p240=mcL[l];
     l=index(i+1,j,k-1,3);
     int p250=mcL[l];
     l=index(i,j+1,k+1,3);
     int p260=mcL[l];
     l=index(i+1,j,k+1,3);
     int p270=mcL[l];
     l=index(i,j+1,k-1,3);
     int p280=mcL[l];
     l=index(i+3,j,k+1,1);
     int p290=mcL[l];
     l=index(i-3,j,k-1,1);
     int p300=mcL[l];
     l=index(i+3,j,k-1,1);
     int p310=mcL[l];
     l=index(i-3,j,k+1,1);
     int p320=mcL[l];
     l=index(i-1,j,k,1);
     int p330=mcL[l];
     l=index(i+1,j,k,1);
     int p340=mcL[l];
     l=index(i+1,j,k,0);
     int p350=mcL[l];
     l=index(i+2,j,k,1);
     int p360=mcL[l];
     l=index(i-1,j,k,0);
     int p370=mcL[l];
     l=index(i-2,j,k,1);
     int p380=mcL[l];
     l=index(i+1,j,k+1,0);
     int p390=mcL[l];
     l=index(i-1,j,k+1,0);
     int p400=mcL[l];
     l=index(i,j,k+1,0);
     int p410=mcL[l];
     l=index(i,j,k,0);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0300303*(p00*p10+p20*p00+p00*p30+p40*p00)-0.00468447*(p50*p00+p60*p00+p00*p70+p00*p80)-0.00568721*(p00*p90+p100*p00+p00*p110+p120*p00)-0.00344484*(p130*p00+p140*p00+p00*p150+p00*p160)-0.0076456*(p170*p00+p180*p00)+0.00666985*(p00*p190+p200*p00)-0.00379933*(p00*p210+p220*p00+p00*p230+p240*p00+p00*p250+p260*p00+p00*p270+p280*p00)+0.000773675*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0417653*(p00*p330*p10+p340*p00*p30+p20*p40*p00)+0.00188433*(p350*p00*p360+p370*p380*p00+p390*p00*p360+p400*p380*p00+p00*p400*p390+p00*p370*p350)-0.00397897*(p410*p390*p00*p340+p400*p410*p330*p00+p420*p350*p00*p340+p370*p420*p330*p00);
     return energy;
  }


  if(b == 2){
     l=index(i,j,k,2);
     int p00=mcL[l];
     l=index(i,j,k,0);
     int p10=mcL[l];
     l=index(i-1,j+1,k,0);
     int p20=mcL[l];
     l=index(i-1,j,k,0);
     int p30=mcL[l];
     l=index(i,j+1,k,0);
     int p40=mcL[l];
     l=index(i,j,k,1);
     int p50=mcL[l];
     l=index(i,j,k-1,1);
     int p60=mcL[l];
     l=index(i-1,j,k-1,1);
     int p70=mcL[l];
     l=index(i-1,j,k,1);
     int p80=mcL[l];
     l=index(i+1,j,k,0);
     int p90=mcL[l];
     l=index(i-2,j+1,k,0);
     int p100=mcL[l];
     l=index(i-2,j,k,0);
     int p110=mcL[l];
     l=index(i+1,j+1,k,0);
     int p120=mcL[l];
     l=index(i+1,j,k,1);
     int p130=mcL[l];
     l=index(i+1,j,k-1,1);
     int p140=mcL[l];
     l=index(i-2,j,k-1,1);
     int p150=mcL[l];
     l=index(i-2,j,k,1);
     int p160=mcL[l];
     l=index(i,j,k,3);
     int p170=mcL[l];
     l=index(i,j,k-1,3);
     int p180=mcL[l];
     l=index(i,j,k+1,2);
     int p190=mcL[l];
     l=index(i,j,k-1,2);
     int p200=mcL[l];
     l=index(i,j,k+1,0);
     int p210=mcL[l];
     l=index(i-1,j+1,k-1,0);
     int p220=mcL[l];
     l=index(i,j,k-1,0);
     int p230=mcL[l];
     l=index(i-1,j+1,k+1,0);
     int p240=mcL[l];
     l=index(i-1,j,k+1,0);
     int p250=mcL[l];
     l=index(i,j+1,k-1,0);
     int p260=mcL[l];
     l=index(i-1,j,k-1,0);
     int p270=mcL[l];
     l=index(i,j+1,k+1,0);
     int p280=mcL[l];
     l=index(i+3,j,k+1,2);
     int p290=mcL[l];
     l=index(i-3,j,k-1,2);
     int p300=mcL[l];
     l=index(i+3,j,k-1,2);
     int p310=mcL[l];
     l=index(i-3,j,k+1,2);
     int p320=mcL[l];
     l=index(i-1,j,k,2);
     int p330=mcL[l];
     l=index(i+1,j,k,2);
     int p340=mcL[l];
     l=index(i-1,j+1,k,3);
     int p350=mcL[l];
     l=index(i+1,j+1,k,3);
     int p360=mcL[l];
     l=index(i-1,j+1,k-1,3);
     int p370=mcL[l];
     l=index(i+1,j+1,k-1,3);
     int p380=mcL[l];
     l=index(i+2,j,k,2);
     int p390=mcL[l];
     l=index(i-2,j,k,2);
     int p400=mcL[l];
     l=index(i,j+1,k-1,3);
     int p410=mcL[l];
     l=index(i,j+1,k,3);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0300303*(p10*p00+p00*p20+p30*p00+p00*p40)-0.00468447*(p00*p50+p00*p60+p70*p00+p80*p00)-0.00568721*(p90*p00+p00*p100+p110*p00+p00*p120)-0.00344484*(p00*p130+p00*p140+p150*p00+p160*p00)-0.0076456*(p00*p170+p00*p180)+0.00666985*(p00*p190+p200*p00)-0.00379933*(p210*p00+p00*p220+p230*p00+p00*p240+p250*p00+p00*p260+p270*p00+p00*p280)+0.000773675*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0417653*(p10*p30*p00+p00*p330*p20+p340*p00*p40)+0.00188433*(p00*p350*p360+p00*p370*p380+p380*p00*p390+p370*p400*p00+p360*p00*p390+p350*p400*p00)-0.00397897*(p00*p340*p410*p380+p330*p00*p370*p410+p00*p340*p420*p360+p330*p00*p350*p420);
     return energy;
  }


  if(b == 3){
     l=index(i,j,k,3);
     int p00=mcL[l];
     l=index(i,j,k,1);
     int p10=mcL[l];
     l=index(i-1,j-1,k,1);
     int p20=mcL[l];
     l=index(i-1,j,k,1);
     int p30=mcL[l];
     l=index(i,j-1,k,1);
     int p40=mcL[l];
     l=index(i-1,j,k,0);
     int p50=mcL[l];
     l=index(i-1,j,k+1,0);
     int p60=mcL[l];
     l=index(i,j,k+1,0);
     int p70=mcL[l];
     l=index(i,j,k,0);
     int p80=mcL[l];
     l=index(i+1,j,k,1);
     int p90=mcL[l];
     l=index(i-2,j-1,k,1);
     int p100=mcL[l];
     l=index(i-2,j,k,1);
     int p110=mcL[l];
     l=index(i+1,j-1,k,1);
     int p120=mcL[l];
     l=index(i-2,j,k,0);
     int p130=mcL[l];
     l=index(i-2,j,k+1,0);
     int p140=mcL[l];
     l=index(i+1,j,k+1,0);
     int p150=mcL[l];
     l=index(i+1,j,k,0);
     int p160=mcL[l];
     l=index(i,j,k,2);
     int p170=mcL[l];
     l=index(i,j,k+1,2);
     int p180=mcL[l];
     l=index(i,j,k+1,3);
     int p190=mcL[l];
     l=index(i,j,k-1,3);
     int p200=mcL[l];
     l=index(i,j,k+1,1);
     int p210=mcL[l];
     l=index(i-1,j-1,k-1,1);
     int p220=mcL[l];
     l=index(i,j,k-1,1);
     int p230=mcL[l];
     l=index(i-1,j-1,k+1,1);
     int p240=mcL[l];
     l=index(i-1,j,k+1,1);
     int p250=mcL[l];
     l=index(i,j-1,k-1,1);
     int p260=mcL[l];
     l=index(i-1,j,k-1,1);
     int p270=mcL[l];
     l=index(i,j-1,k+1,1);
     int p280=mcL[l];
     l=index(i+3,j,k+1,3);
     int p290=mcL[l];
     l=index(i-3,j,k-1,3);
     int p300=mcL[l];
     l=index(i+3,j,k-1,3);
     int p310=mcL[l];
     l=index(i-3,j,k+1,3);
     int p320=mcL[l];
     l=index(i-1,j,k,3);
     int p330=mcL[l];
     l=index(i+1,j,k,3);
     int p340=mcL[l];
     l=index(i+1,j-1,k,2);
     int p350=mcL[l];
     l=index(i+2,j,k,3);
     int p360=mcL[l];
     l=index(i-1,j-1,k,2);
     int p370=mcL[l];
     l=index(i-2,j,k,3);
     int p380=mcL[l];
     l=index(i+1,j-1,k+1,2);
     int p390=mcL[l];
     l=index(i-1,j-1,k+1,2);
     int p400=mcL[l];
     l=index(i,j-1,k+1,2);
     int p410=mcL[l];
     l=index(i,j-1,k,2);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0300303*(p10*p00+p00*p20+p30*p00+p00*p40)-0.00468447*(p50*p00+p60*p00+p00*p70+p00*p80)-0.00568721*(p90*p00+p00*p100+p110*p00+p00*p120)-0.00344484*(p130*p00+p140*p00+p00*p150+p00*p160)-0.0076456*(p170*p00+p180*p00)+0.00666985*(p00*p190+p200*p00)-0.00379933*(p210*p00+p00*p220+p230*p00+p00*p240+p250*p00+p00*p260+p270*p00+p00*p280)+0.000773675*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0417653*(p10*p30*p00+p00*p330*p20+p340*p00*p40)+0.00188433*(p350*p00*p360+p370*p380*p00+p390*p00*p360+p400*p380*p00+p00*p400*p390+p00*p370*p350)-0.00397897*(p410*p390*p00*p340+p400*p410*p330*p00+p420*p350*p00*p340+p370*p420*p330*p00);
     return energy;
  }


}


double Monte_Carlo::normalized_pointenergy(int i, int j, int k, int b){
  double energy = 0.0;
  int l; 
  if(b == 0){
     l=index(i,j,k,0);
     int p00=mcL[l];
     l=index(i,j,k,2);
     int p10=mcL[l];
     l=index(i+1,j-1,k,2);
     int p20=mcL[l];
     l=index(i+1,j,k,2);
     int p30=mcL[l];
     l=index(i,j-1,k,2);
     int p40=mcL[l];
     l=index(i+1,j,k,3);
     int p50=mcL[l];
     l=index(i+1,j,k-1,3);
     int p60=mcL[l];
     l=index(i,j,k-1,3);
     int p70=mcL[l];
     l=index(i,j,k,3);
     int p80=mcL[l];
     l=index(i-1,j,k,2);
     int p90=mcL[l];
     l=index(i+2,j-1,k,2);
     int p100=mcL[l];
     l=index(i+2,j,k,2);
     int p110=mcL[l];
     l=index(i-1,j-1,k,2);
     int p120=mcL[l];
     l=index(i+2,j,k,3);
     int p130=mcL[l];
     l=index(i+2,j,k-1,3);
     int p140=mcL[l];
     l=index(i-1,j,k-1,3);
     int p150=mcL[l];
     l=index(i-1,j,k,3);
     int p160=mcL[l];
     l=index(i,j-1,k,1);
     int p170=mcL[l];
     l=index(i,j-1,k-1,1);
     int p180=mcL[l];
     l=index(i,j,k+1,0);
     int p190=mcL[l];
     l=index(i,j,k-1,0);
     int p200=mcL[l];
     l=index(i,j,k-1,2);
     int p210=mcL[l];
     l=index(i+1,j-1,k+1,2);
     int p220=mcL[l];
     l=index(i,j,k+1,2);
     int p230=mcL[l];
     l=index(i+1,j-1,k-1,2);
     int p240=mcL[l];
     l=index(i+1,j,k-1,2);
     int p250=mcL[l];
     l=index(i,j-1,k+1,2);
     int p260=mcL[l];
     l=index(i+1,j,k+1,2);
     int p270=mcL[l];
     l=index(i,j-1,k-1,2);
     int p280=mcL[l];
     l=index(i+3,j,k+1,0);
     int p290=mcL[l];
     l=index(i-3,j,k-1,0);
     int p300=mcL[l];
     l=index(i+3,j,k-1,0);
     int p310=mcL[l];
     l=index(i-3,j,k+1,0);
     int p320=mcL[l];
     l=index(i-1,j,k,0);
     int p330=mcL[l];
     l=index(i+1,j,k,0);
     int p340=mcL[l];
     l=index(i-1,j,k,1);
     int p350=mcL[l];
     l=index(i+1,j,k,1);
     int p360=mcL[l];
     l=index(i-1,j,k-1,1);
     int p370=mcL[l];
     l=index(i+1,j,k-1,1);
     int p380=mcL[l];
     l=index(i+2,j,k,0);
     int p390=mcL[l];
     l=index(i-2,j,k,0);
     int p400=mcL[l];
     l=index(i,j,k-1,1);
     int p410=mcL[l];
     l=index(i,j,k,1);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0150152*(p00*p10+p20*p00+p00*p30+p40*p00)-0.00234223*(p00*p50+p00*p60+p70*p00+p80*p00)-0.0028436*(p00*p90+p100*p00+p00*p110+p120*p00)-0.00172242*(p00*p130+p00*p140+p150*p00+p160*p00)-0.0038228*(p00*p170+p00*p180)+0.00333493*(p00*p190+p200*p00)-0.00189967*(p00*p210+p220*p00+p00*p230+p240*p00+p00*p250+p260*p00+p00*p270+p280*p00)+0.000386838*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0139218*(p00*p330*p10+p340*p00*p30+p20*p40*p00)+0.000628108*(p00*p350*p360+p00*p370*p380+p380*p00*p390+p370*p400*p00+p360*p00*p390+p350*p400*p00)-0.000994741*(p00*p340*p410*p380+p330*p00*p370*p410+p00*p340*p420*p360+p330*p00*p350*p420);
     return energy;
  }


  if(b == 1){
     l=index(i,j,k,1);
     int p00=mcL[l];
     l=index(i,j,k,3);
     int p10=mcL[l];
     l=index(i+1,j+1,k,3);
     int p20=mcL[l];
     l=index(i+1,j,k,3);
     int p30=mcL[l];
     l=index(i,j+1,k,3);
     int p40=mcL[l];
     l=index(i,j,k,2);
     int p50=mcL[l];
     l=index(i,j,k+1,2);
     int p60=mcL[l];
     l=index(i+1,j,k+1,2);
     int p70=mcL[l];
     l=index(i+1,j,k,2);
     int p80=mcL[l];
     l=index(i-1,j,k,3);
     int p90=mcL[l];
     l=index(i+2,j+1,k,3);
     int p100=mcL[l];
     l=index(i+2,j,k,3);
     int p110=mcL[l];
     l=index(i-1,j+1,k,3);
     int p120=mcL[l];
     l=index(i-1,j,k,2);
     int p130=mcL[l];
     l=index(i-1,j,k+1,2);
     int p140=mcL[l];
     l=index(i+2,j,k+1,2);
     int p150=mcL[l];
     l=index(i+2,j,k,2);
     int p160=mcL[l];
     l=index(i,j+1,k,0);
     int p170=mcL[l];
     l=index(i,j+1,k+1,0);
     int p180=mcL[l];
     l=index(i,j,k+1,1);
     int p190=mcL[l];
     l=index(i,j,k-1,1);
     int p200=mcL[l];
     l=index(i,j,k-1,3);
     int p210=mcL[l];
     l=index(i+1,j+1,k+1,3);
     int p220=mcL[l];
     l=index(i,j,k+1,3);
     int p230=mcL[l];
     l=index(i+1,j+1,k-1,3);
     int p240=mcL[l];
     l=index(i+1,j,k-1,3);
     int p250=mcL[l];
     l=index(i,j+1,k+1,3);
     int p260=mcL[l];
     l=index(i+1,j,k+1,3);
     int p270=mcL[l];
     l=index(i,j+1,k-1,3);
     int p280=mcL[l];
     l=index(i+3,j,k+1,1);
     int p290=mcL[l];
     l=index(i-3,j,k-1,1);
     int p300=mcL[l];
     l=index(i+3,j,k-1,1);
     int p310=mcL[l];
     l=index(i-3,j,k+1,1);
     int p320=mcL[l];
     l=index(i-1,j,k,1);
     int p330=mcL[l];
     l=index(i+1,j,k,1);
     int p340=mcL[l];
     l=index(i+1,j,k,0);
     int p350=mcL[l];
     l=index(i+2,j,k,1);
     int p360=mcL[l];
     l=index(i-1,j,k,0);
     int p370=mcL[l];
     l=index(i-2,j,k,1);
     int p380=mcL[l];
     l=index(i+1,j,k+1,0);
     int p390=mcL[l];
     l=index(i-1,j,k+1,0);
     int p400=mcL[l];
     l=index(i,j,k+1,0);
     int p410=mcL[l];
     l=index(i,j,k,0);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0150152*(p00*p10+p20*p00+p00*p30+p40*p00)-0.00234223*(p50*p00+p60*p00+p00*p70+p00*p80)-0.0028436*(p00*p90+p100*p00+p00*p110+p120*p00)-0.00172242*(p130*p00+p140*p00+p00*p150+p00*p160)-0.0038228*(p170*p00+p180*p00)+0.00333493*(p00*p190+p200*p00)-0.00189967*(p00*p210+p220*p00+p00*p230+p240*p00+p00*p250+p260*p00+p00*p270+p280*p00)+0.000386838*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0139218*(p00*p330*p10+p340*p00*p30+p20*p40*p00)+0.000628108*(p350*p00*p360+p370*p380*p00+p390*p00*p360+p400*p380*p00+p00*p400*p390+p00*p370*p350)-0.000994741*(p410*p390*p00*p340+p400*p410*p330*p00+p420*p350*p00*p340+p370*p420*p330*p00);
     return energy;
  }


  if(b == 2){
     l=index(i,j,k,2);
     int p00=mcL[l];
     l=index(i,j,k,0);
     int p10=mcL[l];
     l=index(i-1,j+1,k,0);
     int p20=mcL[l];
     l=index(i-1,j,k,0);
     int p30=mcL[l];
     l=index(i,j+1,k,0);
     int p40=mcL[l];
     l=index(i,j,k,1);
     int p50=mcL[l];
     l=index(i,j,k-1,1);
     int p60=mcL[l];
     l=index(i-1,j,k-1,1);
     int p70=mcL[l];
     l=index(i-1,j,k,1);
     int p80=mcL[l];
     l=index(i+1,j,k,0);
     int p90=mcL[l];
     l=index(i-2,j+1,k,0);
     int p100=mcL[l];
     l=index(i-2,j,k,0);
     int p110=mcL[l];
     l=index(i+1,j+1,k,0);
     int p120=mcL[l];
     l=index(i+1,j,k,1);
     int p130=mcL[l];
     l=index(i+1,j,k-1,1);
     int p140=mcL[l];
     l=index(i-2,j,k-1,1);
     int p150=mcL[l];
     l=index(i-2,j,k,1);
     int p160=mcL[l];
     l=index(i,j,k,3);
     int p170=mcL[l];
     l=index(i,j,k-1,3);
     int p180=mcL[l];
     l=index(i,j,k+1,2);
     int p190=mcL[l];
     l=index(i,j,k-1,2);
     int p200=mcL[l];
     l=index(i,j,k+1,0);
     int p210=mcL[l];
     l=index(i-1,j+1,k-1,0);
     int p220=mcL[l];
     l=index(i,j,k-1,0);
     int p230=mcL[l];
     l=index(i-1,j+1,k+1,0);
     int p240=mcL[l];
     l=index(i-1,j,k+1,0);
     int p250=mcL[l];
     l=index(i,j+1,k-1,0);
     int p260=mcL[l];
     l=index(i-1,j,k-1,0);
     int p270=mcL[l];
     l=index(i,j+1,k+1,0);
     int p280=mcL[l];
     l=index(i+3,j,k+1,2);
     int p290=mcL[l];
     l=index(i-3,j,k-1,2);
     int p300=mcL[l];
     l=index(i+3,j,k-1,2);
     int p310=mcL[l];
     l=index(i-3,j,k+1,2);
     int p320=mcL[l];
     l=index(i-1,j,k,2);
     int p330=mcL[l];
     l=index(i+1,j,k,2);
     int p340=mcL[l];
     l=index(i-1,j+1,k,3);
     int p350=mcL[l];
     l=index(i+1,j+1,k,3);
     int p360=mcL[l];
     l=index(i-1,j+1,k-1,3);
     int p370=mcL[l];
     l=index(i+1,j+1,k-1,3);
     int p380=mcL[l];
     l=index(i+2,j,k,2);
     int p390=mcL[l];
     l=index(i-2,j,k,2);
     int p400=mcL[l];
     l=index(i,j+1,k-1,3);
     int p410=mcL[l];
     l=index(i,j+1,k,3);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0150152*(p10*p00+p00*p20+p30*p00+p00*p40)-0.00234223*(p00*p50+p00*p60+p70*p00+p80*p00)-0.0028436*(p90*p00+p00*p100+p110*p00+p00*p120)-0.00172242*(p00*p130+p00*p140+p150*p00+p160*p00)-0.0038228*(p00*p170+p00*p180)+0.00333493*(p00*p190+p200*p00)-0.00189967*(p210*p00+p00*p220+p230*p00+p00*p240+p250*p00+p00*p260+p270*p00+p00*p280)+0.000386838*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0139218*(p10*p30*p00+p00*p330*p20+p340*p00*p40)+0.000628108*(p00*p350*p360+p00*p370*p380+p380*p00*p390+p370*p400*p00+p360*p00*p390+p350*p400*p00)-0.000994741*(p00*p340*p410*p380+p330*p00*p370*p410+p00*p340*p420*p360+p330*p00*p350*p420);
     return energy;
  }


  if(b == 3){
     l=index(i,j,k,3);
     int p00=mcL[l];
     l=index(i,j,k,1);
     int p10=mcL[l];
     l=index(i-1,j-1,k,1);
     int p20=mcL[l];
     l=index(i-1,j,k,1);
     int p30=mcL[l];
     l=index(i,j-1,k,1);
     int p40=mcL[l];
     l=index(i-1,j,k,0);
     int p50=mcL[l];
     l=index(i-1,j,k+1,0);
     int p60=mcL[l];
     l=index(i,j,k+1,0);
     int p70=mcL[l];
     l=index(i,j,k,0);
     int p80=mcL[l];
     l=index(i+1,j,k,1);
     int p90=mcL[l];
     l=index(i-2,j-1,k,1);
     int p100=mcL[l];
     l=index(i-2,j,k,1);
     int p110=mcL[l];
     l=index(i+1,j-1,k,1);
     int p120=mcL[l];
     l=index(i-2,j,k,0);
     int p130=mcL[l];
     l=index(i-2,j,k+1,0);
     int p140=mcL[l];
     l=index(i+1,j,k+1,0);
     int p150=mcL[l];
     l=index(i+1,j,k,0);
     int p160=mcL[l];
     l=index(i,j,k,2);
     int p170=mcL[l];
     l=index(i,j,k+1,2);
     int p180=mcL[l];
     l=index(i,j,k+1,3);
     int p190=mcL[l];
     l=index(i,j,k-1,3);
     int p200=mcL[l];
     l=index(i,j,k+1,1);
     int p210=mcL[l];
     l=index(i-1,j-1,k-1,1);
     int p220=mcL[l];
     l=index(i,j,k-1,1);
     int p230=mcL[l];
     l=index(i-1,j-1,k+1,1);
     int p240=mcL[l];
     l=index(i-1,j,k+1,1);
     int p250=mcL[l];
     l=index(i,j-1,k-1,1);
     int p260=mcL[l];
     l=index(i-1,j,k-1,1);
     int p270=mcL[l];
     l=index(i,j-1,k+1,1);
     int p280=mcL[l];
     l=index(i+3,j,k+1,3);
     int p290=mcL[l];
     l=index(i-3,j,k-1,3);
     int p300=mcL[l];
     l=index(i+3,j,k-1,3);
     int p310=mcL[l];
     l=index(i-3,j,k+1,3);
     int p320=mcL[l];
     l=index(i-1,j,k,3);
     int p330=mcL[l];
     l=index(i+1,j,k,3);
     int p340=mcL[l];
     l=index(i+1,j-1,k,2);
     int p350=mcL[l];
     l=index(i+2,j,k,3);
     int p360=mcL[l];
     l=index(i-1,j-1,k,2);
     int p370=mcL[l];
     l=index(i-2,j,k,3);
     int p380=mcL[l];
     l=index(i+1,j-1,k+1,2);
     int p390=mcL[l];
     l=index(i-1,j-1,k+1,2);
     int p400=mcL[l];
     l=index(i,j-1,k+1,2);
     int p410=mcL[l];
     l=index(i,j-1,k,2);
     int p420=mcL[l];

     energy = energy+0.0385761*(p00)-0.0150152*(p10*p00+p00*p20+p30*p00+p00*p40)-0.00234223*(p50*p00+p60*p00+p00*p70+p00*p80)-0.0028436*(p90*p00+p00*p100+p110*p00+p00*p120)-0.00172242*(p130*p00+p140*p00+p00*p150+p00*p160)-0.0038228*(p170*p00+p180*p00)+0.00333493*(p00*p190+p200*p00)-0.00189967*(p210*p00+p00*p220+p230*p00+p00*p240+p250*p00+p00*p260+p270*p00+p00*p280)+0.000386838*(p00*p290+p300*p00+p00*p310+p320*p00)-0.0139218*(p10*p30*p00+p00*p330*p20+p340*p00*p40)+0.000628108*(p350*p00*p360+p370*p380*p00+p390*p00*p360+p400*p380*p00+p00*p400*p390+p00*p370*p350)-0.000994741*(p410*p390*p00*p340+p400*p410*p330*p00+p420*p350*p00*p340+p370*p420*p330*p00);
     return energy;
  }


}



 
//************************************************************ 
 
void Monte_Carlo::pointcorr(int i, int j, int k, int b){
  int l; 
  if(b == 0){
     l=index(i,j,k,0); 
     double p0=mcL[l]; 
     l=index(i,j,k,2); 
     double p1=mcL[l]; 
     l=index(i+1,j-1,k,2); 
     double p2=mcL[l]; 
     l=index(i+1,j,k,2); 
     double p3=mcL[l]; 
     l=index(i,j-1,k,2); 
     double p4=mcL[l]; 
     l=index(i+1,j,k,3); 
     double p5=mcL[l]; 
     l=index(i+1,j,k-1,3); 
     double p6=mcL[l]; 
     l=index(i,j,k-1,3); 
     double p7=mcL[l]; 
     l=index(i,j,k,3); 
     double p8=mcL[l]; 
     l=index(i-1,j,k,2); 
     double p9=mcL[l]; 
     l=index(i+2,j-1,k,2); 
     double p10=mcL[l]; 
     l=index(i+2,j,k,2); 
     double p11=mcL[l]; 
     l=index(i-1,j-1,k,2); 
     double p12=mcL[l]; 
     l=index(i+2,j,k,3); 
     double p13=mcL[l]; 
     l=index(i+2,j,k-1,3); 
     double p14=mcL[l]; 
     l=index(i-1,j,k-1,3); 
     double p15=mcL[l]; 
     l=index(i-1,j,k,3); 
     double p16=mcL[l]; 
     l=index(i,j-1,k,1); 
     double p17=mcL[l]; 
     l=index(i,j-1,k-1,1); 
     double p18=mcL[l]; 
     l=index(i,j,k+1,0); 
     double p19=mcL[l]; 
     l=index(i,j,k-1,0); 
     double p20=mcL[l]; 
     l=index(i,j,k-1,2); 
     double p21=mcL[l]; 
     l=index(i+1,j-1,k+1,2); 
     double p22=mcL[l]; 
     l=index(i,j,k+1,2); 
     double p23=mcL[l]; 
     l=index(i+1,j-1,k-1,2); 
     double p24=mcL[l]; 
     l=index(i+1,j,k-1,2); 
     double p25=mcL[l]; 
     l=index(i,j-1,k+1,2); 
     double p26=mcL[l]; 
     l=index(i+1,j,k+1,2); 
     double p27=mcL[l]; 
     l=index(i,j-1,k-1,2); 
     double p28=mcL[l]; 
     l=index(i+3,j,k+1,0); 
     double p29=mcL[l]; 
     l=index(i-3,j,k-1,0); 
     double p30=mcL[l]; 
     l=index(i+3,j,k-1,0); 
     double p31=mcL[l]; 
     l=index(i-3,j,k+1,0); 
     double p32=mcL[l]; 
     l=index(i-1,j,k,0); 
     double p33=mcL[l]; 
     l=index(i+1,j,k,0); 
     double p34=mcL[l]; 
     l=index(i-1,j,k,1); 
     double p35=mcL[l]; 
     l=index(i+1,j,k,1); 
     double p36=mcL[l]; 
     l=index(i-1,j,k-1,1); 
     double p37=mcL[l]; 
     l=index(i+1,j,k-1,1); 
     double p38=mcL[l]; 
     l=index(i+2,j,k,0); 
     double p39=mcL[l]; 
     l=index(i-2,j,k,0); 
     double p40=mcL[l]; 
     l=index(i,j,k-1,1); 
     double p41=mcL[l]; 
     l=index(i,j,k,1); 
     double p42=mcL[l]; 

     AVcorr[0]+=1.0/4;
     AVcorr[1]+=(p0)/4;
     AVcorr[3]+=(p0*p1+p2*p0+p0*p3+p4*p0)/16;
     AVcorr[5]+=(p0*p5+p0*p6+p7*p0+p8*p0)/16;
     AVcorr[8]+=(p0*p9+p10*p0+p0*p11+p12*p0)/16;
     AVcorr[9]+=(p0*p13+p0*p14+p15*p0+p16*p0)/16;
     AVcorr[12]+=(p0*p17+p0*p18)/8;
     AVcorr[18]+=(p0*p19+p20*p0)/8;
     AVcorr[23]+=(p0*p21+p22*p0+p0*p23+p24*p0+p0*p25+p26*p0+p0*p27+p28*p0)/32;
     AVcorr[43]+=(p0*p29+p30*p0+p0*p31+p32*p0)/16;
     AVcorr[47]+=(p0*p33*p1+p34*p0*p3+p2*p4*p0)/12;
     AVcorr[54]+=(p0*p35*p36+p0*p37*p38+p38*p0*p39+p37*p40*p0+p36*p0*p39+p35*p40*p0)/24;
     AVcorr[82]+=(p0*p34*p41*p38+p33*p0*p37*p41+p0*p34*p42*p36+p33*p0*p35*p42)/16;
     return;
  }


  if(b == 1){
     l=index(i,j,k,1); 
     double p0=mcL[l]; 
     l=index(i,j,k,3); 
     double p1=mcL[l]; 
     l=index(i+1,j+1,k,3); 
     double p2=mcL[l]; 
     l=index(i+1,j,k,3); 
     double p3=mcL[l]; 
     l=index(i,j+1,k,3); 
     double p4=mcL[l]; 
     l=index(i,j,k,2); 
     double p5=mcL[l]; 
     l=index(i,j,k+1,2); 
     double p6=mcL[l]; 
     l=index(i+1,j,k+1,2); 
     double p7=mcL[l]; 
     l=index(i+1,j,k,2); 
     double p8=mcL[l]; 
     l=index(i-1,j,k,3); 
     double p9=mcL[l]; 
     l=index(i+2,j+1,k,3); 
     double p10=mcL[l]; 
     l=index(i+2,j,k,3); 
     double p11=mcL[l]; 
     l=index(i-1,j+1,k,3); 
     double p12=mcL[l]; 
     l=index(i-1,j,k,2); 
     double p13=mcL[l]; 
     l=index(i-1,j,k+1,2); 
     double p14=mcL[l]; 
     l=index(i+2,j,k+1,2); 
     double p15=mcL[l]; 
     l=index(i+2,j,k,2); 
     double p16=mcL[l]; 
     l=index(i,j+1,k,0); 
     double p17=mcL[l]; 
     l=index(i,j+1,k+1,0); 
     double p18=mcL[l]; 
     l=index(i,j,k+1,1); 
     double p19=mcL[l]; 
     l=index(i,j,k-1,1); 
     double p20=mcL[l]; 
     l=index(i,j,k-1,3); 
     double p21=mcL[l]; 
     l=index(i+1,j+1,k+1,3); 
     double p22=mcL[l]; 
     l=index(i,j,k+1,3); 
     double p23=mcL[l]; 
     l=index(i+1,j+1,k-1,3); 
     double p24=mcL[l]; 
     l=index(i+1,j,k-1,3); 
     double p25=mcL[l]; 
     l=index(i,j+1,k+1,3); 
     double p26=mcL[l]; 
     l=index(i+1,j,k+1,3); 
     double p27=mcL[l]; 
     l=index(i,j+1,k-1,3); 
     double p28=mcL[l]; 
     l=index(i+3,j,k+1,1); 
     double p29=mcL[l]; 
     l=index(i-3,j,k-1,1); 
     double p30=mcL[l]; 
     l=index(i+3,j,k-1,1); 
     double p31=mcL[l]; 
     l=index(i-3,j,k+1,1); 
     double p32=mcL[l]; 
     l=index(i-1,j,k,1); 
     double p33=mcL[l]; 
     l=index(i+1,j,k,1); 
     double p34=mcL[l]; 
     l=index(i+1,j,k,0); 
     double p35=mcL[l]; 
     l=index(i+2,j,k,1); 
     double p36=mcL[l]; 
     l=index(i-1,j,k,0); 
     double p37=mcL[l]; 
     l=index(i-2,j,k,1); 
     double p38=mcL[l]; 
     l=index(i+1,j,k+1,0); 
     double p39=mcL[l]; 
     l=index(i-1,j,k+1,0); 
     double p40=mcL[l]; 
     l=index(i,j,k+1,0); 
     double p41=mcL[l]; 
     l=index(i,j,k,0); 
     double p42=mcL[l]; 

     AVcorr[0]+=1.0/4;
     AVcorr[1]+=(p0)/4;
     AVcorr[3]+=(p0*p1+p2*p0+p0*p3+p4*p0)/16;
     AVcorr[5]+=(p5*p0+p6*p0+p0*p7+p0*p8)/16;
     AVcorr[8]+=(p0*p9+p10*p0+p0*p11+p12*p0)/16;
     AVcorr[9]+=(p13*p0+p14*p0+p0*p15+p0*p16)/16;
     AVcorr[12]+=(p17*p0+p18*p0)/8;
     AVcorr[18]+=(p0*p19+p20*p0)/8;
     AVcorr[23]+=(p0*p21+p22*p0+p0*p23+p24*p0+p0*p25+p26*p0+p0*p27+p28*p0)/32;
     AVcorr[43]+=(p0*p29+p30*p0+p0*p31+p32*p0)/16;
     AVcorr[47]+=(p0*p33*p1+p34*p0*p3+p2*p4*p0)/12;
     AVcorr[54]+=(p35*p0*p36+p37*p38*p0+p39*p0*p36+p40*p38*p0+p0*p40*p39+p0*p37*p35)/24;
     AVcorr[82]+=(p41*p39*p0*p34+p40*p41*p33*p0+p42*p35*p0*p34+p37*p42*p33*p0)/16;
     return;
  }


  if(b == 2){
     l=index(i,j,k,2); 
     double p0=mcL[l]; 
     l=index(i,j,k,0); 
     double p1=mcL[l]; 
     l=index(i-1,j+1,k,0); 
     double p2=mcL[l]; 
     l=index(i-1,j,k,0); 
     double p3=mcL[l]; 
     l=index(i,j+1,k,0); 
     double p4=mcL[l]; 
     l=index(i,j,k,1); 
     double p5=mcL[l]; 
     l=index(i,j,k-1,1); 
     double p6=mcL[l]; 
     l=index(i-1,j,k-1,1); 
     double p7=mcL[l]; 
     l=index(i-1,j,k,1); 
     double p8=mcL[l]; 
     l=index(i+1,j,k,0); 
     double p9=mcL[l]; 
     l=index(i-2,j+1,k,0); 
     double p10=mcL[l]; 
     l=index(i-2,j,k,0); 
     double p11=mcL[l]; 
     l=index(i+1,j+1,k,0); 
     double p12=mcL[l]; 
     l=index(i+1,j,k,1); 
     double p13=mcL[l]; 
     l=index(i+1,j,k-1,1); 
     double p14=mcL[l]; 
     l=index(i-2,j,k-1,1); 
     double p15=mcL[l]; 
     l=index(i-2,j,k,1); 
     double p16=mcL[l]; 
     l=index(i,j,k,3); 
     double p17=mcL[l]; 
     l=index(i,j,k-1,3); 
     double p18=mcL[l]; 
     l=index(i,j,k+1,2); 
     double p19=mcL[l]; 
     l=index(i,j,k-1,2); 
     double p20=mcL[l]; 
     l=index(i,j,k+1,0); 
     double p21=mcL[l]; 
     l=index(i-1,j+1,k-1,0); 
     double p22=mcL[l]; 
     l=index(i,j,k-1,0); 
     double p23=mcL[l]; 
     l=index(i-1,j+1,k+1,0); 
     double p24=mcL[l]; 
     l=index(i-1,j,k+1,0); 
     double p25=mcL[l]; 
     l=index(i,j+1,k-1,0); 
     double p26=mcL[l]; 
     l=index(i-1,j,k-1,0); 
     double p27=mcL[l]; 
     l=index(i,j+1,k+1,0); 
     double p28=mcL[l]; 
     l=index(i+3,j,k+1,2); 
     double p29=mcL[l]; 
     l=index(i-3,j,k-1,2); 
     double p30=mcL[l]; 
     l=index(i+3,j,k-1,2); 
     double p31=mcL[l]; 
     l=index(i-3,j,k+1,2); 
     double p32=mcL[l]; 
     l=index(i-1,j,k,2); 
     double p33=mcL[l]; 
     l=index(i+1,j,k,2); 
     double p34=mcL[l]; 
     l=index(i-1,j+1,k,3); 
     double p35=mcL[l]; 
     l=index(i+1,j+1,k,3); 
     double p36=mcL[l]; 
     l=index(i-1,j+1,k-1,3); 
     double p37=mcL[l]; 
     l=index(i+1,j+1,k-1,3); 
     double p38=mcL[l]; 
     l=index(i+2,j,k,2); 
     double p39=mcL[l]; 
     l=index(i-2,j,k,2); 
     double p40=mcL[l]; 
     l=index(i,j+1,k-1,3); 
     double p41=mcL[l]; 
     l=index(i,j+1,k,3); 
     double p42=mcL[l]; 

     AVcorr[0]+=1.0/4;
     AVcorr[1]+=(p0)/4;
     AVcorr[3]+=(p1*p0+p0*p2+p3*p0+p0*p4)/16;
     AVcorr[5]+=(p0*p5+p0*p6+p7*p0+p8*p0)/16;
     AVcorr[8]+=(p9*p0+p0*p10+p11*p0+p0*p12)/16;
     AVcorr[9]+=(p0*p13+p0*p14+p15*p0+p16*p0)/16;
     AVcorr[12]+=(p0*p17+p0*p18)/8;
     AVcorr[18]+=(p0*p19+p20*p0)/8;
     AVcorr[23]+=(p21*p0+p0*p22+p23*p0+p0*p24+p25*p0+p0*p26+p27*p0+p0*p28)/32;
     AVcorr[43]+=(p0*p29+p30*p0+p0*p31+p32*p0)/16;
     AVcorr[47]+=(p1*p3*p0+p0*p33*p2+p34*p0*p4)/12;
     AVcorr[54]+=(p0*p35*p36+p0*p37*p38+p38*p0*p39+p37*p40*p0+p36*p0*p39+p35*p40*p0)/24;
     AVcorr[82]+=(p0*p34*p41*p38+p33*p0*p37*p41+p0*p34*p42*p36+p33*p0*p35*p42)/16;
     return;
  }


  if(b == 3){
     l=index(i,j,k,3); 
     double p0=mcL[l]; 
     l=index(i,j,k,1); 
     double p1=mcL[l]; 
     l=index(i-1,j-1,k,1); 
     double p2=mcL[l]; 
     l=index(i-1,j,k,1); 
     double p3=mcL[l]; 
     l=index(i,j-1,k,1); 
     double p4=mcL[l]; 
     l=index(i-1,j,k,0); 
     double p5=mcL[l]; 
     l=index(i-1,j,k+1,0); 
     double p6=mcL[l]; 
     l=index(i,j,k+1,0); 
     double p7=mcL[l]; 
     l=index(i,j,k,0); 
     double p8=mcL[l]; 
     l=index(i+1,j,k,1); 
     double p9=mcL[l]; 
     l=index(i-2,j-1,k,1); 
     double p10=mcL[l]; 
     l=index(i-2,j,k,1); 
     double p11=mcL[l]; 
     l=index(i+1,j-1,k,1); 
     double p12=mcL[l]; 
     l=index(i-2,j,k,0); 
     double p13=mcL[l]; 
     l=index(i-2,j,k+1,0); 
     double p14=mcL[l]; 
     l=index(i+1,j,k+1,0); 
     double p15=mcL[l]; 
     l=index(i+1,j,k,0); 
     double p16=mcL[l]; 
     l=index(i,j,k,2); 
     double p17=mcL[l]; 
     l=index(i,j,k+1,2); 
     double p18=mcL[l]; 
     l=index(i,j,k+1,3); 
     double p19=mcL[l]; 
     l=index(i,j,k-1,3); 
     double p20=mcL[l]; 
     l=index(i,j,k+1,1); 
     double p21=mcL[l]; 
     l=index(i-1,j-1,k-1,1); 
     double p22=mcL[l]; 
     l=index(i,j,k-1,1); 
     double p23=mcL[l]; 
     l=index(i-1,j-1,k+1,1); 
     double p24=mcL[l]; 
     l=index(i-1,j,k+1,1); 
     double p25=mcL[l]; 
     l=index(i,j-1,k-1,1); 
     double p26=mcL[l]; 
     l=index(i-1,j,k-1,1); 
     double p27=mcL[l]; 
     l=index(i,j-1,k+1,1); 
     double p28=mcL[l]; 
     l=index(i+3,j,k+1,3); 
     double p29=mcL[l]; 
     l=index(i-3,j,k-1,3); 
     double p30=mcL[l]; 
     l=index(i+3,j,k-1,3); 
     double p31=mcL[l]; 
     l=index(i-3,j,k+1,3); 
     double p32=mcL[l]; 
     l=index(i-1,j,k,3); 
     double p33=mcL[l]; 
     l=index(i+1,j,k,3); 
     double p34=mcL[l]; 
     l=index(i+1,j-1,k,2); 
     double p35=mcL[l]; 
     l=index(i+2,j,k,3); 
     double p36=mcL[l]; 
     l=index(i-1,j-1,k,2); 
     double p37=mcL[l]; 
     l=index(i-2,j,k,3); 
     double p38=mcL[l]; 
     l=index(i+1,j-1,k+1,2); 
     double p39=mcL[l]; 
     l=index(i-1,j-1,k+1,2); 
     double p40=mcL[l]; 
     l=index(i,j-1,k+1,2); 
     double p41=mcL[l]; 
     l=index(i,j-1,k,2); 
     double p42=mcL[l]; 

     AVcorr[0]+=1.0/4;
     AVcorr[1]+=(p0)/4;
     AVcorr[3]+=(p1*p0+p0*p2+p3*p0+p0*p4)/16;
     AVcorr[5]+=(p5*p0+p6*p0+p0*p7+p0*p8)/16;
     AVcorr[8]+=(p9*p0+p0*p10+p11*p0+p0*p12)/16;
     AVcorr[9]+=(p13*p0+p14*p0+p0*p15+p0*p16)/16;
     AVcorr[12]+=(p17*p0+p18*p0)/8;
     AVcorr[18]+=(p0*p19+p20*p0)/8;
     AVcorr[23]+=(p21*p0+p0*p22+p23*p0+p0*p24+p25*p0+p0*p26+p27*p0+p0*p28)/32;
     AVcorr[43]+=(p0*p29+p30*p0+p0*p31+p32*p0)/16;
     AVcorr[47]+=(p1*p3*p0+p0*p33*p2+p34*p0*p4)/12;
     AVcorr[54]+=(p35*p0*p36+p37*p38*p0+p39*p0*p36+p40*p38*p0+p0*p40*p39+p0*p37*p35)/24;
     AVcorr[82]+=(p41*p39*p0*p34+p40*p41*p33*p0+p42*p35*p0*p34+p37*p42*p33*p0)/16;
     return;
  }


}
