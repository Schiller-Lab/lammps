/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Ternary Lattice Boltzmann Implementation

   Contributing authors: Ulf D. Schiller <uschiller@mailaps.org>,
                         Fang Wang <fwang8@clemson.edu>
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(lb/multicomponent, FixLbMulticomponent)
// clang-format on
#else

#ifndef LMP_FIX_LB_MULTICOMPONENT_H
#define LMP_FIX_LB_MULTICOMPONENT_H

#include "fix.h"
#include "fix_lb_fluid.h"

namespace LAMMPS_NS {

static constexpr double CS2 = 1. / 3.;        // D3Q19 lattice speed of sound
static constexpr double WG19[19][3][3] = {    // weights for gradient terms in equilibrium
    {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
    {{5. / 36., 0., 0.}, {0., -1. / 9., 0.}, {0., 0., -1. / 9.}},
    {{-1. / 9., 0., 0.}, {0., 5. / 36., 0.}, {0., 0., -1. / 9.}},
    {{5. / 36., 0., 0.}, {0., -1. / 9., 0.}, {0., 0., -1. / 9.}},
    {{-1. / 9., 0., 0.}, {0., 5. / 36., 0.}, {0., 0., -1. / 9.}},
    {{-1. / 9., 0., 0.}, {0., -1. / 9., 0.}, {0., 0., 5. / 36.}},
    {{-1. / 9., 0., 0.}, {0., -1. / 9., 0.}, {0., 0., 5. / 36.}},
    {{-1. / 72., 1. / 12., 0.}, {1. / 12., -1. / 72., 0.}, {0., 0., 1. / 36.}},
    {{-1. / 72., -1. / 12., 0.}, {-1. / 12., -1. / 72., 0.}, {0., 0., 1. / 36.}},
    {{-1. / 72., -1. / 12., 0.}, {-1. / 12., -1. / 72., 0.}, {0., 0., 1. / 36.}},
    {{-1. / 72., 1. / 12., 0.}, {1. / 12., -1. / 72., 0.}, {0., 0., 1. / 36.}},
    {{-1. / 72., 0., 1. / 12.}, {0., 1. / 36., 0.}, {1. / 12., 0., -1. / 72.}},
    {{-1. / 72., 0., -1. / 12.}, {0., 1. / 36., 0.}, {-1. / 12., 0., -1. / 72.}},
    {{-1. / 72., 0., -1. / 12.}, {0., 1. / 36., 0.}, {-1. / 12., 0., -1. / 72.}},
    {{-1. / 72., 0., 1. / 12.}, {0., 1. / 36., 0.}, {1. / 12., 0., -1. / 72.}},
    {{1. / 36., 0., 0.}, {0., -1. / 72., 1. / 12.}, {0., 1. / 12., -1. / 72.}},
    {{1. / 36., 0., 0.}, {0., -1. / 72., -1. / 12.}, {0., -1. / 12., -1. / 72.}},
    {{1. / 36., 0., 0.}, {0., -1. / 72., -1. / 12.}, {0., -1. / 12., -1. / 72.}},
    {{1. / 36., 0., 0.}, {0., -1. / 72., 1. / 12.}, {0., 1. / 12., -1. / 72.}}};

class FixLbMulticomponent : public FixLbFluid {

 public:
  FixLbMulticomponent(class LAMMPS *, int, char **);
  ~FixLbMulticomponent() override;

  int setmask() override;
  void initial_integrate(int) override;
  void end_of_step() override;
  void get_total_momentum(double &, double &, double &);
  void get_total_mass(double &);

 private:
  double tau_r, tau_p, tau_s;
  double gamma_p, gamma_s;
  double kappa1, kappa2, kappa3;
  double kappa_rr, kappa_pp, kappa_ss, kappa_rp, kappa_ps, kappa_rs;
  double alpha;

  int seed;
  double C1, C2, C3;
  int radius;
  double C1_drop, C2_drop;
  double C1_film, C2_film;
  double thickness;

  double ****g_lb;
  double ****gnew;
  double ****geq;
  double ****k_lb;
  double ****knew;
  double ****keq;

  double ***pressure_lb;
  double ***phi_lb;
  double ***psi_lb;

  double ****density_gradient;
  double ****phi_gradient;
  double ****psi_gradient;
  double ***laplace_rho;
  double ***laplace_phi;
  double ***laplace_psi;

  double ***mu_rho;
  double ***mu_phi;
  double ***mu_psi;

  void init_parameters(int, char **);
  void init_lattice();
  void destroy_lattice();

  enum init_type { MIXTURE, DROPLET, LIQUID_LENS, DOUBLE_EMULSION, FILM, MIXED_DROPLET };
  init_type init_method = MIXTURE;

  void init_fluid();
  void init_mixture();
  void init_droplet(double radius);
  void init_liquid_lens(double radius);
  void init_double_emulsion(double radius);
  void init_film(double thickness, double C1, double C2);
  void init_mixed_droplet(double radius, double C1, double C2);

  void lb_update();
  void calc_moments_full();
  void collide_stream(int x, int y, int z);
  void update_cube(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
  void update_slab(int x, int ymin, int ymax, int zmin, int zmax);
  void update_column(int x, int y, int zmin, int zmax);
  void read_slab(int x, int ymin, int ymax, int zmin, int zmax);
  void read_column(int x, int y, int zmin, int zmax);
  void read_site(int x, int y, int z);
  void write_site(int x, int y, int z);

  void calc_moments(int x, int y, int z);
  void calc_chemical_potentials(int x, int y, int z);
  void calc_equilibrium(int x, int y, int z);
  void calc_feq(int x, int y, int z);
  void calc_geq(int x, int y, int z);
  void calc_keq(int x, int y, int z);
  void calc_gradient_laplacian(int x, int y, int z, double ***field, double ****gradient,
                               double ***laplacian);
  void calc_rho_gradients(int x, int y, int z);
  void calc_phi_gradients(int x, int y, int z);
  void calc_psi_gradients(int x, int y, int z);
  double pressure(double rho, double phi, double psi);

  static const int numrequests = 12;
  MPI_Request requests[numrequests];
  MPI_Datatype fluid_scalar_field_mpitype;
  MPI_Datatype fluid_vector_field_mpitype;

  int halo_extent[3];    // extent of halo in x,y,z

  void init_halo();
  void destroy_halo();
  void halo_wait();
  void halo_comm();
  void halo_comm(int dir);

  void init_output();
  void destroy_output();
  void dump_xdmf(int);
};

}    // namespace LAMMPS_NS
#endif
#endif
