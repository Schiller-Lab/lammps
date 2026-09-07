/* ----------------------------------------------------------------------
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

   References:
   [1] Semprebon et al., Phys. Rev. E 93, 033305 (2016)
       https://doi.org/10.1103/PhysRevE.93.033305
   [2] Pooley and Furtado, Phys. Rev. E 77, 046702 (2008)
       https://doi.org/10.1103/PhysRevE.77.046702
   [3] Boyer and Lapuerta, ESAIM Math. Model. Numer. Anal. 40, 653-687 (2006)
       https://doi.org/10.1051/m2an:2006028
   [4] Swift et al., Phys. Rev. E 54, 5041 (1996)
       https://doi.org/10.1103/PhysRevE.54.5041
------------------------------------------------------------------------- */

#include "fix_lb_multicomponent.h"
#include "latboltz_const.h"

#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "random_mars.h"
#include "update.h"

using namespace LAMMPS_NS;

static const char cite_fix_lbmulticomponent[] =
    "fix lb/multicomponent command: doi:10.1016/j.cpc.2023.108898\n\n"
    "@Article{arumugam_kumar_implementation_2024,\n"
    "  author  = {Arumugam Kumar, Gokul Raman and Andrews, James P. and Schiller, Ulf D.},\n"
    "  title   = {Implementation of a Ternary Lattice Boltzmann Model in LAMMPS},\n"
    "  journal = {Comput.~Phys.~Commun.},\n"
    "  year    = {2024},\n"
    "  volume  = {294},\n"
    "  pages   = {108898}\n"
    "}\n\n";

int FixLbMulticomponent::setmask()
{
  return FixConst::INITIAL_INTEGRATE | FixConst::END_OF_STEP;
}

void FixLbMulticomponent::initial_integrate(int vflag)
{
  this->lb_update();
}

void FixLbMulticomponent::end_of_step()
{
  dump_xdmf(update->ntimestep);
}

void FixLbMulticomponent::get_total_momentum(double &jx, double &jy, double &jz)
{
  double local[3] = {0.0, 0.0, 0.0};
  double global[3] = {0.0, 0.0, 0.0};
    
  for (int x = halo_extent[0]; x < subNbx - halo_extent[0]; ++x)
    for (int y = halo_extent[1]; y < subNby - halo_extent[1]; ++y)
      for (int z = halo_extent[2]; z < subNbz - halo_extent[2]; ++z)
        for(int i = 0; i < numvel; ++i){
            local[0] += f_lb[x][y][z][i] * e19[i][0];
            local[1] += f_lb[x][y][z][i] * e19[i][1];
            local[2] += f_lb[x][y][z][i] * e19[i][2];
        }
  
  MPI_Reduce(local, global, 3, MPI_DOUBLE, MPI_SUM, 0, world);
  
  if (me == 0){
    jx = global[0];
    jy = global[1];
    jz = global[2];
  }
}

void FixLbMulticomponent::get_total_mass(double &m)
{
  double local = 0.0;
  double global = 0.0;
    
  for (int x = halo_extent[0]; x < subNbx - halo_extent[0]; ++x)
    for (int y = halo_extent[1]; y < subNby - halo_extent[1]; ++y)
      for (int z = halo_extent[2]; z < subNbz - halo_extent[2]; ++z)
        for(int i = 0; i < numvel; ++i){
            local += f_lb[x][y][z][i];
        }
  
  MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  
  if (me == 0)   
    m = global;
}

void FixLbMulticomponent::lb_update()
{

#if 1
  // no overlap of communication and computation
  halo_comm();
  update_cube(0, subNbx, 0, subNby, 0, subNbz);
#else    // the following is not thoroughly tested as it appeared to be slower [uschille 2022/12/15]
  // communication and computation can be partially overlapped
  // some computations at the box boundary are duplicated which could be optimized

  // communicate in z-direction
  halo_comm(2);
  // update inner cube
  update_cube(2, subNbx - 2, 2, subNby - 2, 2, subNbz - 2);
  // wait for communication of z-slabs
  halo_wait();

  // communicate in y-direction
  halo_comm(1);
  // update z-slabs that are now available
  update_cube(2, subNbx - 2, 2, subNby - 2, 0, 4);
  update_cube(2, subNbx - 2, 2, subNby - 2, subNbz - 4, subNbz);
  // wait for communication of y-slabs
  halo_wait();

  // communicate in x-direction
  halo_comm(0);
  // update y-slabs that are now available
  update_cube(2, subNbx - 2, 0, 4, 0, subNbz);
  update_cube(2, subNbx - 2, subNby - 4, subNby, 0, subNbz);
  // wait for communication of z-slabs
  halo_wait();

  // update x-slabs that are now available
  update_cube(0, 4, 0, subNby, 0, subNbz);
  update_cube(subNbx - 4, subNbx, 0, subNby, 0, subNbz);
#endif

  /* swap the pointers of the lattice copies */
  std::swap(f_lb, fnew);
  std::swap(g_lb, gnew);
  std::swap(k_lb, knew);
}

void FixLbMulticomponent::update_cube(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
  int x;
  read_slab(xmin, ymin, ymax, zmin, zmax);
  read_slab(xmin + 1, ymin, ymax, zmin, zmax);
  for (x = xmin + 2; x < xmax; ++x) { update_slab(x, ymin, ymax, zmin, zmax); }
}

void FixLbMulticomponent::update_slab(int x, int ymin, int ymax, int zmin, int zmax)
{
  int y;
  read_column(x, ymin, zmin, zmax);
  read_column(x, ymin + 1, zmin, zmax);
  for (y = ymin + 2; y < ymax; ++y) { update_column(x, y, zmin, zmax); }
}

void FixLbMulticomponent::update_column(int x, int y, int zmin, int zmax)
{
  int z;
  read_site(x, y, zmin);
  read_site(x, y, zmin + 1);
  for (z = zmin + 2; z < zmax; ++z) {
    read_site(x, y, z);
    write_site(x - 1, y - 1, z - 1);
  }
}

void FixLbMulticomponent::read_slab(int x, int ymin, int ymax, int zmin, int zmax)
{
  int y;
  for (y = ymin; y < ymax; ++y) { read_column(x, y, zmin, zmax); }
}

void FixLbMulticomponent::read_column(int x, int y, int zmin, int zmax)
{
  int z;
  for (z = zmin; z < zmax; ++z) { read_site(x, y, z); }
}

void FixLbMulticomponent::read_site(int x, int y, int z)
{
  calc_moments(x, y, z);
}

void FixLbMulticomponent::write_site(int x, int y, int z)
{
  collide_stream(x, y, z);
}

void FixLbMulticomponent::collide_stream(int x, int y, int z)
{
  int i, xnew, ynew, znew;
  calc_equilibrium(x, y, z);
  for (i = 0; i < numvel; ++i) {
    xnew = x + e19[i][0];
    ynew = y + e19[i][1];
    znew = z + e19[i][2];
    fnew[xnew][ynew][znew][i] = f_lb[x][y][z][i] - (f_lb[x][y][z][i] - feq[x][y][z][i]) / tau_r;
    gnew[xnew][ynew][znew][i] = g_lb[x][y][z][i] - (g_lb[x][y][z][i] - geq[x][y][z][i]) / tau_p;
    knew[xnew][ynew][znew][i] = k_lb[x][y][z][i] - (k_lb[x][y][z][i] - keq[x][y][z][i]) / tau_s;
  }
}

void FixLbMulticomponent::calc_moments(int x, int y, int z)
{
  double rho, phi, psi, j[3], fi, gi, ki;
  int i;
  rho = phi = psi = j[0] = j[1] = j[2] = 0.0;
  for (i = 0; i < numvel; ++i) {
    fi = f_lb[x][y][z][i];
    gi = g_lb[x][y][z][i];
    ki = k_lb[x][y][z][i];
    rho += fi;
    phi += gi;
    psi += ki;
    j[0] += fi * e19[i][0];
    j[1] += fi * e19[i][1];
    j[2] += fi * e19[i][2];
  }
  density_lb[x][y][z] = rho;
  phi_lb[x][y][z] = phi;
  psi_lb[x][y][z] = psi;
  u_lb[x][y][z][0] = j[0] / rho;
  u_lb[x][y][z][1] = j[1] / rho;
  u_lb[x][y][z][2] = j[2] / rho;
  pressure_lb[x][y][z] = pressure(rho, phi, psi);
}

void FixLbMulticomponent::calc_equilibrium(int x, int y, int z)
{
  calc_gradient_laplacian(x, y, z, density_lb, density_gradient, laplace_rho);
  calc_gradient_laplacian(x, y, z, phi_lb, phi_gradient, laplace_phi);
  calc_gradient_laplacian(x, y, z, psi_lb, psi_gradient, laplace_psi);
  calc_chemical_potentials(x, y, z);
  calc_feq(x, y, z);
  calc_geq(x, y, z);
  calc_keq(x, y, z);
}

void FixLbMulticomponent::calc_gradient_laplacian(int x, int y, int z, double ***field,
                                                  double ****gradient, double ***laplacian)
{
  int i, xp, yp, zp, dir;
  laplacian[x][y][z] = 0.0;
  for (dir = 0; dir < 3; dir++) gradient[x][y][z][dir] = 0.0;
  for (i = 0; i < numvel; i++) {
    xp = x + e19[i][0];
    yp = y + e19[i][1];
    zp = z + e19[i][2];
    for (dir = 0; dir < 3; dir++) {
      gradient[x][y][z][dir] += 3. * w_lb19[i] * field[xp][yp][zp] * e19[i][dir];
    }
    laplacian[x][y][z] += 6. * w_lb19[i] * (field[xp][yp][zp] - field[x][y][z]);
  }
}

double FixLbMulticomponent::pressure(double rho, double phi, double psi)
{
  const double rho2 = rho * rho;
  const double rho3 = rho2 * rho;
  const double rho4 = rho3 * rho;
  const double phi2 = phi * phi;
  const double phi3 = phi2 * phi;
  const double phi4 = phi3 * phi;
  const double psi2 = psi * psi;
  const double psi3 = psi2 * psi;
  const double psi4 = psi3 * psi;
  double p0;

#if 0
  /* sympy ccode of '\sum rho mu_\rho - F' */
  p0 = rho*CS2;
  p0 += kappa1*((3.0/32.0)*pow(phi, 4) - 3.0/8.0*pow(phi, 3)*psi + (3.0/8.0)*pow(phi, 3)*rho - 1.0/4.0*pow(phi, 3) + (9.0/16.0)*pow(phi, 2)*pow(psi, 2) - 9.0/8.0*pow(phi, 2)*psi*rho + (3.0/4.0)*pow(phi, 2)*psi + (9.0/16.0)*pow(phi, 2)*pow(rho, 2) - 3.0/4.0*pow(phi, 2)*rho + (1.0/8.0)*pow(phi, 2) - 3.0/8.0*phi*pow(psi, 3) + (9.0/8.0)*phi*pow(psi, 2)*rho - 3.0/4.0*phi*pow(psi, 2) - 9.0/8.0*phi*psi*pow(rho, 2) + (3.0/2.0)*phi*psi*rho - 1.0/4.0*phi*psi + (3.0/8.0)*phi*pow(rho, 3) - 3.0/4.0*phi*pow(rho, 2) + (1.0/4.0)*phi*rho + (3.0/32.0)*pow(psi, 4) - 3.0/8.0*pow(psi, 3)*rho + (1.0/4.0)*pow(psi, 3) + (9.0/16.0)*pow(psi, 2)*pow(rho, 2) - 3.0/4.0*pow(psi, 2)*rho + (1.0/8.0)*pow(psi, 2) - 3.0/8.0*psi*pow(rho, 3) + (3.0/4.0)*psi*pow(rho, 2) - 1.0/4.0*psi*rho + (3.0/32.0)*pow(rho, 4) - 1.0/4.0*pow(rho, 3) + (1.0/8.0)*pow(rho, 2)) + kappa2*((3.0/32.0)*pow(phi, 4) + (3.0/8.0)*pow(phi, 3)*psi - 3.0/8.0*pow(phi, 3)*rho + (1.0/4.0)*pow(phi, 3) + (9.0/16.0)*pow(phi, 2)*pow(psi, 2) - 9.0/8.0*pow(phi, 2)*psi*rho + (3.0/4.0)*pow(phi, 2)*psi + (9.0/16.0)*pow(phi, 2)*pow(rho, 2) - 3.0/4.0*pow(phi, 2)*rho + (1.0/8.0)*pow(phi, 2) + (3.0/8.0)*phi*pow(psi, 3) - 9.0/8.0*phi*pow(psi, 2)*rho + (3.0/4.0)*phi*pow(psi, 2) + (9.0/8.0)*phi*psi*pow(rho, 2) - 3.0/2.0*phi*psi*rho + (1.0/4.0)*phi*psi - 3.0/8.0*phi*pow(rho, 3) + (3.0/4.0)*phi*pow(rho, 2) - 1.0/4.0*phi*rho + (3.0/32.0)*pow(psi, 4) - 3.0/8.0*pow(psi, 3)*rho + (1.0/4.0)*pow(psi, 3) + (9.0/16.0)*pow(psi, 2)*pow(rho, 2) - 3.0/4.0*pow(psi, 2)*rho + (1.0/8.0)*pow(psi, 2) - 3.0/8.0*psi*pow(rho, 3) + (3.0/4.0)*psi*pow(rho, 2) - 1.0/4.0*psi*rho + (3.0/32.0)*pow(rho, 4) - 1.0/4.0*pow(rho, 3) + (1.0/8.0)*pow(rho, 2)) + kappa3*((3.0/2.0)*pow(psi, 4) - 2*pow(psi, 3) + (1.0/2.0)*pow(psi, 2));
#else
  p0 = rho *
          CS2    // Eq. (43) Semprebon et al. (note that there is a typo in this equation in the paper)
      + (kappa1 + kappa2) *
          (3. / 32. * (rho4 + phi4 + psi4) - 1. / 4. * (rho3 + rho * psi - psi3) +
           1. / 8. * (rho2 + phi2 + psi2) - 3. / 8. * (rho3 * psi + psi3 * rho) +
           9. / 16. * (rho2 * phi2 + rho2 * psi2 + phi2 * psi2) +
           3. / 4. * (rho2 * psi - rho * phi2 - rho * psi2 + phi2 * psi) -
           9. / 8. * phi2 * psi * rho) +
      (kappa1 - kappa2) *
          (3. / 8. * (rho3 * phi + rho * phi3 - phi3 * psi - phi * psi3) +
           1. / 4. * (rho * phi - phi * psi - phi3) +
           9. / 8. * (phi * psi2 * rho - phi * psi * rho2) - 3. / 4. * (rho2 * phi + phi * psi2) +
           3. / 2. * phi * psi * rho) +
      kappa3 * (3. / 2. * psi4 - 2. * psi3 + 1. / 2. * psi2);
#endif

  return p0;
}

void FixLbMulticomponent::calc_chemical_potentials(int x, int y, int z)
{
  const double alpha2 = alpha * alpha;
  const double rho = density_lb[x][y][z];
  const double phi = phi_lb[x][y][z];
  const double psi = psi_lb[x][y][z];
  const double D2rho = laplace_rho[x][y][z];
  const double D2phi = laplace_phi[x][y][z];
  const double D2psi = laplace_psi[x][y][z];

#if 0
  /* mu_rho is not needed for the calculations.
     The following expression has not been tested!
     (Eq. (38) in Semprebon et al. may contain typos) */
  mu_rho[x][y][z] = // Eq. (38) in Semprebon et al.
      kappa1/8.*(rho+psi-phi)*(rho+phi-psi-2.)*(rho+phi-psi-1.)
    + kappa2/8.*(rho-phi-psi)*(rho-phi-psi-2.)*(rho-phi-psi-1.)
    - alpha2/4.*((kappa1+kappa2)*(D2rho-D2psi)-(kappa1-kappa2)*D2phi);
#endif

  mu_phi[x][y][z] =    // Eq. (39) in Semprebon et al.
      kappa1 / 8. * (rho + phi - psi) * (rho + phi - psi - 2.) * (rho + phi - psi - 1.) -
      kappa2 / 8. * (rho - phi - psi) * (rho - phi - psi - 2.) * (rho - phi - psi - 1.) -
      alpha2 / 4. * ((kappa1 - kappa2) * (D2rho - D2psi) + (kappa1 + kappa2) * D2phi);

  mu_psi[x][y][z] =    // Eq. (40) in Semprebon et al.
      -kappa1 / 8. * (rho + phi - psi) * (rho + phi - psi - 2.) * (rho + phi - psi - 1.) -
      kappa2 / 8. * (rho - phi - psi) * (rho - phi - psi - 2.) * (rho - phi - psi - 1.) +
      kappa3 * psi * (psi - 1.) * (2. * psi - 1.) +
      alpha2 / 4. *
          ((kappa1 + kappa2) * D2rho + (kappa1 - kappa2) * D2phi -
           (kappa1 + kappa2 + 4. * kappa3) * D2psi);
}

void FixLbMulticomponent::calc_feq(int x, int y, int z)
{
  const double rho = density_lb[x][y][z];
  const double phi = phi_lb[x][y][z];
  const double psi = psi_lb[x][y][z];
  const double p0 = pressure_lb[x][y][z];
  const double *u = u_lb[x][y][z];
  const double *Drho = density_gradient[x][y][z];
  const double *Dphi = phi_gradient[x][y][z];
  const double *Dpsi = psi_gradient[x][y][z];
  const double D2rho = laplace_rho[x][y][z];
  const double D2phi = laplace_phi[x][y][z];
  const double D2psi = laplace_psi[x][y][z];
  double fi, ruu[3][3], G[3][3];
  int i;

  ruu[0][0] = rho * u[0] * u[0];
  ruu[1][1] = rho * u[1] * u[1];
  ruu[2][2] = rho * u[2] * u[2];
  ruu[0][1] = rho * u[0] * u[1];
  ruu[1][2] = rho * u[1] * u[2];
  ruu[2][0] = rho * u[2] * u[0];

  G[0][0] =
      kappa_rr * Drho[0] * Drho[0] + kappa_pp * Dphi[0] * Dphi[0] + kappa_ss * Dpsi[0] * Dpsi[0];
  G[1][1] =
      kappa_rr * Drho[1] * Drho[1] + kappa_pp * Dphi[1] * Dphi[1] + kappa_ss * Dpsi[1] * Dpsi[1];
  G[2][2] =
      kappa_rr * Drho[2] * Drho[2] + kappa_pp * Dphi[2] * Dphi[2] + kappa_ss * Dpsi[2] * Dpsi[2];
  G[0][1] =
      kappa_rr * Drho[0] * Drho[1] + kappa_pp * Dphi[0] * Dphi[1] + kappa_ss * Dpsi[0] * Dpsi[1];
  G[1][2] =
      kappa_rr * Drho[1] * Drho[2] + kappa_pp * Dphi[1] * Dphi[2] + kappa_ss * Dpsi[1] * Dpsi[2];
  G[2][0] =
      kappa_rr * Drho[2] * Drho[0] + kappa_pp * Dphi[2] * Dphi[0] + kappa_ss * Dpsi[2] * Dpsi[0];

  double sumf = 0.0;
  for (i = 1; i < numvel; ++i) {    // Eq. (52) Semprebon et al.
    fi = 3. * w_lb19[i] * p0;
    fi += 3. * w_lb19[i] * rho * (u[0] * e19[i][0] + u[1] * e19[i][1] + u[2] * e19[i][2]);
    fi += 9. / 2. * w_lb19[i] *
        ((ruu[0][0] * e19[i][0] + 2. * ruu[0][1] * e19[i][1]) * e19[i][0] +
         (ruu[1][1] * e19[i][1] + 2. * ruu[1][2] * e19[i][2]) * e19[i][1] +
         (ruu[2][2] * e19[i][2] + 2. * ruu[2][0] * e19[i][0]) * e19[i][2]);
    fi -= 3. / 2. * w_lb19[i] * (ruu[0][0] + ruu[1][1] + ruu[2][2]);
    fi -=
        3. * w_lb19[i] * (kappa_rr * rho * D2rho + kappa_pp * phi * D2phi + kappa_ss * psi * D2psi);
    fi -= 3. * w_lb19[i] *
        (kappa_rp * (rho * D2phi + phi * D2rho) + kappa_rs * (rho * D2psi + psi * D2rho) +
         kappa_ps * (phi * D2psi + psi * D2phi));
    fi += 3. *
        (WG19[i][0][0] * G[0][0] + WG19[i][1][1] * G[1][1] + WG19[i][2][2] * G[2][2] +
         WG19[i][0][1] * G[0][1] + WG19[i][1][2] * G[1][2] + WG19[i][2][0] * G[2][0]);
    fi += 6. * kappa_rp *
        (WG19[i][0][0] * Drho[0] * Dphi[0] + WG19[i][1][1] * Drho[1] * Dphi[1] +
         WG19[i][2][2] * Drho[2] * Dphi[2]);
    fi += 6. * kappa_rs *
        (WG19[i][0][0] * Drho[0] * Dpsi[0] + WG19[i][1][1] * Drho[1] * Dpsi[1] +
         WG19[i][2][2] * Drho[2] * Dpsi[2]);
    fi += 6. * kappa_ps *
        (WG19[i][0][0] * Dphi[0] * Dpsi[0] + WG19[i][1][1] * Dphi[1] * Dpsi[1] +
         WG19[i][2][2] * Dphi[2] * Dpsi[2]);
    fi += 3. * kappa_rp *
        (WG19[i][0][1] * (Drho[0] * Dphi[1] + Drho[1] * Dphi[0]) +
         WG19[i][1][2] * (Drho[1] * Dphi[2] + Drho[2] * Dphi[1]) +
         WG19[i][2][0] * (Drho[2] * Dphi[0] + Drho[0] * Dphi[2]));
    fi += 3. * kappa_rs *
        (WG19[i][0][1] * (Drho[0] * Dpsi[1] + Drho[1] * Dpsi[0]) +
         WG19[i][1][2] * (Drho[1] * Dpsi[2] + Drho[2] * Dpsi[1]) +
         WG19[i][2][0] * (Drho[2] * Dpsi[0] + Drho[0] * Dpsi[2]));
    fi += 3. * kappa_ps *
        (WG19[i][0][1] * (Dphi[0] * Dpsi[1] + Dphi[1] * Dpsi[0]) +
         WG19[i][1][2] * (Dphi[1] * Dpsi[2] + Dphi[2] * Dpsi[1]) +
         WG19[i][2][0] * (Dphi[2] * Dpsi[0] + Dphi[0] * Dpsi[2]));
    feq[x][y][z][i] = fi;
    sumf += fi;
  }
  feq[x][y][z][0] = rho - sumf;
}

void FixLbMulticomponent::calc_geq(int x, int y, int z)
{
  const double phi = phi_lb[x][y][z];
  const double mu_p = mu_phi[x][y][z];
  const double *u = u_lb[x][y][z];
  double gi, puu[3][3];
  int i;

  puu[0][0] = phi * u[0] * u[0];
  puu[1][1] = phi * u[1] * u[1];
  puu[2][2] = phi * u[2] * u[2];
  puu[0][1] = phi * u[0] * u[1];
  puu[1][2] = phi * u[1] * u[2];
  puu[2][0] = phi * u[2] * u[0];

  double sumg = 0.0;
  for (i = 1; i < numvel; ++i) {    // Eq. (53) Semprebon et al.
    gi = 3. * w_lb19[i] * gamma_p * mu_p;
    gi += 3. * w_lb19[i] * phi * (u[0] * e19[i][0] + u[1] * e19[i][1] + u[2] * e19[i][2]);
    gi += 9. / 2. * w_lb19[i] *
        ((puu[0][0] * e19[i][0] + 2. * puu[0][1] * e19[i][1]) * e19[i][0] +
         (puu[1][1] * e19[i][1] + 2. * puu[1][2] * e19[i][2]) * e19[i][1] +
         (puu[2][2] * e19[i][2] + 2. * puu[2][0] * e19[i][0]) * e19[i][2]);
    gi -= 3. / 2. * w_lb19[i] * (puu[0][0] + puu[1][1] + puu[2][2]);
    geq[x][y][z][i] = gi;
    sumg += gi;
  }
  geq[x][y][z][0] = phi - sumg;
}

void FixLbMulticomponent::calc_keq(int x, int y, int z)
{
  const double psi = psi_lb[x][y][z];
  const double mu_s = mu_psi[x][y][z];
  const double *u = u_lb[x][y][z];
  double ki, puu[3][3];
  int i;

  puu[0][0] = psi * u[0] * u[0];
  puu[1][1] = psi * u[1] * u[1];
  puu[2][2] = psi * u[2] * u[2];
  puu[0][1] = psi * u[0] * u[1];
  puu[1][2] = psi * u[1] * u[2];
  puu[2][0] = psi * u[2] * u[0];

  double sumk = 0.0;
  for (i = 1; i < numvel; ++i) {    // Eq. (54) Semprebon et al.
    ki = 3. * w_lb19[i] * gamma_s * mu_s;
    ki += 3. * w_lb19[i] * psi * (u[0] * e19[i][0] + u[1] * e19[i][1] + u[2] * e19[i][2]);
    ki += 9. / 2. * w_lb19[i] *
        ((puu[0][0] * e19[i][0] + 2. * puu[0][1] * e19[i][1]) * e19[i][0] +
         (puu[1][1] * e19[i][1] + 2. * puu[1][2] * e19[i][2]) * e19[i][1] +
         (puu[2][2] * e19[i][2] + 2. * puu[2][0] * e19[i][0]) * e19[i][2]);
    ki -= 3. / 2. * w_lb19[i] * (puu[0][0] + puu[1][1] + puu[2][2]);
    keq[x][y][z][i] = ki;
    sumk += ki;
  }
  keq[x][y][z][0] = psi - sumk;
}

void FixLbMulticomponent::calc_moments_full()
{
  for (int x = halo_extent[0]; x < subNbx - halo_extent[0]; x++) {
    for (int y = halo_extent[1]; y < subNby - halo_extent[1]; y++) {
      for (int z = halo_extent[2]; z < subNbz - halo_extent[2]; z++) { calc_moments(x, y, z); }
    }
  }
}

// homogeneous mixture of C1, C2, and C3 with random concentration fluctuations
void FixLbMulticomponent::init_mixture()
{
  double rho, phi, psi;
  double C1_init, C2_init, C3_init;
  double C1tot = 0., C2tot = 0., C3tot = 0.;
  double C1tot_global = 0., C2tot_global = 0., C3tot_global = 0.;
  int x, y, z, i;

  RanMars *random = new RanMars(lmp, seed + comm->me);

  for (x = halo_extent[0]; x < subNbx - halo_extent[0]; x++) {
    for (y = halo_extent[1]; y < subNby - halo_extent[1]; y++) {
      for (z = halo_extent[2]; z < subNbz - halo_extent[2]; z++) {
        C1_init = C1 + 0.01 * random->gaussian();
        C2_init = C2 + 0.01 * random->gaussian();
        C3_init = 1.0 - C1_init - C2_init;
        rho = densityinit;
        phi = densityinit * (C1_init - C2_init);
        psi = densityinit * C3_init;
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho;
          g_lb[x][y][z][i] = w_lb19[i] * phi;
          k_lb[x][y][z][i] = w_lb19[i] * psi;
        }
        C1tot += C1_init;
        C2tot += C2_init;
        C3tot += C3_init;
      }
    }
  }

  MPI_Reduce(&C1tot, &C1tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C2tot, &C2tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C3tot, &C3tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  double vol = Nbx * Nby * Nbz;
  if (comm->me == 0) {
    error->message(FLERR, "Initialized ternary mixture with <C1> = {:f}, <C2> = {:f}, <C3> = {:f}",
                   C1tot_global / vol, C2tot_global / vol, C3tot_global / vol);
  }

  delete (random);
}

// droplet composed of component C1 and C2 (C3=0)
void FixLbMulticomponent::init_droplet(double radius)
{
  double rho = 1.0, phi, psi = 0.0;
  double pos[3], r2;
  int x, y, z, i;

  for (x = 0; x < subNbx; x++) {
    pos[0] = domain->sublo[0] + (x - halo_extent[0]) * dx_lb;
    for (y = 0; y < subNby; y++) {
      pos[1] = domain->sublo[1] + (y - halo_extent[1]) * dx_lb;
      for (z = 0; z < subNbz; z++) {
        pos[2] = domain->sublo[2] + (z - halo_extent[2]) * dx_lb;
        r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        phi = r2 < radius * radius ? 1.0 - 0.000125 : -1.0 + 0.000125;
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho * densityinit;
          g_lb[x][y][z][i] = w_lb19[i] * phi * densityinit;
          k_lb[x][y][z][i] = w_lb19[i] * psi * densityinit;
        }
      }
    }
  }
}

// liquid lens of component C3 between layers of C1 and C2
void FixLbMulticomponent::init_liquid_lens(double radius)
{
  double rho = 1.0, phi, psi;
  double pos[3], r2;
  int x, y, z, i;

  for (x = 0; x < subNbx; x++) {
    pos[0] = domain->sublo[0] + (x - halo_extent[0]) * dx_lb;
    for (y = 0; y < subNby; y++) {
      pos[1] = domain->sublo[1] + (y - halo_extent[1]) * dx_lb;
      for (z = 0; z < subNbz; z++) {
        pos[2] = domain->sublo[2] + (z - halo_extent[2]) * dx_lb;
        r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        if (r2 < radius * radius) {
          phi = 0.0;
          psi = 1.0;
        } else if (pos[2] > 0) {
          phi = 1.0;
          psi = 0.0;
        } else {
          phi = -1.0;
          psi = 0.0;
        }
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho * densityinit;
          g_lb[x][y][z][i] = w_lb19[i] * phi * densityinit;
          k_lb[x][y][z][i] = w_lb19[i] * psi * densityinit;
        }
      }
    }
  }
}

// double emulsion droplet of C1 and C2 surrounded by C3
void FixLbMulticomponent::init_double_emulsion(double radius)
{
  double rho = 1.0, phi, psi;
  double pos[3], r2;
  int x, y, z, i;

  for (x = 0; x < subNbx; x++) {
    pos[0] = domain->sublo[0] + (x - halo_extent[0]) * dx_lb;
    for (y = 0; y < subNby; y++) {
      pos[1] = domain->sublo[1] + (y - halo_extent[1]) * dx_lb;
      for (z = 0; z < subNbz; z++) {
        pos[2] = domain->sublo[2] + (z - halo_extent[2]) * dx_lb;
        r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
        if (r2 > radius * radius) {
          phi = 0.0;
          psi = 1.0;
        } else if (pos[0] < 0) {
          phi = 1.0;
          psi = 0.0;
        } else {
          phi = -1.0;
          psi = 0.0;
        }
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho * densityinit;
          g_lb[x][y][z][i] = w_lb19[i] * phi * densityinit;
          k_lb[x][y][z][i] = w_lb19[i] * psi * densityinit;
        }
      }
    }
  }
}

void FixLbMulticomponent::init_film(double thickness, double C1_film, double C2_film)
{
  double rho, phi, psi;
  double C1_init, C2_init, C3_init;
  double C1tot = 0., C2tot = 0., C3tot = 0.;
  double C1tot_global = 0., C2tot_global = 0., C3tot_global = 0.;
  double pos[3];
  int x, y, z, i;

  RanMars *random = new RanMars(lmp, seed + comm->me);

  double filmlo = domain->boxlo[1] + 0.5 * (1.0 - thickness) * domain->yprd;
  double filmhi = domain->boxhi[1] - 0.5 * (1.0 - thickness) * domain->yprd;

  for (x = halo_extent[0]; x < subNbx - halo_extent[0]; x++) {
    for (y = halo_extent[1]; y < subNby - halo_extent[1]; y++) {
      pos[1] = domain->sublo[1] + (y - halo_extent[1]) * dx_lb;
      for (z = halo_extent[2]; z < subNbz - halo_extent[2]; z++) {
        if (pos[1] > filmlo && pos[1] < filmhi) {
          C1_init = C1_film + 0.01 * random->gaussian();
          C2_init = C2_film + 0.01 * random->gaussian();
          C3_init = 1.0 - C1_init - C2_init;
        } else {
          C1_init = C1 + 0.01 * random->gaussian();
          C2_init = C2 + 0.01 * random->gaussian();
          C3_init = 1.0 - C1_init - C2_init;
        }
        rho = densityinit;
        phi = densityinit * (C1_init - C2_init);
        psi = densityinit * C3_init;
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho;
          g_lb[x][y][z][i] = w_lb19[i] * phi;
          k_lb[x][y][z][i] = w_lb19[i] * psi;
        }
        C1tot += C1_init;
        C2tot += C2_init;
        C3tot += C3_init;
      }
    }
  }

  MPI_Reduce(&C1tot, &C1tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C2tot, &C2tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C3tot, &C3tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  double vol = Nbx * Nby * Nbz;
  if (comm->me == 0) {
    error->message(FLERR, "Initialized ternary film with <C1> = {:f}, <C2> = {:f}, <C3> = {:f}",
                   C1tot_global / vol, C2tot_global / vol, C3tot_global / vol);
  }

  delete (random);
}

// mixed droplet of component C1 and C2 within pure C3
void FixLbMulticomponent::init_mixed_droplet(double radius, double C1, double C2)
{
  double rho = 1.0, C1_init, C2_init, C3_init, phi, psi;
  double C1tot = 0., C2tot = 0., C3tot = 0.;
  double C1tot_global = 0., C2tot_global = 0., C3tot_global = 0.;
  double pos[3], r2;
  int x, y, z, i;

  RanMars *random = new RanMars(lmp, seed + comm->me);

  for (x = 0; x < subNbx; x++) {
    pos[0] = domain->sublo[0] + (x - halo_extent[0]) * dx_lb;
    for (y = 0; y < subNby; y++) {
      pos[1] = domain->sublo[1] + (y - halo_extent[1]) * dx_lb;
      for (z = 0; z < subNbz; z++) {
        pos[2] = domain->sublo[2] + (z - halo_extent[2]) * dx_lb;
        r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];

        if (r2 < radius * radius) {
          C1_init = C1 + 0.01 * random->gaussian();
          C2_init = 1. - C1_init;
          C3_init = 0.0;
        } else {
          C1_init = 0.0;
          C2_init = 0.0;
          C3_init = 1.0;
        }
        rho = densityinit;
        phi = densityinit * (C1_init - C2_init);
        psi = densityinit * C3_init;
        for (i = 0; i < numvel; i++) {
          f_lb[x][y][z][i] = w_lb19[i] * rho * densityinit;
          g_lb[x][y][z][i] = w_lb19[i] * phi * densityinit;
          k_lb[x][y][z][i] = w_lb19[i] * psi * densityinit;
        }
        C1tot += C1_init;
        C2tot += C2_init;
        C3tot += C3_init;
      }
    }
  }

  MPI_Reduce(&C1tot, &C1tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C2tot, &C2tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(&C3tot, &C3tot_global, 1, MPI_DOUBLE, MPI_SUM, 0, world);
  double vol = Nbx * Nby * Nbz;
  if (comm->me == 0) {
    error->message(FLERR, "Initialized mixed droplet with <C1> = {:f}, <C2> = {:f}, <C3> = {:f}",
                   C1tot_global / vol, C2tot_global / vol, C3tot_global / vol);
  }

  delete (random);
}

void FixLbMulticomponent::init_fluid()
{

  switch (init_method) {
    case MIXTURE:
      init_mixture();
      break;
    case DROPLET:
      init_droplet(radius * dx_lb);
      break;
    case LIQUID_LENS:
      init_liquid_lens(radius * dx_lb);
      break;
    case DOUBLE_EMULSION:
      init_double_emulsion(radius * dx_lb);
      break;
    case FILM:
      init_film(thickness, C1_film, C2_film);
      break;
    case MIXED_DROPLET:
      init_mixed_droplet(radius, C1_drop, C2_drop);
      break;
  }
}

void FixLbMulticomponent::halo_comm(int dir)
{
  int tag_low = 15, tag_high = 25;
  for (int i = 0; i < 12; ++i) requests[i] = MPI_REQUEST_NULL;
  switch (dir) {
    case 2:
      MPI_Isend(&f_lb[2][2][2][0], 2, passzf, comm->procneigh[2][0], tag_low, world, &requests[0]);
      MPI_Irecv(&f_lb[2][2][0][0], 2, passzf, comm->procneigh[2][0], tag_high, world, &requests[1]);
      MPI_Isend(&f_lb[2][2][subNbz - 4][0], 2, passzf, comm->procneigh[2][1], tag_high, world,
                &requests[2]);
      MPI_Irecv(&f_lb[2][2][subNbz - 2][0], 2, passzf, comm->procneigh[2][1], tag_low, world,
                &requests[3]);

      MPI_Isend(&g_lb[2][2][2][0], 2, passzf, comm->procneigh[2][0], tag_low, world, &requests[4]);
      MPI_Irecv(&g_lb[2][2][0][0], 2, passzf, comm->procneigh[2][0], tag_high, world, &requests[5]);
      MPI_Isend(&g_lb[2][2][subNbz - 4][0], 2, passzf, comm->procneigh[2][1], tag_high, world,
                &requests[6]);
      MPI_Irecv(&g_lb[2][2][subNbz - 2][0], 2, passzf, comm->procneigh[2][1], tag_low, world,
                &requests[7]);

      MPI_Isend(&k_lb[2][2][2][0], 2, passzf, comm->procneigh[2][0], tag_low, world, &requests[8]);
      MPI_Irecv(&k_lb[2][2][0][0], 2, passzf, comm->procneigh[2][0], tag_high, world, &requests[9]);
      MPI_Isend(&k_lb[2][2][subNbz - 4][0], 2, passzf, comm->procneigh[2][1], tag_high, world,
                &requests[10]);
      MPI_Irecv(&k_lb[2][2][subNbz - 2][0], 2, passzf, comm->procneigh[2][1], tag_low, world,
                &requests[11]);
      break;
    case 1:
      MPI_Isend(&f_lb[2][2][0][0], 2, passyf, comm->procneigh[1][0], tag_low, world, &requests[0]);
      MPI_Irecv(&f_lb[2][0][0][0], 2, passyf, comm->procneigh[1][0], tag_high, world, &requests[1]);
      MPI_Isend(&f_lb[2][subNby - 4][0][0], 2, passyf, comm->procneigh[1][1], tag_high, world,
                &requests[2]);
      MPI_Irecv(&f_lb[2][subNby - 2][0][0], 2, passyf, comm->procneigh[1][1], tag_low, world,
                &requests[3]);

      MPI_Isend(&g_lb[2][2][0][0], 2, passyf, comm->procneigh[1][0], tag_low, world, &requests[4]);
      MPI_Irecv(&g_lb[2][0][0][0], 2, passyf, comm->procneigh[1][0], tag_high, world, &requests[5]);
      MPI_Isend(&g_lb[2][subNby - 4][0][0], 2, passyf, comm->procneigh[1][1], tag_high, world,
                &requests[6]);
      MPI_Irecv(&g_lb[2][subNby - 2][0][0], 2, passyf, comm->procneigh[1][1], tag_low, world,
                &requests[7]);

      MPI_Isend(&k_lb[2][2][0][0], 2, passyf, comm->procneigh[1][0], tag_low, world, &requests[8]);
      MPI_Irecv(&k_lb[2][0][0][0], 2, passyf, comm->procneigh[1][0], tag_high, world, &requests[9]);
      MPI_Isend(&k_lb[2][subNby - 4][0][0], 2, passyf, comm->procneigh[1][1], tag_high, world,
                &requests[10]);
      MPI_Irecv(&k_lb[2][subNby - 2][0][0], 2, passyf, comm->procneigh[1][1], tag_low, world,
                &requests[11]);
      break;
    case 0:
      MPI_Isend(&f_lb[2][0][0][0], 2, passxf, comm->procneigh[0][0], tag_low, world, &requests[0]);
      MPI_Irecv(&f_lb[0][0][0][0], 2, passxf, comm->procneigh[0][0], tag_high, world, &requests[1]);
      MPI_Isend(&f_lb[subNbx - 4][0][0][0], 2, passxf, comm->procneigh[0][1], tag_high, world,
                &requests[2]);
      MPI_Irecv(&f_lb[subNbx - 2][0][0][0], 2, passxf, comm->procneigh[0][1], tag_low, world,
                &requests[3]);

      MPI_Isend(&g_lb[2][0][0][0], 2, passxf, comm->procneigh[0][0], tag_low, world, &requests[4]);
      MPI_Irecv(&g_lb[0][0][0][0], 2, passxf, comm->procneigh[0][0], tag_high, world, &requests[5]);
      MPI_Isend(&g_lb[subNbx - 4][0][0][0], 2, passxf, comm->procneigh[0][1], tag_high, world,
                &requests[6]);
      MPI_Irecv(&g_lb[subNbx - 2][0][0][0], 2, passxf, comm->procneigh[0][1], tag_low, world,
                &requests[7]);

      MPI_Isend(&k_lb[2][0][0][0], 2, passxf, comm->procneigh[0][0], tag_low, world, &requests[8]);
      MPI_Irecv(&k_lb[0][0][0][0], 2, passxf, comm->procneigh[0][0], tag_high, world, &requests[9]);
      MPI_Isend(&k_lb[subNbx - 4][0][0][0], 2, passxf, comm->procneigh[0][1], tag_high, world,
                &requests[10]);
      MPI_Irecv(&k_lb[subNbx - 2][0][0][0], 2, passxf, comm->procneigh[0][1], tag_low, world,
                &requests[11]);
      break;
  }
}

void FixLbMulticomponent::halo_wait()
{
  MPI_Waitall(numrequests, requests, MPI_STATUS_IGNORE);
}

void FixLbMulticomponent::halo_comm()
{
  halo_comm(2);
  halo_wait();
  halo_comm(1);
  halo_wait();
  halo_comm(0);
  halo_wait();
}

void FixLbMulticomponent::init_halo()
{

  // Create MPI datatypes to pass the f,g,j and feq,geq,keq arrays
  int size;
  MPI_Aint lb, extent;
  MPI_Datatype slice[3];

  MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);

  MPI_Type_free(&passxf);
  MPI_Type_free(&passyf);
  MPI_Type_free(&passzf);

  MPI_Type_vector(subNbz, numvel, numvel, MPI_DOUBLE, &oneslice);
  MPI_Type_create_hvector(subNby, 1, numvel * subNbz * extent, oneslice, &slice[0]);
  MPI_Type_create_resized(slice[0], 0, subNby * subNbz * numvel * extent, &passxf);

  MPI_Type_create_hvector(subNbx - 4, 1, numvel * subNby * subNbz * extent, oneslice, &slice[1]);
  MPI_Type_create_resized(slice[1], 0, subNbz * numvel * extent, &passyf);

  MPI_Type_vector(subNby - 4, numvel, numvel * subNbz, MPI_DOUBLE, &oneslice);
  MPI_Type_create_hvector(subNbx - 4, 1, numvel * subNby * subNbz * extent, oneslice, &slice[2]);
  MPI_Type_create_resized(slice[2], 0, numvel * extent, &passzf);

  MPI_Type_commit(&passxf);
  MPI_Type_commit(&passzf);
  MPI_Type_commit(&passyf);
}

void FixLbMulticomponent::destroy_halo()
{
  // MPI datatypes are freed in parent destructor
}

void FixLbMulticomponent::dump_xdmf(const int step)
{
  if (dump_interval && step % dump_interval == 0) {
    calc_moments_full();
    // Write XDMF grid entry for time step
    if (me == 0) {
      long int block = (long int) fluid_global_n0[0] * fluid_global_n0[1] * fluid_global_n0[2] *
          sizeof(MPI_DOUBLE);
      long int offset = (step / dump_interval) * block *
          (4 +
           3); /* This should be changed to account for dumps actually written.  This offset could malfunction on a restart. */
      double time = update->ntimestep * dt_lb;

      fprintf(dump_file_handle_xdmf,
              "      <Grid Name=\"%d\">\n"
              "        <Time Value=\"%f\"/>\n\n"
              "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n"
              "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
              "          <DataItem Dimensions=\"3\">\n"
              "            %f %f %f\n"
              "          </DataItem>\n"
              "          <DataItem Dimensions=\"3\">\n"
              "            %f %f %f\n"
              "          </DataItem>\n"
              "        </Geometry>\n\n",
              step, time, fluid_global_n0[2], fluid_global_n0[1], fluid_global_n0[0],
              domain->boxlo[2], domain->boxlo[1], domain->boxlo[0], dx_lb, dx_lb, dx_lb);
      fprintf(dump_file_handle_xdmf,
              "        <Attribute Name=\"density\">\n"
              "          <DataItem ItemType=\"Function\" Function=\"$0 * %f\" Dimensions=\"%d %d "
              "%d\">\n"
              "            <DataItem Precision=\"%zd\" Format=\"Binary\" Seek=\"%ld\" "
              "Dimensions=\"%d %d %d\">\n"
              "              %s\n"
              "            </DataItem>\n"
              "          </DataItem>\n"
              "        </Attribute>\n\n",
              dm_lb / (dx_lb * dx_lb * dx_lb), fluid_global_n0[2], fluid_global_n0[1],
              fluid_global_n0[0], sizeof(MPI_DOUBLE), offset, fluid_global_n0[2],
              fluid_global_n0[1], fluid_global_n0[0], dump_file_name_raw.c_str());
      fprintf(dump_file_handle_xdmf,
              "        <Attribute Name=\"phi\">\n"
              "          <DataItem ItemType=\"Function\" Function=\"$0 * %f\" Dimensions=\"%d %d "
              "%d\">\n"
              "            <DataItem Precision=\"%zd\" Format=\"Binary\" Seek=\"%ld\" "
              "Dimensions=\"%d %d %d\">\n"
              "              %s\n"
              "            </DataItem>\n"
              "          </DataItem>\n"
              "        </Attribute>\n\n",
              dm_lb / (dx_lb * dx_lb * dx_lb), fluid_global_n0[2], fluid_global_n0[1],
              fluid_global_n0[0], sizeof(MPI_DOUBLE), offset + block * 1, fluid_global_n0[2],
              fluid_global_n0[1], fluid_global_n0[0], dump_file_name_raw.c_str());
      fprintf(dump_file_handle_xdmf,
              "        <Attribute Name=\"psi\">\n"
              "          <DataItem ItemType=\"Function\" Function=\"$0 * %f\" Dimensions=\"%d %d "
              "%d\">\n"
              "            <DataItem Precision=\"%zd\" Format=\"Binary\" Seek=\"%ld\" "
              "Dimensions=\"%d %d %d\">\n"
              "              %s\n"
              "            </DataItem>\n"
              "          </DataItem>\n"
              "        </Attribute>\n\n",
              dm_lb / (dx_lb * dx_lb * dx_lb), fluid_global_n0[2], fluid_global_n0[1],
              fluid_global_n0[0], sizeof(MPI_DOUBLE), offset + block * 2, fluid_global_n0[2],
              fluid_global_n0[1], fluid_global_n0[0], dump_file_name_raw.c_str());
      fprintf(dump_file_handle_xdmf,
              "        <Attribute Name=\"pressure\">\n"
              "          <DataItem ItemType=\"Function\" Function=\"$0 * %f\" Dimensions=\"%d %d "
              "%d\">\n"
              "            <DataItem Precision=\"%zd\" Format=\"Binary\" Seek=\"%ld\" "
              "Dimensions=\"%d %d %d\">\n"
              "              %s\n"
              "            </DataItem>\n"
              "          </DataItem>\n"
              "        </Attribute>\n\n",
              dm_lb / (dx_lb * dx_lb * dx_lb), fluid_global_n0[2], fluid_global_n0[1],
              fluid_global_n0[0], sizeof(MPI_DOUBLE), offset + block * 3, fluid_global_n0[2],
              fluid_global_n0[1], fluid_global_n0[0], dump_file_name_raw.c_str());
      fprintf(dump_file_handle_xdmf,
              "        <Attribute Name=\"velocity\" AttributeType=\"Vector\">\n"
              "          <DataItem ItemType=\"Function\" Function=\"$0 * %f\" Dimensions=\"%d %d "
              "%d 3\">\n"
              "            <DataItem Precision=\"%zd\" Format=\"Binary\" Seek=\"%ld\" "
              "Dimensions=\"%d %d %d 3\">\n"
              "              %s\n"
              "            </DataItem>\n"
              "          </DataItem>\n"
              "        </Attribute>\n\n",
              dx_lb / dt_lb, fluid_global_n0[2], fluid_global_n0[1], fluid_global_n0[0],
              sizeof(MPI_DOUBLE), offset + block * 4, fluid_global_n0[2], fluid_global_n0[1],
              fluid_global_n0[0], dump_file_name_raw.c_str());
      fprintf(dump_file_handle_xdmf, "      </Grid>\n\n");
    }

    // Write raw data
    {
      int lbox[3];
      lbox[0] = subNbx;
      lbox[1] = subNby;
      lbox[2] = subNbz;

      const size_t lvol = lbox[0] * lbox[1] * lbox[2];

      // Transpose local arrays to fortran-order for paraview output
      std::vector<double> density_2_fort(lvol);
      std::vector<double> phi_2_fort(lvol);
      std::vector<double> psi_2_fort(lvol);
      std::vector<double> pressure_2_fort(lvol);
      std::vector<double> velocity_2_fort(lvol * 3);
      int indexf = 0;
      for (int k = 0; k < lbox[2]; k++) {
        for (int j = 0; j < lbox[1]; j++) {
          for (int i = 0; i < lbox[0]; i++) {
            indexf = i + lbox[0] * (j + lbox[1] * k);
            density_2_fort[indexf] = density_lb[i][j][k];
            phi_2_fort[indexf] = phi_lb[i][j][k];
            psi_2_fort[indexf] = psi_lb[i][j][k];
            pressure_2_fort[indexf] = pressure_lb[i][j][k];
            velocity_2_fort[0 + 3 * indexf] = u_lb[i][j][k][0];
            velocity_2_fort[1 + 3 * indexf] = u_lb[i][j][k][1];
            velocity_2_fort[2 + 3 * indexf] = u_lb[i][j][k][2];
          }
        }
      }

      MPI_File_write_all(dump_file_handle_raw, &density_2_fort[0], 1, fluid_scalar_field_mpitype,
                         MPI_STATUS_IGNORE);
      MPI_File_write_all(dump_file_handle_raw, &phi_2_fort[0], 1, fluid_scalar_field_mpitype,
                         MPI_STATUS_IGNORE);
      MPI_File_write_all(dump_file_handle_raw, &psi_2_fort[0], 1, fluid_scalar_field_mpitype,
                         MPI_STATUS_IGNORE);
      MPI_File_write_all(dump_file_handle_raw, &pressure_2_fort[0], 1, fluid_scalar_field_mpitype,
                         MPI_STATUS_IGNORE);
      MPI_File_write_all(dump_file_handle_raw, &velocity_2_fort[0], 1, fluid_vector_field_mpitype,
                         MPI_STATUS_IGNORE);
    }
  }
}

static MPI_Datatype mpiTypeGlobalWrite(const int local_ghost, const int *local_size,
                                       const int *global_offset, const int global_ghost,
                                       const int *global_size, const MPI_Datatype mpitype)
{
  MPI_Datatype global_mpitype;

  {
    bool endpoint_lower[] = {global_offset[0] == 0, global_offset[1] == 0, global_offset[2] == 0};
    bool endpoint_upper[] = {global_offset[0] + local_size[0] == global_size[0],
                             global_offset[1] + local_size[1] == global_size[1],
                             global_offset[2] + local_size[2] == global_size[2]};

    int sizes[] = {global_ghost + global_size[0] + global_ghost,
                   global_ghost + global_size[1] + global_ghost,
                   global_ghost + global_size[2] + global_ghost};
    int subsizes[] = {
        (global_ghost * endpoint_lower[0] + local_size[0] + global_ghost * endpoint_upper[0]),
        (global_ghost * endpoint_lower[1] + local_size[1] + global_ghost * endpoint_upper[1]),
        (global_ghost * endpoint_lower[2] + local_size[2] + global_ghost * endpoint_upper[2])};
    int starts[] = {global_ghost * !endpoint_lower[0] + global_offset[0],
                    global_ghost * !endpoint_lower[1] + global_offset[1],
                    global_ghost * !endpoint_lower[2] + global_offset[2]};

    // Note Fortran ordering as we switch order for paraview output
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, mpitype,
                             &global_mpitype);
  }

  return global_mpitype;
}

static MPI_Datatype mpiTypeLocalWrite(const int local_ghost, const int *local_size,
                                      const int *global_offset, const int global_ghost,
                                      const int *global_size, const MPI_Datatype mpitype)
{
  MPI_Datatype local_mpitype;

  {
    bool endpoint_lower[] = {global_offset[0] == 0, global_offset[1] == 0, global_offset[2] == 0};
    bool endpoint_upper[] = {global_offset[0] + local_size[0] == global_size[0],
                             global_offset[1] + local_size[1] == global_size[1],
                             global_offset[2] + local_size[2] == global_size[2]};

    int sizes[] = {local_ghost + local_size[0] + local_ghost,
                   local_ghost + local_size[1] + local_ghost,
                   local_ghost + local_size[2] + local_ghost};
    int subsizes[] = {
        (global_ghost * endpoint_lower[0] + local_size[0] + global_ghost * endpoint_upper[0]),
        (global_ghost * endpoint_lower[1] + local_size[1] + global_ghost * endpoint_upper[1]),
        (global_ghost * endpoint_lower[2] + local_size[2] + global_ghost * endpoint_upper[2])};
    int starts[] = {local_ghost - global_ghost * endpoint_lower[0],
                    local_ghost - global_ghost * endpoint_lower[1],
                    local_ghost - global_ghost * endpoint_lower[2]};

    // Fortran ordering for paraview output
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, mpitype,
                             &local_mpitype);
  }

  return local_mpitype;
}

static MPI_Datatype mpiTypeDumpGlobal_ternary(const int *local_size, const int *global_offset,
                                              const int *global_size)
{
  MPI_Datatype dump_ternary;

  // MPI types for local position of the global dump file
  {
    MPI_Datatype real3_mpitype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &real3_mpitype);
    MPI_Type_commit(&real3_mpitype);

    // Scalar and vector types for the local chunk of the global file
    MPI_Datatype scalar_mpitype =
        mpiTypeGlobalWrite(2, local_size, global_offset, 0, global_size, MPI_DOUBLE);
    MPI_Datatype vector_mpitype =
        mpiTypeGlobalWrite(2, local_size, global_offset, 0, global_size, real3_mpitype);

    // rho, phi, psi, pressure, velocity
    {
      MPI_Aint lb, extent;
      MPI_Type_get_extent(scalar_mpitype, &lb, &extent);

      int blocklengths[] = {1, 1, 1, 1, 1};
      MPI_Aint displacements[] = {0, lb + extent, 2 * (lb + extent), 3 * (lb + extent),
                                  4 * (lb + extent)};
      MPI_Datatype datatypes[] = {scalar_mpitype, scalar_mpitype, scalar_mpitype, scalar_mpitype,
                                  vector_mpitype};

      MPI_Type_create_struct(5, blocklengths, displacements, datatypes, &dump_ternary);
    }

    // Free local MPI types
    MPI_Type_free(&real3_mpitype);
    MPI_Type_free(&scalar_mpitype);
    MPI_Type_free(&vector_mpitype);
  }

  return dump_ternary;
}

void FixLbMulticomponent::init_output(void)
{
  fluid_global_n0[0] = Nbx + (domain->periodicity[0] == 0);
  fluid_global_n0[1] = Nby + (domain->periodicity[1] == 0);
  fluid_global_n0[2] = Nbz + (domain->periodicity[2] == 0);

  fluid_local_n0[0] = subNbx - 2 * halo_extent[0] +
      (domain->periodicity[0] == 0 && (comm->myloc[0] == comm->procgrid[0] - 1));
  fluid_local_n0[1] = subNby - 2 * halo_extent[1] +
      (domain->periodicity[1] == 0 && (comm->myloc[1] == comm->procgrid[1] - 1));
  fluid_local_n0[2] = subNbz - 2 * halo_extent[2] +
      (domain->periodicity[2] == 0 && (comm->myloc[2] == comm->procgrid[2] - 1));

  fluid_global_o0[0] = (fluid_local_n0[0]) * comm->myloc[0];
  fluid_global_o0[1] = (fluid_local_n0[1]) * comm->myloc[1];
  fluid_global_o0[2] = (fluid_local_n0[2]) * comm->myloc[2];

  // Local write MPI types for our portion of the global dump file
  fluid_scalar_field_mpitype =
      mpiTypeLocalWrite(2, fluid_local_n0, fluid_global_o0, 0, fluid_global_n0, MPI_DOUBLE);
  fluid_vector_field_mpitype =
      mpiTypeLocalWrite(2, fluid_local_n0, fluid_global_o0, 0, fluid_global_n0, realType3_mpitype);

  // Global write MPI type for our portion of the global dump file
  dump_file_mpitype = mpiTypeDumpGlobal_ternary(fluid_local_n0, fluid_global_o0, fluid_global_n0);

  MPI_Type_commit(&fluid_scalar_field_mpitype);
  MPI_Type_commit(&fluid_vector_field_mpitype);
  MPI_Type_commit(&dump_file_mpitype);

  // Output
  if (dump_interval) {
    if (me == 0) {
      dump_file_handle_xdmf = fopen(dump_file_name_xdmf.c_str(), "w");
      if (!dump_file_handle_xdmf) {
        error->one(FLERR, "Unable to truncate/create \"{}\": {}", dump_file_name_xdmf,
                   utils::getsyserror());
      }
      fprintf(dump_file_handle_xdmf,
              "<?xml version=\"1.0\" ?>\n"
              "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
              "<Xdmf Version=\"2.0\">\n"
              "  <Domain>\n"
              "    <Grid Name=\"fluid\" GridType=\"Collection\" CollectionType=\"Temporal\">\n\n");
    }
    MPI_File_open(world, const_cast<char *>(dump_file_name_raw.c_str()),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &dump_file_handle_raw);
    MPI_File_set_size(dump_file_handle_raw, 0);
    MPI_File_set_view(dump_file_handle_raw, 0, MPI_DOUBLE, dump_file_mpitype, "native",
                      MPI_INFO_NULL);
  }
}

void FixLbMulticomponent::destroy_output()
{

  MPI_Type_free(&fluid_scalar_field_mpitype);
  MPI_Type_free(&fluid_vector_field_mpitype);
}

void FixLbMulticomponent::init_lattice()
{

  // Set halo extent to 2 for gradient calculations
  halo_extent[0] = halo_extent[1] = halo_extent[2] = 2;
  subNbx += 2 * halo_extent[0] - 2;    // -2 comes from prior inintialization in FixLbFluid
  subNby += 2 * halo_extent[1] - 2;
  subNbz += 2 * halo_extent[2] - 2;

  // Destroy memory created in FixLbFluid constructor previously
  memory->destroy(f_lb);
  memory->destroy(fnew);
  memory->destroy(feq);
  memory->destroy(density_lb);
  memory->destroy(u_lb);

  memory->create(feq, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:feq");
  memory->create(f_lb, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:f_lb");
  memory->create(fnew, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:fnew");

  memory->create(geq, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:geq");
  memory->create(g_lb, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:g_lb");
  memory->create(gnew, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:gnew");

  memory->create(keq, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:keq");
  memory->create(k_lb, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:k_lb");
  memory->create(knew, subNbx, subNby, subNbz, numvel, "FixLbMulticomponent:knew");

  memory->create(density_lb, subNbx, subNby, subNbz, "FixLbMulticomponent:rho_lb");
  memory->create(u_lb, subNbx, subNby, subNbz, 3, "FixLBMulticomponent:u_lb");
  memory->create(phi_lb, subNbx, subNby, subNbz, "FixLbMulticomponent:phi_lb");
  memory->create(psi_lb, subNbx, subNby, subNbz, "FixLbMulticomponent:psi_lb");
  memory->create(pressure_lb, subNbx, subNby, subNbz, "FixLbMulticomponent:pressure_lb");
  memory->create(mu_rho, subNbx, subNby, subNbz, "FixLbMulticomponent:mu_rho");
  memory->create(mu_phi, subNbx, subNby, subNbz, "FixLbMulticomponent:mu_phi");
  memory->create(mu_psi, subNbx, subNby, subNbz, "FixLbMulticomponent:mu_psi");

  memory->create(density_gradient, subNbx, subNby, subNbz, 3,
                 "FixLbMulticomponent:density_gradient");
  memory->create(phi_gradient, subNbx, subNby, subNbz, 3, "FixLbMulticomponent:phi_gradient");
  memory->create(psi_gradient, subNbx, subNby, subNbz, 3, "FixLbMulticomponent:psi_gradient");
  memory->create(laplace_rho, subNbx, subNby, subNbz, "FixLbMulticomponent:laplace_rho");
  memory->create(laplace_phi, subNbx, subNby, subNbz, "FixLbMulticomponent:laplace_phi");
  memory->create(laplace_psi, subNbx, subNby, subNbz, "FixLbMulticomponent:laplace_psi");
}

void FixLbMulticomponent::destroy_lattice()
{

  memory->destroy(f_lb);
  memory->destroy(g_lb);
  memory->destroy(k_lb);
  memory->destroy(fnew);
  memory->destroy(gnew);
  memory->destroy(knew);
  memory->destroy(feq);
  memory->destroy(geq);
  memory->destroy(keq);
  memory->destroy(phi_lb);
  memory->destroy(psi_lb);
  memory->destroy(pressure_lb);
  memory->destroy(mu_rho);
  memory->destroy(mu_phi);
  memory->destroy(mu_psi);
  memory->destroy(density_gradient);
  memory->destroy(phi_gradient);
  memory->destroy(psi_gradient);
  memory->destroy(laplace_phi);
  memory->destroy(laplace_rho);
  memory->destroy(laplace_psi);
}

void FixLbMulticomponent::init_parameters(int argc, char **argv)
{

  if (argc < 9)
    error->all(FLERR,
               "Illegal fix lb/multicomponent command: must start with `fix * * lb/multicomponent "
               "* * $rho D3Q19 dx 1`");

  // default parameter values
  seed = 12345;
  alpha = 1.0;
  C1 = 0.333333;
  C2 = 0.333333;
  C3 = 0.333334;    // concentrations
  kappa1 = 0.01;
  kappa2 = 0.01, kappa3 = 0.01;    // surface tensions
  tau_r = 1.0;
  tau_p = 1.0;
  tau_s = 0.666667;    // relaxation times
  gamma_p = 1.0;
  gamma_s = 1.0;            // mobility coefficients
  init_method = MIXTURE;    // initialization

  // parse optional parameters
  int argi = 9;
  while (argi < argc) {
    if (strcmp(argv[argi], "tau_r") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      tau_r = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "tau_p") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      tau_p = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "tau_s") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      tau_s = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "kappa1") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      kappa1 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "kappa2") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      kappa2 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "kappa3") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      kappa3 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "alpha") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      alpha = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "gamma_p") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      gamma_p = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "gamma_s") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      gamma_s = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "C1") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      C1 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "C2") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      C2 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "C3") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      C3 = utils::numeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "dumpxdmf") == 0) {
      if (argi + 3 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      dump_interval = utils::inumeric(FLERR, argv[argi + 1], false, lmp);
      dump_file_name_xdmf = std::string(argv[argi + 2]) + std::string(".xdmf");
      dump_file_name_raw = std::string(argv[argi + 2]) + std::string(".raw");
      argi += 3;
    } else if (strcmp(argv[argi], "seed") == 0) {
      if (argi + 2 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      seed = utils::inumeric(FLERR, argv[argi + 1], false, lmp);
      argi += 2;
    } else if (strcmp(argv[argi], "init") == 0) {
      if (argi + 1 > argc)
        error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
      argi += 1;
      if (strcmp(argv[argi], "mixture") == 0) {
        if (argi + 1 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        init_method = MIXTURE;
        argi += 1;
      } else if (strcmp(argv[argi], "droplet") == 0) {
        if (argi + 2 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        radius = utils::numeric(FLERR, argv[argi + 1], false, lmp);
        init_method = DROPLET;
        argi += 2;
      } else if (strcmp(argv[argi], "liquid_lens") == 0) {
        if (argi + 2 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        radius = utils::numeric(FLERR, argv[argi + 1], false, lmp);
        init_method = LIQUID_LENS;
        argi += 2;
      } else if (strcmp(argv[argi], "double_emulsion") == 0) {
        if (argi + 2 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        radius = utils::numeric(FLERR, argv[argi + 1], false, lmp);
        init_method = DOUBLE_EMULSION;
        argi += 2;
      } else if (strcmp(argv[argi], "film") == 0) {
        if (argi + 4 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        thickness = utils::numeric(FLERR, argv[argi + 1], false, lmp);
        C1_film = utils::numeric(FLERR, argv[argi + 2], false, lmp);
        C2_film = utils::numeric(FLERR, argv[argi + 3], false, lmp);
        init_method = FILM;
        argi += 4;
      } else if (strcmp(argv[argi], "mixed_droplet") == 0) {
        if (argi + 4 > argc)
          error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                     argv[argi]);
        radius = utils::numeric(FLERR, argv[argi + 1], false, lmp);
        C1_drop = utils::numeric(FLERR, argv[argi + 2], false, lmp);
        C2_drop = utils::numeric(FLERR, argv[argi + 3], false, lmp);
        init_method = MIXED_DROPLET;
        argi += 4;
      } else
        error->all(FLERR, "Illegal fix lb/multicomponent command: {} {}", argv[argi - 1],
                   argv[argi]);
    } else
      error->all(FLERR, "Illegal fix lb/multicomponent command: {}", argv[argi]);
  }

  kappa_rr = kappa_pp = (kappa1 + kappa2) / 4.;
  kappa_ss = (kappa1 + kappa2 + 4. * kappa3) / 4.;
  kappa_rp = (kappa1 - kappa2) / 4.;
  kappa_ps = -kappa_rp;
  kappa_rs = -kappa_rr;
}

FixLbMulticomponent::~FixLbMulticomponent()
{

  destroy_output();
  destroy_halo();
  destroy_lattice();
}

FixLbMulticomponent::FixLbMulticomponent(LAMMPS *lmp, int argc, char **argv) :
    FixLbFluid(lmp, 9, argv),    // use only the first 9 arguments to parse in FixLbFluid
    g_lb(nullptr), gnew(nullptr), geq(nullptr), k_lb(nullptr), knew(nullptr), keq(nullptr),
    phi_lb(nullptr), psi_lb(nullptr), pressure_lb(nullptr), mu_phi(nullptr), mu_psi(nullptr),
    density_gradient(nullptr), phi_gradient(nullptr), psi_gradient(nullptr), laplace_rho(nullptr),
    laplace_phi(nullptr), laplace_psi(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_lbmulticomponent);

  init_parameters(argc, argv);
  init_lattice();
  init_halo();
  init_output();
  init_fluid();
  dump_xdmf(update->ntimestep);
}
