/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Originally modified from CG-DNA/fix_nve_dotc_langevin.cpp.

   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include "fix_brownian_anisof.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "random_mars.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBrownianAnisof::FixBrownianAnisof(LAMMPS *lmp, int narg, char **arg) : FixBrownianBase(lmp, narg, arg)
{
  if (dipole_flag || gamma_t_eigen_flag || gamma_r_eigen_flag || gamma_r_flag) {
    error->all(FLERR, "Illegal fix brownian command.");
  }
  if (!gamma_t_flag) { error->all(FLERR, "Illegal fix brownian command."); }
}

/* ---------------------------------------------------------------------- */

void FixBrownianAnisof::init()
{
  FixBrownianBase::init();
  g1 /= gamma_t;
  g2 *= sqrt(gamma_t);
}

/* ---------------------------------------------------------------------- */

void FixBrownianAnisof::initial_integrate(int /*vflag */)
{
  if (domain->dimension == 2) {
    if (!noise_flag) {
      initial_integrate_templated<0, 0, 1>();
    } else if (gaussian_noise_flag) {
      initial_integrate_templated<0, 1, 1>();
    } else {
      initial_integrate_templated<1, 0, 1>();
    }
  } else {
    if (!noise_flag) {
      initial_integrate_templated<0, 0, 0>();
    } else if (gaussian_noise_flag) {
      initial_integrate_templated<0, 1, 0>();
    } else {
      initial_integrate_templated<1, 0, 0>();
    }
  }
  return;
}

/* ---------------------------------------------------------------------- */

template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D> void FixBrownianAnisof::initial_integrate_templated()
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dx, dy, dz;
  double mux, muy, muz;
  double fx, fy, fz;

  double mudotf, mumag;

  for (int i = 0; i < nlocal; i++) {
    mux = mu[i][0];
    muy = mu[i][1];
    muz = mu[i][2];
    fx = f[i][0];
    fy = f[i][1];
    fz = f[i][2];

    mudotf = mux * fx + muy * fy + muz * fz;
    mumag = sqrt(mux * mux + muy * muy + muz * muz);

    if (mask[i] & groupbit) {
      if (Tp_2D) {
        dz = 0;
        if (Tp_UNIFORM) {
          dx = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * mux + \
                          1/(1 - eps) * fx) + g2 * (rng->uniform() - 0.5));
          dy = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * muy + \
                          1/(1 - eps) * fy) + g2 * (rng->uniform() - 0.5));
        } else if (Tp_GAUSS) {
          dx = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * mux + \
                          1/(1 - eps) * fx) + g2 * (rng->gaussian()));
          dy = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * muy + \
                          1/(1 - eps) * fy) + g2 * (rng->gaussian()));
        } else {
          dx = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * mux + \
                          1/(1 - eps) * fx));
          dy = dt * (g1 * ((1 - 1/(1 - eps)) * mudotf / mumag * muy + \
                          1/(1 - eps) * fy));
        }
      } else {
        if (Tp_UNIFORM) {
          dx = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * mux + \
                          fx) + g2 * (rng->uniform() - 0.5));
          dy = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muy + \
                          fy) + g2 * (rng->uniform() - 0.5));
          dz = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muz + \
                          fz) + g2 * (rng->uniform() - 0.5));
        } else if (Tp_GAUSS) {
          dx = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * mux + \
                          fx) + g2 * (rng->gaussian()));
          dy = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muy + \
                          fy) + g2 * (rng->gaussian()));
          dz = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muz + \
                          fz) + g2 * (rng->gaussian()));
        } else {
          dx = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * mux + \
                          fx));
          dy = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muy + \
                          fy));
          dz = dt * (g1 * ((1/(1 - 0.8) - 1) * mudotf / (mumag*mumag) * muz + \
                          fz));
        }
/*     
        dx = dt * g1 * fx;
        dy = dt * g1 * fy;
        dz = dt * g1 * fz;	
*/    }

      x[i][0] += dx;
      v[i][0] = dx / dt;

      x[i][1] += dy;
      v[i][1] = dy / dt;

      x[i][2] += dz;
      v[i][2] = dz / dt;
    }
  }
  return;
}
