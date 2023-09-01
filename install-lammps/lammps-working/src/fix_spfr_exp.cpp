// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Katherine Copenhagen
------------------------------------------------------------------------- */

#include "fix_spfr_exp.h"

#include "comm.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "error.h"
#include "math_const.h"
#include "respa.h"
#include "update.h"
#include "random_mars.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixSPFRexp::FixSPFRexp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11 ) error->all(FLERR,"Illegal fix spfr command");

  af_sig = utils::numeric(FLERR, arg[3], false, lmp);
  
  seed = utils::inumeric(FLERR, arg[4], false, lmp);
  if (seed <= 0) error->all(FLERR, "Fix spfr seed must be > 0.");

  rev_mu = utils::numeric(FLERR, arg[5], false, lmp);
  rev_del = utils::numeric(FLERR, arg[6], false, lmp);

  noise_r = utils::numeric(FLERR, arg[7], false, lmp);
  noise_t = utils::numeric(FLERR, arg[8], false, lmp);

  k_agar = utils::numeric(FLERR, arg[9], false, lmp);
  zc = utils::numeric(FLERR, arg[10], false, lmp);

  rng = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

int FixSPFRexp::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPFRexp::init()
{
  int i, i1, i2, n;
  int nbondlist=neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;
  double nmumag, caf;
  double *ttr = atom->ttr;
  double **x = atom->x;
  double *af = atom->af;
  int *mid = atom->molecule;
  double **nmu = atom->nmu;
  double **mu = atom->mu;
  double **amu = atom->amu;
  int *mask = atom->mask;
  double delx, dely, delz;
  double rsq, r;
  double dt = update->dt;
  int dim = domain->dimension;

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      ttr[i] = rng->uniform()*rev_mu/dt;
      nmu[i][0] = rng->gaussian(0.0,1.0);
      nmu[i][1] = rng->gaussian(0.0,1.0);
      if (dim == 3)
        nmu[i][2] = rng->gaussian(0.0,1.0);
      else
        nmu[i][2] = 0;

      nmumag = sqrt(nmu[i][0]*nmu[i][0] + 
          nmu[i][1]*nmu[i][1] + nmu[i][2]*nmu[i][2]);
      if (nmumag > 0) {
        nmu[i][0] = nmu[i][0]/nmumag;
        nmu[i][1] = nmu[i][1]/nmumag;
        nmu[i][2] = nmu[i][2]/nmumag;
      }
      amu[i][0] = mu[i][0];
      amu[i][1] = mu[i][1];
      amu[i][2] = mu[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPFRexp::pre_force(int /*vflag*/)
{

  int i1, i2, i, j, n;
  double delx, dely, delz;
  double rsq, r; 
  double mux, muy, muz, mumagsq, mumag;
  double amux, amuy, amuz, amumagsq, amumag;
  double n_x, n_y, n_z, n_mag;
  double nmumag;
  double caf;
  double dt = update->dt;

  double **x = atom->x;
  double **f = atom->f;
  double **mu = atom->mu;
  double **nmu = atom->nmu;
  double **amu = atom->amu;
  double *af = atom->af;
  double *ttr = atom->ttr;
  tagint *tag = atom->tag;

  int *mid = atom->molecule;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int *num_bond = atom->num_bond;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int dim = domain->dimension;

  double rev_p = (1/rev_mu)*dt;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      // Check for reversals.
      if (atom->map(mid[i]) == i) {
        ttr[i]+=1;
        if (ttr[i] >= rev_del/dt) {
          if (rng->uniform() < rev_p) {
            if (af[i] < 0)
              af[i] = rng->gaussian(1.0, af_sig);
            else
              af[i] = -rng->gaussian(1.0, af_sig);
	    ttr[i] = 0;
	  }
        }
      }
    // Clear dipole.
      mu[i][0] = 0;
      mu[i][1] = 0;
      mu[i][2] = 0;
      amu[i][0] = 0.0;
      amu[i][1] = 0.0;
      amu[i][2] = 0.0;

      //Move noise vector.
      n_x = rng->gaussian(0.0,1.0);
      n_y = rng->gaussian(0.0,1.0);
      if (dim == 3)
        n_z = rng->gaussian(0.0,1.0);
      else
        n_z = 0;

      n_mag = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);
      if (n_mag > 0) {
        n_x = noise_t*dt*n_x/n_mag;
        n_y = noise_t*dt*n_y/n_mag;
        n_z = noise_t*dt*n_z/n_mag;
      } else {
        n_x = 0;
        n_y = 0;
        n_z = 0;
      }
      nmu[i][0] = nmu[i][0]+n_x;
      nmu[i][1] = nmu[i][1]+n_y;
      nmu[i][2] = nmu[i][2]+n_z;

      nmumag = sqrt(nmu[i][0]*nmu[i][0] + 
          nmu[i][1]*nmu[i][1] + nmu[i][2]*nmu[i][2]);
      if (nmumag > 0) {
        nmu[i][0] = nmu[i][0]/nmumag;
        nmu[i][1] = nmu[i][1]/nmumag;
        nmu[i][2] = nmu[i][2]/nmumag;
      }
    }
  }
  for (n = 0; n < nbondlist; n++) {
    // Point dipoles at previous atom.
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    delx = x[i2][0] - x[i1][0];
    dely = x[i2][1] - x[i1][1];
    delz = x[i2][2] - x[i1][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    if (i1 < nlocal) {
      if (af[atom->map(mid[i1])] >= 0) {
        mu[i1][0] = delx / r;
        mu[i1][1] = dely / r;
        mu[i1][2] = delz / r;
      } else {
        amu[i1][0] = delx / r;
        amu[i1][1] = dely / r;
        amu[i1][2] = delz / r;
      }
    }
    if (i2 < nlocal) {
      if (af[atom->map(mid[i2])] < 0) {
        mu[i2][0] = delx / r;
        mu[i2][1] = dely / r;
        mu[i2][2] = delz / r;
      } else {
        amu[i2][0] = delx / r;
        amu[i2][1] = dely / r;
        amu[i2][2] = delz / r;
      }
    }
  }

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      caf = af[atom->map(mid[i])];
      if (caf == 0) {
        af[atom->map(mid[i])] = rng->gaussian(1.0,0.1);
	caf = af[atom->map(mid[i])];
      }
      mux = mu[i][0];
      muy = mu[i][1];
      muz = mu[i][2];
      mumagsq = mux*mux + muy*muy + muz*muz;
      mumag = sqrt(mumagsq);
      if (mumag > 0) {
        mu[i][0] = caf*mux/mumag;
        mu[i][1] = caf*muy/mumag;
        mu[i][2] = caf*muz/mumag;
      } else {
        amux = amu[i][0] + noise_r*nmu[i][0];
        amuy = amu[i][1] + noise_r*nmu[i][1];
        amuz = amu[i][2] + noise_r*nmu[i][2];
        amumagsq = amux*amux + amuy*amuy + amuz*amuz;
        amumag = sqrt(amumagsq);
        if (amumag > 0) {
          mu[i][0] = caf*amux/amumag;
          mu[i][1] = caf*amuy/amumag;
          mu[i][2] = caf*amuz/amumag;
        } else {
          mu[i][0] = caf;
	  mu[i][1] = 0;
	  mu[i][2] = 0;
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPFRexp::post_force(int /*vflag*/)
{
  double **f = atom->f;
  double **x = atom->x;
  double nlocal = atom->nlocal;
  double dz;
  int i;
  for (i = 0; i < nlocal; i++) {
    dz = x[i][2] - zc;
    if (x[i][2] < zc)
      f[i][2] += -k_agar * dz;
  }
}
