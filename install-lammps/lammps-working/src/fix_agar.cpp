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

#include "fix_agar.h"

#include "atom.h"
#include "error.h"
#include "group.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAgar::FixAgar(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix agar command");

  k_agar = utils::numeric(FLERR,arg[3],false,lmp);
  zc = utils::numeric(FLERR,arg[4],false,lmp);
  fmt::print(screen,"Agar setup\n");
}

/* ---------------------------------------------------------------------- */

int FixAgar::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAgar::init(){
  fmt::print(screen,"Agar init\n");
}

/* ---------------------------------------------------------------------- */

void FixAgar::pre_force()
{
  double dz;
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  double nlocal = atom->nlocal;

  int i;
  fmt::print(screen,"Agar pre force\n");
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dz = x[i][2] - zc;
      f[i][2] += -10;// .0*x[i][2];
    }
  }
}
