/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(spfrexp,FixSPFRexp);
// clang-format on
#else

#ifndef LMP_FIX_SPFR_EXP_H
#define LMP_FIX_SPFR_EXP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSPFRexp : public Fix {
 public:
  FixSPFRexp(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void pre_force(int);
  void post_force(int);

 protected:
  int seed;
  double af_sig, rev_mu, rev_del, noise_r, noise_t;
  double k_agar, zc;

  class RanMars *rng;
	

};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
