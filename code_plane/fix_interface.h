#ifndef LMP_FIX_INTERFACE_H
#define LMP_FIX_INTERFACE_H

//#define PI 3.141592653589793238462643383279

#include "fix.h"

namespace LAMMPS_NS {

class FixInterface : public Fix {
 public:
  int nwall;
  int wallwhich[6];
  double coord0[6];
  int xflag;           // 1 if any wall position is a variable
  int xstyle[6];
  int xindex[6];
  char *xstr[6];

  FixInterface(class LAMMPS *, int, char **);
  virtual ~FixInterface();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  void pre_force(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

  virtual void precompute(int) = 0;
  virtual void wall_particle(int, int, double) = 0;

 protected:
  /*
    gamma0 is the surface tension of air and vapor, depending on the 
    coarsed grained particle, temperatrue...
    sigma is the radius of the particle
    costheta is the cosine of the contact angle, depending on the 
    interaction between nanoparticle and LJ particle.
   */
  double gamma0[6],sigma[6],costheta[6];
  double ewall[7],ewall_all[7];
  double xscale,yscale,zscale;
  int estyle[6],sstyle[6],wstyle[6];
  int eindex[6],sindex[6];
  char *estr[6],*sstr[6];
  int varflag;                // 1 if any wall position,epsilon,sigma is a var
  int eflag;                  // per-wall flag for energy summation
  int ilevel_respa;
  int fldflag;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Wall defined twice in fix wall command

Self-explanatory.

E: Fix wall cutoff <= 0.0

Self-explanatory.

E: Cannot use fix wall zlo/zhi for a 2d simulation

Self-explanatory.

E: Cannot use fix wall in periodic dimension

Self-explanatory.

E: Variable name for fix wall does not exist

Self-explanatory.

E: Variable for fix wall is invalid style

Only equal-style variables can be used.

E: Variable evaluation in fix wall gave bad value

The returned value for epsilon or sigma < 0.0.

*/
