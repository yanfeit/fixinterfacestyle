#include <math.h>
#include "fix_interface_pieranski.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixInterfacePieranski::FixInterfacePieranski(LAMMPS *lmp, int narg, char **arg) :
  FixInterface(lmp, narg, arg){}


/* 
   surface tension of air/water: gamma0
   contact angle : costheta
   radius of particle: sigma
   distance of center of particle and interface : delta
   
 */
void FixInterfacePieranski::precompute(int m)
{
  double my_pi = 3.141592653589793238462643383279;
  coeff1[m] = my_pi * gamma0[m];
  coeff2[m] = - 2 * costheta[m] * my_pi * gamma0[m] * sigma[m];
  coeff3[m] = - 2 * my_pi * gamma0[m];
  coeff4[m] = 2 * my_pi * gamma0[m] * costheta[m] * sigma[m];
}

void FixInterfacePieranski::wall_particle(int m, int which, double coord)
{
  double delta, fwall;

  double **x = atom->x;
  double **f = atom->f;
  int * mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which/2;
  int side = which%2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i ++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= sigma[m] || delta <= -sigma[m]) continue;

      fwall = side * (coeff3[m] * delta + coeff4[m]);
      f[i][dim] -= fwall;
      ewall[0] += coeff1[m] * delta * delta + coeff2[m] * delta;
      ewall[m+1] += fwall;
	
    }
 
}
