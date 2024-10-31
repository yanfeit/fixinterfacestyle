#include <math.h>
#include "fix_interface_pieranski_sphere.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixInterfacePieranskiSphere::FixInterfacePieranskiSphere(LAMMPS *lmp, int narg, char **arg) :
  FixInterface(lmp, narg, arg){}


/* 
   surface tension of air/water: gamma0
   contact angle : costheta
   radius of particle: sigma
   distance of center of particle and interface : delta
   
 */
void FixInterfacePieranskiSphere::precompute(int m)
{
  double my_pi = 3.141592653589793238462643383279;
  coeff1[m] = my_pi * gamma0[m];
  coeff2[m] = - 2 * costheta[m] * my_pi * gamma0[m] * sigma[m];
  coeff3[m] = - 2 * my_pi * gamma0[m];
  coeff4[m] = 2 * my_pi * gamma0[m] * costheta[m] * sigma[m];
}

/*
  Notice that here coord is the plane's location. We always
  need a plane to determine the sphere.
 */

void FixInterfacePieranskiSphere::wall_particle(int m, int which, double coord)
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

      // Maybe not necessary...
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      //if (delta <= 0) continue;

      double xprime, yprime, zprime;
      xprime = x[i][0] - center_x[m];
      yprime = x[i][1] - center_y[m];
      zprime = x[i][2] - center_z[m];

      double dis = sqrt(xprime*xprime + yprime*yprime + zprime*zprime);

      double z0 = radius_s[m] - dis; // reletive distance with particle and surface


      if (z0 <= -sigma[m] || z0 >= sigma[m] ) continue;

      fwall =  (coeff3[m] * z0 + coeff4[m]);
      f[i][0] -= fwall * xprime/dis;
      f[i][1] -= fwall * yprime/dis;
      f[i][2] -= fwall * zprime/dis;
      ewall[0] += coeff1[m] * z0 * z0 + coeff2[m] * z0;
      ewall[m+1] += fwall;
	
    }
 
}
