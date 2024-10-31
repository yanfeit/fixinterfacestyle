#ifdef FIX_CLASS

FixStyle(interface/pieranski/sphere, FixInterfacePieranskiSphere)

#else

#ifndef LMP_FIX_INTERFACE_PIERANSKI_SPHERE_H
#define LMP_FIX_INTERFACE_PIERANSKI_SPHERE_H
//#define PI 3.141592653589793238462643383279

#include "fix_interface.h"

namespace LAMMPS_NS{

  class FixInterfacePieranskiSphere : public FixInterface {
  public:
    FixInterfacePieranskiSphere(class LAMMPS *, int, char **);
    void precompute(int);
    void wall_particle(int, int, double);

  private:
    double coeff1[6], coeff2[6], coeff3[6], coeff4[6];
  };
  
}

#endif
#endif
