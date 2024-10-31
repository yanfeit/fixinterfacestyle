#ifdef FIX_CLASS

FixStyle(interface/pieranski, FixInterfacePieranski)

#else

#ifndef LMP_FIX_INTERFACE_PIERANSKI_H
#define LMP_FIX_INTERFACE_PIERANSKI_H


#include "fix_interface.h"

namespace LAMMPS_NS{

  class FixInterfacePieranski : public FixInterface {
  public:
    FixInterfacePieranski(class LAMMPS *, int, char **);
    void precompute(int);
    void wall_particle(int, int, double);

  private:
    double coeff1[6], coeff2[6], coeff3[6], coeff4[6];
  };
  
}

#endif
#endif
