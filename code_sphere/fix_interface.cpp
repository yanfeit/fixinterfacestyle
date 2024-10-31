#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_interface.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixInterface::FixInterface(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  nwall(0)
{
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // parse args

  int scaleflag = 1;
  fldflag = 0;
  int pbcflag = 0;

  for (int i = 0; i < 6; i++) xstr[i] = estr[i] = sstr[i] = NULL;

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+9 > narg) error->all(FLERR,"Illegal fix wall command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        xstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        xstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        xstr[nwall] = new char[n];
        strcpy(xstr[nwall],&arg[iarg+1][2]);
      } else {
        xstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }

      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        cxstr[nwall] = new char[n];
        strcpy(cxstr[nwall],&arg[iarg+2][2]);
        cxstyle[nwall] = VARIABLE;
      } else {
        center_x[nwall] = force->numeric(FLERR,arg[iarg+2]);
        cxstyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        cystr[nwall] = new char[n];
        strcpy(cystr[nwall],&arg[iarg+3][2]);
        cystyle[nwall] = VARIABLE;
      } else {
        center_y[nwall] = force->numeric(FLERR,arg[iarg+3]);
        cystyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        czstr[nwall] = new char[n];
        strcpy(czstr[nwall],&arg[iarg+4][2]);
        czstyle[nwall] = VARIABLE;
      } else {
        center_z[nwall] = force->numeric(FLERR,arg[iarg+4]);
        czstyle[nwall] = CONSTANT;
      }
      
      if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
        int n = strlen(&arg[iarg+5][2]) + 1;
        rsstr[nwall] = new char[n];
        strcpy(rsstr[nwall],&arg[iarg+5][2]);
        rsstyle[nwall] = VARIABLE;
      } else {
        radius_s[nwall] = force->numeric(FLERR,arg[iarg+5]);
        rsstyle[nwall] = CONSTANT;
      }
      
      if (strstr(arg[iarg+6],"v_") == arg[iarg+6]) {
        int n = strlen(&arg[iarg+6][2]) + 1;
        estr[nwall] = new char[n];
        strcpy(estr[nwall],&arg[iarg+6][2]);
        estyle[nwall] = VARIABLE;
      } else {
        gamma0[nwall] = force->numeric(FLERR,arg[iarg+6]);
        estyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+7],"v_") == arg[iarg+7]) {
        int n = strlen(&arg[iarg+7][2]) + 1;
        sstr[nwall] = new char[n];
        strcpy(sstr[nwall],&arg[iarg+7][2]);
        sstyle[nwall] = VARIABLE;
      } else {
        sigma[nwall] = force->numeric(FLERR,arg[iarg+7]);
        sstyle[nwall] = CONSTANT;
      }

      if (strstr(arg[iarg+8],"v_") == arg[iarg+8]) {
        int n = strlen(&arg[iarg+8][2]) + 1;
        cstr[nwall] = new char[n];
        strcpy(cstr[nwall],&arg[iarg+8][2]);
        cstyle[nwall] = VARIABLE;
      } else {
        costheta[nwall] = force->numeric(FLERR,arg[iarg+8]);
        cstyle[nwall] = CONSTANT;
      }

      nwall++;
      iarg += 9;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"fld") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"no") == 0) fldflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) fldflag = 1;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pbc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall command");
      if (strcmp(arg[iarg+1],"yes") == 0) pbcflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pbcflag = 0;
      else error->all(FLERR,"Illegal fix wall command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall command");
  }

  size_vector = nwall;

  // error checks

  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");
  for (int m = 0; m < nwall; m++)
    if ( (costheta[m] > 1.0) || (costheta[m] < -1.0) )
      error->all(FLERR,"Invalid Contact angle!, Set it to (-1, 1)");

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,"Cannot use fix wall zlo/zhi for a 2d simulation");

  if (!pbcflag) {
    for (int m = 0; m < nwall; m++) {
      if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
      if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
      if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
        error->all(FLERR,"Cannot use fix wall in periodic dimension");
    }
  }
  
  // scale factors for wall position for CONSTANT and VARIABLE walls

  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (xstyle[m] != EDGE) flag = 1;

  if (flag) {
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    for (int m = 0; m < nwall; m++) {
      if (xstyle[m] != CONSTANT) continue;
      if (wallwhich[m] < YLO) coord0[m] *= xscale;
      else if (wallwhich[m] < ZLO) coord0[m] *= yscale;
      else coord0[m] *= zscale;
    }
  }

  // set xflag if any wall positions are variable
  // set varflag if any wall positions or parameters are variable
  // set wstyle to VARIABLE if any other parameters (except wall position) is a variable

  varflag = xflag = 0;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) xflag = 1;
    if (xflag || estyle[m] == VARIABLE || sstyle[m] == VARIABLE || cstyle[m] == VARIABLE || cxstyle[m] == VARIABLE || cystyle[m] == VARIABLE || czstyle[m] == VARIABLE || rsstyle[m] == VARIABLE) varflag = 1;
    if (estyle[m] == VARIABLE || sstyle[m] == VARIABLE || cstyle[m] == VARIABLE || cxstyle[m] == VARIABLE || cystyle[m] == VARIABLE || czstyle[m] == VARIABLE || rsstyle[m] == VARIABLE) wstyle[m] = VARIABLE;
    else wstyle[m] = CONSTANT;
  }

  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixInterface::~FixInterface()
{
  for (int m = 0; m < nwall; m++) {
    delete [] xstr[m];
    delete [] estr[m];
    delete [] sstr[m];
  }
}

/* ---------------------------------------------------------------------- */

int FixInterface::setmask()
{
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInterface::init()
{
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      xindex[m] = input->variable->find(xstr[m]);
      if (xindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(xindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (cxstyle[m] == VARIABLE) {
      cxindex[m] = input->variable->find(cxstr[m]);
      if (cxindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(cxindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (cystyle[m] == VARIABLE) {
      cyindex[m] = input->variable->find(cystr[m]);
      if (cyindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(cyindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (czstyle[m] == VARIABLE) {
      czindex[m] = input->variable->find(czstr[m]);
      if (czindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(czindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (rsstyle[m] == VARIABLE) {
      rsindex[m] = input->variable->find(rsstr[m]);
      if (rsindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(rsindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (estyle[m] == VARIABLE) {
      eindex[m] = input->variable->find(estr[m]);
      if (eindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(eindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
    if (sstyle[m] == VARIABLE) {
      sindex[m] = input->variable->find(sstr[m]);
      if (sindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(sindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
     if (cstyle[m] == VARIABLE) {
      cindex[m] = input->variable->find(cstr[m]);
      if (cindex[m] < 0)
        error->all(FLERR,"Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(cindex[m]))
        error->all(FLERR,"Variable for fix wall is invalid style");
    }
  }

  // setup coefficients

  for (int m = 0; m < nwall; m++) precompute(m);

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixInterface::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    if (!fldflag) post_force(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixInterface::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   only called if fldflag set, in place of post_force
------------------------------------------------------------------------- */

void FixInterface::pre_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixInterface::post_force(int vflag)
{
  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;
 
  // coord = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for epsilon/sigma variables need to re-invoke precompute()

  if (varflag) modify->clearstep_compute();

  double coord;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(xindex[m]);
      if (wallwhich[m] < YLO) coord *= xscale;
      else if (wallwhich[m] < ZLO) coord *= yscale;
      else coord *= zscale;
    } else coord = coord0[m];
    if (wstyle[m] == VARIABLE) {
      if (cxstyle[m] == VARIABLE) {
        center_x[m] = input->variable->compute_equal(cxindex[m]);
      }
      if (cystyle[m] == VARIABLE) {
        center_y[m] = input->variable->compute_equal(cyindex[m]);
      }
      if (czstyle[m] == VARIABLE) {
        center_z[m] = input->variable->compute_equal(czindex[m]);
      }
      if (rsstyle[m] == VARIABLE) {
        radius_s[m] = input->variable->compute_equal(rsindex[m]);
      }
      if (estyle[m] == VARIABLE) {
        gamma0[m] = input->variable->compute_equal(eindex[m]);
        if (gamma0[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall gave bad value");
      }
      if (sstyle[m] == VARIABLE) {
        sigma[m] = input->variable->compute_equal(sindex[m]);
        if (sigma[m] < 0.0)
          error->all(FLERR,"Variable evaluation in fix wall gave bad value");
      }
      if (cstyle[m] == VARIABLE) {
        costheta[m] = input->variable->compute_equal(cindex[m]);
      }
      precompute(m);
    }

    wall_particle(m,wallwhich[m],coord);
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixInterface::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixInterface::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixInterface::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixInterface::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall,ewall_all,nwall+1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return ewall_all[n+1];
}

