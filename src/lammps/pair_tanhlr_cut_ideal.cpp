/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Arben Jusufi, Axel Kohlmeyer (Temple U.)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_tanhlr_cut_ideal.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairTanhlrCutIdeal::PairTanhlrCutIdeal(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairTanhlrCutIdeal::~PairTanhlrCutIdeal()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(htanh);
    memory->destroy(sigmah);
    memory->destroy(rmh);
    memory->destroy(rmh4);
    memory->destroy(ptanh);
    memory->destroy(offset);
    memory->destroy(ideal_potential);
  }
}

/* ---------------------------------------------------------------------- */

void PairTanhlrCutIdeal::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,imol,jmol;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,rexp,utanh,factor_lj, tanr;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    // 
    imol = atom->molecule[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jmol = atom->molecule[j];

      // i and j are index that start with zeor, so need to use floor
      // i and j are only local index; use tag to map them back to global
      int iglobal, jglobal, jmi;

      iglobal = atom->tag[i];
      jglobal = atom->tag[j];
      jmi = abs(jglobal-iglobal);

      // not the same molecule
      if ( (imol != jmol) ) continue;
      // not in the convergent orientation
      if ( (type[i] == 1) && (iglobal > jglobal) ) continue;
      if ( (type[j] == 1) && (jglobal > iglobal) ) continue;

      if ( (type[i] == 2) && (iglobal < jglobal) ) continue;
      if ( (type[j] == 2) && (jglobal < iglobal) ) continue;

      int ctcfL_l, ctcfL_r, ctcfR_l, ctcfR_r, dL, dR;
      if (jglobal > iglobal) {
          ctcfL_l = ctcfList[iglobal-1][0];
          ctcfL_r = ctcfList[iglobal-1][1];
          ctcfR_l = ctcfList[jglobal-1][0];
          ctcfR_r = ctcfList[jglobal-1][1];
      } else {
          ctcfL_l = ctcfList[jglobal-1][0];
          ctcfL_r = ctcfList[jglobal-1][1];
          ctcfR_l = ctcfList[iglobal-1][0];
          ctcfR_r = ctcfList[iglobal-1][1];
      }
      if ( (ctcfL_l==-1) || (ctcfR_r==-1) ) {
          continue;
      } else if ( (ctcfR_r - ctcfL_r > maxCTCFSep) || (ctcfR_l - ctcfL_l > maxCTCFSep) ) {
          // finite extrusion effect
          continue;
      }

      if ( (rsq < cutsq[itype][jtype]) && (jmi<nideal) ){
        r = sqrt(rsq);

        double ideal_alpha=0.0;
        if ( (type[i]==3) && (type[j]==3) ) {
            ideal_alpha = ideal_potential[jmi][0];
        } else if ( (type[i]!=3) && (type[j]!=3) ) {
            ideal_alpha = ideal_potential[jmi][2];
        } else {
            ideal_alpha = ideal_potential[jmi][1];
        }

        if (r <= rmh[itype][jtype]) { 
            rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
            tanr = tanh(rexp);
            utanh = ideal_alpha*(1.0+ tanr);
            // p is 0.5 * h, the extra negative sign is taken care in f
            fpair = factor_lj/r*ideal_alpha*sigmah[itype][jtype]*(1-tanr*tanr);
        } else {
            utanh = ideal_alpha*rmh4[itype][jtype] / rsq / rsq;
            // there is an extra r here which is to normalize position vector
            fpair = factor_lj/r*ideal_alpha*rmh4[itype][jtype]* (4.0) /rsq/rsq/r;
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = utanh - offset[itype][jtype] * ideal_alpha;
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(htanh,n+1,n+1,"pair:htanh");
  memory->create(sigmah,n+1,n+1,"pair:sigmah");
  memory->create(rmh,n+1,n+1,"pair:rmh");
  memory->create(rmh4,n+1,n+1,"pair:rmh4");
  memory->create(ptanh,n+1,n+1,"pair:ptanh");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::settings(int narg, char **arg)
{

// format for the ideal potential
// start with 0, 
// 0    0.1
// 1    0.1
// 2    0.1 

// the first argument is the cutoff
// the second argument is the input file for ideal chromosome potential
// the third argument is the exponent coefficient of the loop forming probability

  if (narg != 4) error->all(FLERR,"Illegal pair_style ideal command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
 
  // pass through the potential file 
  FILE *fp = fopen(arg[1],"r");
  char line[1024];
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ideal chromosome potential file");

  // get the number of entries in the ideal chromosome
  nideal = 0;
  while(fgets(line,1024,fp)) ++nideal;
  rewind(fp);

  memory->create(ideal_potential,nideal,3,"pair:ideal");

  char *ptr;
  int idx;

  // Loop through the rest of the lines
  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    idx = atoi(ptr);
    ptr = strtok(NULL," \t\n\r\f");

    // The second site
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted ideal chromosome potential file");
    ideal_potential[idx][0] = force->numeric(FLERR, ptr) * 0.5; // note that 0.5 is added here. This is because the contact function
                                                                // has 0.5. It is been added here to the ideal potential so that
                                                                // there is no need to do this for force calculation. If look 
                                                                // at function init_one, this is also why ptanh is half of htanh

    ptr = strtok(NULL," \t\n\r\f");
    ideal_potential[idx][1] = force->numeric(FLERR, ptr) * 0.5;

    ptr = strtok(NULL," \t\n\r\f");
    ideal_potential[idx][2] = force->numeric(FLERR, ptr) * 0.5;
  }
  fclose(fp); 

  //
  // pass through the ctcf id file
  fp = fopen(arg[2],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open the ctcf id file");

  int ntotal;
  ntotal = atom->natoms;

  memory->create(ctcfList,ntotal,2, "pair:ctcfList");

  // Loop through the rest of the lines
  idx = 0;
  while(fgets(line,102400,fp)) {
    ptr = strtok(line," \t\n\r\f");
    ctcfList[idx][0] = atoi(ptr); 
    ptr = strtok(NULL," \t\n\r\f");
    ctcfList[idx][1] = atoi(ptr); 
    idx++;
  }
  fclose(fp); 

  // maximum CTCF pair separation
  maxCTCFSep = atoi(arg[3]);
  
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR, arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR, arg[1],atom->ntypes,jlo,jhi);

  double htanh_one = force->numeric(FLERR,arg[2]);
  double rmh_one = force->numeric(FLERR,arg[3]);
  double sigmah_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      htanh[i][j] = htanh_one;
      sigmah[i][j] = sigmah_one;
      rmh[i][j] = rmh_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTanhlrCutIdeal::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
  ptanh[i][j] = htanh[i][j] * 0.5;

  if (offset_flag) {
    double rexp = (rmh[i][j]-cut[i][j])*sigmah[i][j];
    //offset[i][j] = (1.0 + tanh(rexp)); // note that in the offset here, no prefactor is calculated.
    offset[i][j] = pow(rmh[i][j],4.0)/pow(cut[i][j],4.0); // note that in the offset here, no prefactor is calculated.
                                       // This is because the prefactor depends on j-i, and is added in
                                       // a case by case fashion.
  } else offset[i][j] = 0.0;

  htanh[j][i] = htanh[i][j];
  sigmah[j][i] = sigmah[i][j];
  rmh[j][i] = rmh[i][j];
  rmh4[i][j] = pow(rmh[i][j],4.0);
  rmh4[j][i] = rmh4[i][j];
  ptanh[j][i] = ptanh[i][j];
  offset[j][i] = offset[i][j];
  cut[j][i] = cut[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&htanh[i][j],sizeof(double),1,fp);
        fwrite(&rmh[i][j],sizeof(double),1,fp);
        fwrite(&sigmah[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&htanh[i][j],sizeof(double),1,fp);
          fread(&rmh[i][j],sizeof(double),1,fp);
          fread(&sigmah[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&htanh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rmh[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigmah[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTanhlrCutIdeal::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairTanhlrCutIdeal::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r, rexp,utanh,phitanh, tanr, ideal_alpha, scale_factor;

  r=sqrt(rsq);

  int imol, jmol;
  imol = atom->molecule[i];
  jmol = atom->molecule[j];
  int iglobal, jglobal, jmi;

  iglobal = atom->tag[i];
  jglobal = atom->tag[j];
  jmi = abs(jglobal-iglobal);

  if ( (itype == 1) && (iglobal > jglobal) ) return 0.0 ;
  if ( (jtype == 1) && (jglobal > iglobal) ) return 0.0 ;

  if ( (itype == 2) && (iglobal < jglobal) ) return 0.0 ;
  if ( (jtype == 2) && (jglobal < iglobal) ) return 0.0 ;

  scale_factor = 1.0;
  int ctcfL, ctcfR, dL, dR;
  if (jglobal > iglobal) {
      ctcfL = ctcfList[iglobal-1][0];
      dL = iglobal - ctcfL;
      ctcfR = ctcfList[jglobal-1][1];
      dR = ctcfR - jglobal;
  } else {
      ctcfL = ctcfList[jglobal-1][0];
      dL = jglobal - ctcfL;
      ctcfR = ctcfList[iglobal-1][1];
      dR = ctcfR - iglobal;
  }
  if ( (ctcfL==-1) || (ctcfR==-1) ) {
      scale_factor = 0.0;
  } else {
      if ( dL>20 && dR>20 ) {
          scale_factor = 0.0;
      } else {
          scale_factor = scale_matrix[ctcfR-ctcfL] ;

      }
      //printf("ever here i? %d, %d, %d, %d, %d, %12.6f %12.6f\n", iglobal, jglobal, ctcfL, ctcfR, ctcfR-ctcfL, scale_factor, qratio);
  }

  ideal_alpha = 0.0;
  if ( (imol==jmol) && (jmi<nideal) ) {
      if ( (itype==3) && (jtype==3) ) {
        ideal_alpha = ideal_potential[jmi][0]* scale_factor;
      //} else if ( ( (jglobal>iglobal) && (itype==1) && (jtype==2) ) || ( (jglobal<iglobal) && (itype==2) && (jtype==1) ) ) {
      } else if ( ( (itype==1) && (jtype==2) ) || ( (itype==2) && (jtype==1) ) || ( (itype==1) && (jtype==4) ) || ( (itype==4) && (jtype==1) ) || ( (itype==2) && (jtype==4) ) || ( (itype==4) && (jtype==2) ) ) {
            ideal_alpha = ideal_potential[jmi][2] * scale_factor;
      } else {
        ideal_alpha = ideal_potential[jmi][1] * scale_factor;
      }
  } // else ideal_alpha = 0.0;
  //ideal_alpha = ideal_alpha * scale_factor;

  //printf("ever here i? %d %d %d %12.6f \n", i, j, jmi, ideal_alpha);

  if (r <= rmh[itype][jtype]) {
      rexp = (rmh[itype][jtype]-r)*sigmah[itype][jtype];
      tanr = tanh(rexp);
      utanh = ideal_alpha*(1.0 + tanr);
      fforce = factor_lj/r*ideal_alpha*(1.0-tanr*tanr)*sigmah[itype][jtype];
  } else {
      utanh = ideal_alpha*pow(rmh[itype][jtype],4.0) / pow(r,4.0);
      // there is an extra r here which is to normalize position vector
      //  negative sign is removed 
      fforce = factor_lj/r*ideal_alpha*pow(rmh[itype][jtype],4.0) * (4.0) /pow(r,5.0);
  }

  phitanh = utanh - offset[itype][jtype] * ideal_alpha;
  return factor_lj*phitanh;
}

/* ---------------------------------------------------------------------- */
double PairTanhlrCutIdeal::memory_usage()
{
  const int n=atom->ntypes;

  double bytes = Pair::memory_usage();

  bytes += 7*((n+1)*(n+1) * sizeof(double) + (n+1)*sizeof(double *));
  bytes += 1*((n+1)*(n+1) * sizeof(int) + (n+1)*sizeof(int *));

  return bytes;
}
