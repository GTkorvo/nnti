/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */
/* It was modified by Courtenay T. Vaughan for use in Zoltan */

#include <stdio.h>
#include <math.h>
#include "lb_const.h"
#include "rib_const.h"

/* macros for routines */
#define max(a, b) ((a) < (b) ? (b) : (a))
#define min(a, b) ((a) > (b) ? (b) : (a))
#define sign(x)   ((x) >= 0 ? 1.0 : -1.0)

static void evals2(double[2][2], double *, double *);
static void eigenvec2(double[2][2], double, double *, double *);

int LB_inertial2d(
     struct Dot_Struct *dotpt,  /* graph data structure for weights */
     int              dotnum,   /* number of vtxs in graph */
     int              wgtflag,  /* are vertex weights being used? */
     double           cm[3],    /* center of mass in each direction */
     double           evec[3],  /* eigenvector */
     double           *value,   /* array for value to sort on */
     MPI_Comm         comm      /* communicator for partition */
)
{
     double tensor[2][2];       /* inertial tensor */
     double xcm, ycm;           /* center of mass in each direction */
     double xx, yy, xy;         /* elements of inertial tensor */
     double xdif, ydif;         /* deviation from center of mass */
     double eval, res;          /* eigenvalue and error in eval calculation */
     double wgt_sum;            /* sum of all the vertex weights */
     int    i;                  /* loop counter */
     double xcmt, ycmt, wgtt;   /* temp for center of mass */
     double xxt, yyt, xyt;      /* temp for tensor */

     /* Compute center of mass and total mass. */

     xcm = ycm = 0.0;
     if (wgtflag) {
        wgt_sum = 0.0;
        for (i = 0; i < dotnum; i++) {
           wgt_sum += dotpt[i].Weight;
           xcm += dotpt[i].Weight*dotpt[i].X[0];
           ycm += dotpt[i].Weight*dotpt[i].X[1];
        }
     }
     else {
        wgt_sum = dotnum;
        for (i = 0; i < dotnum; i++) {
           xcm += dotpt[i].X[0];
           ycm += dotpt[i].X[1];
        }
     }

     /* Sum weights across processors */

     MPI_Allreduce(&xcm,&xcmt,1,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&ycm,&ycmt,1,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&wgt_sum,&wgtt,1,MPI_DOUBLE,MPI_SUM,comm);

     xcm = xcmt/wgtt;
     ycm = ycmt/wgtt;

     /* Generate 3 elements of Inertial tensor. */
     xx = yy = xy = 0.0;
     if (wgtflag)
        for (i = 0; i < dotnum; i++) {
           xdif = dotpt[i].X[0] - xcm;
           ydif = dotpt[i].X[1] - ycm;
           xx += dotpt[i].Weight*xdif*xdif;
           yy += dotpt[i].Weight*ydif*ydif;
           xy += dotpt[i].Weight*xdif*ydif;
        }
     else
        for (i = 0; i < dotnum; i++) {
           xdif = dotpt[i].X[0] - xcm;
           ydif = dotpt[i].X[1] - ycm;
           xx += xdif*xdif;
           yy += ydif*ydif;
           xy += xdif*ydif;
        }

     /* Sum tensor across processors */

     MPI_Allreduce(&xx,&xxt,1,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&yy,&yyt,1,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&xy,&xyt,1,MPI_DOUBLE,MPI_SUM,comm);

     /* Compute eigenvector with maximum eigenvalue. */

     tensor[0][0] = xxt;
     tensor[1][1] = yyt;
     tensor[1][0] = tensor[0][1] = xyt;
     evals2(tensor, &res, &eval);
     eigenvec2(tensor, eval, evec, &res);

     /* Calculate value to sort/split on for each cell. */
     /* This is inner product with eigenvector. */
     for (i = 0; i < dotnum; i++)
        value[i] = (dotpt[i].X[0] - xcm)*evec[0] +
                   (dotpt[i].X[1] - ycm)*evec[1];

     /* KDDKDD -- Do we need to set cm[0] and cm[1]? */
     cm[0] = xcm;
     cm[1] = ycm;

     /* zero unused third dimension */
     cm[2] = evec[2] = 0.0;

     return(LB_OK);
}


/* Find eigenvalues of 2x2 symmetric system by solving quadratic. */
static void evals2(
     double H[2][2],            /* symmetric matrix for eigenvalues */
     double *eval1,             /* smallest eigenvalue */
     double *eval2              /* middle eigenvalue */
)
{
     double M[2][2];            /* normalized version of matrix */
     double b, c;               /* coefficents of cubic equation */
     double root1, root2;       /* roots of quadratic */
     double xmax;               /* largest matrix element */
     int    i, j;               /* loop counters */

     xmax = 0.0;
     for (i = 0; i < 2; i++)
        for (j = i; j < 2; j++)
           if (fabs(H[i][j]) > xmax)
              xmax = fabs(H[i][j]);
     if (xmax != 0)
        for (i = 0; i < 2; i++)
           for (j = 0; j < 2; j++)
              M[i][j] = H[i][j] / xmax;

     b = -M[0][0] - M[1][1];
     c = M[0][0] * M[1][1] - M[1][0] * M[1][0];
     root1 = -.5 * (b + sign(b) * sqrt(max(0.0, b * b - 4 * c)));
     root2 = c / root1;

     root1 *= xmax;
     root2 *= xmax;
     *eval1 = min(root1, root2);
     *eval2 = max(root1, root2);
}


/* Solve for eigenvector of SPD 2x2 matrix, with given eigenvalue. */
static void eigenvec2(
     double A[2][2],            /* matrix */
     double eval,               /* eigenvalue */
     double evec[2],            /* eigenvector returned */
     double *res                /* normalized residual */
)
{
     double norm;               /* norm of eigenvector */
     double res1, res2;         /* components of residual vector */
     int    i;                  /* loop counter */

     if (fabs(A[0][0] - eval) > fabs(A[1][1] - eval)) {
        evec[0] = -A[1][0];
        evec[1] = A[0][0] - eval;
     }
     else {
        evec[0] = A[1][1] - eval;
        evec[1] = -A[1][0];
     }

     /* Normalize eigenvector and calculate a normalized eigen-residual. */
     norm = sqrt(evec[0] * evec[0] + evec[1] * evec[1]);
     if (norm == 0) {
        evec[0] = 1;
        evec[1] = 0;
        norm = 1;
     }
     for (i = 0; i < 2; i++)
        evec[i] /= norm;
     res1 = (A[0][0] - eval) * evec[0] + A[1][0] * evec[1];
     res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1];
     *res = sqrt(res1 * res1 + res2 * res2);

     res1 = fabs(A[0][0]) + fabs(A[1][0]);
     res2 = fabs(A[1][1]) + fabs(A[1][0]);
     *res /= max(res1, res2);
}
