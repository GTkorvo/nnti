/*! \file phys.h
    \brief Physics constants
*/
#ifndef PHYS_H
#define PHYS_H 1

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

/* Define some physical constants (SI) */

/* From CODATA... See Rev Mod Phys V72 #2, 2000 */
#define CMM_MAGNETIC_CONST     (12.566370614e-7)  /* Magnetic constant (vacuum magn susc)           */
#define CMM_ELECTRIC_CONST     (8.854187817e-12)  /* Electric constant (vacuum dielectric constant) */

/*! \def CMM_LIGHT
Speed of light \f$ c \f$ in vacuum. Units are m/sec.
*/
#define CMM_LIGHT         (2.99792458e8)

/*! \def CMM_HBAR
\f$ \hbar \f$, Planck's constant divided by \f$2\pi\f$. Units are J-sec.
*/
#define CMM_HBAR          (1.054571596e-34)

/*! \def CMM_FINESTRUCT
The fine structure constant, \f$ \alpha \f$. Dimensionless.
*/
#define CMM_FINESTRUCT    (1./137.03599976)

/*! \def CMM_GRAV
Newton's gravitational constant. Units are N m^2/kg^2.
*/
#define CMM_GRAV          (6.673e-11)

/*! \def CMM_CHARGE
Charge on the electron. Units are Coulombs.
*/
#define CMM_CHARGE        (1.602176462e-19)

/*! \def CMM_BOHR
Bohr radius. Units are m.
*/
#define CMM_BOHR          (0.5291772083e-10)

/*! \def CMM_ELECTRONVOLT
Electron volt. Units are J.
*/
#define CMM_ELECTRONVOLT  (CMM_CHARGE)

/*! \def CMM_KCAL
Kilocalorie. Units are J.
*/
#define CMM_KCAL          (4186.8)

/*! \def CMM_BOLTZMAN
Boltzman's constant. Units are J/K.
*/
#define CMM_BOLTZMAN      (1.3806503e-23)

/*! \def CMM_AVOGADRO
Avogadro's number (1 mole). Dimensionless.
*/
#define CMM_AVOGADRO      (6.02214199e23)

/*! \def CMM_KCALPERMOL
Kcal/mole. Units are J.
*/
#define CMM_KCALPERMOL    (CMM_KCAL/CMM_AVOGADRO)

/*! \def CMM_EMASS
Electron mass. Units are kg.
*/
#define CMM_EMASS         (9.10938188e-31)

/*! \def CMM_PMASS
Proton mass. Units are kg.
*/
#define CMM_PMASS         (1.67262158e-27)

/*! \def CMM_AMU
Atomic mass unit. Units are kg.
*/
#define CMM_AMU           (1.66053873e-27)

/*! \def CMM_HARTREE
Hartree energy. Units are kg.
// A Hartree is approx 4.3597e-18 J, or 27.21 eV
*/
#define CMM_HARTREE       (CMM_CHARGE*CMM_CHARGE/(4.*CMM_PI*CMM_ELECTRICMM_CONST*CMM_BOHR))

/*! \def CMM_RYD
Rydberg energy. Units are kg.
*/
#define CMM_RYD           (.5*CMM_HARTREE)

/*! \def CMM_G0
Gravitational acceleration at Earth surface. Units are m/sec^2.
*/
#define CMM_G0            (9.81)


/* To convert temperature in deg K to kcal/mol

   t_kcalmol = t_kelvin * (CMM_BOLTZMAN / CMM_KCAL_PERMOL);

*/

/*

   Force in Newtons between two charges in vacuum is

               q_1 q_2
   f_12 =   --------------         (r_1 - r_2)
            4*pi * eps_0 |r_1 - r_2|^3

   where charges are given in Coulombs, dist in meters, eps_0 is CMM_ELECTRIC_CONST

*/

/* Define some conversion factors */

/* misc */
/*! \def CMM_DEG_RADIAN
Convert degrees to radians.
*/
#define CMM_DEG_RADIAN    (CMM_PI/180.)     /* convert from degrees to radians */

/*! \def CMM_KCALPERMOL_EV
Convert Kcal/mole to ev.
*/
#define CMM_KCALPERMOL_EV (CMM_KCALPERMOL/CMM_ELECTRONVOLT)   /* convert from kcal/mol to ev */

/*! \def CMM_RYD_EV
Convert Rydberg to ev.
*/
#define CMM_RYD_EV        (.5*CMM_RYD/CMM_ELECTRONVOLT)       /* convert from Rydbergs to ev */

/* stupid US units */

/*! \def CMM_INCH_CM
Convert inches to cm.
*/
#define CMM_INCH_CM       (2.54)  /* convert from inches to cm */

/*! \def CMM_LB_KG
Convert lb (mass) to kg.
*/
#define CMM_LB_KG         (0.4535924)     /* convert from lb mass to kg */

/*! \def CMM_GAL_LITER
Convert gal to L.
*/
#define CMM_GAL_LITER     (3.785412)      /* convert from gallons to liters */

#ifdef __cplusplus
}
#endif                          /* __cplusplus */


#endif
