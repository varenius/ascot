      PROGRAM XHF2002_IERS
*+
*  - - - - - - - - - - - - -
*   X H F 2 0 0 2 _ I E R S
*  - - - - - - - - - - - - -
*
*  This program is part of the International Earth Rotation and
*  Reference Systems Service (IERS) Conventions software collection.
*
*  This program is a driver for the function HF2002_IERS.
*  This program computes TCB-TCG at the geocenter using an approximation of the
*  time ephemeris TE405 (Irwin and Fukushima) at a given TT value.
*
*  In general, Class 1, 2, and 3 models represent physical effects that
*  act on geodetic parameters while canonical models provide lower-level
*  representations or basic computations that are used by Class 1, 2, or
*  3 models.
*
*  Status: Class 1
*
*     Class 1 models are those recommended to be used a priori in the
*     reduction of raw space geodetic data in order to determine
*     geodetic parameter estimates.
*     Class 2 models are those that eliminate an observational
*     singularity and are purely conventional in nature.
*     Class 3 models are those that are not required as either Class
*     1 or 2.
*     Canonical models are accepted as is and cannot be classified as
*     a Class 1, 2, or 3 model.
*
*  Notes:
*
*  1) See the function HF2002_IERS for special notes.
*
*  Called:
*     HF2002_IERS Computes TCB-TCG at the geocenter using an approximation of the time ephemeris TE405
*
*  Test cases:
*     expected output:
*  Year      JD(TT)      TCB-TCG / s
*  1600 2305445.0   -176.1773555610856D0
*  1700 2341970.0   -129.4460764716997D0
*  1800 2378495.0    -82.7147442912425D0
*  1900 2415020.0    -35.9834220078961D0
*  2000 2451545.0     10.7478546790489D0
*  2100 2488070.0     57.4792210328145D0
*  2200 2524595.0    104.2104805049757D0
*
*  References:
*
*     Harada, W. and Fukushima, T., 2003, ``Harmonic Decomposition of
*     Time Ephemeris TE405", Astron. J., 126, pp. 2557-2561

*     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
*     IERS Technical Note No. 36, BKG (2010)
*
*  Revisions:
*  2009 August 28 B.E.Stetzler  Initial standardization of function
*                               and provided a test case
*  2009 August 28 B.E.Stetzler  Capitalized all variables for FORTRAN
*                               77 compatibility
*  2009 August 31 B.E.Stetzler  Provided a test case
*  2010 May    25 G.  Petit     Xhf2002_IERS
*  2010 June   22 B.E.Stetzler  Added linear correction to function
*                               hf2002 and updated test case
*  2010 June   23 B.E.Stetzler  Changed function call to hf2002_iers
*  2010 August 13 B.E.Stetzler  Included platform independent comments
*                               from P. Wallace and used proper code 
*                               indentation and syntax
*  2010 August 26 B.E.Stetzler  Negated use of external data file by
*                               placing all data in DATA statements
*                               and updated test case
*-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER JYEAR
      DOUBLE PRECISION T0, TJD, HF2002_IERS, W

      T0 = 2451545D0
      WRITE ( *, '(T2,A,T8,A,T20,A)' ) 'YEAR','JD(TT)','TCB-TCG / s'

      DO JYEAR=1600,2200,100
         TJD = T0+DBLE(JYEAR-2000)*365.25D0
         W   = HF2002_IERS(TJD)
         WRITE ( *, '(I5,F10.1,F21.13)' ) JYEAR, TJD, W
      END DO

*+----------------------------------------------------------------------
*
*  Copyright (C) 2008
*  IERS Conventions Center
*
*  ==================================
*  IERS Conventions Software License
*  ==================================
*
*  NOTICE TO USER:
*
*  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
*  WHICH APPLY TO ITS USE.
*
*  1. The Software is provided by the IERS Conventions Center ("the
*     Center").
*
*  2. Permission is granted to anyone to use the Software for any
*     purpose, including commercial applications, free of charge,
*     subject to the conditions and restrictions listed below.
*
*  3. You (the user) may adapt the Software and its algorithms for your
*     own purposes and you may distribute the resulting "derived work"
*     to others, provided that the derived work complies with the
*     following requirements:
*
*     a) Your work shall be clearly identified so that it cannot be
*        mistaken for IERS Conventions software and that it has been
*        neither distributed by nor endorsed by the Center.
*
*     b) Your work (including source code) must contain descriptions of
*        how the derived work is based upon and/or differs from the
*        original Software.
*
*     c) The name(s) of all modified routine(s) that you distribute
*        shall be changed.
*
*     d) The origin of the IERS Conventions components of your derived
*        work must not be misrepresented; you must not claim that you
*        wrote the original Software.
*
*     e) The source code must be included for all routine(s) that you
*        distribute.  This notice must be reproduced intact in any
*        source distribution.
*
*  4. In any published work produced by the user and which includes
*     results achieved by using the Software, you shall acknowledge
*     that the Software was used in obtaining those results.
*
*  5. The Software is provided to the user "as is" and the Center makes
*     no warranty as to its use or performance.   The Center does not
*     and cannot warrant the performance or results which the user may
*     obtain by using the Software.  The Center makes no warranties,
*     express or implied, as to non-infringement of third party rights,
*     merchantability, or fitness for any particular purpose.  In no
*     event will the Center be liable to the user for any consequential,
*     incidental, or special damages, including any lost profits or lost
*     savings, even if a Center representative has been advised of such
*     damages, or for any claim by any third party.
*
*  Correspondence concerning IERS Conventions software should be
*  addressed as follows:
*
*                     Gerard Petit
*     Internet email: gpetit[at]bipm.org
*     Postal address: IERS Conventions Center
*                     Time, frequency and gravimetry section, BIPM
*                     Pavillon de Breteuil
*                     92312 Sevres  FRANCE
*
*     or
*
*                     Brian Luzum
*     Internet email: brian.luzum[at]usno.navy.mil
*     Postal address: IERS Conventions Center
*                     Earth Orientation Department
*                     3450 Massachusetts Ave, NW
*                     Washington, DC 20392
*
*
*-----------------------------------------------------------------------
      END
