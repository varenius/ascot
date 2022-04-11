
      subroutine parametric (modn,dtq,a1,t1,a2,t2,d)
C
C Author: Zuheir Altamimi (zuheir.altamimi@ign.fr), IGN France
C
C Last updated: December 02, 2021
C
C Compute the post-seismic deformation/correction "d" using parametric models:
C #    Model
C 0    PWL (Piece-Wise Linear Function)
C 1    Logarithmic Function
C 2    Exponential Function
C 3    Logarithmic + Exponential
C 4    Two Exponential Functions
C 5    Two Logarithmic Functions
C
C IN:
C modn: model #
C dtq : time difference (t - t_Earthquake) in decimal year (but see note below)
C a1: amplitude 1 of the parametric model, if modn = 1 or 2 (or 3 or 4, if a2 & t2 are supplied)
C a2: amplitude 2 of the parametric model, if modn = 3 or 4
C t1: relaxation time 1, if modn = 1 or 2 (or 3 or 4, if a2 & t2 are supplied)
C t2: relaxtaion time 2, if modn = 3 or 4
C
C OUT:
C d: post-seismic correction
C
C Units: - mm for a1, a2, d
C        - year for t1, t2
C
C Note: Time unit is decimal year. It is advised to compute "dtq" by:
C (MJD - MJD_Earthquake)/365.25 where MJD is the modified julian day.
C

      implicit none

      doubleprecision dtq,a1,t1,a2,t2,te1,te2,d
      integer modn

      d = 0.d0

      if (modn.eq.0) return

      if (modn.eq.1) then
       d = a1*dlog( 1+ dtq/t1)
       return
      end if

      if (modn.eq.2) then
       te1 = dtq/t1
       d =  a1*(1-dexp(-te1))
       return
      end if

      if (modn.eq.3) then
       te2 = dtq/t2
       d = a1*dlog( 1+ dtq/t1) + a2*(1-dexp(-te2))
       return
      end if

      if (modn.eq.4) then
       te1 = dtq/t1
       te2 = dtq/t2
       d = a1*(1-dexp(-te1)) + a2*(1-dexp(-te2))
       return
      end if

      if (modn.eq.5) then
       d = a1*dlog( 1+ dtq/t1) + a2*dlog( 1+ dtq/t2)
       return
      end if

      end
