SUBROUTINE CALC_HFEOP ( TIME, LHFEOP_FILE, DELTA_T, EOPOUT )

      use hfeop_xyu
      implicit none
      character*(80) lhfeop_file
      real*8 EOP(4,2)
      real*8 EOPOUT(3,1)
      real*8 time
      real*8 delta_t
!      lhfeop_file="./iers2010_xyu.txt"
      call import_tides_xyu(lhfeop_file)
      call calc_hf_eop_xyu(time, delta_t, eop)
      eopout(1:3,1)=eop(1:3,1);
      RETURN
      END
