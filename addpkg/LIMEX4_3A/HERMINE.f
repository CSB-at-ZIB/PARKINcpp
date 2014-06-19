      subroutine HERMINE ( n, k, Dense, t1, t2, tEval, yEval )
c
      implicit none
c
c-----------------------------------------------------------------------
c
cdi      integer            i, j, k, Max_Row_Tab, n
      integer            i, k, Max_Row_Tab, n
c
      double precision   t1, t2, tEval
c
      double precision   yEval ( * )
c
c-----------------------------------------------------------------------
c
c     Define size of problem via the include file
c
c        'LIMEX4_3_Size_Definitions.h'.
c
c     See the installation notes for a detailed description.
c
c-----------------------------------------------------------------------
c
      integer            Max_Nr_of_Equations, Max_Non_Zeros_Jacobian,
     2                   Max_Non_Zeros_B,
     3                   Max_Lower_Diagonals, Max_Upper_Diagonals,
     4                   Max_It_Vectors
c
c-----------------------------------------------------------------------
c
      include 'LIMEX4_3_Size_Definitions.h'
c

      parameter ( Max_Row_Tab =   7 )
c
c-----------------------------------------------------------------------
c
      integer            ipt   ( Max_Row_Tab ), nj ( Max_Row_Tab )
c
      double precision   Dense ( Max_Nr_of_Equations, 30 )
c
      double precision   Work ( Max_Nr_of_Equations , Max_Row_Tab + 1 )
c
c-----------------------------------------------------------------------
c
c     Define stepsize sequence and ipt pointers
c
c-----------------------------------------------------------------------
c
      do i = 1, Max_Row_Tab
         nj(i) = i
      end do
c
      ipt(1) = 3
c
      do i = 2, Max_Row_Tab
         ipt(i) = ipt(i-1) + nj(i)
      end do
c
c-----------------------------------------------------------------------
c
      if ( tEval .le. t2 .and. tEval .gt. t1 ) then
c
         call Comp_Herm ( n, Max_Nr_of_Equations, Dense, k, ipt, nj,
     2                    Work )
c
         call Eval_Herm ( n, Max_Nr_of_Equations, Dense, k,
     2                    ( tEval - t1 ) / ( t2 - t1 ), yEval )
c
      else
c
         write ( *, '(a,f12.5,a,2f12.5)' )
     2      ' Evalution point ', tEval,
     3      ' in HERMINE outside of current interval ', t1, t2
c
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
