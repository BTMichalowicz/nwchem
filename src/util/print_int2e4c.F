      subroutine print_int2e4c(basis,tol2e)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer nsh,basis,max1e,max2e,mem1,mem2
      integer ish,jsh,ksh,lsh,ilo,ihi,jlo,jhi,klo,khi,llo,lhi
      integer g_ilo,g_ihi,g_jlo,g_jhi
      integer nfi,nfj,nfk,nfl,lnijkl
      integer p,q,a,b,c,i,j,k
      double precision tol2e,schwarz1,schwarz2
      double precision, allocatable :: eri1(:),eri2(:),scr(:)
      integer me, nproc
      me = ga_nodeid()
      nproc = ga_nnodes()
      call ga_sync()
      if (.not. bas_numcont(basis,nsh)) then
          call errquit('print: bas_numcont',ish, BASIS_ERR)
      endif
      call int_mem_1e(max1e, mem1)
      call int_mem_2e4c(max2e, mem2)
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi)) then
           call errquit('print: bas_cn2bfr',ish, BASIS_ERR)
        endif
        nfi=ihi-ilo+1
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi)) then
             call errquit('print: bas_cn2bfr',jsh, BASIS_ERR)
          endif
          nfj=jhi-jlo+1
          do ksh=1,nsh
            if (.not. bas_cn2bfr(basis,ksh,klo,khi)) then
               call errquit('print: bas_cn2bfr',ksh, BASIS_ERR)
            endif
            nfk=khi-klo+1
            do lsh=1,nsh
              if (.not. bas_cn2bfr(basis,lsh,llo,lhi)) then
                 call errquit('print: bas_cn2bfr',lsh, BASIS_ERR)
              endif
              nfl=lhi-llo+1
              lnijkl=nfi*nfj*nfk*nfl
              schwarz1=schwarz_shell(ish,jsh)*schwarz_shell(ksh,lsh)
              schwarz2=schwarz_shell(ish,lsh)*schwarz_shell(ksh,jsh)
              if (schwarz1.ge.tol2e)then
                  call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                          mem2, scr, max2e, eri1 )
              else
                call dfill(lnijkl,0.0d00,eri1,1)
              endif
              if (schwarz2.ge.tol2e)then
                call int_2e4c(basis, ish, lsh, basis, ksh, jsh,
     $                        mem2, scr, max2e, eri2 )
              else
                call dfill(lnijkl,0.0d00,eri2,1)
              endif
            enddo
          enddo
        enddo
      enddo
      end
