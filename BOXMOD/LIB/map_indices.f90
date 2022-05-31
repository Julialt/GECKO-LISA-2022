      SUBROUTINE map_indices

      USE module_data_gecko_main

      IMPLICIT NONE

!----------------------------------------------------------------------
! create index map from species to reactions
!----------------------------------------------------------------------

      WRITE(6,*) 'starting index map'
      lpmap%nloss = 0
      lpmap%nprod = 0
      DO ire=1,numre
        DO i=1,numstoi(ire,1)
          ispec = idrestoi(ire,i)
          lpmap(ispec)%nloss = lpmap(ispec)%nloss + 1
        ENDDO
        DO i=1,numstoi(ire,2)
          ispec = idpdstoi(ire,i)
          lpmap(ispec)%nprod = lpmap(ispec)%nprod + 1
        ENDDO
      ENDDO

      WRITE(6,*) ' ... index map allocating'
      DO ispec = 1, numsp
       ALLOCATE(lpmap(ispec)%idl(lpmap(ispec)%nloss), &
                lpmap(ispec)%idp(lpmap(ispec)%nprod), &
                lpmap(ispec)%stpd(lpmap(ispec)%nprod))
      ENDDO

      lpmap%nloss = 0
      lpmap%nprod = 0

      WRITE(6,*) ' ... index mapping'
      DO ire=1,numre
        DO i=1,numstoi(ire,1)
          ispec = idrestoi(ire,i)
          nloss = lpmap(ispec)%nloss + 1
          lpmap(ispec)%nloss = nloss
          lpmap(ispec)%idl(nloss) = ire
        ENDDO
        DO i=1,numstoi(ire,2)
          ispec = idpdstoi(ire,i)
          nprod = lpmap(ispec)%nprod + 1
          lpmap(ispec)%nprod = nprod
          lpmap(ispec)%idp(nprod) = ire
          lpmap(ispec)%stpd(nprod) = pdstoicf(ire,i)
        ENDDO
      ENDDO
      WRITE(6,*) 'end of index map'

!==================================
      END SUBROUTINE map_indices
!==================================
