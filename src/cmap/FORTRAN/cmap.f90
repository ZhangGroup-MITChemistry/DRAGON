Program main
!
! part of the code for reading DCD file is copied from
! http://www.ks.uiuc.edu/Research/namd/wiki/index.cgi?ReadingDCDinFortran
!
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    ! variable definitions
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !

    implicit none
    character*500      :: dcdFileName, instr, stype

    ! variables to read in dcd file
    double precision  :: d
    real              :: t,dummyr
    real,allocatable  :: x(:,:),y(:,:),z(:,:)
    integer           :: nset,natom,dummyi,i,j,nframes
    character*4       :: dummyc

    ! parameters for contact map
    double precision, parameter :: r_cut=1.76, sigma=3.72
    integer, parameter:: nCS=15                         ! number of chromatin types
    integer           :: nCG, cgfactor, atomSta, atomStap1,&
                         atomEnd, first_frame

    ! optimization variables
    double precision :: dx, dy, dz, rij, r2, pij, cgmean, r_cut4
    double precision, allocatable :: cmap(:,:),       &   ! per bead
                                     cmapCG(:,:),     &   ! per chromatin state
                                     cvInd(:),        &   ! cmap from HiC
                                     cvCorr(:,:),     &
                                     counter(:)

    ! general variables
    integer           :: ib, jb, ifr, kc, lc
    integer           :: reclength, natomh
    
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    ! function starts here
    ! =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~ !
    
    ! getarg retrieves command line arguments, and sets first to "filename"
    call getarg(1,dcdFileName)
    call getarg(2,instr)
    read (instr, '(i10)') cgfactor
    call getarg(3,instr)
    read (instr, '(i10)') atomSta
    call getarg(4,instr)
    read (instr, '(i10)') atomEnd
    call getarg(5,instr)
    read (instr, '(i10)') first_frame 

    !
    ! read dcd file
    open(10,file=trim(dcdFileName),status='old',form='unformatted')
    read(10) dummyc, nframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) natom

    allocate(x(natom,nframes))
    allocate(y(natom,nframes))
    allocate(z(natom,nframes))

    do i = 1, nframes
       read(10) (d, j=1, 6)
       read(10) (x(j,i),j=1,natom)
       read(10) (y(j,i),j=1,natom)
       read(10) (z(j,i),j=1,natom)
    end do

    print*
    print *, 'The dcd contains ',nframes,' frames and ',natom,' atoms'
    print*
    close(10)

    ! 
    ! allocate variable
    natomh = atomEnd-atomSta

    allocate(cmap(natomh,natomh))
    allocate(cvInd(nAtomh))
    allocate(counter(nAtomh))

    atomStap1 = atomSta + 1
    ! 
    ! calculate the contact probability
    cmap = 0.0
    r_cut4 = r_cut**4
    do ifr = first_frame, nframes

        do ib = atomStap1, atomEnd
            do jb = ib+1, atomEnd

                dx = x(ib,ifr)-x(jb,ifr) 
                dy = y(ib,ifr)-y(jb,ifr) 
                dz = z(ib,ifr)-z(jb,ifr) 
                r2 = dx**2 + dy**2 + dz**2
                rij = sqrt(r2)

                if (rij <= r_cut) then
                    pij = 0.5*( 1.0+tanh(sigma*(r_cut-rij)) )
                else 
                    pij = 0.5*r_cut4 / r2 / r2
                endif 
                
                cmap(ib-atomSta,jb-atomSta) = cmap(ib-atomSta,jb-atomSta) + pij  

            enddo
        enddo

    enddo

    cmap = cmap / (nframes - first_frame + 1)

    ! 
    ! initialize the contact map
    cvInd = 0.0
    counter = 0

    do ib = atomStap1, atomEnd
        do jb = ib+1, atomEnd

            pij = cmap(ib-atomSta, jb-atomSta)

            cvind(jb-ib) = cvind(jb-ib) + pij
            counter(jb-ib)    = counter(jb-ib) + 1

        enddo
    enddo

    cvind(1:(nAtomh-1)) =  cvind(1:(nAtomh-1)) / counter(1:(nAtomh-1))

    !
    ! save the contact matrix
    open(unit=10, file='cp_scaling.txt')
    do ib = 1, natomh-1
        write(10, '(E14.7)') cvInd(ib)
    enddo
    close(10)


    !
    ! coarse grain the matrix
    NCG      = natomh / cgfactor
    allocate(cmapCG(NCG, NCG))

    cmapCG = 0.0
    do kc = 1, NCG
        do lc = kc, NCG
            do ib = (kc-1)*cgfactor+1, kc*cgfactor
                do jb = (lc-1)*cgfactor+1, lc*cgfactor
                    cmapCG(kc, lc) = cmapCG(kc, lc) + cmap(ib,jb)
                enddo
            enddo

            cmapCG(kc,lc) = cmapCG(kc,lc) / cgfactor / cgfactor
            cmapCG(lc,kc) = cmapCG(kc,lc)
        enddo
    enddo

    cgmean = 0.0
    do kc = 1, nCG-1
       cgmean = cgmean + cmapCG(kc,kc+1)
    enddo
    cgmean = cgmean / (nCG-1)
    cmapCG = cmapCG / cgmean

    !
    ! save the contact matrix
    open(unit=10, file='contact_map_CG.txt')
    do ib = 1, nCG
        do jb = 1, nCG-1
            write(10, '(f12.8)', advance='no') cmapCG(ib, jb)
        enddo
        write(10, '(f12.8)') cmapCG(ib,nCG)
    enddo
    close(10)

    open(11, file='nframes.txt')
    write(11, *) nframes - first_frame + 1
    close(11)

    !
    ! release the memory used to store coordinates. 
    deallocate(x)
    deallocate(y)
    deallocate(z)

    deallocate(cmap)
    deallocate(cmapCG)

end program main
