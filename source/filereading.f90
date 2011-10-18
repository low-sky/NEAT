module mod_abundIO
use mod_abundtypes
implicit none!

contains

subroutine read_ilines(ILs, Iint)        
        TYPE(line), DIMENSION(62) :: ILs
        INTEGER :: Iint 

        Iint = 1

        301 FORMAT(A11, 1X, A6, 1X, F7.2, 1X, A20,1X,A4)
        OPEN(201, file="source/Ilines_levs", status='old')
                DO WHILE (Iint < 62)!(.true.)
                        READ(201,301,end=401) ILs(Iint)%name, ILs(Iint)%ion, ILs(Iint)%wavelength, ILs(Iint)%transition ,ILs(Iint)%zone!end condition breaks loop.  
                        Iint = Iint + 1
                END DO
                401 PRINT*, "done reading important lines, Iint = ", Iint 
        CLOSE(201)
end subroutine        

end module

module mod_abundmaths
use mod_abundtypes
implicit none!

contains

!this fantastically ugly function gets the location of certain ions in the important ions array using their name as a key.

integer function get_ion(ionname, iontable, Iint)
        CHARACTER*11 :: ionname
        TYPE(line), DIMENSION(62) :: iontable 
        INTEGER :: i
        INTEGER, INTENT(IN) :: Iint

        do i = 1, Iint

                !PRINT*, trim(iontable(i)%name), trim(ionname)

                if(trim(iontable(i)%name) == trim(ionname))then
                        get_ion = i
                        return
                endif
        end do

        PRINT*, "Nudge Nudge, wink, wink error. Ion not found, say no more.", ionname

end function        


subroutine element_assign(ILs, linelist, Iint, listlength)
        TYPE(line), DIMENSION(62), INTENT(OUT) :: ILs
        TYPE(line), DIMENSION(:) :: linelist 
        INTEGER, INTENT(IN) :: Iint, listlength
        INTEGER :: i, j

        do i = 1, Iint
                do j = 1, listlength
                        if(linelist(j)%wavelength == ILs(i)%wavelength)then 
                                ILs(i)%intensity = linelist(j)%intensity
                                ILs(i)%int_err   = linelist(j)%int_err
                                cycle
                        endif        
                end do 
        end do

end subroutine

subroutine get_H(H_BS, linelist, listlength)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: H_BS
        TYPE(line), DIMENSION(:) :: linelist 
        INTEGER :: i, j, listlength
        REAL*8 :: HW = 0.00000000
        CHARACTER*10 :: blank 
        !another ugly kludge, but it works.

        do i = 1, 4
                if(i == 1)then
                        blank = "Halpha     "
                        HW = 6562.77D0 
                elseif(i == 2)then
                        blank = "Hbeta      "
                        HW = 4861.33D0
                elseif(i == 3)then
                        blank = "Hgamma     "
                        HW = 4340.47D0
                elseif(i == 4)then
                        blank = "Hdelta     "
                        HW = 4101.74D0
                else
                        PRINT*, "This is an EX-PARROT!!"
                endif        

                do j = 1, listlength
                         if (linelist(j)%wavelength-HW==0) then
                                H_BS(i)%name = blank
                                H_BS(i)%wavelength = linelist(j)%wavelength
                                H_BS(i)%intensity = linelist(j)%intensity
                                H_BS(i)%int_err = linelist(j)%int_err
                        endif
                end do
        end do

end subroutine

subroutine get_He(He_lines, linelist,listlength)
        TYPE(line), DIMENSION(4), INTENT(OUT) :: He_lines
        TYPE(line), DIMENSION(:), INTENT(IN) :: linelist
        INTEGER :: i, j, listlength
        REAL*8 :: HW
        CHARACTER*10 :: blank
        !another ugly kludge, but it works.  
        do i = 1, 4
                if(i == 1)then
                        blank = "HeII4686   " 
                        HW = 4685.68D0 
                elseif(i == 2)then
                        blank = "HeI4471    "
                        HW = 4471.50D0
                elseif(i == 3)then
                        blank = "HeI5876    "
                        HW = 5875.66D0
                elseif(i == 4)then
                        blank = "HeI6678    "
                        HW = 6678.16D0
                else
                        PRINT*, "This is an EX-PARROT!!"
                endif
                He_lines(i)%name = blank
                He_lines(i)%wavelength = HW
                He_lines(i)%intensity = 0.0
                He_lines(i)%int_err = 0.0
                do j = 1, listlength
                        if(linelist(j)%wavelength == HW) then 
                                He_lines(i)%intensity = linelist(j)%intensity
                                He_lines(i)%int_err = linelist(j)%int_err
                        endif
                end do
        end do

end subroutine

!extinction laws now in mod_extinction

end module mod_abundmaths 

module mod_atomicdata !all data types and subroutines related to reading and storing atomic data, required by all other modules that use atomic data
use mod_recombination_lines
implicit none
!data types first, then subroutines
TYPE atomic_data  !data type for storing atomic data relevent to forbidden lines
        CHARACTER*10 ion
	CHARACTER*20 label(150)
        REAL*8 T(35)
        REAL*8 ROOTT(35)
        REAL*8 wave_num(150)
        REAL*8 A_Coefficients(150,150)
        REAL*8 Collision_strengths(35,150,150)
        INTEGER IRATS
        INTEGER G(150)
        INTEGER NTEMPS, NLEVS

END TYPE

!uses data types definte in rec_lines.f90 for recombination line data

contains
!Atomic data is required by lines 300-400 of diagnostic_equib and lines 247-347 of abundance_equib

!***ATOMIC DATA FORMAT***!

!Integer NCOMS, number of subsequent lines of comments
!NCOMS lines of comments
!NLEVS (number of levels), NTEMPS (number of temperatures included in data), 0 (input format, always 0)
!NLEV labels for levels
!NTEMP temperatures
!0 (switch) [if 0 then collision strengths are tabulated]
!table of data
! First collision strengths, can ignore temps as each section of table is in same order as temperatures listed previously
! Then Transition probabilities
! Finally statistical weights of levels


subroutine read_data(ion,label,wave_num,A_coeffs,T,ROOTT,col_str,switch,G,NTEMPS,NLEVS) !performs actual reading of data from files
	CHARACTER*10  :: ion
	CHARACTER*20, DIMENSION(:) :: label
	CHARACTER*1 COMMENTS(78)
	integer NCOMS, NLEVS, NLEV1, NTEMPS, N, GX, I, J, K, KP1, L, ID(2), JD(2), switch, temp
        integer, DIMENSION(:) :: G
        REAL*8 WN, AX, QX
        REAL*8, DIMENSION(:) :: T, ROOTT, wave_num
        REAL*8, DIMENSION(:,:) :: A_coeffs
        REAL*8, DIMENSION(:,:,:) :: col_str

OPEN(Unit=1,STATUS='OLD',FILE='atomic_data06/'//ion(1:(INDEX(ION,' ') - 1))//'.dat',ACTION='READ')
print*, ion
READ(1,*) NCOMS
do I=1,NCOMS
	READ(1,1003) COMMENTS
enddo

READ(1,*)NLEVS,NTEMPS

do I=1,NLEVS
	READ(1,1002) label(I) !read in level labels
enddo

do I=1,NTEMPS
	READ(1,*) T(I) !read in temperatures
            T(I) = LOG10 (T(I))
            ROOTT(I) = SQRT(T(I))
enddo

READ(1,*) switch

!read collision strengths (diagnostic_equib ln154-171)
QX=1
K = 1
        DO WHILE (QX .gt. 0)
                READ(1,*) ID(2), JD(2), QX
                IF (QX.eq.0.D0) exit
                if (ID(2) .eq. 0) then
                   ID(2) = ID(1)
                   K = K + 1
                else
                   ID(1) = ID(2)
                   K = 1
                endif
                if (JD(2) .eq. 0) then
                   JD(2) = JD(1)
                else
                   JD(1) = JD(2)
                endif
                if (QX .ne. 0.D0) then
                I = ID(2) 
                J = JD(2) 
                col_str(K,I,J) = QX
                endif
        enddo
        NLEV1 = NLEVS-1
      DO K = 1,NLEV1
        KP1 = K + 1 
          DO L = KP1, NLEVS
            READ (1,*) I, J, AX  !read transition probabilities
            A_coeffs(J,I) = AX 
          ENDDO 
      ENDDO 

DO I=1,NLEVS
	READ(1,*) N, GX, WN !read wavenumbers
        G(N) = GX
        wave_num(N) = WN
enddo
	CLOSE(UNIT=1)


 1002 FORMAT(A20)
 1003 FORMAT(78A1)
end subroutine

subroutine get_RL_data(oiiRLs,niiRLs,ciiRLs,neiiRLs,xiiiRLs)
        INTEGER I
        TYPE(oiiRL), DIMENSION(415) :: oiiRLs
        TYPE(niiRL), DIMENSION(99) :: niiRLs
        TYPE(ciiRL), DIMENSION(57) :: ciiRLs
        TYPE(neiiRL), DIMENSION(38) :: neiiRLs
        TYPE(xiiiRL), DIMENSION(6) :: xiiiRLs

!Format statements
            301 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, F7.4,     &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4) !oii format

            302 FORMAT (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, 1X, F7.4, &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4) !nii format

       303 FORMAT (F7.2, 1X, F6.4, 1X, F7.4, 1X, F7.4, 1X, F7.4, 1X, F7.4) !cii format

       304 FORMAT (F7.2, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F7.4, 1X, F6.3) !neii format

       305 FORMAT (A3,1X,F7.2, 1X, F5.3, 1X, F6.3, 1X, F5.3, 1X, F5.3, 1X, F5.4) !xiii format

!Filereading

            OPEN(201, file="atomic_data06/Roii.dat", status='old') !oii
            DO i = 1,415
            READ(201,301) oiiRLs(i)%ION, oiiRLs(i)%Wave, oiiRLs(i)%Hyb, &
     &oiiRLs(i)%Rem1, oiiRLs(i)%Rem2, oiiRLs(i)%Rem3, oiiRLs(i)%Rem4,   &
     &oiiRLs(i)%gf1, oiiRLs(i)%q_gf1, oiiRLs(i)%gf2, oiiRLs(i)%q_gf2,   &
     &oiiRLs(i)%Mult, oiiRLs(i)%E1, oiiRLs(i)%n_E1, oiiRLs(i)%n_E1GA,   &
     &oiiRLs(i)%g1, oiiRLs(i)%n_g1, oiiRLs(i)%Term1, oiiRLs(i)%E2,      &
     &oiiRLs(i)%n_E2, oiiRLs(i)%n_E2GA, oiiRLs(i)%g2, oiiRLs(i)%n_g2,   &
     &oiiRLs(i)%Term2, oiiRLs(i)%Br_A, oiiRLs(i)%Br_B, oiiRLs(i)%Br_C 

            oiiRLs(i)%EM = 0
            oiiRLs(i)%Int = 0
            oiiRLs(i)%Obs = 0
            oiiRLs(i)%abundance = 0
            END DO 
      CLOSE(201)

            OPEN(201, file="atomic_data06/Rnii.dat", status='old') !nii
            DO i = 1,99
            READ(201,302) niiRLs(i)%ION, niiRLs(i)%Wave, niiRLs(i)%Hyb, &
     &niiRLs(i)%Rem1, niiRLs(i)%Rem2, niiRLs(i)%Rem3, niiRLs(i)%Rem4,   &
     &niiRLs(i)%gf1, niiRLs(i)%q_gf1, niiRLs(i)%gf2, niiRLs(i)%q_gf2,   &
     &niiRLs(i)%Mult, niiRLs(i)%E1, niiRLs(i)%n_E1, niiRLs(i)%n_E1GA,   &
     &niiRLs(i)%g1, niiRLs(i)%n_g1, niiRLs(i)%Term1, niiRLs(i)%E2,      &
     &niiRLs(i)%n_E2, niiRLs(i)%n_E2GA, niiRLs(i)%g2, niiRLs(i)%n_g2,   &
     &niiRLs(i)%Term2, niiRLs(i)%Br_LS

            niiRLs(i)%Int = 0
            niiRLs(i)%Obs = 0
            niiRLs(i)%abundance = 0
            END DO
      CLOSE(201)

       OPEN(201, file="atomic_data06/Rcii.dat", status='old') !cii
       DO i = 1,57
         READ(201,303) ciiRLs(i)%Wave, ciiRLs(i)%a, ciiRLs(i)%b, &
         & ciiRLs(i)%c, ciiRLs(i)%d, ciiRLs(i)%f

            ciiRLs(i)%Int = 0
            ciiRLs(i)%Obs = 0
            ciiRLs(i)%abundance = 0
       END DO
       CLOSE(201)

       OPEN(201, file="atomic_data06/Rneii.dat", status='old') !neii
       DO i = 1,38
         READ(201,304) neiiRLs(i)%Wave, neiiRLs(i)%a, neiiRLs(i)%b, &
         & neiiRLs(i)%c, neiiRLs(i)%d, neiiRLs(i)%f, neiiRLs(i)%Br

            neiiRLs(i)%Int = 0
            neiiRLs(i)%Obs = 0
            neiiRLs(i)%abundance = 0
       END DO
       CLOSE(201)

       OPEN(201, file="atomic_data06/Rxiii.dat", status='old') !xiii
       DO i = 1,6
         READ(201,305) xiiiRLs(i)%ion, xiiiRLs(i)%Wave, xiiiRLs(i)%a, &
         & xiiiRLs(i)%b, xiiiRLs(i)%c, xiiiRLs(i)%d, xiiiRLs(i)%Br

            xiiiRLs(i)%Int = 0
            xiiiRLs(i)%Obs = 0
            xiiiRLs(i)%abundance = 0
       END DO
       CLOSE(201)



end subroutine

subroutine get_atomicdata(atomicdata_array,No_ions,oiiRLs,niiRLs,ciiRLs,neiiRLs,xiiiRLs) !organises data into a single object

        integer :: No_ions, I
        TYPE(atomic_data), DIMENSION(No_ions) :: atomicdata_array
        TYPE(oiiRL), DIMENSION(415) :: oiiRLs
        TYPE(niiRL), DIMENSION(99) :: niiRLs
        TYPE(ciiRL), DIMENSION(57) :: ciiRLs
        TYPE(neiiRL), DIMENSION(38) :: neiiRLs
        TYPE(xiiiRL), DIMENSION(6) :: xiiiRLs

atomicdata_array(1)%ion = "ariii     "
atomicdata_array(2)%ion = "ariv      "
atomicdata_array(3)%ion = "arv       "
atomicdata_array(4)%ion = "ci        "
atomicdata_array(5)%ion = "ciii      "
atomicdata_array(6)%ion = "cliii     "
atomicdata_array(7)%ion = "neii      "
atomicdata_array(8)%ion = "neiii     "
atomicdata_array(9)%ion = "neiv      "
atomicdata_array(10)%ion = "nev       "
atomicdata_array(11)%ion = "nii       "
atomicdata_array(12)%ion = "niii      "
atomicdata_array(13)%ion = "oi        "
atomicdata_array(14)%ion = "oii       "
atomicdata_array(15)%ion = "oiii      "
atomicdata_array(16)%ion = "sii       "
atomicdata_array(17)%ion = "siii      "
atomicdata_array(18)%ion = "siv       "

DO I = 1,No_ions
        CALL read_data(atomicdata_array(I)%ion,atomicdata_array(I)%label,atomicdata_array(I)%wave_num,atomicdata_array(I)%A_Coefficients,atomicdata_array(I)%T,atomicdata_array(I)%ROOTT,atomicdata_array(I)%Collision_strengths,atomicdata_array(I)%IRATS,atomicdata_array(I)%G,atomicdata_array(I)%NTEMPS,atomicdata_array(I)%NLEVS)
enddo

        CALL get_RL_data(oiiRLs,niiRLs,ciiRLs,neiiRLs,xiiiRLs)

end subroutine

end module mod_atomicdata
