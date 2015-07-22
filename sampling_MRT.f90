!======================================================================
!     Version:  June 30
!     Sctering algorithm
!     Matrix operations in diferent file
!     Commented print statements are for degguging and testing the 
!     final version won't have those
!======================================================================



!!!!!!The program test is just for testing my sampling and 
!!!!!!my matrices !!!!!!!!!!!!!!!!!!!!!!!!!!!
!program test
!
!    real:: j
!    real, dimension(4):: vect
!    vect(1) = 0.00
!    vect(2) = 0.00
!    vect(3) = 1.00
!    vect(4) = 0.00
!    j = 0
!
!    !!!!!!!!!! test cases for rorating from photon frame to lab frame!!!!!
!    !call test_rotate_in_lab_frame    
!    CALL RANDOM_SEED()
!do while(j <=3)
!    call scat( vect(1), vect(2), vect(3) )
!    print*,"new_vect =", "(", vect,")"
!    j = j+1
!end do
!end program test

!!!!!!!!!!!!!!!!!!Main prgram to test the different subroutine!!!!!!!!!!!!!!!!!!!
subroutine scat(val_u, val_v, val_w)
    real val_u, val_v, val_w
    !!!! g, pl, pc, s are the for the dust scattering matrix later !!!!
    !real g, pl, pc, s
    !we need to cal a subroutine to caculte Pmax for complicated function.
    real, dimension(4):: vect_in_lab_frame, S_initial_ray, S_initial_dust, S_final_dust, S_final_ray
    real, dimension(4):: new_vect_in_lab_frame1, new_vect_in_lab_frame2, vect_in_photon_frame1, vect_in_photon_frame2
    real u, v, w, e
    u = val_u
    v = val_v
    w = val_w
    e = 0.00
    print*, "entering scat"
    !Initial phonton direction in lab frame
    vect_in_lab_frame = (/u, v, w, e/)
    call normalize(vect_in_lab_frame)
    S_initial_ray  = (/1.00, 0.00, 0.00, 0.00/)
    S_initial_dust  = (/1.00, 0.00, 0.00, 0.00/)
    print*, "The initial direction of the particule is (", vect_in_lab_frame, ")"
    print*, "chose l to be  1"
    !read(*,*) l
    l =1
    do while(l == 1)

        ! n = 0 is Raylieght scattering and n= 1 is dust
        n = 0
        if(n == 0) then
            call monteCarloRay(vect_in_lab_frame, vect_in_photon_frame1, new_vect_in_lab_frame1, S_initial_ray, S_final_ray)
            print*, "new  direction of the particule is (", new_vect_in_lab_frame1, ")"
            val_u = new_vect_in_lab_frame1(1)
            val_v = new_vect_in_lab_frame1(2)
            val_w = new_vect_in_lab_frame1(3)
            print*, "S_initial_ray = (",S_initial_ray, ")"
            print*, "S_final_ray = (",S_final_ray,")"
            l = 2
        else if(n == 1) then
            call monteCarloDust(vect_in_lab_frame, vect_in_photon_frame2, new_vect_in_lab_frame2, S_initial_dust, S_final_dust)
            print*, "new scat direction=(",new_vect_in_lab_frame2,")"
            val_u = new_vect_in_lab_frame2(1)
            val_v = new_vect_in_lab_frame2(2)
            val_w = new_vect_in_lab_frame2(3)
            print*, "S_initial_dust = (",S_initial_dust, ")"
            print*, "S_final_dust = (",S_final_dust,")"
            l = 2
        else
           print*, "you did not enter 0 or 1 try again"
           l =1
        end if
    end do 
    print*, "leaving scat"
        
end subroutine scat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!Monte Dust!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine monteCarloDust(vect_in_lab_frame, vect_in_photon_frame2, new_vect_in_lab_frame2, S_initial_dust, S_final_dust)
    real i2, M2, Pmax2
    !!!! g, pl, pc, s are the for the dust scattering matrix later !!!!
    real g, pl, pc, s
    !we need to cal a subroutine to caculte Pmax for complicated function.
    real, dimension(4):: vect_in_lab_frame, S_initial_dust, S_final_dust
    real, dimension(4):: new_vect_in_lab_frame2, vect_in_photon_frame2

    Pmax2 = 2.00
    !print*, "g is the scattering asymetry parameter ranging from 0 to 1"
    !print*, "set g ="
    !read(*,*) g
    g = 0.50
    !print*, "pl is the maximum linear polarization"
    !print*, "set pl ="
    !read(*,*) pl
    pl = .60
    !print*, "pc is the peak circular polarization"
    !print*, "pc ="
    !read(*,*)pc
    pc = 0.00
    !print*, "we are using 1 for s, the skew factor"
    s = 1.00
    !Initial phonton direction in lab frame

    
    ! random direction sampling in photon frame using the rejection method
    !print*, "started sampling"
    call sample_dust_scattering(i2, M2, Pmax2, g)
    
    !            print*, "i M =(", i2, ",", M2, ")"
    !print*, "done sampling"
    call Stoke_Vect_Dust(M2, g, pl, pc, s, S_initial_dust, vect_in_lab_frame, S_final_dust)
    vect_in_photon_frame2 = (/sin ( acos(M2) ) * cos(i2), sin( acos(M2)) * sin(i2), M2, 0.00/)
    call normalize(vect_in_photon_frame2)
    print*, "vect in photon frame = (",vect_in_photon_frame2,")"
    
    ! new direction after scattering 
    call rotate_in_lab_frame(vect_in_lab_frame, vect_in_photon_frame2, new_vect_in_lab_frame2)
    print*, "new scat direction=(",new_vect_in_lab_frame2,")"
end subroutine monteCarloDust
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!Monte rayleihgth!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine monteCarloRay(vect_in_lab_frame, vect_in_photon_frame1, new_vect_in_lab_frame1, S_initial_ray, S_final_ray)
    real Pmax1, i1, M1
    !we need to cal a subroutine to caculte Pmax for complicated function.
    real, dimension(4):: vect_in_lab_frame, S_initial_ray, S_final_ray
    real, dimension(4):: new_vect_in_lab_frame1, vect_in_photon_frame1
    real u, v, w, e
    
    Pmax1 = 2.00

    print*, "I am in monte carlo"
    ! random direction sampling in photon frame using the rejection method
    call sample_rayleight_scattering(i1, M1, Pmax1)
    print*, "M and i1 = (", M1, ",", i1,")"
    call Stoke_Vect_Ray(M1, S_initial_ray, vect_in_lab_frame, S_final_ray)

    vect_in_photon_frame1 = (/sin ( acos(M1) ) * cos(i1), sin( acos(M1)) * sin(i1), M1, 0.00/)
    call normalize(vect_in_photon_frame1)
    print*, "vect in photon =(", vect_in_photon_frame1, ")"
    print*, "vect in in lab frame=(", vect_in_lab_frame, ")"

    ! new direction after scattering 
    call rotate_in_lab_frame(vect_in_lab_frame, vect_in_photon_frame1, new_vect_in_lab_frame1)
    call normalize(new_vect_in_lab_frame1)
    print*, "new vect in lab =(", new_vect_in_lab_frame1, ")"

end subroutine monteCarloRay

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Rayleight sampling!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sample_rayleight_scattering(i, M, Pmax)
    ! to hold the random number genrated!
    ! i (azymuthal angle, M cos^2(theta) !)
    ! initializing variables !
    real :: epsilon1 = 0, epsilon2 = 0
    real :: i, M, discriminant
    real :: pi = acos( -1.00)

    ! set polorized to 0 for non-polorized and 1 for polorized
    polorized = 0
    j = 1
    do while( j > 0 )
        call random_number(epsilon1)
        call random_number(epsilon2)
        !PRINT *, epsilon1
        !print *,"epsilon1 =", epsilon2
        !print *,"epsilon2 =", epsilon1

        ! This for non-polorized case 
        if(polorized == 0) then
            i = 2 * pi * epsilon1
            discriminant = 32*epsilon2 - 7.00 
            !print*, "discriminat =", discriminant
            if(discriminant > 0.00) then
                M =  0.50*(sqrt(discriminant) -3.00)
                if(M <=1 .and. M>=-1) then
                    j = 0
                else
                    i = 1
                end if
            else
                j = 1
            end if
        ! for polorized case
        else
            i = 2 * pi * epsilon1
            ! brute force sampling method below bu slower
            P = epsilon2 * Pmax
            M = 2 * epsilon2 - 1
            P_i_M = 1 + M * M
            if( P < P_i_M ) then
                j = 0
            end if
        end if
    end do
!    print *, 'i is the azymuthal angle'
!    print *,'i =', i
!    print *, 'M is the direction cos(theta)'
!    print *, 'M =', M
!    print *, 'i is the randomly sample using epsilon1'
!    print *, 'P =', P
!    print *, 'P_i_M is the vavule of the probability distribution calculated from i and M '
!    print *, 'P_i_M=', P_i_M
end subroutine sample_rayleight_scattering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Dust sampling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sample_dust_scattering(i, M, Pmax, g)
    real :: epsilon1, epsilon2, P, i, P_i_M, Pmax, M, g
    real :: pi
    pi  = acos(-1.00)
    !print*, "enter 1 for unpolorized scattering"
    !read(*,*) option 
    polorized = 0
    k = 1
    do while(k > 0)
        call random_number(epsilon1)
        call random_number(epsilon2)
    
        if(polorized  == 0) then

            if((1 - g + 2 * g * epsilon2) <= 0.0000001 .and.(1 - g + 2 * g * epsilon2)>= -.0000001) then
                k = 1
            else
                M = ( 1 + g*g - ((1 - g * g) / (1 - g + 2 * g * epsilon2)) * ((1 - g * g) / (1 - g + 2 * g * epsilon2)) ) / (2 * g)
                i = 2 * pi * epsilon1
                !print*, "i M =(", i, ",", M, ")"
                k = -1
            end if

        ! polorized case
        else
            !PRINT *, epsilon1
            !print *, epsilon2
            i = 2 * pi * epsilon1
            P = epsilon2 * Pmax
            M = 2 * epsilon2 - 1
            P_i_M = 1 + M * M
            if(P < P_i_M) then
                k = -1
            end if
        end if
    end do
    !print*, 'Dust Scatering **************'
    !print *, 'i is the azymuthal angle'
    !print *,'i =', i
    !print *, 'M is the direction cos(theta)'
    !print *, 'M =', M
    !print *, 'i is the randomly sample using epsilon2'
    !print *, 'P =', P
    !print *, 'P_i_M is the vavule of the probability distribution calculated from i and M '
    !print *, 'P_i_M=', P_i_M
end subroutine sample_dust_scattering

!!!!!!! this subroutine calculates the Stoke Vectord to Rayleigh Scattering!!!!!!!!!!!!!
subroutine Stoke_Vect_Ray(M, S_initial, vect, S_final)
real:: M
!!S_initial and S_final represnet the stoke vectors
real, dimension(4)::vect, S_initial, S_final
real, dimension(4, 4):: scatMatrix, rot_inMatrix, rot_outMatrix, LR, LRL
call Rayleight_Scattering_Matrix(M, scatMatrix)
!print*, "M =", M
!print*, scatMatrix
call RotateInMatrix(vect, rot_inMatrix)
call RotateOutMatrix(rot_inMatrix, rot_outMatrix)
call MatrixMult_4X4_4X4(rot_inMatrix, scatMatrix, LR)
call MatrixMult_4X4_4X4(LR, rot_outMatrix, LRL)
call MatrixMult_4X4_4X1(LRL, S_initial, S_final)
!print*, S_final
end subroutine Stoke_Vect_Ray


!!!!!!!! This subroutine calsultes te Stike Vectors for Dust Scattering !!!!!!!!!!!!!
subroutine Stoke_vect_Dust(M, g, pl, pc, s, S_initial, vect, S_final)
real:: M, g, pl, pc, s
!!S_initial and S_final represnet the stoke vectors
real, dimension(4):: vect, S_initial, S_final
real, dimension(4, 4):: scatMatrix, rot_inMatrix, rot_outMatrix, LR, LRL
call Dust_Scattering_Matrix(M, g, pl, pc, s, scatMatrix)
call RotateInMatrix(vect, rot_inMatrix)
call RotateOutMatrix(rot_inMatrix, rot_outMatrix)
call MatrixMult_4X4_4X4(rot_inMatrix, scatMatrix, LR)
call MatrixMult_4X4_4X4(LR, rot_outMatrix, LRL)
call MatrixMult_4X4_4X1(LRL, S_initial, S_final)
end subroutine Stoke_vect_Dust

subroutine test_rotate_in_lab_frame

    real, dimension(4):: new_vect_in_lab_frame, vect_in_photon_frame, vect_in_lab_frame
    vect_in_photon_frame(1) = 0
    vect_in_photon_frame(2) = 0
    vect_in_photon_frame(3) = 1
    vect_in_photon_frame(4) = 0

    vect_in_lab_frame(1) = 0
    vect_in_lab_frame(2) = 0
    vect_in_lab_frame(3) = 1
    vect_in_lab_frame(4) = 0
    call rotate_in_lab_frame(vect_in_lab_frame, vect_in_photon_frame, new_vect_in_lab_frame)
    !print*, "lab = (", vect_in_lab_frame,")"
    !print*, "phton frame = (", vect_in_photon_frame,")"
    !print*, "new lab = (", new_vect_in_lab_frame,")"

    vect_in_photon_frame(1) = 0
    vect_in_photon_frame(2) = 0
    vect_in_photon_frame(3) = 1
    vect_in_photon_frame(4) = 0

    vect_in_lab_frame(1) = 0
    vect_in_lab_frame(2) = 1
    vect_in_lab_frame(3) = 0
    vect_in_lab_frame(4) = 0
    call rotate_in_lab_frame(vect_in_lab_frame, vect_in_photon_frame, new_vect_in_lab_frame)
    !print*, "lab = (", vect_in_lab_frame,")"
    !print*, "phton frame = (", vect_in_photon_frame,")"
    !print*, "new lab = (", new_vect_in_lab_frame,")"

end subroutine test_rotate_in_lab_frame



