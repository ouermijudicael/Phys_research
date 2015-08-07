!====================================================================
!	Version: 2 August 6 2015
!       Sactering algorithm
!       getting subroutine from different file
!       Comment are for duggugin purposes
!=====================================================================

!!!!!!!!!!!!!!!!!!!
! scat subroutine !
subroutine scat2(u, v, w, S_old)
    real:: u, v, w
    real, dimension(4):: S_old, S_new
    real:: i1, theta1, phi1, big_theta, i2, theta2, phi2
    call getThetaPhi(u, v, w, phi1, theta1)
    call samplei1(i1)
    call sampleBigTheta(i1, big_theta)
    call calculatei2ThetaPhi(i1, theta1, phi1, big_theta, i2, theta2, phi2)
    call calculateS(i1, i2, big_theta, S_old, S_new)
    call getUVW(phi, theta, u, v, w)
    S_old = S_new 

end subroutine scat2

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to sample i1 !
subroutine samplei1(i1)
    real :: i1, epsilon1 = 0
    real :: pi = acos(-1.00)
    
    call random_number(epsilon1)
    i1 = 2 * pi * epsilon1 
end subroutine samplei1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to sample big Theta from from the sactering matrix R !
subroutine sampleBigTheta(i1, big_theta)
    ! initializing variables !
    real :: epsilon2 = 0
    real :: i1, M, discriminant, big_theta
    real :: pi = acos( -1.00)

    ! set polorized to 0 for non-polorized and 1 for polorized
    polorized = 0
    j = 1
    do while( j > 0 )
        call random_number(epsilon2)

        ! This for non-polorized case 
        if(polorized == 0) then
            discriminant = 32*epsilon2 - 7.00
            !print*, "discriminat =", discriminant
            if(discriminant > 0.00) then
                M =  0.50*(sqrt(discriminant) -3.00)
                if(M <=1 .and. M>=-1) then
                    j = 0
                end if
            end if
        ! for polorized case
        !else
        !    i = 2 * pi * epsilon1
        !    ! brute force sampling method below bu slower
        !    P = epsilon2 * Pmax
        !    M = 2 * epsilon2 - 1
        !    P_i_M = 1 + M * M
        !    if( P < P_i_M ) then
        !        j = 0
        !    end if
        end if
    end do
    big_theta = acos(M)
end subroutine sampleBigTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the new i2, theta and phi      !
! Th order in which the calculations are done is important !
subroutine calculatei2ThetaPhi(i1, theta1, phi1, big_theta, i2, theta2, phi2)
    real:: i1, theta1, phi1, big_theta, i2, theta2, phi2
    
    ! calculation of theta2 !
    theta2 = acos( cos(big_theta)*cos(theta1) + sin(big_theta)*sin(theta1)*cos(i1) )
    ! calculation of i2 !
    i2 = acos( ( cos(theta1)-cos(theta2)*cos(big_theta) ) / ( sin(theta2)*sin(big_theta) ) )

    ! calculate phi2 !
    ! what happens when big_theta is positive or negative !
    phi2 = -acos( (cos(big_theta)-cos(theta1)*cos(theta2) ) / ( sin(theta1)*sin(theta2) ) ) + phi1

end subroutine calculatei2ThetaPhi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine calculates the new Stoke parameter !
! using S_new = L(pi - i2)R L(-i1)                   !
subroutine calculateS(i1, i2, big_theta, S_old, S_new)
    real :: i1, i2, big_theta
    real, dimension(4):: Sold, S_new
    real, dimension(4, 4)::L1, L2, R, LR, LRL 
    call MullerMatrix(i2, L2)
    call MullerMatrix(i1, L1)

    ! cos(big_theta) = M
    call Rayleight_Scattering_Matrix(cos(big_theta), R)
    call MatrixMult_4x4_4X4(L2, R, LR)
    call MatrixMult_4x4_4X4(LR, L2, LRL)
    call MatrixMult_4x4_4X1(LRL, S_old, S_new)
end subroutine calculateS
