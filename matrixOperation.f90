!======================================================================
!
!     Version:  August 6 2015
!     This file contains tools for our scattering to operate 
!     properly (matix operation)
!
!======================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! This subroutine rotates into the phton frame !!!!!!!!!
subroutine rotate_in_photon_frame(vect, output_vect)
    real, dimension(4):: vect, output_vect
   
    !photon direction in lab frame
    real u, v, w
    
    real, dimension(4,4):: xymatrix, zymatrix, totalmatrix
    real, dimension(4,4):: totalmatrix1
    !!! initializing direction of photon in lab frame !!!
    u = vect(1)
    v = vect(2)
    w = vect(3)
    e = 0.00

    !print*, "vector in lab frame :", vect

    !!!!!!!!! Rotation matrix of y-x about the z axis !!!!!!!!
    xymatrix(1,1) = u / sqrt(1- w*w)
    xymatrix(1,2) = v / sqrt(1- w*w)
    xymatrix(1,3) = 0.00
    xymatrix(1,4) = 0.00
    xymatrix(2,1) = - v / sqrt(1- w*w)
    xymatrix(2,2) = u / sqrt(1- w*w)
    xymatrix(2,3) = 0.00
    xymatrix(2,4) = 0.00
    xymatrix(3,1) = 0.00
    xymatrix(3,2) = 0.00
    xymatrix(3,3) = 1.00
    xymatrix(3,4) = 4.00
    xymatrix(4,1) = 0.00
    xymatrix(4,2) = 0.00
    xymatrix(4,3) = 0.00
    xymatrix(4,4) = 1.00

    !print*,"xymatrix \n", xymatrix

    !!!!!!!! Roatation matrix of z-y about the x axis!!!!!!!!!!!
    zymatrix(1,1) = w
    zymatrix(1,2) = 0.00
    zymatrix(1,3) = - sqrt(1 - w*w)
    zymatrix(1,4) = 0.00
    zymatrix(2,1) = 0.00
    zymatrix(2,2) = 1.00
    zymatrix(2,3) = 0.00
    zymatrix(2,4) = 0.00
    zymatrix(3,1) = sqrt(1 - w*w)
    zymatrix(3,2) = 0.00
    zymatrix(3,3) = w
    zymatrix(3,4) = 0.00
    zymatrix(4,1) = 0.00
    zymatrix(4,2) = 0.00
    zymatrix(4,3) = 0.00
    zymatrix(4,4) = 1.00

    !print*, "zymatrix", zymatrix

    !!!!!! total rotation into photon frame !!!!!!!!!!!!!!
    call RotateInMatrix(vect, totalmatrix)
    !!!!!!!!!!!!This is just to test my xy and zy matrices !!!!!!
    !call MatrixMult_4X4_4X4(zymatrix, xymatrix, totalmatrix1)
    !print*, "total matrix", totalmatrix
    !print*, "total matrix1", totalmatrix1

    !!!!!!!!!! vector in the photon frame !!!!!!!!!!!!!!!!!!!
    call MatrixMult_4X4_4X1(totalmatrix, vect, output_vect)
    !print*, "in photon frame", output_vect
end subroutine rotate_in_photon_frame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! This subroutine rotates out of the photon frame !!!!!!!!
subroutine rotate_in_lab_frame(input_vect1, input_vect2, output_vect)
    real, dimension(4,4):: totalmatrix, invertedM
    real, dimension(4):: input_vect1, input_vect2, output_vect

    call RotateInMatrix(input_vect1, totalmatrix)
    !!!!!!!!!! invert total matrix !!!!!!!!!!!!!!!!!!!!!
    call InvertMatrix_4X4(totalmatrix, invertedM)

    !!!!!!!!!! vector in the photon frame !!!!!!!!!!!!!!!!!!!
    call MatrixMult_4X4_4X1(invertedM, input_vect2, output_vect)
end subroutine rotate_in_lab_frame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!Rotate in photon frame!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RotateInMatrix(input_vect1, totalmatrix)
    real, dimension(4, 4):: totalmatrix
    real, dimension(4):: input_vect1
    real u, v, w, e
    u = input_vect1(1)
    v = input_vect1(2)
    w = input_vect1(3)
    e = 0.00
    !when w = 1 no need to rotate!
    if(w == 1) then
        totalmatrix(1,1) = 1.00  
        totalmatrix(1,2) = 0.00
        totalmatrix(1,3) = 0.00  
        totalmatrix(1,4) = 0.00 
        totalmatrix(2,1) = 0.00  
        totalmatrix(2,2) = 1.00 
        totalmatrix(2,3) = 0.00  
        totalmatrix(2,4) = 0.00  
        totalmatrix(3,1) = 0.00
        totalmatrix(3,2) = 0.00
        totalmatrix(3,3) = 1.00  
        totalmatrix(3,4) = 0.00  
        totalmatrix(4,1) = 0.00  
        totalmatrix(4,2) = 0.00 
        totalmatrix(4,3) = 0.00  
        totalmatrix(4,4) = 1.00  
    else
        !!!!!! total rotation into photon frame !!!!!!!!!!!!!!
        totalmatrix(1,1) = w*u/sqrt(1-w*w)  
        totalmatrix(1,2) = u*v/sqrt(1-w*w)  
        totalmatrix(1,3) = - sqrt(1-w*w)  
        totalmatrix(1,4) = 0.00 
        totalmatrix(2,1) = - v/sqrt(1-w*w)  
        totalmatrix(2,2) = u/sqrt(1-w*w)  
        totalmatrix(2,3) = 0.00  
        totalmatrix(2,4) = 0.00  
        totalmatrix(3,1) = u  
        totalmatrix(3,2) = v 
        totalmatrix(3,3) = w  
        totalmatrix(3,4) = 0.00  
        totalmatrix(4,1) = 0.00  
        totalmatrix(4,2) = 0.00 
        totalmatrix(4,3) = 0.00  
        totalmatrix(4,4) = 1.00
    end if  
end subroutine RotateInMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!Rotate out of photon frame in lab frame!!!!!!!!!!!
subroutine RotateOutMatrix(A, B)
    real, dimension(4, 4) :: A, B
    call InvertMatrix_4X4(A, B)
end subroutine RotateOutMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine multiply nxn matric with an nxn matrix !!!!
subroutine MatrixMult_4X4_4X4(A, B, C)
    real, dimension(4, 4):: A, B, C
    do i=1, 4, 1
           C(i,1) = A(i,1) * B(1,1) + A(i, 2) * B(2,1) + A(i,3) * B(3, 1) + A(i, 4) * B(4, 1)
           C(i,2) = A(i,1) * B(1,2) + A(i, 2) * B(2,2) + A(i,3) * B(3, 2) + A(i, 4) * B(4, 2)
           C(i,3) = A(i,1) * B(1,3) + A(i, 2) * B(2,3) + A(i,3) * B(3, 3) + A(i, 4) * B(4, 3)
           C(i,4) = A(i,1) * B(1,4) + A(i, 2) * B(2,4) + A(i,3) * B(3, 4) + A(i, 4) * B(4, 4)
    end do
    !print*,"I am printing C", C
end subroutine MatrixMult_4X4_4X4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This subroutine multi[plies nxn with nxm matrix !!!!!!!!!!!!!!!!
subroutine MatrixMult_4X4_4X1(A, vect, final_vect)
    real, dimension(4, 4):: A
    real, dimension(4):: vect, final_vect 
    do i=1, 4, 1
           final_vect(i) = A(i,1) * vect(1) + A(i, 2) * vect(2) + A(i,3) * vect(3) + A(i, 4) * vect(4)
    end do
    !print*, "I am printing final vect", final_vect
end subroutine MatrixMult_4X4_4X1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! This subrouinte takes the inverse of a 3x3 matrix !!!!!!!!!!!!!
subroutine InvertMatrix_4X4(A, C)
    real, dimension(4,4):: A, C
    C(1,1) = A(1,1)
    C(1,2) = A(2,1)
    C(1,3) = A(3,1)
    C(1,4) = A(4,1)
    C(2,1) = A(1,2)
    C(2,2) = A(2,2)
    C(2,3) = A(3,2)
    C(2,4) = A(4,2)
    C(3,1) = A(1,3)
    C(3,2) = A(2,3)
    C(3,3) = A(3,3)
    C(3,4) = A(4,3)
    C(4,1) = A(1,4)
    C(4,2) = A(2,4)
    C(4,3) = A(3,4)
    C(4,4) = A(4,4)
end subroutine InvertMatrix_4X4
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!Muller Matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MullerMatrix(alpha, L)
    real, dimension(4, 4):: L
    real alpha
    L(1, 1) = 1.00
    L(1, 2) = 0.00
    L(1, 3) = 0.00
    L(1, 4) = 0.00
    L(2, 1) = 0.00
    L(2, 2) = cos(2 * alpha)
    L(2, 3) = sin(2 * alpha) 
    L(3, 1) = 0.00
    L(3, 2) = -sin(2 * alpha)
    L(3, 3) = cos(2 * alpha)
    L(3, 4) = 0.00
    L(4, 1) = 0.00
    L(4, 2) = 0.00
    L(4, 3) = 0.00
    L(4, 4) = 1.00
end subroutine MullerMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Rauleight scatering Matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Rayleight_Scattering_Matrix(M, R)
    real, dimension(4, 4):: R
    real M
    R(1, 1) = M * M + 1.00
    R(1, 2) = M * M - 1.00
    R(1, 3) = 0.00
    R(1, 4) = 0.00
    R(2, 1) = M * M - 1.00
    R(2, 2) = M * M + 1.00
    R(2, 3) = 0.00
    R(3, 1) = 0.00
    R(3, 2) = 0.00
    R(3, 3) = 2.00 * M
    R(3, 4) = 0.00
    R(4, 1) = 0.00
    R(4, 2) = 0.00
    R(4, 3) = 0.00
    R(4, 4) = 2.00 * M
end subroutine Rayleight_Scattering_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! Dust Scattering Matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Dust_Scattering_Matrix(M, g, pl, pc, s, R)
    real, dimension(4, 4):: R
    real :: Mf, pi, M
    !!!! test this later!!!
    pi = acos(-1.00)
    Mf = cos( acos(M) * (1 + 3.13 * s * exp(-7 * acos(M) / pi)) )
    
    !!!! test application of the power here
    R(1, 1) = ( 1 - g* g)/ (1 + g * g - 2 * g * M)**(3/2)
    R(1, 2) = -pl * R(1, 1) * (1 - M * M) / (1 + M * M)
    R(1, 3) = 0.00
    R(1, 4) = 0.00
    R(2, 1) = -pl * R(1, 1) * (1 - M * M) / (1 + M * M)
    R(2, 2) = (1 - g* g)/ (1 + g * g - 2 * g * M)**(3/2) 
    R(2, 3) = 0.00
    R(3, 1) = 0.00
    R(3, 2) = 0.00
    R(3, 3) = R(1, 1) * (2 * M) / ( 1 + M * M)
    R(3, 4) = 0.00
    R(4, 1) = 0.00
    R(4, 2) = 0.00
    R(4, 3) = -pc * R(1, 1) * (1 - Mf * Mf) / (1 + Mf * Mf)  
    R(4, 4) =  R(1, 1) * (2 * M) / ( 1 + M * M)

end subroutine Dust_Scattering_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!This subroutine normalizes vectors!!!!!!!!!!!!!!!!
subroutine normalize(vect)
    real, dimension(4):: vect
    real:: temp

    temp = sqrt( vect(1)*vect(1) + vect(2)*vect(2) + vect(3)*vect(3) + vect(4)*vect(4) )

    vect(1) = vect(1) / temp
    vect(2) = vect(2) / temp
    vect(3) = vect(3) / temp
    vect(4) = vect(4) / temp
end subroutine normalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine uses vetor componants to !
! find angle theta and phi                 !
subroutine getThetaPhi(u, v, w, phi1, theta1)
    real :: u, v, w
    real :: theta1, phi1
    if(u == 0) then
        phi1 = 0
    else
        phi1 = atan(v/u)
    endif
    if(w==0) then
        theta1 = 0
    else
        theta1 = atan( sqrt(u*u + v*v) / w)
    endif

end subroutine getThetaPhi


subroutine getUVW(phi, theta, u, v, w)
    u = cos(phi)*sin(theta)
    v = sin(phi)*sin(theta)
    w = cos(theta)
end subroutine getUVW
