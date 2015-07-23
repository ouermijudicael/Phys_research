program testing
    call testRotationMatrices
    call samplingTestRayleight
    call monteTestRayleight
end program testing

subroutine testRotationMatrices
    real, dimension(4):: vect_in_lab_frame, vect_in_photon_frame
    real, dimension(4):: vect_in_lab_frame_result, vect_in_photon_frame1
    real, dimension(4):: vect0, vect1, vect2, vect3, vect4, vect5, vect6, vect7, vect8


    real u, v, w, e
    u = 1/sqrt(3.00)
    v = 1/sqrt(3.00)
    w = 1/sqrt(3.00)
    e = 0.00
    vect_in_lab_frame = (/u, v, w, e/)
    vect_in_photon_frame1 = (/0, 0, 1, 0/)
    vect0 = (/0, 0, 1, 0/) ! in the z direction
    vect3 = (/0, 1, 0, 0/) ! in y direction
    vect6 = (/1, 0, 0, 0/) ! in x direction
    
    !!!!!!!!!!! call for rotating in the photon frame !!!!!!!!!!
    print*, "#################  testing ################################"
    print*, "initial photon direction in lab frame:  ", vect_in_lab_frame 
    call rotate_in_photon_frame(vect_in_lab_frame, vect_in_photon_frame)
    print*,"photon direction in photon fram after rotation to photon frame:", vect_in_photon_frame
    call rotate_in_lab_frame(vect_in_lab_frame, vect_in_photon_frame1, vect_in_lab_frame_result)
    print*, "photon direction in lab frame after rotation from photon frame to lab frame:", vect_in_lab_frame_result
    
    print*, "#################  testing ################################"
    print*, "initial photon direction in lab frame:  ", vect0 
    call rotate_in_photon_frame(vect0, vect1)
    print*,"photon direction in photon fram after rotation to photon frame:", vect1
    call rotate_in_lab_frame(vect0, vect1, vect2)
    print*, "photon direction in lab frame after rotation from photon frame to lab frame:", vect2

    print*, "#################  testing ################################"
    print*, "initial photon direction in lab frame:  ", vect3 
    call rotate_in_photon_frame(vect3, vect4)
    print*,"photon direction in photon fram after rotation to photon frame:", vect1
    call rotate_in_lab_frame(vect3, vect4, vect5)
    print*, "photon direction in lab frame after rotation from photon frame to lab frame:", vect5

    print*, "#################  testing ################################"
    print*, "initial photon direction in lab frame:  ", vect6 
    call rotate_in_photon_frame(vect6, vect7)
    print*,"photon direction in photon fram after rotation to photon frame:", vect7
    call rotate_in_lab_frame(vect6, vect7, vect8)
    print*, "photon direction in lab frame after rotation from photon frame to lab frame:", vect8

end subroutine testRotationMatrices


subroutine samplingTestRayleight
    real:: i, M, Pmax, j
    Pmax = 2
    j = 0
    do while(j < 10)    
        call sample_rayleight_scattering(i, M, Pmax)
        print*, "####################testing sampling Rayleight #################################"
        if(i < 6.3 .and. M < 2) then
            print*, "i, M = (", i,", ", M, ")"
        else
           print*, " THERE IS A PROBLEM FIX IT"
           print*, "i, M = (", i,", ", M, ")"
        end if
        j = j + 1
    end do

end subroutine samplingTestRayleight


subroutine  monteTestRayleight
real, dimension(4):: vect0, vect1, vect2, S0, S1
    vect0 = (/0, 0, 1, 0/)
    S0 = (/1, 0, 0, 0/)
    call monteCarloRay(vect0, vect1, vect3, S0, S1)
    print*, "photon direction in lab =", vect0
    print*, "photon direction in photon frame =", vect1
    print*, "new photon direction in lab =", vect2
    print*, "S initial =", S0
    print*, "S final =", S1
end subroutine monteTestRayleight
