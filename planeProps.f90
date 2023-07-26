program main
        use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
        use :: compMat
        use :: csv_file
        implicit none
        
        real(dp) :: Ef, Gf, vf, Poisson_f, Em, Gm, vm, Poisson_m, E1, E2, G12, Poisson_12, t   ! Non-array  inputs
        real(dp) :: Ex, Ey, Gxy, Poisson_xy, Ebx, Eby, Gbxy, Poisson_bxy
        real(dp), dimension(3,3) :: Q, Qbar   ! Q matrix and its rotation Qbar
        real(dp), dimension(6,6) :: mainABD, ABD   ! MainABD is the sum of each iterating ABD
        real(dp), allocatable :: theta(:), tk(:), y(:)   ! Array holding the direction of each layer
        integer :: i, j, k, n   ! Loop variables, i and j are used to display matrices, k iterates over the layer loop and n is the
                                ! number of layers
        open(unit = 1, file = "compMat.csv", status = "unknown")

        ! Define laminate fiber and matrix properties
        Ef = 230e9   ! All relevant units in pascals
        Gf = 15e9
        vf = 0.3345
        Poisson_f = 0.2

        Em = 4e9
        Gm = 1.481e9
        vm = 1 - vf
        Poisson_m = 0.35
        
        ! Write material properties to csv
        call csv_write(1, "Ef", .false.)
        call csv_write(1, Ef, .true.)
        call csv_write(1, "Gf", .false.)
        call csv_write(1, Gf, .true.)
        call csv_write(1, "vf", .false.)
        call csv_write(1, vf, .true.)
        call csv_write(1, "Poisson_f", .false.)
        call csv_write(1, Poisson_f, .true.)
        call csv_write(1, "Em", .false.)
        call csv_write(1, Em, .true.)
        call csv_write(1, "Gm", .false.)
        call csv_write(1, Gm, .true.)
        call csv_write(1, "vm", .false.)
        call csv_write(1, vm, .true.)
        call csv_write(1, "Poisson_m", .false.)
        call csv_write(1, Poisson_m, .true.)
        call csv_write(1, "", .true.)

        ! Obtain plane properties
        call planeProperties(Ef, Gf, vf, Poisson_f, Em, Gm, vm, Poisson_m, E1, E2, G12, Poisson_12)

        print *, "E1 ="
        print "(ES12.6)", E1
        print *, "E2 ="
        print "(ES12.6)", E2
        print *, "G12 ="
        print "(ES12.6)", G12
        print *, "Poisson_12 =", Poisson_12
        
        ! Write plane properties to csv
        call csv_write(1, "E1", .false.)
        call csv_write(1, E1, .true.)
        call csv_write(1, "E2", .false.)
        call csv_write(1, E2, .true.)
        call csv_write(1, "G12", .false.)
        call csv_write(1, G12, .true.)
        call csv_write(1, "Poisson_12", .false.)
        call csv_write(1, Poisson_12, .true.)
        call csv_write(1, "", .true.)
        
        ! Main layer loop, finds the ABD matrix for each layer and adds them all up on the mainABD matrix to find the complete
        ! matrix for the laminate.

        n = 8            ! Define number of layers, must coincide with the amount of entries in theta and z
        mainABD(:,:) = 0 ! Main ABD matrix is initialized as a 0 matrix
        theta = [0.0, 0.0, -45.0, 45.0, 45.0, -45.0, 0.0, 0.0]  ! Define each layer's angle and thickness prior to the loop
        tk = [0.0635, 0.0635, 0.0635, 0.0635, 0.0635, 0.0635, 0.0635, 0.0635]/1000
        t = sum(tk)
        y = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        if(MOD(n, 2) == 0) then
            do k = 1,n
                    if(k < 4) then
                            y(k) = -(tk(4) + sum(tk(k:3))) + tk(k)/2
                    elseif(k == 5) then
                            y(k) = -(tk(4) + sum(tk(k:3))) + tk(k)/2
                    else
                            y(k) = (tk(5) + sum(tk(6:k))) - tk(k)/2
                    endif
            enddo
        else
            do k = 1,n
                    if(k < 4) then
                            y(k) = -(tk(4)/1 + sum(tk(k:3))) + tk(k)/2
                    elseif(k == 4) then
                            y(k) = 0 
                    else
                            y(k) = (tk(4)/2 + sum(tk(5:k))) - tk(k)/2
                    endif
            enddo
        endif

        !print *, "y =", y   ! For debugging purposes.

        do k = 1,n 
                print *, "--------------------"
                print *, "Angle = ", theta(k)

                call QMatrix(E1, E2, G12, Poisson_12, Q)
                print *, "Q # ", k
                do i=1, 3
                        write(*, '(3(ES12.6, "  "))') (Q(i,j), j=1,3)
                end do
                call csv_write(1, "Q @ theta = ", .false.)
                call csv_write(1, theta(k), .true.)
                call csv_write(1, Q)
                call csv_write(1, "", .true.)

                call getQbarMatrix(Q, theta(k), Qbar)
                print *, "Qbar # ", k
                do i=1, 3
                        write(*, '(3(ES12.6, "  "))') (Qbar(i,j), j=1,3)
                end do
                call csv_write(1, "Qbar @ theta = ", .false.)
                call csv_write(1, theta(k), .true.)
                call csv_write(1, Qbar)
                call csv_write(1, "", .true.)
                
                ! A matrix in N/m
                ! B matrix in N
                ! D matrix in Nm

                call getABDMatrix(Qbar, tk(k), y(k), ABD)
                print *, "ABD # ", k
                do i=1, 6
                        write(*, '(6(ES16.6, "  "))') (ABD(i,j), j=1,6)
                end do
                call csv_write(1, "Component ABD @ theta = ", .false.)
                call csv_write(1, theta(k), .true.)
                call csv_write(1, ABD)
                call csv_write(1, "", .true.)

                mainABD = mainABD + ABD

                print *, "--------------------"
        end do
        
        ! Print out the complete ABD matrix, mainABD
        print *, "--------------------"
                print *, "Full ABD matrix"
                do i=1, 6
                        write(*, '(6(ES16.6, "  "))') (mainABD(i,j), j=1,6)
                end do       
        print *, "--------------------"

        call csv_write(1, "ABD matrix", .true.)
        call csv_write(1, mainABD)
        call csv_write(1, "", .true.)
        
        ! Find the laminate's equivalent moduli for both tension and beding stresses
        call getEquivalentMod(mainABD, t, Ex, Ey, Gxy, Poisson_xy, Ebx, Eby, Gbxy, Poisson_bxy)
        print *, "--------------------"
        print *, "Equivalent moduli (Units in pascals)"
        print *, "Ex ="
        print "(ES12.6)", Ex
        print *, "Ey ="
        print "(ES12.6)", Ey
        print *, "Gxy ="
        print "(ES12.6)", Gxy
        print *, "Poisson_xy =", Poisson_xy
        print *, "Ebx ="
        print "(ES12.6)", Ebx
        print *, "Eby ="
        print "(ES12.6)", Eby
        print *, "Gbxy ="
        print "(ES12.6)", Gbxy
        print *, "Poisson_bxy =", Poisson_bxy

        ! Write equivalent moduli to csv
        call csv_write(1, "Ex", .false.)
        call csv_write(1, Ex, .true.)
        call csv_write(1, "Ey", .false.)
        call csv_write(1, Ey, .true.)
        call csv_write(1, "Gxy", .false.)
        call csv_write(1, Gxy, .true.)
        call csv_write(1, "Poisson_xy", .false.)
        call csv_write(1, Poisson_xy, .true.)
        call csv_write(1, "Ebx", .false.)
        call csv_write(1, Ebx, .true.)
        call csv_write(1, "Eby", .false.)
        call csv_write(1, Eby, .true.)
        call csv_write(1, "Gbxy", .false.)
        call csv_write(1, Gbxy, .true.)
        call csv_write(1, "Poisson_bxy", .false.)
        call csv_write(1, Poisson_bxy, .true.)
end program main
