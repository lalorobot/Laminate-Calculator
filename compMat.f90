module compMat
contains
    subroutine planeProperties(Ef, Gf, vf, Poisson_f, Em, Gm, vm, Poisson_m, E1, E2, G12, Poisson_12)
        use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
        implicit none

        real(dp), intent(in) :: Ef, Gf, vf, Poisson_f, Em, Gm, vm, Poisson_m
        real(dp), intent(out) :: E1, E2, G12, Poisson_12

        E1 = vm*Em + vf*Ef
        E2 = 1/(vm/Em + vf/Ef)
        G12 = 1/(vm/Gm + vf/Gf)
        Poisson_12 = Poisson_m*vm + Poisson_f*vf

    end subroutine planeProperties

    subroutine QMatrix(E1, E2, G12, Poisson_12, Q)
        use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
        implicit none
        
        real(dp), intent(in) :: E1, E2, G12, Poisson_12
        real(dp), dimension(3,3), intent(out) :: Q

        real(dp) :: d

        Q(:,:) = 0 ! Q matrix initialization
        d = 1 - Poisson_12**2 * E2/E1 ! Delta factor

        Q(1,1) = E1/d
        Q(1,2) = Poisson_12 * E2/d
        Q(2,1) = Q(1,2)
        Q(2,2) = E2/d
        Q(3,3) = G12
        
    end subroutine QMatrix

    subroutine getQbarMatrix(Q, thetaDeg, Qbar)
        use, intrinsic :: iso_fortran_env, only: sp =>real32, dp=>real64
        implicit none

        real(dp), dimension(3,3), intent(in) :: Q
        real(dp), intent(in) :: thetaDeg
        real(dp), dimension(3,3), intent(out) :: Qbar
        real(dp) :: c, s, theta, pi

        Qbar(:,:) = 0 ! Qbar matrix initialization
        
        pi = 3.1415927
        theta = thetaDeg*pi/180

        c = cos(theta)
        s = sin(theta)

        Qbar(1,1) = Q(1,1)*c**4 + 2*(Q(1,2) + 2*Q(3,3))*s**2 * c**2 + Q(2,2)*s**4
        Qbar(1,2) = (Q(1,1) + Q(2,2) - 4*Q(3,3))*c**2 * s**2 + Q(1,2)*(c**4 + s**4)
        Qbar(2,1) = Qbar(1,2)
        Qbar(2,2) = Q(1,1)*s**4 + 2*(Q(1,2) + 2*Q(3,3))*s**2*c**2 + Q(2,2)*c**4
        Qbar(1,3) = (Q(1,1) - Q(1,2) - 2*Q(3,3))*c**3 * s + (Q(1,2) - Q(2,2) + 2*Q(3,3))*c*s**3
        Qbar(3,1) = Qbar(1,3)
        Qbar(2,3) = (Q(1,1) - Q(1,2) - 2*Q(3,3))*s**3 * c + (Q(1,2) - Q(2,2) + 2*Q(3,3))*s*c**3
        Qbar(3,2) = Qbar(2,3)
        Qbar(3,3) = (Q(1,1) + Q(2,2) - 2*Q(1,2) - 2*Q(3,3))*c**2 * s**2 + Q(3,3)*(c**4 + s**4)

    end subroutine getQbarMatrix

    subroutine getABDMatrix(Qbar, t, y, ABD)
        use, intrinsic :: iso_fortran_env, only: sp =>real32, dp=>real64
        implicit none

        real(dp), dimension(3,3), intent(in) :: Qbar
        real(dp), intent(in) :: t, y
        real(dp), dimension(6,6), intent(out) :: ABD

        real(dp), dimension(3,3) :: A, B, D
        
        A = Qbar*t
        B = Qbar*t*y
        D = Qbar*(t**3/12 + t*y**2)
        
        !print *, "D = ", D   ! Unnecessary on normal operation, utilized for debugging.
    
        ABD(1:3, 1:3) = A
        ABD(1:3, 4:6) = B
        ABD(4:6, 1:3) = B
        ABD(4:6, 4:6) = D

    end subroutine getABDMatrix

    subroutine getEquivalentMod(ABD, t, Ex, Ey, Gxy, Poisson_xy, Ebx, Eby, Gbxy, Poisson_bxy)
        use, intrinsic :: iso_fortran_env, only: sp =>real32, dp=>real64
        implicit none

        real(dp), dimension(6,6), intent(in) :: ABD
        real(dp), intent(in) :: t
        real(dp), intent(out) :: Ex, Ey, Gxy, Poisson_xy, Ebx, Eby, Gbxy, Poisson_bxy

        real(dp), dimension(3,3) :: A, D

        A = ABD(1:3,1:3)
        D = ABD(4:6,4:6)

        ! Find equivalent moduli for tension and compression. These are not valid for bending, use b moduli instead.
        Ex = (A(1,1)*A(2,2) - A(1,2)**2)/(t*A(2,2))
        Ey = (A(1,1)*A(2,2) - A(1,2)**2)/(t*A(1,1))
        Gxy = A(3,3)/t
        Poisson_xy = A(1,2)/A(2,2)

        ! Find equivalent moduli for bending.
        Ebx = 12*(D(1,1) * D(2,2) - D(1,2)**2)/(t**3 * D(2,2))
        Eby = 12*(D(1,1) * D(2,2) - D(1,2)**2)/(t**3 * D(1,1))
        Gbxy = 12*D(3,3)/t**3
        Poisson_bxy = D(1,2)/D(2,2)

    end subroutine getEquivalentMod

end module 
