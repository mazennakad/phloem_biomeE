!gfortran -o my_program my_program.f90 -llapack -lblas


!!!! This subroutine solves steady state with constant viscosity

program phloem_biomee
    implicit none

    type(cohort_type),pointer :: cc
    ! Declare constants
    real(8), parameter :: pi_val = 3.141592653589793   ! pi constant
    real(8), parameter :: Rg = 8.3145                  ! ideal gas constant J/(mol K)
    real(8), parameter :: Mw = 342.2965                ! Molecular mass of sucrose g/mol
    real(8), parameter :: rho = 1000                   ! Density of water kg/m3
    real(8), parameter :: g_const = 9.81               ! Gravitational constant m/(s2)
    real(8), parameter :: D = 4.0e-10                  ! Sucrose Diffusivity in water m2/s
    integer, parameter :: n = 100                      ! number of node cells
    real(8), parameter :: T     = 293.15               ! phloem temperature K
    real(8), parameter :: nu0  = 1.5e-3                 ! sap viscosity Pa s
    real(8), parameter :: k   = 5e-14                  ! xylem-phloem membrane permeability m/(Pa s)
    real(8) :: dz = 1/n                                ! grid size
    ! Declare variables
    real(8), allocatable :: c(:)    ! Sucrose concentration at the center of the cell
    real(8), allocatable :: u(:)    ! Longitudinal velocity at the face of the cell
    real(8), allocatable :: p(:)    ! Dynamic fluid pressure at the center of the cell
    real(8), allocatable :: v(:)    ! Transverse velocity at the center of the cell
    real(8), allocatable :: nu(:)    ! Dynamic viscosity at the center of the cell
    real(8), allocatable :: Psi(:)    ! Xylem water potential at the center of the cell
    real(8) :: dx, L  ! diameter at breast height and tree height
    real(8), allocatable :: Resp_s, GResp_s, Growth_s  ! stem demand sucrose g/(ms)
    real(8), allocatable :: Resp_r, GResp_r, Growth_r  ! root demand sucrose g/(ms)
    real(8), allocatable :: Psi_L, psi0             ! leaf water potential and its scaling
    real(8) :: Sstem, Sstemr, Sstemg, Sroot, Sleaf
    ! scaling parameters
    real(8) :: c0, cw, os, es, v0, u0, p0, t0, G, X0, Mu, Pe, Ss, Sl, Sr


    L = cc%height
    dx = cc%dbh
    a  = cc%p_thicnkness
    Resp_s =   cc%DR_stem*Mw/(12*mol_C)
    Resp_r =   cc%DR_root*Mw/(12*mol_C)
    GResp_s =   cc%DGR_stem*Mw/(12*mol_C)
    GResp_r =   cc%DGR_root*Mw/(12*mol_C)
    Growth_s =   cc%DG_stem*Mw/(12*mol_C)
    Growth_r =   cc%DG_root*Mw/(12*mol_C)
    Sstem = (Resp_s + GResp_s + Growth_s)/(3600*L*pi*dx)
    Sroot = (Resp_r + GResp_r + Growth_r)/(3600*a*pi*dx)
    Sleaf = (Sstem*L + Sroot*a)/a
    Psi_L    = cc%Psi_L
    psi0    = - Psi_L

    ! Non-dimensional and scaling quantities
    c0 = psi0 * Mw * pi_val * dx / (Rg * T)
    cw = c0 / (Mw * pi_val * dx)
    os = Rg * T * c0 / (Mw * pi_val * dx)
    es = a / L
    v0 = k * os
    u0 = v0 / es
    p0 = L * nu0 * u0 / (a**2)
    t0 = (a**2) / D
    G = rho * g_const * L / os
    X0 = psi0 / os
    Mu = k * nu0 * (L**2) / (a**3)
    Pe = a * v0 / D
    Ss = Sstem / (u0 * c0)
    Sl = Sleaf / (u0 * c0)
    Sr = Sroot / (u0 * c0)

    ! Allocate variables to be solved
    allocate(c(n),u(n-1),p(n),v(n),nu(n),Psi(n))
    ! Initial condition variables
    allocate(co(n), uo(n-1), vo(n), po(n),nuo(n))

    do i=1,n-1
      co(i) = cc%c(i)
      uo(i) = cc%u(i)
      vo(i) = cc%v(i)
      po(i) = cc%p(i)
    end do
    co(n) = cc%c(n)
    vo(n) = cc%v(n)
    po(n) = cc%p(n)

    ! xylem water potential
    call linspace(data%Psi_L(j) - rho*g_const*L, data%Psi_W(j), n, Psi)
    Psi = Psi / psi0

    call iterate(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, co, uo, vo, po, nuo, cw, &
     c, u, v, p, nu)

    ! Deallocate arrays
    deallocate(c, u, p, v, nu, co, uo, vo, po, nuo,Psi)

end program phloem_biomee

subroutine linspace(start, end_val, n, array)
    integer, intent(in) :: n
    real(8), intent(in) :: start, end_val
    real(8), intent(out) :: array(n)
    integer :: i
    do i = 1, n
        array(i) = start + (end_val - start) * (i - 1) / real(n - 1, 8)
    end do
end subroutine linspace

subroutine iterate(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, co, uo, vo, po, nuo, &
                  cw, ck, uk, vk, pk, nuk)

    ! Input parameters
    integer, intent(in) :: n
    real(8), intent(in) :: dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, cw
    real(8), dimension(n), intent(in) :: co, uo, vo, po, nuo
    real(8), dimension(n-1), intent(in) :: uo

    ! Output parameters
    real(8), dimension(n), intent(out) :: ck, uk, vk, pk, nuk
    real(8), dimension(n-1), intent(out) :: uk


    ! Local variables
    logical :: check_C
    real(8), allocatable :: vector_S(:), Sk(:)
    integer :: vector_length
    integer :: iteration

    ! Initialize variables
    check_C = .true.
    iteration = 1

    ! Allocate arrays for newton results
    vector_length = 5*n - 1
    allocate(vector_S(vector_length))
    allocate(Sk(vector_length))

    ! Main iteration loop
    do while (check_C)
        ! Set the new variables
        if (iteration /= 1) then
            co = ck
            uo = uk
            vo = vk
            po = pk
            nuo = nuk
        end if

        ! Call newton subroutine
        call newton(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, &
                   co, uo, vo, po, nuo, cw, vector_S, Sk)

        ! Check convergence
        call check_convergence(vector_S, Sk, vector_length, check_C)

        ! Update solution vectors
        pk = Sk(1:n)
        uk = Sk(n+1:2*n-1)
        vk = Sk(2*n:3*n-1)
        ck = Sk(3*n:4*n-1)
        nuk = Sk(4*n:5*n-1)

        ! Check maximum iterations
        if (iteration > 50) then
            exit
        end if

        iteration = iteration + 1
    end do

    ! Deallocate temporary arrays
    deallocate(vector_S)
    deallocate(Sk)

end subroutine iterate

subroutine check_convergence(vector_S, Sk, size, check_C)
    implicit none
    ! Input parameters
    integer, intent(in) :: size
    real(8), intent(in), dimension(size) :: vector_S, Sk
    ! Output parameter
    logical, intent(out) :: check_C
    ! Local variable
    real(8) :: rms

    ! Compute RMS error
    rms = sqrt(sum((vector_S - Sk)**2) / dble(size))

    ! Check convergence
    if (rms <= 1.0d-4) then
        check_C = .false.
    else
        check_C = .true.
    end if

end subroutine check_convergence

subroutine newton(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, &
    c, u, v, p, nu, cw, vector_S, Sk)


    ! Input parameters
    integer, intent(in) :: n
    real, intent(in) :: dz
    real, intent(in) :: Mu, G, X0, Pe, Sl, Ss, Sr, es, cw
    real, intent(in) :: c(n), u(n-1), v(n), p(n), nu(n), Psi(n)

    ! Output parameters
    real, allocatable, intent(out) :: vector_S(:), Sk(:)

    ! Internal variables
    real, allocatable :: PP(:,:), PU(:,:), PV(:,:), PC(:,:), PN(:,:), F1(:)
    real, allocatable :: UP(:,:), UU(:,:), UV(:,:), UC(:,:), UN(:,:), F2(:)
    real, allocatable :: VP(:,:), VU(:,:), VV(:,:), VC(:,:), VN(:,:), F3(:)
    real, allocatable :: CP(:,:), CU(:,:), CV(:,:), CC(:,:), CN(:,:), F4(:)
    real, allocatable :: NP(:,:), NU(:,:), NV(:,:), NC(:,:), NN(:,:), F5(:)
    real, allocatable :: m(:,:), F(:), sol(:), Fk(:)
    real(8) :: resk
    integer :: i

    ! Call BuildP subroutine
    call BuildP(n, dz, Mu, G, X0, Psi, c, p, nu, PP, PU, PV, PC, PN, F1)

    ! Call BuildU subroutine
    call BuildU(n, dz, u, p, nu, UP, UU, UV, UC, UN, F2)

    ! Call BuildV subroutine
    call BuildV(n, dz, v, p, nu, c, G, X0, Psi, VP, VU, VV, VC, VN, F3)

    ! Call BuildC subroutine
    call BuildC(n, dz, Pe, Sl, Ss, Sr, es, c, u, v, &
        CP, CU, CV, CC, CN, F4)

    allocate(NP(n), NU(n, n-1),NV(n),NC(n),NN(n),F5(n))
    NP = 0.0d0
    NU = 0.0d0
    NV = 0.0d0
    NC = 0.0d0
    NN = 0.0d0
    do i = 1, n
        NN(i) = 1.0
    end do
    F5 = 0.0d0


    ! Construct matrix m
    allocate(m(5*n, 5*n))
    m = 0.0d0
    ! Populate m with sparse blocks (using dense representation)
    ! Upper blocks
    m(1:n,1:n)       = PP
    m(1:n,n+1:2*n-1) = PU
    m(1:n,2*n:3*n-1) = PV
    m(1:n,3*n:4*n-1) = PC
    m(1:n,4*n+1:5*n) = PN

    m(n+1:2*n-1,1:n)       = UP
    m(n+1:2*n-1,n+1:2*n-1) = UU
    m(n+1:2*n-1,2*n:3*n-1) = UV
    m(n+1:2*n-1,3*n:4*n-1) = UC
    m(n+1:2*n-1,4*n+1:5*n) = UN

    m(2*n:3*n-1,1:n)       = VP
    m(2*n:3*n-1,n+1:2*n-1) = VU
    m(2*n:3*n-1,2*n:3*n-1) = VV
    m(2*n:3*n-1,3*n:4*n-1) = VC
    m(2*n:3*n-1,4*n+1:5*n) = VN

    m(3*n:4*n-1,1:n)       = CP
    m(3*n:4*n-1,n+1:2*n-1) = CU
    m(3*n:4*n-1,2*n:3*n-1) = CV
    m(3*n:4*n-1,3*n:4*n-1) = CC
    m(3*n:4*n-1,4*n+1:5*n) = CN

    m(4*n+1:5*n,1:n)       = NP
    m(4*n+1:5*n,n+1:2*n-1) = NU
    m(4*n+1:5*n,2*n:3*n-1) = NV
    m(4*n+1:5*n,3*n:4*n-1) = NC
    m(4*n+1:5*n,4*n+1:5*n) = NN

    ! Construct F vector
    allocate(F(5*n))
    F = 0.0d0
    F(1:n)         = -F1
    F(n+1:2*n-1)   = -F2
    F(2*n:3*n-1)   = -F3
    F(3*n:4*n-1)   = -F4
    F(4*n+1:5*n)   = -F5

    ! Solve the linear system: (m + 1e-6 * I) * sol = F
    allocate(sol(5*n))
    ! Add regularization term to the diagonal
    do i = 1, 5*n
        m(i,i) = m(i,i) + 1e-6
    end do
    ! Call a linear solver (e.g., LAPACK) here. Placeholder:
    call linear_solver(m, F, sol)

    ! Construct vector_S by concatenating p, u, v, c, nu
    allocate(vector_S(5*n))
    vector_S(1:n)          = p
    vector_S(n+1:2*n-1)    = u
    vector_S(2*n:3*n-1)    = v
    vector_S(3*n:4*n-1)    = c
    vector_S(4*n+1:5*n)    = nu

    ! Compute Sk = vector_S + sol
    allocate(Sk(5*n))
    Sk = vector_S
    Sk = Sk + sol

    ! Extract updated variables from Sk
    ! pk  = Sk(1:n)
    ! uk  = Sk(n+1:2*n-1)
    ! vk  = Sk(2*n:3*n-1)
    ! ck  = Sk(3*n:4*n-1)
    ! nuk = Sk(4*n+1:5*n)

!    ! Recompute Fk with updated variables
!    call BuildP(n, dz, Mu, G, X0, Psi, Sk(1:n), Sk(n+1:2*n-1), Sk(4*n+1:5*n), F1)
!    call BuildU(n, dz, Sk(n+1:2*n-1), Sk(1:n), Sk(4*n+1:5*n), F2)
!    call BuildV(n, dz, Sk(2*n:3*n-1), Sk(1:n), Sk(4*n+1:5*n), Sk(3*n:4*n-1), G, X0, Psi, F3)
!    call BuildC(n, dz, Pe, Sl, Ss, Sr, es, Sk(3*n:4*n-1), Sk(n+1:2*n-1), &
!               Sk(2*n:3*n-1), F4)

!    ! Construct Fk vector
!    allocate(Fk(5*n))
!    Fk = 0.0
!    Fk(1:n)         = F1
!    Fk(n+1:2*n-1)   = F2
!    Fk(2*n:3*n-1)   = F3
!    Fk(3*n:4*n-1)   = F4
!    Fk(4*n+1:5*n)   = F5
!
!    ! Compute the norm of Fk
!    resk = 0.0
!    do i = 1, 5*n
!        resk = resk + Fk(i)**2
!    end do
!    resk = sqrt(resk)

    ! Deallocate allocated arrays
    deallocate(PP, PU, PV, PC, PN, F1)
    deallocate(UP, UU, UV, UC, UN, F2)
    deallocate(VP, VU, VV, VC, VN, F3)
    deallocate(CP, CU, CV, CC, CN, F4)
    deallocate(NP, NU, NV, NC, NN, F5)
    deallocate(m, F, sol, Fk, vector_S, Sk)

end subroutine newton

subroutine linear_solver(A, b, x)

    ! Input parameters
    real, intent(in) :: A(:,:), b(:)
    real, intent(out) :: x(size(b))

    ! Local variables
    integer :: n, lda, info
    integer, allocatable :: ipiv(:)

    n = size(b)
    lda = n

    ! Allocate pivot array
    allocate(ipiv(n))

    ! Call LAPACK routine to solve Ax = b
    call dgesv(n, 1, A, lda, ipiv, b, n, info)

    ! Check for successful exit
    if (info /= 0) then
        print *, "Error: The linear solver failed with info =", info
        stop
    end if

    ! Copy the solution to x
    x = b

    ! Deallocate pivot array
    deallocate(ipiv)

end subroutine linear_solver



subroutine BuildP(n, dz, Mu, G, X0, Psi, c, p, nu, PP, PU, PV, PC, PN, F1)

    ! Input arguments
    integer, intent(in) :: n
    real(8), intent(in) :: dz, Mu, G, X0
    real(8), dimension(n), intent(in) :: Psi, c, p, nu

    ! Output arguments
    real(8), dimension(n, n), intent(out) :: PP, PV, PC, PN
    real(8), dimension(n, n-1), intent(out) :: PU
    real(8), dimension(n), intent(out) :: F1

    ! Local variables
    real(8), dimension(n) :: diagP, diagN
    real(8), dimension(n-1) :: uppP, lowP, uppN, lowN
    integer :: i

    ! Initialize matrices to zero
    PP = 0.0d0
    PU = 0.0d0
    PV = 0.0d0
    PC = 0.0d0
    PN = 0.0d0

    ! -----------------------------------------------------
    ! P variable (PP matrix)
    ! -----------------------------------------------------
    diagP = 0.0d0
    uppP = 0.0d0
    lowP = 0.0d0
    diagP(1) = -Mu
    diagP(n) = -Mu
    diagN = 0.0d0
    uppN = 0.0d0
    lowN = 0.0d0
    F1 = 0.0d0
    do i = 2, n-1
      diagP(i) = -nu(i) * 2.0d0 / (3.0d0 * dz**2) - Mu
      uppP(i) = (nu(i+1) - nu(i-1)) / (12.0d0 * dz**2) + nu(i) / (3.0d0 * dz**2)
      lowP(i-1) = -(nu(i+1) - nu(i-1)) / (12.0d0 * dz**2) + nu(i) / (3.0d0 * dz**2)
      diagN(i) = (p(i+1) - 2.0d0 * p(i) + p(i-1)) / (3.0d0 * dz**2)
      uppN(i) = (p(i+1) - p(i-1)) / (12.0d0 * dz**2)
      lowN(i-1) = -uppN(i)
      F1(i) = (nu(i+1) - nu(i-1)) * (p(i+1) - p(i-1)) / (6.0d0 * dz**2) + &
              nu(i) * (p(i+1) - 2.0d0 * p(i) + p(i-1)) / (3.0d0 * dz**2) - &
              Mu * p(i) + c(i) + X0 * Psi(i) - G * dz * (i - 0.5d0)
    end do
    F1(1) = -Mu * p(1) + c(1) + X0 * Psi(1) - G * dz * 0.5d0
    F1(n) = -Mu * p(n) + c(n) + X0 * Psi(n) - G * dz * (n - 0.5d0)


    do i = 1, n
        PC(i, i) = 1.0d0
        PN(i, i) = diagN(i)
        PP(i, i) = diagP(i)
        if (i < n) then
          PN(i, i+1) = uppN(i)
          PP(i, i+1) = uppP(i)
        end if
        if (i > 1) then
          PN(i, i-1) = lowN(i-1)
          PP(i, i-1) = lowP(i-1)
        end if
    end do


end subroutine BuildP

subroutine BuildU(n, dz, u, p, nu, UP, UU, UV, UC, UN, F2)

    ! Input arguments
    integer, intent(in) :: n
    real(8), intent(in) :: dz
    real(8), dimension(n-1), intent(in) :: u
    real(8), dimension(n), intent(in) :: p, nu

    ! Output arguments
    real(8), dimension(n-1, n), intent(out) :: UP, UV, UC, UN
    real(8), dimension(n-1, n-1), intent(out) :: UU
    real(8), dimension(n-1), intent(out) :: F2

    ! Local variables
    integer :: i

    UV = 0.0d0
    UC = 0.0d0
    UP = 0.0d0
    UN = 0.0d0
    UU = 0.0d0
    F2 = 0.0d0
    do i = 1, n-1
        UP(i, i) = -(nu(i+1) + nu(i)) / (6.0d0 * dz)
        UN(i, i) = (p(i+1) - p(i)) / (6.0d0 * dz)
        UU(i, i) = 1.0d0
        F2(i) = u(i) + (nu(i+1) + nu(i)) * (p(i+1) - p(i)) / (6.0d0 * dz)
        if (i < n-1) then
            UP(i, i+1) = (nu(i+1) + nu(i)) / (6.0d0 * dz)
            UN(i, i+1) = (p(i+1) - p(i)) / (6.0d0 * dz)
        end if
    end do
    UP(n-1, n) = (nu(n) + nu(n-1)) / (6.0d0 * dz)
    UN(n-1, n) = (p(n) - p(n-1)) / (6.0d0 * dz)

end subroutine BuildU

subroutine BuildV(n, dz, v, p, nu, VP, VU, VV, VC, VN, F3)

    ! Input arguments
    integer, intent(in) :: n
    real(8), intent(in) :: dz
    real(8), dimension(n), intent(in) :: v, p, nu,

    ! Output arguments
    real(8), dimension(n, n), intent(out) :: VP, VV, VC, VN
    real(8), dimension(n, n-1), intent(out) :: VU
    real(8), dimension(n), intent(out) :: F3

    ! Local variables
    real(8), dimension(n) :: diagP, diagN
    real(8), dimension(n-1) :: uppP, lowP, uppN, lowN
    integer :: i


    diagP = 0.0d0
    uppP = 0.0d0
    lowP = 0.0d0
    diagN = 0.0d0
    uppN = 0.0d0
    lowN = 0.0d0
    F3 = 0.0d0
    do i = 2, n-1
        diagP(i) = -nu(i) / (4.0d0 * dz**2)
        uppP(i) = (nu(i+1) - nu(i-1)) / (32.0d0 * dz**2) + nu(i) / (8.0d0 * dz**2)
        lowP(i-1) = -(nu(i+1) - nu(i-1)) / (32.0d0 * dz**2) + nu(i) / (8.0d0 * dz**2)
        diagN(i) = (p(i+1) - 2.0d0 * p(i) + p(i-1)) / (8.0d0 * dz**2)
        uppN(i) = (p(i+1) - p(i-1)) / (32.0d0 * dz**2)
        lowN(i-1) = -uppN(i)
        F3(i) = v(i) + (nu(i+1) - nu(i-1)) * (p(i+1) - p(i-1)) / (32.0d0 * dz**2) + &
                nu(i) * (p(i+1) - 2.0d0 * p(i) + p(i-1)) / (8.0d0 * dz**2)
    end do

    VC = 0.0d0
    VU = 0.0d0
    VP = 0.0d0
    VV = 0.0d0
    VN = 0.0d0
    do i = 1, n
        VP(i, i) = diagP(i)
        VV(i, i) = 1.0d0
        VN(i, i) = diagN(i)
        if (i < n) then
          VP(i, i+1) = uppP(i)
          VN(i, i+1) = uppN(i)
        end if
        if (i > 1) then
          VP(i, i-1) = lowP(i-1)
          VN(i, i-1) = lowN(i-1)
        end if
    end do


end subroutine BuildV


subroutine BuildC(n, dz, Pe, Sl, Ss, Sr, es, c, u, v, CP, CU, CV, CC, CN, F4)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: dz, Pe, Sl, Ss, Sr, es
    real, intent(in) :: c(n), u(n-1), v(n)
    real, allocatable, intent(out) :: CP(:,:), CU(:,:), CV(:,:), CC(:,:), CN(:,:), F4(:)

    ! Internal variables
    integer :: i
    real, allocatable :: a1(:), a2(:), a3(:), a4(:), diagC(:), uppC(:), lowC(:)

    ! Allocate matrices
    allocate(CP(n, n), CU(n, n-1), CV(n, n), CC(n, n), CN(n, n), F4(n))

    ! Initialize matrices to zero
    CP = 0.0d0
    CN = 0.0d0

    ! Allocate helper arrays
    allocate(a1(n-1), a2(n-2), a3(n-1), a4(n-2), diagC(n), uppC(n-1), lowC(n-1))

    ! Construct CU matrix (U variable)
    do i = 1, n-1
        a1(i) = ( Pe - (7.0/45.0)*(Pe**2)*v(i) )*c(i)/dz &
              - (4.0/105.0)*( (Pe/dz)**2 )* u(i) * (c(i+1) - c(i)) &
              - 2.0*(7.0/120.0)*(Ss/es)*(Pe**2)*dz*(real(i))/dz
        a3(i) = - (7.0/45.0)*(Pe**2) * u(i) * c(i) / dz
    end do

    do i = 1, n-2
        a2(i) = - ( Pe - (7.0/45.0)*(Pe**2)*v(i) )*c(i)/dz &
               + (4.0/105.0)*( (Pe/dz)**2 )* u(i) * (c(i+1) - c(i)) &
               + 2.0*(7.0/120.0)*(Ss/es)*(Pe**2)*dz*(real(i))/dz
        a4(i) = (7.0/45.0)*(Pe**2) * u(i) * c(i) / dz
    end do

    CU = 0.0d0
    CV = 0.0d0
    do i = 1, n-1
        CU(i, i) = a1(i)
        CV(i, i) = a3(i)
        if (i > 1) then
          CU(i, i-1) = a2(i-1)
          CV(i, i-1) = a4(i-1)
        end if
    end do
    CU(n, n-1) = - ( Pe - (7.0/45.0)*(Pe**2)*v(n-1) )*c(n-1)/dz &
                 + (4.0/10.5.0)*( (Pe/dz)**2 )*u(n-1)*( c(n) - c(n-1) ) &
                 + 2.0*(7.0/120.0)*(Ss/es)*(Pe**2)*dz*(n-1)/dz

    CV(n, n-1) = (7.0/45.0)*(Pe**2) * u(n-1) * c(n-1) / dz
    CV(1,1) = 0

    allocate(F4(n))
    F4 = 0.0d0

    do i = 2, n-1
        diagC(i) = ( Pe - (7.0/45.0)*(Pe**2)*v(i) )*u(i) / dz &
                 + ((2.0/105.0)*(Pe**2)*(u(i)**2) + es**2) / (dz**2) &
                 + ((2.0/105.0)*(Pe**2)*(u(i-1)**2) + es**2) / (dz**2)
        F4(i) = (Pe - (7.0/45.0)*(Pe**2)*v(i))*u(i)*c(i)/dz &
              - ((2.0/105.0)*(Pe**2)*(u(i)**2) + es**2) * (c(i+1) - c(i)) / (dz**2) &
              + 2.0 * ((Pe*(dz**2)/(2.0*es))*(real(i)**2) - (7.0/120.0)*(Pe**2)*dz*u(i)*real(i)/es) * Ss / dz &
              - (Pe - (7.0/45.0)*(Pe**2)*v(i-1))*u(i-1)*c(i-1)/dz &
              + ((2.0/105.0)*(Pe**2)*(u(i-1)**2) + es**2) * (c(i) - c(i-1)) / (dz**2) &
              - 2.0 * ((Pe*(dz**2)/(2.0*es))*(real(i-1)**2) - (7.0/120.0)*(Pe**2)*dz*u(i-1)*real(i-1)/es) * Ss / dz
    end do

    diagC(1) = Pe*u(1)/dz + ( (2.0/105.0)*(Pe**2)*(u(1)**2) + es**2 ) / (dz**2)
    diagC(n) = ((2.0/105.0)*(Pe**2)*(u(n-1)**2) + es**2) / (dz**2)

    ! Boundary conditions for F4
    F4(1) = (Pe - (7.0/45.0)*(Pe**2)*v(1))*u(1)*c(1)/dz &
          - ((2.0/105.0)*(Pe**2)*(u(1)**2) + es**2) * (c(2) - c(1)) / (dz**2) &
          + 2.0 * ((Pe*(dz**2)/(2.0*es)) - (7.0/120.0)*(Pe**2)*dz*u(1)/es) * Ss / dz &
          - Pe*Sl/dz

    F4(n) = (Pe*Sr + (Pe/es)*Ss)/dz &
          - (Pe - (7.0/45.0)*(Pe**2)*v(n-1))*u(n-1)*c(n-1)/dz &
          + ((2.0/105.0)*(Pe**2)*(u(n-1)**2) + es**2) * (c(n) - c(n-1)) / (dz**2) &
          - 2.0 * ((Pe*(dz**2)/(2.0*es))*(real(n-1)**2) - (7.0/120.0)*(Pe**2)*dz*u(n-1)*real(n-1)/es) * Ss / dz

    do i = 1, n-1
        uppC(i) = - ((2.0/105.0)*(Pe**2)*(u(i)**2) + es**2) / (dz**2)
        lowC(i) = - ( Pe - (7.0/45.0)*(Pe**2)*v(i) )*u(i)/dz &
                - ((2.0/105.0)*(Pe**2)*(u(i)**2) + es**2) / (dz**2)
    end do

    CC = 0.0
    do i = 1, n
        CC(i, i) = diagC(i)
        if (i>1) CC(i,i-1) = lowC(i-1)
        if (i<n) CC(i,i+1) = uppC(i)
    end do


    ! Deallocate helper arrays
    deallocate(a1, a2, a3, a4, diagC, uppC, lowC)

end subroutine BuildC


!!!!!!!!!!!!!!!!!!! for variable viscosity

!subroutine BuildN(n, cw, bb, c, nu, NP, NU, NV, NC, NN, F5)
!
!    ! Input arguments
!    integer, intent(in) :: n
!    real(8), intent(in) :: cw
!    real(8), dimension(5), intent(in) :: bb
!    real(8), dimension(n), intent(in) :: c, nu
!
!    ! Output arguments
!    real(8), dimension(n, n), intent(out) :: NP, NV, NC, NN
!    real(8), dimension(n, n-1), intent(out) :: NU
!    real(8), dimension(n), intent(out) :: F5
!
!    ! Local variables
!    real(8), dimension(n) :: n1
!    integer :: i
!
!    ! Initialize NP, NU, NV, and NN to zeros
!    NP = 0.0d0
!    NV = 0.0d0
!    NU = 0.0d0
!    NN = 0.0d0
!
!    ! Calculate diagonal elements of NC
!    n1 = ( bb(2) * cw - 2.0d0 * bb(3) * (cw**2) * c - 3.0d0 * bb(4) * (cw**3) * c**2 - &
!           4.0d0 * bb(5) * (cw**4) * c**3 ) * &
!         exp( - bb(2) * cw * (c - 1.0d0) - bb(3) * (cw**2) * (c**2 - 1.0d0) - &
!              bb(4) * (cw**3) * (c**3 - 1.0d0) - bb(5) * (cw**4) * (c**4 - 1.0d0) )
!
!    NC = 0.0d0
!    do i = 1, n
!        NC(i, i) = n1(i)
!!    end do
!
!    ! Set NN as the identity matrix
!    do i = 1, n
!        NN(i, i) = 1.0d0
!    end do
!
!    ! Calculate F5
!    F5 = nu - exp( - bb(2) * cw * (c - 1.0d0) - bb(3) * (cw**2) * (c**2 - 1.0d0) - &
!                   bb(4) * (cw**3) * (c**3 - 1.0d0) - bb(5) * (cw**4) * (c**4 - 1.0d0) )
!
!end subroutine BuildN





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for the initial state



subroutine initial(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, cw, c, u, nu, v, p)

    ! Input parameters
    integer, intent(in) :: n, m, nts
    real(8), intent(in) :: dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi(n), cw

    ! Output arrays
    real(8), dimension(n), intent(out) :: c, nu, v, p
    real(8), dimension(n-1), intent(out) :: u

    ! Local variables
    real(8), dimension(n) :: co, vo, po, nuo
    real(8), dimension(n-1) :: uo

    call guess(n, dz, G, Mu, X0, Psi, cw, co, uo, vo, po, nuo)
    call iterate(n, dz, Mu, G, X0, Pe, Sl, Ss, Sr, es, Psi, &
    co, uo, vo, po, nuo, cw, c, u, v, p, nu)


end subroutine initial

subroutine guess(n, dz, G, Mu, X0, Psi, cw, ci, ui, vi, pi, nui)

    ! Input parameters
    integer, intent(in) :: n
    real(8), intent(in) :: dz, G, Mu, X0, Psi(n), cw

    ! Output arrays
    real(8), dimension(n), intent(out) :: co, vo, po, nuo
    real(8), dimension(n-1), intent(out) :: uo


    integer :: i

    ! Calculate concentration
    do i = 1, n
        co(i) = 3.0d0 * (G * dz * (dble(i) - 0.5d0) - X0 * Psi(i))
    end do

    nuo = 1.0d0

    ! Calculate velocity components
    call VelocityS(co, nuo, Mu, Psi, X0, G, dz, n, uo, vo, po)

end subroutine guess

subroutine VelocityS(co, nuo, Mu, Psi, X0, G, dz, n, uo, vo, po)
    ! Arguments
    integer, intent(in) :: n
    real(8), intent(in) :: co(n), nuo(n), Mu, Psi(n), X0, G, dz
    real(8), intent(out) :: uo(n-1), vo(n), vo(n)

    ! Local variables
    real(8) :: diagP(n), uppP(n), lowP(n), eqP(n)
    integer :: i

    ! Initialize arrays
    diagP = 0.0d0
    uppP = 0.0d0
    lowP = 0.0d0
    eqP = 0.0d0

    diagP(1) = -Mu
    diagP(n) = -Mu
    eqP(1) = -co(1) - X0*Psi(1) + G*dz*0.5d0
    eqP(n) = -co(n) - X0*Psi(n) + G*dz*(dble(n) - 0.5d0)
    do i = 2, n-1
        diagP(i) = -nuo(i)*2.0d0/(3.0d0*(dz**2)) - Mu
        uppP(i) = (nuo(i+1) - nuo(i-1))/(12.0d0*(dz**2)) + nuo(i)/(3.0d0*(dz**2))
        lowP(i) = -(nuo(i+1) - nuo(i-1))/(12.0d0*(dz**2)) + nuo(i)/(3.0d0*(dz**2))
        eqP(i) = -co(i) - X0*Psi(i) + G*dz*(dble(i) - 0.5d0)
    end do


    ! Solve tridiagonal system
    call Thomas(lowP, diagP, uppP, eqP, po, n)

    ! Calculate axial velocity
    do i = 1, n-1
        uo(i) = (nuo(i+1) + nuo(i))*(po(i) - po(i+1))/(6.0d0*dz)
    end do

    ! Calculate radial velocity
    vo = 0.0d0
    do i = 2, n-1
        vo(i) = (3.0d0/8.0d0)*dz*(uo(i) - uo(i-1))
    end do

end subroutine VelocityS

subroutine thomas(aa, bb, cc, dd, q, n)
    ! Arguments
    integer, intent(in) :: n
    real(8), intent(in) :: aa(n), bb(n), cc(n), dd(n)
    real(8), intent(out) :: q(n)

    ! Local variables
    real(8) :: bet(n), gam(n)
    integer :: i

    ! Initialize first elements
    bet(1) = bb(1)
    gam(1) = dd(1)/bb(1)

    ! Forward elimination
    do i = 2, n
        bet(i) = bb(i) - (aa(i)*cc(i-1)/bet(i-1))
        gam(i) = (dd(i) - aa(i)*gam(i-1))/bet(i)
    end do

    ! Back substitution
    q(n) = gam(n)

    do i = n-1, 1, -1
        q(i) = gam(i) - (cc(i)*q(i+1)/bet(i))
    end do

end subroutine thomas





end program phloem_biomee
