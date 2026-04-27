

!==============================================================
! BuP module — d(z) et G_eff pour CAMB
! À inclure dans equations.f90 ou sources.f90
!
! Auteur : BuP Paper 3
! Usage  : ajouter "use bup_geff" en tête de module CAMB
!==============================================================

module bup_geff
    implicit none

    ! Constantes BuP
    real(dl), parameter :: DC        = 3.059842935509462_dl
    real(dl), parameter :: AC        = 0.7252150458197096_dl
    real(dl), parameter :: BETA_BUP  = 1.4134453781512604_dl
    real(dl), parameter :: ALPHA_MIC = 1.78_dl

    ! Paramètres Paper 2 best-fit (à mettre à jour si besoin)
    real(dl) :: bup_X0   = 0.537_dl
    real(dl) :: bup_beta = 2.000_dl
    real(dl) :: bup_dd   = 0.878_dl

contains

    !----------------------------------------------------------
    ! Dimension effective d(z)
    !----------------------------------------------------------
    function bup_d_of_z(z) result(d)
        real(dl), intent(in) :: z
        real(dl) :: d, X

        X = bup_X0 * (1.0_dl + z)**bup_beta
        d = DC - bup_dd / (1.0_dl + X**ALPHA_MIC)
        ! Clamp physique
        d = max(2.0_dl, min(d, DC))

    end function bup_d_of_z

    !----------------------------------------------------------
    ! G_eff(z) / G = 2/(d(z)-1)
    !----------------------------------------------------------
    function bup_geff(z) result(geff)
        real(dl), intent(in) :: z
        real(dl) :: geff, d

        d    = bup_d_of_z(z)
        geff = 2.0_dl / (d - 1.0_dl)

    end function bup_geff

    !----------------------------------------------------------
    ! G_eff avec coupure d'échelle (Version B, Paper 3)
    ! geff(z,k) = geff(z) / (1 + (k/k_star)^2)
    !----------------------------------------------------------
    function bup_geff_scale(z, k, k_star) result(geff)
        real(dl), intent(in) :: z, k, k_star
        real(dl) :: geff

        geff = bup_geff(z) / (1.0_dl + (k/k_star)**2)

    end function bup_geff_scale

end module bup_geff
