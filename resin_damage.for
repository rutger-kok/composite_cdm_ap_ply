! This is a VUMAT subroutine implementing a continuum damage mechanics
! framework for composite materials in Abaqus.
! Copyright (C) 2022 Rutger Wouter Kok

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
! 02110-1301 USA


      subroutine cdm_resin(k,nblock,nstatev,strain,stateOld,C_init,
     1                     e,GMPlus,GMMinus,YT,YC,lch,stress,
     2                     stateNew)
        ! Computes the degraded stress in resin regions using an
        ! isotropic CDM model
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, dimension(6,6), intent(inout) :: C_init  ! initial stiffness matrix
        real*8, dimension(6), intent(in) :: strain  ! strain array
        real*8, intent(in) :: e  ! isotropic stiffness
        real*8, intent(in) :: GMPlus,GMMinus  ! fracture toughnesses
        real*8, intent(in) :: YT,YC  ! strengths
        real*8, intent(in) :: lch  ! characteristic length
        ! local variables
        real*8, dimension(6,6) :: C_deg  ! degraded stiffness matrix
        real*8, dimension(6) :: trialStress  ! trial stress array
        real*8 :: dM,dMState,dMStateOld  ! scalar damage variables
        ! output variables
        real*8, dimension(6), intent(out) :: stress  ! final stress
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array

        dM = 0.0d0  ! initialize matrix damage variable
        trialStress = matmul(C_init, strain)  ! calculate trial stress
        call resin_damage(k,nblock,nstatev,stateOld,stateNew,
     1                    trialStress,YC,YT,lch,e,GMPlus,GMMinus,dM)
        dMStateOld = stateOld(k,50)  ! load old dM variable from state array
        dMState = min(0.95d0, max(dMStateOld,dM))
        stateNew(k,50) = dMState  ! assign new dM variable to state array
        C_deg = C_init * (1.0d0 - dMState)  ! calculate degraded stiffness matrix
        stress = matmul(C_deg, strain)  ! calculate degraded stresses
      end subroutine cdm_resin
      
      
      subroutine resin_damage(k,nblock,nstatev,stateOld,stateNew,
     1                        stress,YC,YT,lch,E_m,GMPlus,GMMinus,dM)
        ! Computes scalar damage variable dM for isotropic material
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, dimension(6), intent(in) :: stress  ! stress array
        real*8, intent(in) :: E_m  ! isotropic stiffness
        real*8, intent(in) :: GMPlus,GMMinus  ! fracture toughnesses
        real*8, intent(in) :: YT,YC  ! strengths
        real*8, intent(in) :: lch  ! characteristic length
        ! local variables
        real*8 :: s_11,s_22,s_33,s_12,s_23,s_13  ! shorthand stresses
        real*8 :: I1,J2  ! invariants of (deviatoric for J2) stress tensor
        real*8 :: L_M  ! characteristic length at damage onset
        real*8 :: FM_T,FM_T_old,FM_T_max,FM_C,FM_C_old,FM_C_max  ! failure criteria
        real*8 :: rMPlus,rMMinus  ! elastic domain thresholds
        real*8 :: AMPlus,AMMinus  ! exponential damage dissipation parameters
        real*8 :: dMPlus,dMMinus  ! tensile/compressive scalar dmaage variables
        ! output variables
        real*8, intent(out) :: dM  ! scalar dmaage variable
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array
        
        ! reassign stresses for brevity
        s_11 = stress(1)
        s_22 = stress(2)
        s_33 = stress(3)
        s_12 = stress(4)
        s_23 = stress(5)
        s_13 = stress(6)

        ! compute first and second ibvariants of deviatoric stress tensor
        I1 = s_11 + s_22 + s_33
        J2 = (1.0d0 / 6.0d0) * (((s_11 - s_22)**2.0d0 +
     1       (s_22 - s_33)**2.0d0 + (s_33 - s_11)**2.0d0))

        ! compute failure indices
        if (I1 >= 0.0d0) then
            FM_T = (3.0d0 * J2 + I1 * (YC - YT)) / (YT * YC)
        else
            FM_C = -((3.0d0 * J2 + I1 * (YC - YT)) / (YT * YC))
        end if

        ! Load old state variables
        FM_T_old = stateOld(k,47)
        FM_C_old = stateOld(k,48)
        ! Compute max failure indices
        FM_T_max = max(FM_T,FM_T_old)
        FM_C_max = max(FM_C,FM_C_old)
        ! Record failure indices in state array
        stateNew(k,47) = FM_T_max
        stateNew(k,48) = FM_C_max

        if ((FM_T_max >= 1.0d0).or.(FM_C_max >= 1.0d0)) then
          ! Check if damage has already initiated
          if (stateOld(k, 49) == 0.0d0) then
            stateNew(k, 49) = lch ! record char. length at onset
            L_M = lch
          else
            ! load fiber direction characteristic length
            L_M = stateOld(k,49)
            stateNew(k,49) = L_M
          end if
          ! calculate linear-exponential softening response
          rMPlus = max(FM_T_max,FM_C_max)  ! elastic domain threshold
          AMPlus = (2.0d0*L_M*YT**2.0d0)/
     1             (2.0d0*E_m*GMPlus - L_M*YT**2.0d0)
          dMPlus = 1.0d0 - (1.0d0/rMPlus)*exp(AMPlus*(1.0d0-rMPlus))
        else
          dMPlus = 0.0d0
        end if

        ! Longitudinal compressive damage
        if (FM_C_max >= 1.0d0) then
          ! Check if damage has already initiated
          if (stateOld(k,49) == 0.0d0) then
            stateNew(k,49) = lch ! record char. length at onset
            L_M = lch
          else
            L_M = stateOld(k,49) ! load char. length
            stateNew(k,49) = L_M
          end if
          rMMinus = FM_C_max
          AMMinus = (2.0d0*L_M*YC**2.0d0)/
     1              (2.0d0*E_m*GMMinus - L_M*YC**2.0d0)
          dMMinus = 1.0d0 - (1.0d0/rMMinus)*exp(AMMinus*(1.0d0-rMMinus))
        else
          dMMinus = 0.0d0
        end if
        
        dM = 1.0d0 - (1.0d0 - dMPlus)*(1.0d0 - dMMinus)

        ! catch NaN dM values
        if (dM /= dM) then
          dM = 0.0d0
        end if

      end subroutine resin_damage