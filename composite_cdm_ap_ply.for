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

      include 'resin_damage.for'  ! import pure resin cdm model
      include 'rotation_matrix.for'  ! import rotation matrix function
      
      subroutine vumat(
        ! Read only - input variables from Abaqus
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
        ! Write only - outputs Abaqus needs from subroutine
     7  stressNew, stateNew, enerInternNew, enerInelasNew)

        include 'vaba_param.inc'
        
        ! Dimension the Abaqus input variables
        dimension props(nprops), density(nblock), coordMp(nblock,*),
     1    charLength(nblock), strainInc(nblock, ndir+nshr),
     2    relSpinInc(nblock, nshr), tempOld(nblock),
     3    stretchOld(nblock, ndir+nshr), fieldOld(nblock, nfieldv),
     4    defgradOld(nblock,ndir+nshr+nshr), stateOld(nblock, nstatev),
     5    stressOld(nblock, ndir+nshr), enerInternOld(nblock),
     6    enerInelasOld(nblock), tempNew(nblock),
     7    stretchNew(nblock, ndir+nshr), fieldNew(nblock, nfieldv),
     8    defgradNew(nblock,ndir+nshr+nshr), enerInelasNew(nblock),
     9    stressNew(nblock,ndir+nshr), stateNew(nblock, nstatev),
     1    enerInternNew(nblock)     
        character*80 cmname

        ! Declare variables used in main part of script
        integer :: k,i  ! loop counters
        real*8, dimension(6) :: stress_total,strain  ! total stress & strain in element
        real*8, dimension(6) :: stress_1,stress_2  ! stress in element constituents
        real*8, dimension(6,6) :: C_tow,C_resin  ! iniital stiffness matrices tow and resin
        real*8, dimension(6,6) :: C_1_global,C_2_global  ! stiffness matrices constituents
        real*8, dimension(6,6) :: T1,T2  ! stress/strain transformation matrices
        real*8 :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23 ! elastic constants
        real*8 :: nu21,nu31,nu32  ! secondary poisson's ratios
        real*8 :: XT,XC,YT,YC,ZT,ZC,SL,SR,ST  ! strengths
        real*8 :: G1Plus,G1Minus,G2Plus,G2Minus,G6  ! fracture toughnesses 
        real*8 :: int_angle  ! interface angle between constituents
        real*8 :: alpha_mean  ! mean out-of-plane angle
        real*8 :: lch,lch_1,lch_2  ! characteristic lengths of element and consituents
        real*8 :: vfrac_1,vfrac_2  ! volume fractions of constituents
        real*8 :: psi,psi_resin  ! stiffness matrix calculation terms
        real*8 :: und_flag,und_angle_1,und_angle_2  ! undulation flag and in-plane angles of constituents
        real*8 :: cpt  ! cured ply thickness
        real*8 :: u_length  ! undulation length
        real*8 :: pi  ! pi constant
        real*8 :: stress_power  ! element internal energy
        real*8 :: F  ! determinant of deformation gradient
        
        ! declare external functions
        real*8, external :: det_3x3

        ! Load elastic constants orthotropic ply
        e11 = props(1) ! stiffness fiber 
        e22 = props(2) ! stiffness transverse  (in-plane)
        e33 = props(3) ! stiffness transverse  (out-of-plane)
        nu12 = props(4) ! Poisson's ratio 12 
        nu13 = props(5) ! Poisson's ratio 13 
        nu23 = props(6) ! Poisson's ratio 23 
        ! shear moduli multiplied by two (tensorial-engineering strain)
        g12 = 2.0d0*props(7) ! shear modulus 12  
        g13 = 2.0d0*props(8) ! shear modulus 13 
        g23 = 2.0d0*props(9) ! shear modulus 23 

        ! Ply strengths
        XT = props(10) ! tensile strength fiber 
        XC = props(11) ! compressive strength fiber 
        YT = props(12) ! in-situ tensile strength transverse 
        YC = props(13) ! in-situ compressive strength transverse
        ZT = props(14) ! in-situ tensile strength through thickness
        ZC = props(15) ! in-situ tensile strength through thickness
        SL = props(16) ! in-situ longitudinal shear strength
        SR = props(17) ! in-situ through thickness shear strength
        ST = props(18) ! in-situ transverse shear strength

        ! Fracture toughnesses (Gc)
        G1Plus = props(19) ! tensile Gc fiber (exponential)
        G1Minus = props(20) ! comp. Gc fiber
        G2Plus = props(21) ! tensile Gc trans.
        G6 = props(22) ! shear Gc
        G2Minus = G6/cos(0.92502450355699d0) ! comp. Gc trans.

        ! Define undulation properties
        pi = 4.0d0*atan(1.0d0)  ! calculate pi
        und_flag = props(23)  ! 1 if undulation, 0 if straight tow
        und_angle_1 = props(24)  ! in-plane angle constituent 1
        und_angle_2 = props(25)  ! in-plane angle constituent 2
        ! handle negative constituent angles
        if (und_angle_1 < 0.0d0) then
            und_angle_1 = 180.0d0 + und_angle_1
        end if
        if (und_angle_2 < 0.0d0) then
            und_angle_2 = 180.0d0 + und_angle_2
        end if
        ! calculate interface angle between constituents
        int_angle = (und_angle_2 - und_angle_1) * (pi/180.0d0)
        cpt = props(26)  ! load cured ply thickness
        u_length = props(27)  ! load undulation length

        ! Calculate remaining Poisson's ratios
        nu21 = nu12*(e22/e11) ! Poisson's ratio 21
        nu31 = nu13*(e33/e11) ! Poisson's ratio 31
        nu32 = nu23*(e33/e22) ! Poisson's ratio 32

        ! Calculate undamaged stiffness matrix for a tow
        psi = 1.0d0/(-e11*e22+e22**2*nu12**2+e11*e33*nu23**2+e22*e33*
     1        nu13**2+e22*e33*nu12*nu13*nu23*2.0d0)
        C_tow = 0.0d0
        C_tow(1,1) = -e11**2*(e22-e33*nu23**2)*psi
        C_tow(1,2) = -e11*e22*(e22*nu12+e33*nu13*nu23)*psi
        C_tow(1,3) = -e11*e22*e33*(nu13+nu12*nu23)*psi
        C_tow(2,1) = -e11*e22*(e22*nu12+e33*nu13*nu23)*psi
        C_tow(2,2) = -e22**2*(e11-e33*nu13**2)*psi
        C_tow(2,3) = -e22*e33*(e11*nu23+e22*nu12*nu13)*psi
        C_tow(3,1) = -e11*e22*e33*(nu13+nu12*nu23)*psi
        C_tow(3,2) = -e22*e33*(e11*nu23+e22*nu12*nu13)*psi
        C_tow(3,3) = -e22*e33*(e11-e22*nu12**2)*psi
        C_tow(4,4) = g12
        C_tow(5,5) = g23
        C_tow(6,6) = g13

        ! Calculate undamaged stiffness for pure resin
        psi_resin = e22/((1.0d0 + nu12) * (1.0d0 - 2.0d0 * nu12))
        C_resin = 0.0d0
        C_resin(1,1) = (1.0d0 - nu12) * psi_resin
        C_resin(1,2) = nu12 * psi_resin
        C_resin(1,3) = nu12 * psi_resin
        C_resin(2,1) = nu12 * psi_resin
        C_resin(2,2) = (1.0d0 - nu12) * psi_resin
        C_resin(2,3) = nu12 * psi_resin
        C_resin(3,1) = nu12 * psi_resin
        C_resin(3,2) = nu12 * psi_resin
        C_resin(3,3) = (1.0d0 - nu12) * psi_resin
        C_resin(4,4) = g12
        C_resin(5,5) = g12
        C_resin(6,6) = g12

        ! Initial elastic step, for Abaqus tests 
        if (stepTime == 0) then
          do k = 1, nblock

            ! calculate characteristic lengths
            lch = charlength(k)
            vfrac_1 = 0.818308233d0  ! volume fraction of constituent 1
            vfrac_2 = 0.181691767d0  ! volume fraction of constituent 1
            lch_1 = lch * vfrac_1**(1.0d0/3.0d0)  ! char. length constituent 1
            lch_2 = lch * vfrac_2**(1.0d0/3.0d0)  ! char. length constituent 2

            ! Initialize state variables
            do i = 1,nstatev
              stateNew(k,i) = 0.d0
            end do
            ! Initial strain
            do i = 1,6
              strain(i) = strainInc(k,i)
            end do
            
            ! Calculate initial elastic stress (multiply stiffness and strain)
            if (und_flag == 1.0d0) then
              ! Calculate average out-of-plane angle
              call undulation_geometry(cpt,u_length,alpha_mean)
              ! Calculate stresses in undulation regions where
              ! constituent 2 is pure resin
              if (und_angle_2 == 2.0d0*pi) then
                stress_2 = matmul(C_resin, strain)
                ! rotate strain to material coordinate system
                call rot_matrix(int_angle, alpha_mean, T1)
                C_1_global = matmul(matmul(T1, C_tow), transpose(T1))
                stress_1 = matmul(C_1_global, strain)
              else
                ! Calculate stresses undulation regions where
                ! both constituents are tows
                ! rotate strain to material coordinate system
                call rot_matrix(int_angle, alpha_mean, T1)
                C_1_global = matmul(matmul(T1, C_tow), transpose(T1))
                stress_1 = matmul(C_1_global, strain)
                ! rotate strain to material coordinate system
                call rot_matrix(0.0d0, alpha_mean, T2)
                C_2_global = matmul(matmul(T2, C_tow), transpose(T2))
                stress_2 = matmul(C_2_global, strain)
              end if
              ! calculate total stress by volume averaging
              stress_total = stress_1 * vfrac_1 + stress_2 * vfrac_2
            else if (und_flag == 2.0d0) then
              ! Calculate stresses in resin rich regions
              stress_2 = matmul(C_resin, strain)
              ! rotate strain to material coordinate system
              call rot_matrix(int_angle, 0.0d0, T1)
              C_1_global = matmul(matmul(T1, C_tow), transpose(T1))
              stress_1 = matmul(C_1_global, strain)
              ! Calculate total stress by volume averaging
              stress_total = stress_1 * vfrac_1 + stress_2 * vfrac_2
            else
              ! Calculate stress straight tow regions
              stress_total = matmul(C_tow, strain)
            end if
            
            ! Assign calculated stresses to stressNew array
            do i = 1,6
              stressNew(k,i) = stress_total(i)
            end do

          end do
        else

          ! If not initial step, calc. stresses according to CDM model
          do k = 1,nblock
            ! calculate characteristic lengths
            lch = charlength(k)
            vfrac_1 = 0.818308233d0
            vfrac_2 = 0.181691767d0
            lch_1 = lch * vfrac_1**(1.0d0/3.0d0)
            lch_2 = lch * vfrac_2**(1.0d0/3.0d0)

            ! Calc. new strain and assign to state variable array (1-6)
            do i = 1,6
              stateNew(k,i) = stateOld(k,i) + strainInc(k,i)
            end do

            ! Assign strain values to strain array
            do i = 1,6
              strain(i) = stateNew(k,i)
            end do

            if (und_flag == 1.0d0) then
              ! if element is in an undulation region calculate
              ! degraded stress for each consistuent independently
              ! first calculate mean oop angle
              call undulation_geometry(cpt,u_length,alpha_mean)
              if (und_angle_2 == 2.0d0*pi) then
                ! calculate stresses in pure resin consituent
                call cdm_resin(k,nblock,nstatev,strain,stateOld,C_resin,
     1                         e11,G2Plus,G2Minus,YT,YC,lch_2,stress_2,
     2                         stateNew)
                ! calculate stresses in tow constituent
                call cdm(k,nblock,nstatev,strain,stateOld,
     1                   C_tow,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2                   XT,XC,YT,YC,ZT,ZC,SL,SR,ST,int_angle,
     3                   alpha_mean,25,G1Plus,G1Minus,G2Plus,G2Minus,
     4                   lch_1,stress_1,stateNew)
              else
                ! calculate stresses in tow 1 (if tow 1 \= pure resin)
                call cdm(k,nblock,nstatev,strain,stateOld,
     1                   C_tow,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2                   XT,XC,YT,YC,ZT,ZC,SL,SR,ST,int_angle,
     3                   0.0d0,6,G1Plus,G1Minus,G2Plus,G2Minus,
     4                   lch_1,stress_1,stateNew)
                ! claculate stresses in tow 2
                call cdm(k,nblock,nstatev,strain,stateOld,
     1                   C_tow,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2                   XT,XC,YT,YC,ZT,ZC,SL,SR,ST,0.0d0,
     3                   alpha_mean,25,G1Plus,G1Minus,G2Plus,G2Minus,
     4                   lch_2,stress_2,stateNew)
              end if
              ! calculate homogenized stresses by volume-averaging
              stress_total = stress_1 * vfrac_1 + stress_2 * vfrac_2
            else if (und_flag == 2.0d0) then
              ! calculate stresses in resin rich regions
              call cdm_resin(k,nblock,nstatev,strain,stateOld,C_resin,
     1                       e11,G2Plus,G2Minus,YT,YC,lch_2,stress_2,
     2                       stateNew)
              ! calculate stresses in tow 1 (if tow 1 \= pure resin)
              call cdm(k,nblock,nstatev,strain,stateOld,
     1                 C_tow,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2                 XT,XC,YT,YC,ZT,ZC,SL,SR,ST,int_angle,
     3                 0.0d0,25,G1Plus,G1Minus,G2Plus,G2Minus,
     4                 lch_1,stress_1,stateNew)
              ! calculate homogenized stresses by volume averaging
              stress_total = stress_1 * vfrac_1 + stress_2 * vfrac_2
            else
              ! if element is not in undulation region calculate
              ! stress as normal
              call cdm(k,nblock,nstatev,strain,stateOld,
     1                 C_tow,e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2                 XT,XC,YT,YC,ZT,ZC,SL,SR,ST,0.0d0,0.0d0,6,
     3                 G1Plus,G1Minus,G2Plus,G2Minus,lch,
     4                 stress_total,stateNew)
            end if

            ! Element deletion criteria
            ! calculate the determinant of the deformation gradient
            F = det_3x3(defGradNew(k,:))
            stateNew(k,45) = F  ! store as state variable
            if (((F > 0.0d0).and.(F < 0.5d0)).or.(F >= 2.5d0)) then
                  stateNew(k,46) = 0.0d0  ! flag element for deletion
            end if

            ! assign stresses to stressNew array
            do i = 1,6
              stressNew(k,i) = stress_total(i)
            end do
            ! calculate internal energy
            stress_power = 0.5d0*((stressNew(k,1)+stressOld(k,1))*
     1                    strainInc(k,1)+(stressNew(k,2)+
     2                    stressOld(k,2))*strainInc(k,2)+
     3                    (stressNew(k,3)+stressOld(k,3))*
     4                    strainInc(k,3)+2.0d0*(stressNew(k,4)+
     5                    stressOld(k,4))*strainInc(k,4)+
     6                    2.0d0*(stressNew(k,5)+stressOld(k,5))*
     7                    strainInc(k,5)+2.0d0*(stressNew(k,6)+
     8                    stressOld(k,6))*strainInc(k,6))
            ! assign energies to enerInternNew array
            enerInternNew(k) = enerInternOld(k)+stress_power/density(k)

          end do
        end if
      end subroutine vumat

      subroutine cdm(k,nblock,nstatev,strain,stateOld,C_tow,
     1               e11,e22,e33,nu12,nu13,nu23,g12,g13,g23,
     2               XT,XC,YT,YC,ZT,ZC,SL,SR,ST,ip_angle,oop_angle,idx,
     3               G1Plus,G1Minus,G2Plus,G2Minus,lch,
     4               stress_global,stateNew)
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        integer, intent(in) :: idx  ! starting index for state variables
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, dimension(6,6), intent(inout) :: C_tow  ! stiffness matrix tow regions
        real*8, dimension(6), intent(in) :: strain  ! total strain
        real*8, intent(in) :: e11,e22,e33,nu12,nu13,nu23,g12,g13,g23  ! elastic constants
        real*8, intent(in) :: XT,XC,YT,YC,ZT,ZC,SL,SR,ST  ! strengths
        real*8, intent(in) :: G1Plus,G1Minus,G2Plus,G2Minus  ! fracture toughnesses
        real*8, intent(in) :: lch  ! characteristic length
        real*8, intent(in) :: ip_angle,oop_angle  ! in-plane and out-of-plane angles
        ! local variables
        real*8, dimension(6,6) :: C_deg  ! degraded stiffness matrix
        real*8, dimension(6,6) :: T_gm,T_mg  ! transformation matrices (global to material and vice versa)
        real*8, dimension(6) :: stress_mat,strain_mat  ! stress/strain in material CSYS
        real*8 :: nu21,nu31,nu32  ! secondary Poisson's ratios
        real*8 :: d1,d2,d3  ! scalar damage variables
        real*8 :: psi,pi  ! stiffness matrix calculation parameter, and pi
        real*8 :: d1State,d1StateOld,d2State,d2StateOld  ! damage state vars.
        real*8 :: d3State,d3StateOld  ! damage state vars.
        ! output variables
        real*8, dimension(6), intent(out) :: stress_global  ! total stress in global CSYS
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array
        
        ! Recalculate additional Poisson's ratio terms
        nu21 = nu12*(e22/e11) ! Poisson's ratio 21
        nu31 = nu13*(e33/e11) ! Poisson's ratio 31
        nu32 = nu23*(e33/e22) ! Poisson's ratio 32

        pi = 4.0d0*atan(1.0d0)  ! calculate pi
        ! rotate strains to material coordinate system
        if ((ip_angle == 0.0d0).or.(ip_angle == 2.0d0*pi)) then
            call rot_matrix(0.0d0, oop_angle, T_gm)
            call rot_matrix(0.0d0, -oop_angle, T_mg)
        else
          call rot_matrix(ip_angle, oop_angle, T_gm)
          call rot_matrix(-ip_angle, -oop_angle, T_mg)
        end if
        strain_mat = matmul(T_gm, strain)
        stress_mat = matmul(C_tow, strain_mat)

        ! Calculate damage in fiber direction
        call long_damage(k,nblock,nstatev,stateOld,stateNew,
     1                   stress_mat,strain_mat,e11,e22,XT,XC,
     2                   G1Plus,G1Minus,lch,idx,d1)
        ! Calculate damage in transverse and through-thickness matrix directions
        call matrix_damage(k,nblock,nstatev,stress_mat,strain_mat,
     1                     stateNew,idx,stateOld,YT,YC,ZT,ZC,SL,SR,ST,
     2                     G2Plus,G2Minus,lch,d2,d3)

        ! Update state variables
        d1StateOld = stateOld(k,idx + 16)
        d1State = min(0.999d0, max(d1StateOld,d1))
        stateNew(k,idx + 16) = d1State
        d2StateOld = stateOld(k,idx + 17)
        d2State = min(0.95d0, max(d2StateOld,d2))
        stateNew(k,idx + 17) = d2State
        d3StateOld = stateOld(k,idx + 18)
        d3State = min(0.95d0, max(d3StateOld,d3))
        stateNew(k,idx + 18) = d3State

        ! Calculate degraded stiffness matrix
        psi = 1.0d0/(e11*e22-e22**2*nu12**2+d1State*e22**2*nu12**2+
     1    d2State*e22**2*nu12**2-e11*e33*nu23**2-e22*e33*nu13**2+
     2    d1State*e22*e33*nu13**2+d2State*e11*e33*nu23**2+d3State*
     3    e11*e33*nu23**2+d3State*e22*e33*nu13**2-d1State*d2State*
     4    e22**2*nu12**2-e22*e33*nu12*nu13*nu23*2.0d0-d1State*d3State*
     5    e22*e33*nu13**2-d2State*d3State*e11*e33*nu23**2+d1State*
     6    e22*e33*nu12*nu13*nu23*2.0d0+d2State*e22*e33*nu12*nu13*nu23*
     7    2.0d0+d3State*e22*e33*nu12*nu13*nu23*2.0d0-d1State*d2State*
     8    e22*e33*nu12*nu13*nu23*2.0d0-d1State*d3State*e22*e33*nu12*
     9    nu13*nu23*2.0d0-d2State*d3State*e22*e33*nu12*nu13*nu23*
     1    2.0d0+d1State*d2State*d3State*e22*e33*nu12*nu13*nu23*2.0d0)

        C_deg = 0.0d0
        C_deg(1,1) = -e11**2*(d1State-1.0d0)*(e22-e33*nu23**2+d2State*
     1    e33*nu23**2+d3State*e33*nu23**2-d2State*d3State*e33*nu23**2)
     2    *psi
        C_deg(1,2) = e11*e22*(d1State-1.0d0)*(d2State-1.0d0)*
     1    (e22*nu12+e33*nu13*nu23-d3State*e33*nu13*nu23)*psi
        C_deg(1,3) = e11*e22*e33*(d1State-1.0d0)*(d3State-1.0d0)*
     1    (nu13+nu12*nu23-d2State*nu12*nu23)*psi
        C_deg(2,1) = e11*e22*(d1State-1.0d0)*(d2State-1.0d0)*
     1    (e22*nu12+e33*nu13*nu23-d3State*e33*nu13*nu23)*psi
        C_deg(2,2) = -e22**2*(d2State-1.0d0)*(e11-e33*nu13**2+d1State*
     1    e33*nu13**2+d3State*e33*nu13**2-d1State*d3State*e33*nu13**2)*
     2    psi
        C_deg(2,3) = e22*e33*(d2State-1.0d0)*(d3State-1.0d0)*
     1    (e11*nu23+e22*nu12*nu13-d1State*e22*nu12*nu13)*psi
        C_deg(3,1) = e11*e22*e33*(d1State-1.0d0)*(d3State-1.0d0)*
     1    (nu13+nu12*nu23-d2State*nu12*nu23)*psi
        C_deg(3,2) = e22*e33*(d2State-1.0d0)*(d3State-1.0d0)*
     1    (e11*nu23+e22*nu12*nu13-d1State*e22*nu12*nu13)*psi
        C_deg(3,3) = -e22*e33*(d3State-1.0d0)*(e11-e22*nu12**2+d1State*
     1    e22*nu12**2+d2State*e22*nu12**2-d1State*d2State*e22*nu12**2)*
     2    psi
        C_deg(4,4) = g12 * (1.0d0 - d1State) * (1.0d0 - d2State)
        C_deg(5,5) = g23 * (1.0d0 - d2State) * (1.0d0 - d3State)
        C_deg(6,6) = g13 * (1.0d0 - d1State) * (1.0d0 - d3State)

        ! calculate stresses in material CSYS
        stress_mat = matmul(C_deg, strain_mat)
        ! rotate stresses to global CSYS
        stress_global = matmul(T_mg, stress_mat)
      end subroutine cdm

      function det_3x3(A) result(det)
        ! Computes the determinant of a 3x3 matrix
        ! matrix should be input as a 9 element vector
        implicit none
        ! input variables
        real*8, dimension(9), intent(in)  :: A ! 3x3 input matrix 
        real*8 :: det  ! determinant of A
        det = A(1)*(A(2)*A(3) - A(5)*A(8)) -
     1        A(4)*(A(7)*A(3) - A(5)*A(6)) +
     2        A(9)*(A(7)*A(8) - A(2)*A(6))
      end function det_3x3

      subroutine long_damage(k,nblock,nstatev,stateOld,stateNew,
     1                       stress,strain,e11,e22,XT,XC,
     2                       G1Plus,G1Minus,lch,idx,d1)
        ! Computes scalar damage variable in fiber direction, d1
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        integer, intent(in) :: idx  ! starting index for state variables
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, dimension(6), intent(in) :: stress, strain  ! input stress and strain
        real*8, intent(in) :: XT,XC  ! longitudinal strengths
        real*8, intent(in) :: G1Plus,G1Minus,lch  ! fracture toughnesses and char. length
        real*8, intent(in) :: e11,e22  ! moduli
        ! local variables
        real*8 :: d1Plus,d1Minus,d1MinusStar  ! sclar damage variables tensile/compressive loading
        real*8 :: r1Plus,r1Minus  ! eleastic domain thresholds
        real*8 :: A1Plus,A1Minus,A1PlusMinus  ! exponentential energy dissipation parameter
        real*8 :: F1_C,F1_T,F1_C_max,F1_T_max,F1_T_old,F1_C_old  ! failure criteria
        real*8 :: L_1  ! characteristic length at damage initiation
        ! output variables
        real*8, intent(out) :: d1  ! longitudinal scalar damage variable
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array
        
        ! Initialize damage activation functions (failure indices)
        F1_T = 0.0d0
        F1_C = 0.0d0

        ! Calculate longitudinal tensile failure indices
        if (stress(1) >= 0.d0) then
          F1_T = strain(1)/(XT/e11)  ! max. strain criteria
        else if (stress(1) < 0.d0) then
          F1_C = (-strain(1))/(XC/e11)  ! max. strain criteria
        end if

        ! Load old state variables
        F1_T_old = stateOld(k,idx + 1)
        F1_C_old = stateOld(k,idx + 2)
        ! Compute max failure indices
        F1_T_max = max(F1_T,F1_T_old)
        F1_C_max = max(F1_C,F1_C_old)
        ! Record failure indices in state array
        stateNew(k,idx + 1) = F1_T_max
        stateNew(k,idx + 2) = F1_C_max
        
        d1Plus = 0.0d0
        if ((F1_T_max >= 1.0d0).or.(F1_C_max >= 1.0d0)) then
          ! Check if damage has already initiated
          if (stateOld(k, idx + 3) == 0.0d0) then
            stateNew(k, idx + 3) = lch ! record char. length at onset
            L_1 = lch
          else
            ! if dmaage has already been triggered, load fiber
            ! direction characteristic length
            L_1 = stateOld(k,idx + 3)
            stateNew(k,idx + 3) = L_1
          end if
          ! calculate linear-exponential softening response
          r1Plus = max(F1_T_max,F1_C_max)  ! elastic domain threshold
          A1Plus = (2.0d0*L_1*XT**2.0d0)/
     1              (2.0d0*e11*G1Plus - L_1*XT**2.0d0)
          d1Plus = 1.0d0 - (1.0d0/r1Plus)*exp(A1Plus*(1.0d0-r1Plus))
        else
          d1Plus = 0.0d0
        end if

        ! Longitudinal compressive damage
        if (F1_C_max >= 1.0d0) then
          ! Check if damage has already initiated
          if (stateOld(k,idx + 3) == 0.0d0) then
            stateNew(k,idx + 3) = lch ! record char. length at onset
            L_1 = lch
          else
            L_1 = stateOld(k,idx + 3) ! load char. length
            stateNew(k,idx + 3) = L_1
          end if
          r1Minus = F1_C_max
          A1Minus = 0.0620116303502072d0 !(2.0d0*L_1*XC**2.0d0)/
!     1              (2.0d0*e11*G1Minus - L_1*XC**2.0d0)
          d1MinusStar = 1.0d0 - (1.0d0/r1Minus)*
     1                  exp(A1Minus*(1.0d0-r1Minus))
          A1PlusMinus = 0.5d0 * ((e11 - e22) / e11)
          d1Minus = 1.0d0 - (1.0d0 - d1MinusStar) *
     1              (1.0d0 - A1PlusMinus*d1Plus)
        else
          d1Minus = 0.0d0
        end if
        
        ! crack closure under load reversal
        if (stress(1) >= 0.0d0) then
            d1 = d1Plus
        else
            d1 = d1Minus
        end if

        ! prevent NaN damage variables
        if (d1 /= d1) then
          d1 = 0.0d0
        end if

      end subroutine long_damage

      subroutine matrix_damage(k,nblock,nstatev,stress,strain,stateNew,
     1                         idx,stateOld,YT,YC,ZT,ZC,SL,SR,ST,GTPlus,
     2                         GTMinus,lch,d2,d3)
        ! Computes transverse and through thickness damage variables
        ! d2 and d3
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        integer, intent(in) :: idx  ! starting index for state variables
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, dimension(6), intent(in) :: stress, strain  ! input stress and strain
        real*8, intent(in) :: YT,YC,ZT,ZC,SL,SR,ST  ! strengths
        real*8, intent(in) :: GTPlus,GTMinus  ! fracture toughnesses
        real*8, intent(in) :: lch  ! characteristic length
        ! local variables
        real*8 :: F2_T,F2_C,F3_T,F3_C  ! failure criteria
        real*8 :: s22,s23,s12,s33,s13  ! shorthand notation stresses
        real*8 :: d2Plus,d2Minus,d3Plus,d3Minus  ! damage variables
        integer :: id_2T,id_2C,id_3T,id_3C  ! indices for failure modes in state variable array
        ! output variables
        real*8, intent(out) :: d2,d3  ! scalar damage variables
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array
        ! external functions
        real*8, external :: mc  ! implements macaulay brackets

        ! reassign stresses for brevity
        s22 = stress(2)
        s33 = stress(3)
        s12 = stress(4)
        s23 = stress(5)
        s13 = stress(6)

        ! calculate failure indices
        F2_T = (mc(s22)/YT)**2.0d0 + (s12/SL)**2.0d0 + (s23/ST)**2.0d0
        F2_C = (mc(-s22)/(2.0d0 * ST))**2.0d0 + 
     1         ((YC/(2.0d0 * ST))**2.0d0) * (s22/YC) +
     2         (s12/SL)**2.0d0
        F3_T = (mc(s33)/ZT)**2.0d0 + (s13/SR)**2.0d0 + (s23/ST)**2.0d0
        F3_C = (mc(-s33)/(2.0d0 * ST))**2.0d0 + 
     1         ((ZC/(2.0d0 * ST))**2.0d0) * (s33/ZC) +
     2         (s13/SR)**2.0d0

        ! Define index for each failure mode in state variable array
        id_2T = idx + 4
        id_2C = idx + 7
        id_3T = idx + 10
        id_3C = idx + 13

        ! initialize damage variables
        d2Plus = 0.0d0
        d2Minus = 0.0d0
        d3Plus = 0.0d0
        d3Minus = 0.0d0

        if (stress(2) >= 0.0d0) then
          ! if failure index exceeds unity or damage has already initiated
          ! calculate value of d2
          if ((F2_T > 1.0d0).or.(stateOld(k,id_2T) /= 0.0d0)) then
            call calculate_d(k,nblock,nstatev,stateNew,stateOld,
     1                     strain(2),stress(2),lch,id_2T,GTPlus,d2Plus)
          else
            ! if failure index <1 d2 = 0
            d2Plus = 0.0d0
          end if
        else
          if ((F2_C > 1.0d0).or.(stateOld(k,id_2C) /= 0.0d0)) then
            call calculate_d(k,nblock,nstatev,stateNew,stateOld,
     1                       strain(2),stress(2),lch,id_2C,GTMinus,
     2                       d2Minus)
          else
            d2Minus = 0.0d0
          end if
        end if

        if (stress(3) >= 0.0d0) then
          if ((F3_T > 1.0d0).or.(stateOld(k,id_3T) /= 0.0d0)) then
            call calculate_d(k,nblock,nstatev,stateNew,stateOld,
     1                       strain(3),stress(3),lch,id_3T,GTPlus,
     2                       d3Plus)
          else
            d3Plus = 0.0d0
          end if
        else
          if ((F3_C > 1.0d0).or.(stateOld(k,id_3C) /= 0.0d0)) then
            call calculate_d(k,nblock,nstatev,stateNew,stateOld,
     1                       strain(3),stress(3),lch,id_3C,GTMinus,
     2                       d3Minus)
          else
            d3Minus = 0.0d0
          end if
        end if
        
        d2 = 1.0d0 - (1.0d0 - d2Minus) * (1.0d0 - d2Plus)
        d3 = 1.0d0 - (1.0d0 - d3Minus) * (1.0d0 - d3Plus)
        ! catch NaN d2 and d3 values
        if (d2 /= d2) then
          d2 = 0.0d0
        end if
        if (d3 /= d3) then
          d3 = 0.0d0
        end if

      end subroutine matrix_damage

      subroutine calculate_d(k,nblock,nstatev,stateNew,stateOld,
     1                       delta,s,lch,stateID,G,d)
        ! Calculate scalar damage variable values
        implicit none
        ! input variables
        integer, intent(in) :: nblock,nstatev  ! num of integration points and num of state variables
        integer, intent(in) :: k  ! current integration point
        integer, intent(in) :: stateID  ! starting index state variable array
        real*8, dimension(nblock,nstatev), intent(in) :: stateOld  ! old state variable array
        real*8, intent(in) :: delta  ! strain (scalar)
        real*8, intent(in) :: s  ! stress (scalar)
        real*8, intent(in) :: lch  ! characteristic length
        real*8, intent(in) :: G  ! fracture toughness
        ! local variables
        real*8 :: delta_0,sigma_0  ! strain/stress at damage onset
        real*8 :: delta_f  ! final failure strain
        integer :: i  ! loop counter
        ! output variables
        real*8, intent(out) :: d  ! scalar damage variable
        real*8, dimension(nblock,nstatev), intent(inout) :: stateNew  ! new state variable array

        ! if damage has not yet been initiated
        if (stateOld(k,stateID) == 0.0d0) then
          delta_0 = delta  ! set strain at damage onset
          sigma_0 = s  ! set stress at damage onset
          delta_f = (2.0d0 * G) / (sigma_0 * lch)  ! calculate final failure strain
          if (abs(delta_f) <= abs(delta_0)) then  ! check final failure strain > strain at damage onset
            print *, 'Error: delta_f > delta_0. Consider reducing the element size.'
            call xplb_exit
          end if
          ! assign variables to state array
          stateNew(k, stateID) = lch
          stateNew(k, stateID + 1) = delta_0
          stateNew(k, stateID + 2) = delta_f
          d = 0.0d0
        else
          ! load state variables
          delta_0 = stateOld(k, stateID + 1)
          delta_f = stateOld(k, stateID + 2)
          do i = stateID, stateID + 2  ! ensure state retention
            stateNew(k, i) = stateOld(k, i)
          end do
          ! calculat d value (linear stress degradation)
          d = (delta_f * (delta - delta_0)) /
     1        (delta * (delta_f - delta_0))
        end if
      end subroutine calculate_d

      function mc(input) result(output)
        ! Macaulay operator
        implicit none
        real*8, intent(in) :: input  ! real input
        real*8 :: output  ! input if input > 0 else 0
        output = (input + abs(input)) / 2.0d0
      end function mc

      subroutine undulation_geometry(cpt, L_u, alpha_mean)
        ! Computes mean out-of-plane angle given undulation geometry
        implicit none
        ! input variables
        real*8, intent(in) :: cpt,L_u  ! cured ply thickness, undulation length
        ! local variables
        real*8, dimension(:), allocatable :: alpha  ! out-of-plane angle
        real*8 :: dx  ! step size
        real*8 :: pi  ! pi constant
        real*8, dimension(:), allocatable :: x  ! length discretized
        real*8, dimension(:), allocatable :: h  ! undulation profile
        real*8, dimension(:), allocatable :: dz  ! change in height
        integer :: i  ! loop counter
        integer :: n  ! number of segments
        ! output variables
        real*8, intent(out) :: alpha_mean  ! average out-of-plane angle

        ! Define discretization for undulation homogenization
        n = 15

        dx = (L_u / 2.0d0) / float(n)  ! define x step size
        pi = 4.0d0*atan(1.0d0)  ! define pi
        
        ! Allocate variable size arrays
        allocate(alpha(n))
        allocate(x(n))
        allocate(h(n))
        allocate(dz(n))

        do i = 1,n
            x(i) = dx * i  ! calculate x at step i
            h(i) = 0.5d0 * cpt + 0.5d0 * cpt * sin((pi / L_u) *
     1             (x(i) - L_u / 2.0d0))  ! calculate height in segment i
            if (i == 1) then  ! assign 0 to first element of dz array
              dz(i) = 0.0d0
            else
              dz(i) = h(i) - h(i - 1)
            end if
            alpha(i) = atan(dz(i) / dx)  ! calculate oop angle in segment i
        end do
        alpha_mean = sum(alpha) / float(n)  ! undulation angle
      end subroutine undulation_geometry