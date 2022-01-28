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
      
      subroutine rot_matrix(ip_angle, oop_angle, R)
        ! Computes rotation matrix for roation about two axes
        implicit none
        ! input variables
        real*8, intent(in) :: ip_angle, oop_angle  ! in-plane & out-of-plane angles
        ! local variables
        real*8 l1,l2,l3,m1,m2,m3,n1,n2,n3  ! components of direction cosines
        ! output variables
        real*8, dimension(6,6), intent(out) :: R  ! rotation matrix

        ! define components of direction cosines l,m,n
        l1 = cos(ip_angle) * cos(oop_angle)
        l2 = cos(oop_angle) * sin(ip_angle)
        l3 = -sin(oop_angle)
        m1 = -sin(ip_angle)
        m2 = cos(ip_angle)
        m3 = 0.0d0
        n1 = cos(ip_angle)*sin(oop_angle)
        n2 = sin(oop_angle)*sin(ip_angle)
        n3 = cos(oop_angle)

        ! create transformation matrix
        R = 0.0d0
        R(1,1) = l1**2.0d0
        R(1,2) = m1**2.0d0
        R(1,3) = n1**2.0d0
        R(1,4) = 2.0d0*l1*m1
        R(1,5) = 2.0d0*m1*n1
        R(1,6) = 2.0d0*n1*l1
        R(2,1) = l2**2.0d0
        R(2,2) = m2**2.0d0
        R(2,3) = n2**2.0d0
        R(2,4) = 2.0d0*l2*m2
        R(2,5) = 2.0d0*m2*n2
        R(2,6) = 2.0d0*n2*l2
        R(3,1) = l3**2.0d0
        R(3,2) = m3**2.0d0
        R(3,3) = n3**2.0d0
        R(3,4) = 2.0d0*l3*m3
        R(3,5) = 2.0d0*m3*n3
        R(3,6) = 2.0d0*n3*l3
        R(4,1) = l1*l2
        R(4,2) = m1*m2
        R(4,3) = n1*n2
        R(4,4) = l1*m2+l2*m1
        R(4,5) = m1*n2+m2*n1
        R(4,6) = n1*l2+n2*l1
        R(5,1) = l2*l3
        R(5,2) = m2*m3
        R(5,3) = n2*n3
        R(5,4) = l2*m3+l3*m2
        R(5,5) = m2*n3+m3*n2
        R(5,6) = l2*n3+l3*n2
        R(6,1) = l1*l3
        R(6,2) = m1*m3
        R(6,3) = n1*n3
        R(6,4) = l1*m3+l3*m1
        R(6,5) = m1*n3+m3*n1
        R(6,6) = l1*n3+l3*n1
      end subroutine rot_matrix