!*******************************************************************************
! Double precision (D) subroutines
!
! Copyright (C) 2022 Jon May, Peter Bird, Michele Carafa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!*******************************************************************************

module DSphere

                ! Basic tools for operations on a sphere
                ! and/or the surface of a spherical planet.
                ! Note that these tools are algebraic, and
                ! do not provide any graphical support.
                ! Exactly parallel to MODULE Sphere.f90,
                ! except that all REAL variables have been
                ! replaced with DOUBLE PRECISION variables,
                ! and subprogram names have "D" added as a
                ! sign that they expect DOUBLE PRECISION
                ! arguments and return DOUBLE PRECISION answers.
    !
    ! By Peter Bird, UCLA, May 1997 - April 2003, & December 2014;
    ! update to REAL*8 (DOUBLE PRECISION) in July 2015.
    ! Copyright (c) 1997, 1998, 1999, 2003, 2014, 2015 by
    ! Peter Bird and the Regents of the University of California.
    !-----------------------------------------------------------------
    !
    !                   CONTENTS OF THIS MODULE:
    !                  --------------------------
    ! INTENDED FOR THE USER TO CALL:
    !           DOUBLE PRECISION FUNCTION  DArc
    !                          SUBROUTINE  DCircles_Intersect
    !                    LOGICAL FUNCTION  DClockways
    !           DOUBLE PRECISION FUNCTION  DCompass
    !                          SUBROUTINE  DCross
    !           DOUBLE PRECISION FUNCTION  DDot
    !           DOUBLE PRECISION FUNCTION  DEasting
    !           DOUBLE PRECISION FUNCTION  DLength
    !                          SUBROUTINE  DLocal_Phi
    !                          SUBROUTINE  DLocal_Theta
    !                          SUBROUTINE  DLonLat_2_ThetaPhi
    !                          SUBROUTINE  DLonLat_2_Uvec
    !           DOUBLE PRECISION FUNCTION  DMagnitude
    !                          SUBROUTINE  DMake_Uvec
    !                          SUBROUTINE  DNorthEast_Convention
    !           DOUBLE PRECISION FUNCTION  DRelative_Compass
    !                          SUBROUTINE  DSpherical_Area
    !                          SUBROUTINE  DThetaPhi_2_LonLat
    !                          SUBROUTINE  DThetaPhi_2_Uvec
    !                          SUBROUTINE  DTurn_To
    !                          SUBROUTINE  DUvec_2_LonLat
    !                          SUBROUTINE  DUvec_2_ThetaPhi
    !                          SUBROUTINE  DUvec_2_PlungeAzimuth
    !
    ! UTILITY ROUTINE FOR DEBUGGING:
    !                          SUBROUTINE  DTraceback
    !-----------------------------------------------------------------

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: Pi        = 3.14159265358979D0 ! These values were confirmed
    DOUBLE PRECISION, PARAMETER :: Pi_over_2 = 1.57079632679490D0 ! with PROGRAM Check_Pi.
    DOUBLE PRECISION, PARAMETER :: Two_Pi    = 6.28318530717959D0
    DOUBLE PRECISION, PARAMETER :: degrees_per_radian = 57.2957795130823D0
    DOUBLE PRECISION, PARAMETER :: radians_per_degree = 0.0174532925199433D0

    ! ---------------------------------------------------------
    ! |       General Note on Unit Vectors in a Sphere        |
    ! | The Cartesian coordinate system in which these unit   |
    ! | vectors are expressed has its origin at the center    |
    ! | of the sphere. 1st axis outcrops at ( 0 E,  0 N).     |
    ! |                2nd axis outcrops at (90 E,  0 N).     |
    ! |                3rd axis outcrops at (?? E, 90 N).     |
    ! ---------------------------------------------------------

    CONTAINS

    DOUBLE PRECISION FUNCTION DArc (from_uvec, to_uvec)
      ! Returns length of great-circle arc in radians.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: from_uvec, to_uvec
        DOUBLE PRECISION :: crossed, dotted
        DOUBLE PRECISION, DIMENSION(3) :: t_vec
        CALL DCross (from_uvec, to_uvec, t_vec)
        crossed = DLength(t_vec) ! >= 0.
        dotted = DDot(from_uvec, to_uvec)
        DArc = DATAN2(crossed, dotted) ! 0.0D0 to Pi
    END FUNCTION DArc

    SUBROUTINE DCircles_Intersect (pole_a_uvec, dot_a, first_a_uvec, last_a_uvec, &
                                 & pole_b_uvec, dot_b, first_b_uvec, last_b_uvec, & ! input
                                 & overlap, number, point1_uvec, point2_uvec)   ! output
      ! Finds all the points of intersection (0, 1, or 2) between
      !   arcs of two small circles on a unit sphere.
      ! Note: The set of small circles includes great circles
      !   as a special case.  The set of arcs includes complete-
      !   circle arcs as a special case.
      ! The two arcs are generically named "a" and "b".
      ! The pole of each is given by a "pole_x_uvec" (unit vector,
      !   from center of sphere; see module data for definitions).
      ! The plane containing the small circle is specified by "dot_x",
      !   its distance from the origin.  (If either distance is
      !   greater than 1.00, there can be no points of intersection.)
      !"First" and "last" points on each arc are specified by uvecs
      !   (unit vectors) "first_x_uvec" and "last_x_uvec",
      !   where one goes counterclockwise around the pole
      !   first to last (assuming a viewpoint outside the
      !   unit sphere).  If "first" = "last", the arc is a full circle.
      ! NOTE: Although "first_x_uvec" and "last_x_uvec" will usually
      !   be points on the corresponding arcs, they need not be.
      !   Only their azimuths from "pole_x_uvec" are used.  In fact,
      !   if "first_x_uvec == last_x_uvec", then even the azimuths are
      !   not used, so ANY vector (even a zero vector) can be sent in
      !   these two positions to signal a complete small circle.
      ! CAUTION: Results will be strange and unpredictable if the
      !   the "unit vectors" supplied are not actually 1.00D0 long!
      ! Values returned are "number" (the count of intersection points)
      !   and unit vectors for as many points as were found.
      ! If only one intersection is found, it is always placed in point1_uvec.
      ! In the special case where the two small circles are the same
      !   circle, and share some common arc, the logical flag "overlap"
      !   is set, and the first and last points of the common arc are
      !   reported, and number = 2 on output.  If both of these overlapped
      !   circles are complete circles, these first and last points are
      !   the same point.
        IMPLICIT NONE
        REAL*8, DIMENSION(3), INTENT(IN)  :: pole_a_uvec, first_a_uvec, last_a_uvec, &
                                          &  pole_b_uvec, first_b_uvec, last_b_uvec
       !{Following DOUBLE PRECISION was present in the original MODULE Map_Projections.}
        DOUBLE PRECISION, INTENT(IN)    :: dot_a, dot_b
        LOGICAL, INTENT(OUT)            :: overlap
        INTEGER, INTENT(OUT)            :: number
        REAL*8, DIMENSION(3), INTENT(OUT) :: point1_uvec, point2_uvec

       !{Following DOUBLE PRECISIONs were present in the original MODULE Map_Projections.}
        DOUBLE PRECISION :: a_dot_b, alpha, beta, deterMinant, t
        DOUBLE PRECISION, DIMENSION(2,2) :: inverse
        LOGICAL :: ring_a, ring_b ! are the arcs actually complete circles?
        REAL*8 :: along, antipole_gap, &
                & first_a_radians, first_b_radians, last_a_radians, &
                & last_b_radians, min_radius, point_wrt_a_radians, point_wrt_b_radians, &
                & pole_gap
        REAL*8, PARAMETER :: tolerance = 0.0000000014D0
        REAL*8, DIMENSION(3) :: line_uvec, offset, t_vec

      ! Check for an easy answer (no intersections):
        number = 0
        overlap = .FALSE.
        IF (dot_a >=  1.0D0) RETURN
        IF (dot_a <= -1.0D0) RETURN
        IF (dot_b >=  1.0D0) RETURN
        IF (dot_b <= -1.0D0) RETURN

      ! Decide whether arcs are complete circles; otherwise, deterMine
      !   (relative) azimuths of endpoints:
        ring_a = (first_a_uvec(1) == last_a_uvec(1)).AND. &
               & (first_a_uvec(2) == last_a_uvec(2)).AND. &
               & (first_a_uvec(3) == last_a_uvec(3))
        ring_b = (first_b_uvec(1) == last_b_uvec(1)).AND. &
               & (first_b_uvec(2) == last_b_uvec(2)).AND. &
               & (first_b_uvec(3) == last_b_uvec(3))
        IF (.NOT.ring_a) THEN
            first_a_radians = -DRelative_Compass(pole_a_uvec, first_a_uvec)
            last_a_radians  = -DRelative_Compass(pole_a_uvec, last_a_uvec)
            IF (first_a_radians > last_a_radians) last_a_radians = last_a_radians + Two_Pi
        END IF
        IF (.NOT.ring_b) THEN
            first_b_radians = -DRelative_Compass(pole_b_uvec, first_b_uvec)
            last_b_radians  = -DRelative_Compass(pole_b_uvec, last_b_uvec)
            IF (first_b_radians > last_b_radians) last_b_radians = last_b_radians + Two_Pi
        END IF

      ! Test for special cases of parallel and antipodal poles:
        t_vec = pole_a_uvec - pole_b_uvec
        pole_gap = DLength(t_vec)
        t_vec = pole_a_uvec + pole_b_uvec
        antipole_gap = DLength(t_vec)
        IF (pole_gap <= tolerance) THEN ! poles virtually identical
            IF (dot_a /= dot_b) RETURN ! no intersection
            ! From here on, assume dot_a == dot_b:
            IF (ring_a) THEN
                overlap = .TRUE.
                number = 2
                point1_uvec = first_b_uvec
                point2_uvec = last_b_uvec
            ELSE IF (ring_b) THEN
                overlap = .TRUE.
                number = 2
                point1_uvec = first_a_uvec
                point2_uvec = last_a_uvec
            ELSE ! neither circle is complete
                IF ((first_b_radians >= first_a_radians).AND. &
                   &(first_b_radians <= last_a_radians)) THEN ! 1st of b is in a.
                    point1_uvec = first_b_uvec
                    IF (last_b_radians <= last_a_radians) THEN
                        point2_uvec = last_b_uvec
                    ELSE
                        point2_uvec = last_a_uvec
                    END IF
                    IF ((point1_uvec(1) == point2_uvec(1)).AND. &
                       &(point1_uvec(2) == point2_uvec(2)).AND. &
                       &(point1_uvec(3) == point2_uvec(3))) THEN ! only one common point
                        overlap = .FALSE.
                        number = 1
                    ELSE ! a common arc
                        overlap = .TRUE.
                        number = 2
                    END IF
                ELSE IF ((first_a_radians >= first_b_radians).AND. &
                        &(first_a_radians <= last_b_radians)) THEN ! 1st of a is in b.
                    point1_uvec = first_a_uvec
                    IF (last_a_radians <= last_b_radians) THEN
                        point2_uvec = last_a_uvec
                    ELSE
                        point2_uvec = last_b_uvec
                    END IF
                    IF ((point1_uvec(1) == point2_uvec(1)).AND. &
                       &(point1_uvec(2) == point2_uvec(2)).AND. &
                       &(point1_uvec(3) == point2_uvec(3))) THEN ! only one common point
                        overlap = .FALSE.
                        number = 1
                    ELSE ! a common arc
                        overlap = .TRUE.
                        number = 2
                    END IF
                END IF ! separate arcs of same circle [no ELSE; we RETURN, with number = 0]
            END IF ! either circle is complete, or neither is
        ELSE IF (antipole_gap <= tolerance) THEN ! antipodal
            IF (dot_a /= -dot_b) RETURN ! no intersection
            ! From here on, assume dot_a == -dot_b:
            IF (ring_a) THEN
                overlap = .TRUE.
                number = 2
                point1_uvec = first_b_uvec
                point2_uvec = last_b_uvec
            ELSE IF (ring_b) THEN
                overlap = .TRUE.
                number = 2
                point1_uvec = first_a_uvec
                point2_uvec = last_a_uvec
            ELSE ! neither circle is complete
                ! Because of antipodal relationship, azimuths are confusing.
                ! Redefine azimuths of arc b in terms of pole a,
                !  while reversing the direction along arc b to make it
                !  counterclockwise about a:
                first_b_radians = -DRelative_Compass(pole_a_uvec, last_b_uvec)
                last_b_radians  = -DRelative_Compass(pole_a_uvec, first_b_uvec)
                IF (last_b_radians < first_b_radians) last_b_radians = last_b_radians + Two_Pi
                IF ((first_b_radians >= first_a_radians).AND. &
                   &(first_b_radians <= last_a_radians)) THEN ! 1st of b is in a.
                    point1_uvec = first_b_uvec
                    IF (last_b_radians <= last_a_radians) THEN
                        point2_uvec = last_b_uvec
                    ELSE
                        point2_uvec = last_a_uvec
                    END IF
                    IF ((point1_uvec(1) == point2_uvec(1)).AND. &
                       &(point1_uvec(2) == point2_uvec(2)).AND. &
                       &(point1_uvec(3) == point2_uvec(3))) THEN ! only one common point
                        overlap = .FALSE.
                        number = 1
                    ELSE ! a common arc
                        overlap = .TRUE.
                        number = 2
                    END IF
                ELSE IF ((first_a_radians >= first_b_radians).AND. &
                        &(first_a_radians <= last_b_radians)) THEN ! 1st of a is in b.
                    point1_uvec = first_a_uvec
                    IF (last_a_radians <= last_b_radians) THEN
                        point2_uvec = last_a_uvec
                    ELSE
                        point2_uvec = last_b_uvec
                    END IF
                    IF ((point1_uvec(1) == point2_uvec(1)).AND. &
                       &(point1_uvec(2) == point2_uvec(2)).AND. &
                       &(point1_uvec(3) == point2_uvec(3))) THEN ! only one common point
                        overlap = .FALSE.
                        number = 1
                    ELSE ! a common arc
                        overlap = .TRUE.
                        number = 2
                    END IF
                END IF ! separate arcs on same circle [no ELSE; we RETURN, with number = 0]
            END IF ! either circle is complete, or neither is

        ELSE ! **** normal case; unrelated poles *****
            ! Each small circle lies in its own plane.
            ! These two planes intersect in a line in 3-D.
            ! Find "offset" = (non-unit) vector to point on this line
            !   which is closest to origin:
            a_dot_b = DDot(pole_a_uvec, pole_b_uvec)
            deterMinant = 1.0D0 - a_dot_b**2
            IF (deterMinant == 0.0D0) THEN
                WRITE (*,"(' ERROR: DeterMinant = 0.0D0 in DCircles_Intersect')")
                CALL DTraceback
            END IF
            t = 1.0D0 / deterMinant
            inverse(1,1) =  t           !  t * matrix(2,2)
            inverse(1,2) = -t * a_dot_b ! -t * matrix(1,2)
            inverse(2,1) = -t * a_dot_b ! -t * matrix(2,1)
            inverse(2,2) =  t           !  t * matrix(1,1)
            alpha = inverse(1,1)*dot_a + inverse(1,2)*dot_b
            beta  = inverse(2,1)*dot_a + inverse(2,2)*dot_b
            offset(1) = alpha*pole_a_uvec(1) + beta*pole_b_uvec(1)
            offset(2) = alpha*pole_a_uvec(2) + beta*pole_b_uvec(2)
            offset(3) = alpha*pole_a_uvec(3) + beta*pole_b_uvec(3)
            min_radius = DLength(offset)
            IF (min_radius == 1.0D0) THEN ! circles osculate at one point
                ! Check longitudes about poles to see if point is in arcs:
                IF (.NOT.ring_a) THEN
                    point_wrt_a_radians = -DRelative_Compass(pole_a_uvec, offset)
                    IF (point_wrt_a_radians < first_a_radians) point_wrt_a_radians = point_wrt_a_radians + Two_Pi
                END IF
                IF (ring_a .OR. &
                   &((point_wrt_a_radians >= first_a_radians).AND. &
                   & (point_wrt_a_radians <= last_a_radians))) THEN
                    IF (.NOT.ring_b) THEN
                        point_wrt_b_radians = -DRelative_Compass(pole_b_uvec, offset)
                        IF (point_wrt_b_radians < first_b_radians) point_wrt_b_radians = point_wrt_b_radians + Two_Pi
                    END IF
                    IF (ring_b .OR. &
                       &((point_wrt_b_radians >= first_b_radians).AND. &
                       & (point_wrt_b_radians <= last_b_radians))) THEN
                        number = 1
                        point1_uvec = offset
                    END IF ! in arc b [no ELSE; we RETURN with number = 0]
                END IF ! in arc a [no ELSE; we RETURN with number = 0]
            ELSE IF (min_radius < 1.) THEN ! two intersection points between circles
                ! Find vector parallel to line of intersection of the two
                !    planes which contain the two small circles:
                CALL DCross (pole_a_uvec, pole_b_uvec, t_vec)
                CALL DMake_Uvec (t_vec, line_uvec)
                ! Now, points on this line are expressed by
                !   vector "offset" plus any multiple of vector "line_uvec".
                ! Call the multiple "along".
                ! Solve the hyperbolic equation that says that the
                !   radius of a point on this line is 1.00
                !  (i.e., it is also a point on the unit sphere):
                !         min_radius**2 + along**2 = 1.00**2
                along = DSQRT ( 1.0D0 - min_radius**2 )
                point1_uvec = offset + along * line_uvec
                point2_uvec = offset - along * line_uvec
                ! Remember, these are only tentative solutions.

                ! Check longitudes about poles to see if point1 is in arcs:
                IF (.NOT.ring_a) THEN
                    point_wrt_a_radians = -DRelative_Compass(pole_a_uvec, point1_uvec)
                    IF (point_wrt_a_radians < first_a_radians) point_wrt_a_radians = point_wrt_a_radians + Two_Pi
                END IF
                IF (ring_a .OR. &
                   &((point_wrt_a_radians >= first_a_radians).AND. &
                   & (point_wrt_a_radians <= last_a_radians))) THEN
                    IF (.NOT.ring_b) THEN
                        point_wrt_b_radians = -DRelative_Compass(pole_b_uvec, point1_uvec)
                        IF (point_wrt_b_radians < first_b_radians) point_wrt_b_radians = point_wrt_b_radians + Two_Pi
                    END IF
                    IF (ring_b .OR. &
                       &((point_wrt_b_radians >= first_b_radians).AND. &
                       & (point_wrt_b_radians <= last_b_radians))) THEN
                        number = 1
                    END IF ! point1 also in arc b: a solution.
                END IF ! point1 is in arc a

                ! Check longitudes about poles to see if point2 is in arcs:
                IF (.NOT.ring_a) THEN
                    point_wrt_a_radians = -DRelative_Compass(pole_a_uvec, point2_uvec)
                    IF (point_wrt_a_radians < first_a_radians) point_wrt_a_radians = point_wrt_a_radians + Two_Pi
                END IF
                IF (ring_a .OR. &
                   &((point_wrt_a_radians >= first_a_radians).AND. &
                   & (point_wrt_a_radians <= last_a_radians))) THEN
                    IF (.NOT.ring_b) THEN
                        point_wrt_b_radians = -DRelative_Compass(pole_b_uvec, point2_uvec)
                        IF (point_wrt_b_radians < first_b_radians) point_wrt_b_radians = point_wrt_b_radians + Two_Pi
                    END IF
                    IF (ring_b .OR. &
                       &((point_wrt_b_radians >= first_b_radians).AND. &
                       & (point_wrt_b_radians <= last_b_radians))) THEN
                        IF (number == 1) THEN ! add a second solution
                            number = 2
                        ELSE ! point2 is the first and only solution
                            number = 1
                            point1_uvec = point2_uvec
                        END IF ! is this the first or the second solution?
                    END IF ! point2 in arc b: a solution.
                END IF ! point2 is in arc a

            END IF ! min_radius = 1.0D0, or < 1.0D0 [no ELSE; we RETURN with number = 0]
        END IF ! poles parallel, or antipodal, or unrelated
    END SUBROUTINE DCircles_Intersect

    DOUBLE PRECISION FUNCTION DCompass (from_uvec, to_uvec)
      ! Returns azimuth (in radians, clockwise from North)
      !   of the great-circle arc from "from_uvec" to "to_uvec",
      !   measured at location "from_uvec".
      ! Does NOT work at North or South pole!  (See DRelative_Compass.)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: from_uvec, to_uvec
        DOUBLE PRECISION :: ve, vn
        DOUBLE PRECISION, DIMENSION(3) :: omega_vec, omega_uvec, &
                                        & Phi, Theta, &
                                        & v_vec
        IF ((from_uvec(1) == 0.0D0).AND.(from_uvec(2) == 0.0D0)) THEN
            WRITE (*,"(' ERROR: Compass undefined at N or S pole.  Use Relative_Compass.')")
            CALL DTraceback
        ELSE IF ((from_uvec(1) == to_uvec(1)).AND. &
           &(from_uvec(2) == to_uvec(2)).AND. &
           &(from_uvec(3) == to_uvec(3))) THEN
            WRITE (*,"(' ERROR: Compass bearing from point to itself undefined.')")
            CALL DTraceback
        ELSE IF ((from_uvec(1) == -to_uvec(1)).AND. &
                &(from_uvec(2) == -to_uvec(2)).AND. &
                &(from_uvec(3) == -to_uvec(3))) THEN
            WRITE (*,"(' ERROR: Compass bearing from point to antipode undefined.')")
            CALL DTraceback
        ELSE
          ! Normal case:
            CALL DCross (from_uvec, to_uvec, omega_vec)
            CALL DMake_Uvec(omega_vec, omega_uvec)
            CALL DCross (omega_uvec, from_uvec, v_vec)
            CALL DLocal_Theta(from_uvec, Theta)
            CALL DLocal_Phi  (from_uvec, Phi)
            vn = -DDot(v_vec, Theta)
            ve =  DDot(v_vec, Phi)
            DCompass = DATAN2(ve, vn)
        END IF
    END FUNCTION DCompass

    SUBROUTINE DCross (a_vec, b_vec, c_vec)
      ! vector cross product: a x b = c
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: a_vec, b_vec
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: c_vec
        c_vec(1) = a_vec(2)*b_vec(3) - a_vec(3)*b_vec(2)
        c_vec(2) = a_vec(3)*b_vec(1) - a_vec(1)*b_vec(3)
        c_vec(3) = a_vec(1)*b_vec(2) - a_vec(2)*b_vec(1)
    END SUBROUTINE DCross

    DOUBLE PRECISION FUNCTION DDot (a_vec, b_vec)
      ! returns scalar (dot) product of two 3-component vectors
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: a_vec, b_vec
        DDot = a_vec(1)*b_vec(1) + a_vec(2)*b_vec(2) + a_vec(3)*b_vec(3)
    END FUNCTION DDot

    DOUBLE PRECISION FUNCTION DEasting(delta_lon_degrees)
        !returns positive result, 0.0D0 ~ 359.999...D0
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: delta_lon_degrees
        DEasting = DMOD((delta_lon_degrees + 720.0D0),360.0D0)
    END FUNCTION DEasting

    INTEGER FUNCTION IAbove(x)
      ! returns first integer >= x, unlike INT() or NINT():
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: x
        INTEGER answer
        DOUBLE PRECISION :: t
        answer = INT(x)
        IF (x > 0.0D0) THEN
            t = 1.00D0 * answer
            IF (x > t) THEN
                answer = answer + 1
            END IF
        END IF
        IAbove = answer
    END FUNCTION IAbove

    DOUBLE PRECISION FUNCTION DLength(a_vec)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: a_vec
        DOUBLE PRECISION :: t
        t = a_vec(1)**2 + &
          & a_vec(2)**2 + &
          & a_vec(3)**2
        IF (t == 0.0D0) THEN
            DLength = 0.0D0
        ELSE
            DLength = DSQRT(t)
        END IF
    END FUNCTION DLength

    SUBROUTINE DLocal_Phi (b_, Phi)
        ! returns local East-pointing unit vector in Cartesian coordinates
        ! for location b_; not intended to work at the poles!
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: b_
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: Phi
        DOUBLE PRECISION, DIMENSION(3)              :: temp
        IF (b_(1) == 0.0D0) THEN
            IF (b_(2) == 0.0D0) THEN
                WRITE (*,"(' ERROR: DLocal_Phi was requested for N or S pole.')")
                CALL DTraceback
             END IF
        END IF
        temp(1) = -b_(2)
        temp(2) = b_(1)
        temp(3) = 0.0D0
        CALL DMake_Uvec(temp, Phi)
    END SUBROUTINE DLocal_Phi

    SUBROUTINE DLocal_Theta (b_, Theta)
        ! returns local South-pointing unit vector in Cartesian coordinates
        ! for location b_; not intended to work at the poles!
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: b_
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: Theta
        DOUBLE PRECISION, DIMENSION(3) :: temp
        DOUBLE PRECISION  :: equat, new_equat
        equat = DSQRT(b_(1)**2 + b_(2)**2)   ! equatorial component
        IF (equat == 0.0D0) THEN
            WRITE (*,"(' ERROR: DLocal_Theta was requested for N or S pole.')")
            CALL DTraceback
        END IF
        new_equat = b_(3)   ! swap components in a meridional plane
        temp(3) = - equat       ! "
        temp(1) = new_equat * b_(1) / equat ! partition new equatorial component
        temp(2) = new_equat * b_(2) / equat ! "
        CALL DMake_Uvec (temp, Theta)
    END SUBROUTINE DLocal_Theta

    SUBROUTINE DLonLat_2_ThetaPhi (lon, lat, theta, phi)
      ! "lon" is East longitude in degrees.
      ! "lat" is North latitude in degrees.
      ! "theta" is co-latitude, from N pole, in radians
      ! "phi" is East longitude in radians
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)  :: lon, lat
        DOUBLE PRECISION, INTENT(OUT) :: theta, phi
        IF ((lat > 90.0D0).OR.(lat < -90.0D0)) THEN
            WRITE (*,"(' ERROR: Latitude outside legal range.')")
            CALL DTraceback
        ELSE
            theta = radians_per_degree * (90.0D0 - lat)
            phi = radians_per_degree * lon
        END IF
    END SUBROUTINE DLonLat_2_ThetaPhi

    SUBROUTINE DLonLat_2_Uvec (lon, lat, uvec)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: lon, lat
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: uvec
        DOUBLE PRECISION :: theta, phi
        CALL DLonLat_2_ThetaPhi (lon, lat, theta, phi)
        CALL DThetaPhi_2_Uvec (theta, phi, uvec)
    END SUBROUTINE DLonLat_2_Uvec

    DOUBLE PRECISION FUNCTION DMagnitude (b_)
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: b_
        DMagnitude = DSQRT(b_(1)**2 +b_(2)**2 + b_(3)**2)
    END FUNCTION DMagnitude

    !SUBROUTINE MATHERRQQ( name, length, info, retcode)
    !    !Provided so that domain errors, underflows, etc.
    !    !  can be trapped and debugged; otherwise program crashes!
    !    USE DFLIB
    !    INTEGER(2) length, retcode
    !    CHARACTER(length) name
    !    RECORD /MTH$E_INFO/ info
    !    PRINT *, "Entered MATHERRQQ"
    !    PRINT *, "Failing function is: ", name
    !    PRINT *, "Error type is: ", info.errcode
    !    IF ((info.ftype == TY$REAL4 ).OR.(info.ftype == TY$REAL8)) THEN
    !        PRINT *, "Type: REAL"
    !        PRINT *, "Enter the desired function result: "
    !        READ(*,*) info.r8res
    !        retcode = 1
    !    ELSE IF ((info.ftype == TY$CMPLX8 ).OR.(info.ftype == TY$CMPLX16)) THEN
    !        PRINT *, "Type: COMPLEX"
    !        PRINT *, "Enter the desired function result: "
    !        READ(*,*) info.c16res
    !        retcode = 1
    !    END IF
    !END SUBROUTINE MATHERRQQ

    SUBROUTINE DMake_Uvec (vector, uvec)
      ! Shortens or lengthens a three-component vector to a unit vector.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: vector
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: uvec
        DOUBLE PRECISION :: factor, size
        size = DLength(vector)
        IF (size > 0.0D0) THEN
           factor = 1.0D0 / size
           uvec = vector * factor
        ELSE
           WRITE (*,"(' ERROR: Cannot DMake_Uvec of (0.0D0, 0.0D0, 0.0D0).')")
           CALL DTraceback
        END IF
    END SUBROUTINE DMake_Uvec

    SUBROUTINE DNorthEast_Convention (location_uvec, north_uvec, east_uvec)
      ! At most positions ("location_uvec") returns "north_uvec"
      ! (same as -DLocal_Theta) and "east_uvec" (same as +DLocal_Phi).
      ! However, within a small distance from either the North or South
      ! pole, it adopts arbitrary conventional directions based on
      ! the limiting directions along the 0E meridian as the latitude
      ! approaches +90N or -90N.  Since both FUNCTION  DRelative_Compass
      ! and SUBROUTINE  DTurn_To call this routine, they can work
      ! together even when location_uvec is at, or near, a pole.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: location_uvec
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: north_uvec, east_uvec
        DOUBLE PRECISION, PARAMETER :: pole_area = 7.615D-7 ! (0.05 degrees -> radians -> squared)
        DOUBLE PRECISION :: equat2
        DOUBLE PRECISION, DIMENSION(3) :: t_vec
        equat2 = location_uvec(1)**2 + location_uvec(2)**2
        IF (equat2 > pole_area) THEN
          ! normal case:
            t_vec(1) = -location_uvec(2)
            t_vec(2) = +location_uvec(1)
            t_vec(3) = 0.0D0
            CALL DMake_Uvec (t_vec, east_uvec)
        ELSE
          ! very close to N or S pole; act as if on 0E meridian:
            east_uvec  = (/ 0.0D0, 1.0D0, 0.0D0 /)
        END IF
        CALL DCross (location_uvec, east_uvec, north_uvec)
    END SUBROUTINE DNorthEast_Convention

    DOUBLE PRECISION FUNCTION DRelative_Compass (from_uvec, to_uvec)
      ! At most points, works exactly like DCompass:
      !   returns azimuth (in radians, clockwise from North)
      !   of the great-circle arc from "from_uvec" to "to_uvec",
      !   measured at location from_uvec.
      ! However, unlike DCompass, it does not crash at the N and S poles,
      !   where Theta and Phi are undefined.  Instead, it uses
      !   SUBROUTINE  DNorthEast_Convention to make an
      !   arbitrary choice of axes, so that RELATIVE azimuths can
      !   be measured from pivot points at the poles, by multiple calls.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: from_uvec, to_uvec
        DOUBLE PRECISION :: ve, vn
        DOUBLE PRECISION, DIMENSION(3) :: omega_vec, omega_uvec, &
                                        & north_uvec, east_uvec, &
                                        & v_vec
        IF ((from_uvec(1) == to_uvec(1)).AND. &
           &(from_uvec(2) == to_uvec(2)).AND. &
           &(from_uvec(3) == to_uvec(3))) THEN
            WRITE (*,"(' ERROR: Compass bearing from point to itself undefined.')")
            CALL DTraceback
        ELSE IF ((from_uvec(1) == -to_uvec(1)).AND. &
                &(from_uvec(2) == -to_uvec(2)).AND. &
                &(from_uvec(3) == -to_uvec(3))) THEN
            WRITE (*,"(' ERROR: Compass bearing from point to antipode undefined.')")
            CALL DTraceback
        ELSE
          ! Normal case:
            CALL DCross (from_uvec, to_uvec, omega_vec)
            CALL DMake_Uvec(omega_vec, omega_uvec)
            CALL DCross (omega_uvec, from_uvec, v_vec)
            CALL DNorthEast_Convention (from_uvec, north_uvec, east_uvec)
            vn = DDot(v_vec, north_uvec)
            ve = DDot(v_vec, east_uvec)
            DRelative_Compass = DATAN2(ve, vn)
        END IF
    END FUNCTION DRelative_Compass

    SUBROUTINE DSpherical_Area(vertex1, vertex2, vertex3, steradians)
      ! Given 3 Cartesian (x, y, z) unit-vectors,
      ! each from the center of a unit-sphere to its surface,
      ! defining 3 surface points (vertex1, vertex2, vertex3),
      ! returns the area of the spherical triangle in steradians.
      !{For dimensional area, multiply this by the square of the radius.}
      ! Reference is:
      !     I. Todhunter [1886]
      !     Spherical Trigonometry
      !     McMillan & Co., London {5th edition};
      !     http://www.gutenberg.org/ebooks/19770 ,
      !     especially sections 97 (large triangles) and
      !     section 109 (small triangles).
      ! Additional feature added in this code (only):
      ! If vertices are specified in counterclockwise order
      !(as seen from outside the sphere), then the area is positive;
      ! if they are specified in clockwise order, then the area is
      ! returned as negative.  If the 3 unit-vectors are coplanar,
      ! or if at least 2 of them are identical, then the area is zero.
      ! For testing of the ideal cross-over point between answer1 and answer2,
      ! see my files Test_DSpherical_Area.f90, ...txt, ...xlsx .
        IMPLICIT NONE
        REAL*8, DIMENSION(3), INTENT(IN) :: vertex1, vertex2, vertex3
        REAL*8, INTENT(OUT) :: steradians
        LOGICAL :: positive
        REAL*8 :: A, alpha, answer1, answer2, azim12, azim13, azim21, azim23, azim31, azim32, &
                & B, beta, C, dotted, E, fraction, gamma, plane_area
        REAL*8, DIMENSION(3) :: crossed, v2_minus_v1, v3_minus_v1
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !Test for sign of answer, including possibility that answer is zero:
        v2_minus_v1(1:3) = vertex2(1:3) - vertex1(1:3)  ! not a uvec
        v3_minus_v1(1:3) = vertex3(1:3) - vertex1(1:3)  ! not a uvec
        CALL DCross (v2_minus_v1, v3_minus_v1, crossed) ! not a uvec
        dotted = DDot(crossed, vertex1)
        IF (dotted == 0.0D0) THEN
            steradians = 0.0D0
            RETURN
        ELSE
            positive = (dotted > 0.0D0)
        END IF
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !Method for large triangles:
        ! (All azimuths are in radians, clockwise from North.)
        azim12 = DRelative_Compass (from_uvec = vertex1, to_uvec = vertex2)
        azim13 = DRelative_Compass (from_uvec = vertex1, to_uvec = vertex3)
        azim21 = DRelative_Compass (from_uvec = vertex2, to_uvec = vertex1)
        azim23 = DRelative_Compass (from_uvec = vertex2, to_uvec = vertex3)
        azim31 = DRelative_Compass (from_uvec = vertex3, to_uvec = vertex1)
        azim32 = DRelative_Compass (from_uvec = vertex3, to_uvec = vertex2)
        A = ABS(azim12 - azim13)
        IF (A > Pi) A = Two_Pi - A ! where Two_Pi is a predefined global
        B = ABS(azim21 - azim23)
        IF (B > Pi) B = Two_Pi - B
        C = ABS(azim31 - azim32)
        IF (C > Pi) C = Two_Pi - C
        E = A + B + C - Pi         ! where Pi is a predefined global
        IF (positive) THEN
            answer1 = E
        ELSE
            answer1 = -E
        END IF
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !Method for small triangles:
        alpha = DArc(vertex1, vertex2) ! N.B. The naming of these 3 side-arcs may not
        beta =  DArc(vertex2, vertex3) !      necessarily agree with naming in Todhunter[1886](?),
        gamma = DArc(vertex3, vertex1) !      but this will make no difference to our answer.
        fraction = (alpha**2 + beta**2 + gamma**2) / 24.0D0
        plane_area = DMagnitude(crossed) / 2.0D0
        IF (positive) THEN
            answer2 = plane_area * (1.0D0 + fraction)
        ELSE
            answer2 = -plane_area * (1.0D0 + fraction)
        END IF
        !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !Selection of final answer:
        IF (E > 1.0D-8) THEN ! "large" triangle formula is more accurate:
            steradians = answer1
        ELSE ! "small" triangle formula is more accurate:
            steradians = answer2
        END IF
    END SUBROUTINE DSpherical_Area

    SUBROUTINE DThetaPhi_2_LonLat (theta, phi, lon, lat)
      ! Converts from theta (co-latitude, from N pole, in radians)
      ! and phi (longitude, E from Greenwich, in radians)
      ! to "lon" (East longitude, in degrees; West is negative)
      ! and "lat" (North latitude, in degrees; South is negative).
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN)  :: theta, phi
        DOUBLE PRECISION, INTENT(OUT) :: lon, lat
        lat = 90.0D0 - degrees_per_radian * DABS(theta)
        lat = MAX (lat, -90.0D0)
        lon = degrees_per_radian * phi
        IF (lon > 180.0D0) lon = lon - 360.0D0
        IF (lon <= -180.0D0) lon = lon + 360.0D0
    END SUBROUTINE DThetaPhi_2_LonLat

    SUBROUTINE DThetaPhi_2_Uvec (theta, phi, uvec)
      ! Converts from theta (co-latitude, from N pole) and
      ! phi (longitude, E from Greenwich) [both in radians]
      ! to a 3-component Cartesian unit vector, which points
      ! from the center of the unit sphere to a surface point.
      ! Its first axis outcrops at (0E, 0N).
      ! Its second axis outcrops at (90E, 0N).
      ! Its third axis outcrops at 90N.
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: theta, phi
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: uvec
        DOUBLE PRECISION :: equat
        uvec(3) = DCOS(theta)
        equat = DSIN(theta)
        uvec(1) = equat * DCOS(phi)
        uvec(2) = equat * DSIN(phi)
    END SUBROUTINE DThetaPhi_2_Uvec

    SUBROUTINE DTraceback ()
	use SharedVars
	use ShellSetSubs
        ! The sole function of this unit is to cause a traceable error,
        !   so that the route into the unit that called it is also traced.
        ! This unit is a good place to put a breakpoint while debugging!
        ! The intentional error must NOT be detected during compilation,
        !   but MUST cause a traceable error at run-time.
        ! If this routine has any error detected during compilation,
        !   then change its code to cause a different intentional error.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3) :: y


	    write(ErrorMsg,'(A/,A)') "Traceback was called to execute an intentional error:", &
        &                        "An array subscript will be intentionally out-of-range."
	    call FatalError(ErrorMsg,ThID)
        DO i = 1, 4
            y(i) = 1.0D0 * i
        END DO
    END SUBROUTINE DTraceback

    SUBROUTINE DTurn_To (azimuth_radians, base_uvec, far_radians, & ! inputs
                      &  omega_uvec, result_uvec)
      ! Computes uvec "result_uvec" (a 3-component Cartesian unit
      ! vector from the center of the planet) which results from
      ! rotating along a great circle beginning at "base_uvec"
      ! for an angle of "far_radians", in the initial direction
      ! given by azimuth_radians" (clockwise, from North).
      ! Also returned is "omega_uvec", the pole of rotation.
      ! NOTE: At the poles, azimuth is undefined.  Near the
      ! poles, it is defined but numerically unstable.  Therefore,
      ! Turn_To uses the same SUBROUTINE  DNorthEast_Convention as
      ! FUNCTION  DRelative_Compass does, so they can work together
      ! to find internal points on a small circle (as in DProcess_L4_
      ! Paths).
        IMPLICIT NONE
        DOUBLE PRECISION,                 INTENT(IN)  :: azimuth_radians, far_radians
        DOUBLE PRECISION, DIMENSION(3),   INTENT(IN)  :: base_uvec
        DOUBLE PRECISION, DIMENSION(3),   INTENT(OUT) :: omega_uvec, result_uvec
        DOUBLE PRECISION, PARAMETER :: pole_width = 1.745D-4 ! 0.01 degrees; must match DNortheast_Convention!
        DOUBLE PRECISION :: complement, cos_size, e_part, n_part, sin_size
        DOUBLE PRECISION, DIMENSION(3)   :: east_uvec, north_uvec, t_uvec
        DOUBLE PRECISION, DIMENSION(3,3) :: rotation_matrix
        CALL DNorthEast_Convention (base_uvec, north_uvec, east_uvec)
        e_part = -DCOS(azimuth_radians)
        n_part =  DSIN(azimuth_radians)
        omega_uvec(1) = e_part*east_uvec(1) + n_part*north_uvec(1)
        omega_uvec(2) = e_part*east_uvec(2) + n_part*north_uvec(2)
        omega_uvec(3) = e_part*east_uvec(3) + n_part*north_uvec(3)
        cos_size = DCOS(far_radians)
        sin_size = DSIN(far_radians)
        complement = 1.00D0 - cos_size
        rotation_matrix(1,1) = cos_size + complement*omega_uvec(1)*omega_uvec(1)
        rotation_matrix(1,2) = complement*omega_uvec(1)*omega_uvec(2) - sin_size*omega_uvec(3)
        rotation_matrix(1,3) = complement*omega_uvec(1)*omega_uvec(3) + sin_size*omega_uvec(2)
        rotation_matrix(2,1) = complement*omega_uvec(2)*omega_uvec(1) + sin_size*omega_uvec(3)
        rotation_matrix(2,2) = cos_size + complement*omega_uvec(2)*omega_uvec(2)
        rotation_matrix(2,3) = complement*omega_uvec(2)*omega_uvec(3) - sin_size*omega_uvec(1)
        rotation_matrix(3,1) = complement*omega_uvec(3)*omega_uvec(1) - sin_size*omega_uvec(2)
        rotation_matrix(3,2) = complement*omega_uvec(3)*omega_uvec(2) + sin_size*omega_uvec(1)
        rotation_matrix(3,3) = cos_size + complement*omega_uvec(3)*omega_uvec(3)
       !Copy base_uvec in case user of this routine plans to change it:
        t_uvec = base_uvec
        result_uvec(1) = rotation_matrix(1,1)*t_uvec(1) + &
                       & rotation_matrix(1,2)*t_uvec(2) + &
                       & rotation_matrix(1,3)*t_uvec(3)
        result_uvec(2) = rotation_matrix(2,1)*t_uvec(1) + &
                       & rotation_matrix(2,2)*t_uvec(2) + &
                       & rotation_matrix(2,3)*t_uvec(3)
        result_uvec(3) = rotation_matrix(3,1)*t_uvec(1) + &
                       & rotation_matrix(3,2)*t_uvec(2) + &
                       & rotation_matrix(3,3)*t_uvec(3)
    END SUBROUTINE DTurn_To

    SUBROUTINE DUvec_2_LonLat (uvec, lon, lat)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: uvec
        DOUBLE PRECISION, INTENT(OUT) :: lon, lat
        DOUBLE PRECISION :: theta, phi
        CALL DUvec_2_ThetaPhi (uvec, theta, phi)
        CALL DThetaPhi_2_LonLat (theta, phi, lon, lat)
    END SUBROUTINE DUvec_2_LonLat

    SUBROUTINE DUvec_2_ThetaPhi (uvec, theta, phi)
      ! converts from Cartesian unit vector to theta (colatitude)
      ! and phi (longitude), both in radians
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: uvec
        DOUBLE PRECISION, INTENT(OUT) :: theta, phi
        DOUBLE PRECISION :: equat, equat2
        equat2 = uvec(1)*uvec(1) + uvec(2)*uvec(2)
        IF (equat2 == 0.0D0) THEN
            phi = 0.0D0 ! actually undefined; provide default 0.0D0
            IF (uvec(3) > 0.0D0) THEN
                theta = 0.0D0 ! N pole
            ELSE
                theta = Pi    ! S pole
            END IF
        ELSE
            equat = DSQRT(equat2)
            theta = DATAN2(equat, uvec(3))
            phi = DATAN2(uvec(2), uvec(1))
        END IF
    END SUBROUTINE DUvec_2_ThetaPhi

    SUBROUTINE DUvec_2_PlungeAzimuth (location_uvec, lineation_uvec, & ! inputs
                                   &  plunge_degrees, azimuth_degrees) ! outputs
      ! converts Cartesian unit vector lineation_uvec
      ! which is at geographic location location_uvec
      ! into plunge toward azimuth, in degrees.
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: location_uvec, lineation_uvec
        INTEGER, INTENT(OUT) :: plunge_degrees, azimuth_degrees
        DOUBLE PRECISION :: equat, horizontal_component, phi_component, theta_component, up_component
        DOUBLE PRECISION, DIMENSION(3) :: phi_uvec, theta_uvec
        equat = DSQRT(location_uvec(1)*location_uvec(1) + location_uvec(2)*location_uvec(2))
        IF (equat == 0.0D0) THEN
            WRITE (*,"(' ERROR: location_uvec sent to DUvec_2_PlungeAzimuth may not be N or S pole.')")
            CALL DTraceback
        END IF
        CALL DLocal_Theta(location_uvec, theta_uvec)
        CALL DLocal_Phi  (location_uvec, phi_uvec)
        up_component    = DDot(lineation_uvec, location_uvec)
        theta_component = DDot(lineation_uvec, theta_uvec)
        phi_component   = DDot(lineation_uvec, phi_uvec)
        horizontal_component = DSQRT(theta_component**2 + phi_component**2)
        IF (horizontal_component <= 0.0D0) THEN
            plunge_degrees = 90
            azimuth_degrees = 0
        ELSE IF (up_component <= 0.0D0) THEN
            plunge_degrees = NINT(degrees_per_radian * DATAN2(-up_component, horizontal_component))
            azimuth_degrees = NINT(degrees_per_radian * DATAN2(phi_component, -theta_component))
            IF (azimuth_degrees < 0) azimuth_degrees = azimuth_degrees + 360
        ELSE ! up_component is positive
            plunge_degrees = NINT(degrees_per_radian * DATAN2(up_component, horizontal_component))
            azimuth_degrees = NINT(degrees_per_radian * DATAN2(-phi_component, theta_component))
            IF (azimuth_degrees < 0) azimuth_degrees = azimuth_degrees + 360
        END IF
    END SUBROUTINE DUvec_2_PlungeAzimuth

    DOUBLE PRECISION FUNCTION DDot_3D (a_, b_)
        ! Dot product of 3-component vectors.
        DOUBLE PRECISION, DIMENSION(3)  :: a_, b_
        DDot_3D = a_(1)*b_(1) + a_(2)*b_(2) + a_(3)*b_(3)
    END FUNCTION DDot_3D

end module DSphere !===============================================


module DDislocation
!
! Computes benchmark offsets/velocities due to shear dislocations
! from earthquakes/frictional slip in the brittle part of each fault,
! at the rate predicted by a finite-element model.
!
! Originally, exactly parallel to MODULE Dislocation.f90,
! except that all REAL variables were
! replaced with DOUBLE PRECISION (REAL*8) variables,
! and subprogram names had "D" added as a
! sign that they expect DOUBLE PRECISION
! arguments and return DOUBLE PRECISION answers.

! By Peter Bird, UCLA, gradually developed 1980-2015;
!                      major logic revision on 2002.08.13;
!                      REAL*8 added (for all REALs) 2015.02.12;
!                      minor numerical improvement on 2015.11.03.
! THEN, in 2020.09 there was a major upgrade to routines DChange, Halo, Aura,
!       so that fault-displacment (Burgers) vectors could have an opening component,
!       whose remote surface displacments are computed from the solution of
!       Okada [1985, BSSA].
! Note that this upgrade was NOT retroactively added to Dislocation.f90.
!
! Usage from a F-E program like Shells/OrbScore, with fault elements:
! make a single call for the entire F-E grid:
!
!     CALL DCoseis (farg,fdip,geothe,geophi,mxnode,   &
!    &              mxfel,numgeo,nfl,nodef,radius,    &
!    &              v,wedge,xnode,ynode,zmnode,ztranf,& ! inputs
!    &              geouth,geouph)                      ! outputs
!
!    For a typical year with no earthquakes,
!    correct the predicted benchmark velocities
!    by subtracting the coseismic part:
!
!     DO i = 1, numgeo
!          predic(1, i) = predic(1, i) - geouth(i)
!          predic(2, i) = predic(2, i) - geouph(i)
!     END DO

! - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

! Usage from a program (like NeoKinema) where the faults are not
! represented by elements: make many calls, one for each combination
! of fault segment and benchmark.  Call DChange directly:
!
!   benchmark_unlocked_vw = benchmark_vw ! whole array; initialize before sum
!   DO ..... fault segments ........
!       ....translate fault data into Change variables...
!       DO i = 1, internal_benchmarks
!           btheta = ...
!           bphi   = ...
!           CALL DChange (argume, &
!         &               btheta, bphi, &
!         &               dipf, lf, &
!         &               ftheta, fphi, &
!         &               radius, &
!         &               slip, &
!         &               wedge, &
!         &               ztop, zbot, &  ! inputs
!         &               duthet, duphi) ! output
!           For a typical year with no earthquakes,
!               correct the observed benchmark velocities
!               by adding the (estimated) missing coseismic part:
!           benchmark_unlocked_vw(2 * i - 1) = benchmark_unlocked_vw(2 * i - 1) + duthet
!           benchmark_unlocked_vw(2 * i    ) = benchmark_unlocked_vw(2 * i    ) + duphi
!       END DO ! i = 1, internal_benchmarks
!   END DO ! fault segments

!NOTE that this typical usage will NOT make use of the new [2020] mode-1 Okada [1985] sections.

!==================================================================

CONTAINS

      SUBROUTINE DCoseis (farg, fdip, geothe, geophi, &
     &                    mxnode, mxfel, numgeo, nfl, nodef, radius, &
     &                    v, wedge, xnode, ynode, zmnode, ztranf, & ! inputs
     &                    uTheta, uPhi) ! outputs

!    Computes the horizontal-plane velocity components
!   (uTheta = South, uPhi = East) for each of NUMGEO benchmark sites,
!    due to earthquake dislocations on the brittle parts of faults.

!   (The usual reason for doing this is so the coseismic part
!    of the benchmark velocities can be subtracted from the
!    velocities predicted by the F-E model before comparing to
!    geodetic data in an aseismic year.)

!    GEOTHE and GEOPHI should be the theta and phi coordinates of
!    the benchmarks, in radians.  Theta is colatitude, measured
!    Southward from the North pole.  Phi is latitude, measured
!    Eastward from the prime meridian.  The same convention
!    should be used for node positions in XNODE (= theta) and
!    YNODE (= phi); the use of x for theta and y for phi is
!    a memnonic convention adopted when flat-planet Plates
!    was rewritten as spherical-planet Shells
!    Similarly, the 1st (x, or theta) component of V should be
!    the Southward velocity, while the 2nd (y, or phi) component
!    of V should be the Eastward velocity at the nodes.

!    Most of the actual work is done by the CALL DChange
!    statement; this routine is only reponsible for cutting
!    each fault into several segments (for greater accuracy)
!    and summing the corrections due to each segment of
!    each fault.

!    The following parameter NSEGME determines how many segments
!    each fault element is divided into; higher values may
!    give more accuracy for benchmarks close to faults with
!    spatially-varying slip rates.

       IMPLICIT NONE
       INTEGER, PARAMETER :: nsegme = 3

!    Arguments (see text above):
       INTEGER,                    INTENT(IN)  :: mxfel, mxnode, nfl, numgeo
       INTEGER, DIMENSION(4, nfl), INTENT(IN)  :: nodef
       DOUBLE PRECISION, DIMENSION(2, mxnode), INTENT(IN) :: v
       REAL*8,                     INTENT(IN)  :: radius, wedge
       REAL*8, DIMENSION(2,mxfel), INTENT(IN)  :: farg, fdip, ztranf
       REAL*8, DIMENSION(mxnode),  INTENT(IN)  :: xnode, ynode, zmnode
       REAL*8, DIMENSION(numgeo),  INTENT(IN)  :: geothe, geophi
       REAL*8, DIMENSION(numgeo),  INTENT(OUT) :: uTheta, uPhi

       DOUBLE PRECISION :: lf
       INTEGER :: i, ilayer, j, m, n1, n2, n3, n4
       REAL*8 :: al, argume, &
     &           crossx, crossy, &
     &           delvx, delvy, dip, duthet, duphi, &
     &           fhival, height, moho, &
     &           s, smid, unitx, unity, &
     &           vx1, vx2, vx3, vx4, vy1, vy2, vy3, vy4, &
     &           x1, x2, x3, x4, xbegin, xend, xmid, &
     &           y1, y2, y3, y4, ybegin, yend, ymid, &
     &           zbot, ztop
       REAL*8, DIMENSION(3) :: slip
!---------------------------------------------------------------------

!                    STATEMENT FUNCTIONS:

!      Interpolation within one fault element:
       fhival(s, x1, x2) = x1 + s * (x2 - x1)

!---------------------------------------------------------------------

!     Initialize to zero before beginning sums:
       DO j = 1, numgeo
            uTheta(j) = 0.0D0
            uPhi(j) = 0.0D0
       END DO

!     Sum contributions from fault elements:
       DO 100 i = 1, nfl
            n1 = nodef(1, i)
            n2 = nodef(2, i)
            n3 = nodef(3, i)
            n4 = nodef(4, i)
            x1 = xnode(n1)
            x2 = xnode(n2)
            x3 = xnode(n3)
            x4 = xnode(n4)
            y1 = ynode(n1)
            y2 = ynode(n2)
            y3 = ynode(n3)
            y4 = ynode(n4)
            vx1 = v(1, n1)
            vx2 = v(1, n2)
            vx3 = v(1, n3)
            vx4 = v(1, n4)
            vy1 = v(2, n1)
            vy2 = v(2, n2)
            vy3 = v(2, n3)
            vy4 = v(2, n4)
            xbegin = x1
            ybegin = y1

!         Sum over nseg segments within each fault:
            DO 90 m = 1, nsegme

                 smid = ((m * 1.0D0) - 0.5D0) / (nsegme * 1.0D0)
                 CALL DOnarc(smid, x1, y1, x2, y2, & ! input
     &                       xmid, ymid) ! output

                 s = (m * 1.0D0) / (nsegme * 1.0D0)
                 IF(m == nsegme) THEN
                      xend = x2
                      yend = y2
                 ELSE
                      CALL DOnarc(s, x1, y1, x2, y2, & ! input
     &                            xend, yend) ! output
                 END IF
                 al = DFltLen(ybegin, yend, radius, xbegin, xend)
                 lf = 0.5D0 * al

!              FARG is the direction from node 1 to node 2,
!              in radians counterclockwise from the
!              +theta (South) direction, and so is ARGUME:
                 argume = DChord(farg(1, i), smid * 1.0D0, farg(2, i))

                 dip = fhival(smid, fdip(1, i), fdip(2, i))

!              We consider that the node-1,2 side of the fault
!              moves while the node-3,4 side is fixed:
                 delvx = fhival(smid, vx1, vx2) - fhival(smid, vx4, vx3)
                 delvy = fhival(smid, vy1, vy2) - fhival(smid, vy4, vy3)

!              UNIT is a horizontal unit vector pointing
!              along the fault from the node-1 end to the node-2 end;
!              this will be the same as the Z1 direction of HALO if
!              the fault is vertical, and equal or opposite to the
!              X1 direction of AURA if the fault is dipping:
                 unitx = DCOS(argume)
                 unity = DSIN(argume)

!              CROSS is a horizontal unit vector pointing
!              across the fault, 90 degrees clockwise from UNIT
!             (when viewed from outside the planet);
!              this will be the same as the Z2 direction of HALO
!              if the fault is exactly vertical or
!              equal or opposite to the X2 direction of AURA
!              if the fault is dipping.
                 crossx = + unity
                 crossy = -unitx

                 slip(1) = delvx * unitx + delvy * unity
!               Note: positive for sinistral slip

                 slip(2) = delvx * crossx + delvy * crossy
!               Note: positive for opening across trace

!               Vertical component is positive when down:
                 IF (ABS(1.570796D0 - dip) > wedge) THEN
                      slip(3) = slip(2) * DTAN(dip)
                 ELSE
                      slip(3) = 0.0D0
                 END IF

!                Consider both crustal and mantle patches:
                 moho = fhival(smid, zmnode(n1), zmnode(n2))
                 DO 80 ilayer = 1, 2

                      height = ztranf(ilayer, i)
                      IF (height > 0.0D0) THEN
                           IF (ilayer == 1) THEN
                               ztop = 0.0D0
                           ELSE
                               ztop = moho
                           END IF
                           zbot = ztop + height

!                          Compute effects of patch at all benchmarks:
                           DO 70 j = 1, numgeo
!                            -------------------------------------------
                                CALL DChange(argume, &
     &                                       geothe(j), geophi(j), &
     &                                       dip, lf, &
     &                                       xmid, ymid, &
     &                                       radius, &
     &                                       slip, &
     &                                       wedge, &
     &                                       ztop, zbot, &  ! inputs
     &                                       duthet, duphi) ! output
!                            -------------------------------------------
                                uTheta(j) = uTheta(j) + duthet
                                uPhi(j) = uPhi(j) + duphi
   70                      CONTINUE
                      END IF

   80            CONTINUE

!                Prepare to loop to next segment:
                 xbegin = xend
                 ybegin = yend
   90       CONTINUE
  100  CONTINUE
       END SUBROUTINE DCoseis
!-----------------------------------------------------------------------

       SUBROUTINE DChange (argume, &
     &                     btheta, bphi, &
     &                     dipf, lf, &
     &                     ftheta, fphi, &
     &                     radius, &
     &                     slip, &
     &                     wedge, &
     &                     ztop, zbot, &  ! inputs
     &                     duthet, duphi, dur) ! output

!   This routine is a "driver" or "wrapper" for the Cartesian (flat-Earth)
!   *Mansinha & Smylie routine Halo (vertical shear-dislocation patch)
!                           or Aura (dipping  shear-dislocation patch)
!   *and tensile (mode-1; opening) crack routine of Okada [1985, BSSA],
!   which permits them to be called from spherical-planet
!   programs because it converts coordinates from
!   (theta, phi) = (colatitude, longitude) in radians
!   to the (z1,z2,z3) Cartesian units of a locally flat planet.
!   Also, output (horizontal) vector (DUTHET, DUPHI) is expressed
!   in the spherical-planet (theta, phi) system, so components are
!  (South, East).  Units are the same as input vector SLIP,
!   which may be either a relative displacement or
!   a relative velocity of the two sides of the fault.

!   The third component of the output displacement[-rate] vector,
!   DUR (radial, or up) is OPTIONAL.

!   FTHETA and FPHI should give (theta, phi) coordinates of the
!   midpoint of the trace of the plane containing the
!   rectangular dislocation patch.

!   BTHETA and BPHI should give (theta, phi) coordinates of the
!   benchmark or test point at which the displacement (-rate)
!   is desired.

!   LF is the half-length (from center to end) of the dislocation
!   patch, measured along a horizontal strike line.
!   Units are the same as RADIUS, ZTOP, ZBOT (below).

!   ARGUME is the argument of the trace of the dislocation patch,
!   measured in radians counterclockwise from +theta (from South).

!   DIPF is fault dip in radians measured clockwise
!   (initially, down) from horizontal on the right side of the
!   fault (when looking in direction ARGUME).

!   The fault/crack dislocation vector (Burgers vector)
!   SLIP must be already rotated into fault-trace-centered coordinates:
!   SLIP(1) is the component parallel to the fault,
!      and it is positive for sinistral sense of slip.
!   SLIP(2) is the horizontal component in the direction
!      perpendicular to the trace of the fault,
!      and it is positive for divergence across the trace.
!   SLIP(3) is the relative vertical component, and it
!      is positive when the right side of the fault
!     (when looking along direction ARGUME) is down.
!   Note that SLIP(2) > 0 with SLIP(3) > 0 is associated with normal faulting;
!       while SLIP(2) < 0 with SLIP(3) < 0 is associated with thrust faulting;
!         and SLIP(2) > 0 with SLIP(3) < 0 is associated with tensile (mode-1) opening.

!   WEDGE is a tolerance (in radians) for dip; if DIPF is
!   within WEDGE of Pi/2, then fault is considered vertical
!   and routine Halo is called; otherwise, Aura is called.

!   ZTOP and ZBOT are (positive) depths to top and bottom of
!   the dislocation patch, respectively.  Units are the same
!   as RADIUS, which gives the size of the planet.
!   Note that a 0.0 input value will be accepted for ZTOP,
!   but in the actual calculation a small positive depth
!   will be substituted for this zero.

!   NOTE: If distance from "F" point to "B" point exceeds
!         one planetary radius, result is approximated as zero and
!         Halo/Aura are never called.

       IMPLICIT NONE
!     -------------------------------------------------------
!   Arguments (see text above):
       REAL*8, INTENT(IN) :: argume, btheta, bphi, dipf, &
                           & ftheta, fphi, radius, wedge, zbot, ztop

       REAL*8, DIMENSION(3), INTENT(IN) :: slip

       REAL*8, INTENT(OUT) :: duthet, duphi ! S and E components
       REAL*8, INTENT(OUT), OPTIONAL :: dur ! radial or up component
!     -------------------------------------------------------

!      Following internal variables must be DP to agree with Halo, Aura:
       DOUBLE PRECISION :: dbot, dtop, lf, theta
       DOUBLE PRECISION, DIMENSION(3) :: u, ucap, x, z

!      Following are DP for precision in determination of
!      variables FAR and ALPHA.  (At short distances,
!      unit vectors must be DP or location precision is
!      lost and benchmark can end up on wrong side of fault!)
       DOUBLE PRECISION :: crossp, dot1, dot2, dotpro
       DOUBLE PRECISION, DIMENSION(3) :: buvec, dvec, fuvec, tvec, uPhi, uTheta

       INTEGER :: i
       LOGICAL :: vert
       REAL*8 :: alpha, cosaz, far, sinaz, tazimf

!     Compute position relative to fault trace;
!     this is where we convert from a spherical-Earth
!     to a flat-Earth model, implicitly using a gnomonic
!     projection (which preserves distance and azimuth
!     from the center of the trace of the dislocation plane).

       fuvec(1) = DSIN(ftheta) * DCOS(fphi)
       fuvec(2) = DSIN(ftheta) * DSIN(fphi)
       fuvec(3) = DCOS(ftheta)
       buvec(1) = DSIN(btheta) * DCOS(bphi)
       buvec(2) = DSIN(btheta) * DSIN(bphi)
       buvec(3) = DCOS(btheta)
       dotpro = fuvec(1) * buvec(1) + fuvec(2) * buvec(2) + fuvec(3) * buvec(3)
       tvec(1) = fuvec(2) * buvec(3) - fuvec(3) * buvec(2)
       tvec(2) = fuvec(3) * buvec(1) - fuvec(1) * buvec(3)
       tvec(3) = fuvec(1) * buvec(2) - fuvec(2) * buvec(1)
       crossp = DSQRT(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)
       far = radius * DATAN2(crossp, dotpro)
       !
       !Notice of important change:
       !Prior to October 2015 the test condition was:
       !IF (far > (40.0D0 * zbot)) THEN
       !However, while plotting some very high-resolution (0.01 degree x 0.01 degree)
       !models of coseismic strains, created with StrainRates2015, I noticed that
       !this gives an objectionable "cutoff" artifact along a circular arc.
       !Thus, the outer radius for application of these solutions needs to be much bigger.
       !On the other hand, these flat-Earth solutions should certainly NOT be trusted
       !in the opposite hemisphere!  As a compromise, I now cut them off at an arc
       !distance of one radian from the center of the trace of the dislocation patch:
       !The result will NOT be accurate in the far-field, but at least it will be
       !smoothly going to zero, so that it does not contaminate computed strain-rates.
       IF (far > radius) THEN
            duthet = 0.0D0
            duphi = 0.0D0
       ELSE
            dvec(1) = buvec(1) - fuvec(1)
            dvec(2) = buvec(2) - fuvec(2)
            dvec(3) = buvec(3) - fuvec(3)
            uTheta(1) = DCOS(ftheta) * DCOS(fphi)
            uTheta(2) = DCOS(ftheta) * DSIN(fphi)
            uTheta(3) = -DSIN(ftheta)
            dot1 = dvec(1) * uTheta(1) + dvec(2) * uTheta(2) + dvec(3) * uTheta(3)
            uPhi(1) = -DSIN(fphi)
            uPhi(2) = DCOS(fphi)
            uPhi(3) = 0.0D0
            dot2 = dvec(1) * uPhi(1) + dvec(2) * uPhi(2) + dvec(3) * uPhi(3)
            alpha = DATAN2(dot2, dot1)
            z(1) = far * DCOS(alpha - argume)
            z(2) = -far * DSIN(alpha - argume)

!           Note: Z(3) = 0 because we believe it is more important
!           to convey that the benchmark is on the free surface
!           than it is to convey the precise angular relationship
!           to the dislocation.  However, some day one might test
!           a positive Z(3), in proportion to square of "far".
            z(3) = 0.0D0

            vert = ABS(dipf - 1.5708) <= wedge
            IF (vert) THEN ! vertical fault or crack, so use Halo:

                 tazimf = argume
                 theta = dipf
!                Z values computed above are still good
                 ucap(1:3) = slip(1:3)
                 dtop = ztop
                 dbot = zbot
!               Kludge to insure that input parameter DTOP is not zero,
!               because, in that case, HALO may give numerically-unstable
!               ("bad") results for near-field test points:
                 dtop = MAX(dtop, 0.0001D0 * dbot)
!              - - - - - - - - - - - - - - - - - - - - - - -
                 CALL Halo (ucap, z, lf, dbot, dtop, & ! inputs
     &                      u) ! output
!              - - - - - - - - - - - - - - - - - - - - - - -

            ELSE  ! non-vertical, dipping fault, so use Aura:

                 IF (dipf <= 1.5708D0) THEN
                      tazimf = argume
                      theta = dipf
                      DO 10 i = 1, 3
                           x(i) = z(i)
                           ucap(i) = slip(i)
   10                 CONTINUE
                 ELSE
!                 Look at this fault from the other end/side,
!                 so that dip THETA will be less than Pi/2:
                      tazimf = argume + 3.14159265358979D0
                      theta = 3.14159265358979D0 - dipf
                      x(1) = -z(1)
                      x(2) = -z(2)
                      x(3) = z(3)
                      ucap(1) = slip(1)
                      ucap(2) = slip(2)
                      ucap(3) = -slip(3)
!                     Note: Reversal of Z1, Z2 axes (but not Z3)
!                     during name-change to X1,X2,X3 axes (for AURA)
!                     combined with change of "moving" block
!                     leaves "sinistral" and "opening" components
!                     unchanged, while "relative vertical" component
!                     reverses.
                 END IF
                 dtop = ztop / DSIN(theta)
                 dbot = zbot / DSIN(theta)
!               Kludge to insure that input parameter DTOP is not zero,
!               because, in that case, AURA may give numerically-unstable
!               ("bad") results for near-field test points:
                 dtop = MAX(dtop, 0.0001D0 * dbot)
!              - - - - - - - - - - - - - - - - - - - - - - -
                 CALL Aura (ucap, theta, x, lf, dbot, dtop, & ! inputs
     &                      u) ! output
!              - - - - - - - - - - - - - - - - - - - - - - -
            END IF

            sinaz = DSIN(tazimf)
            cosaz = DCOS(tazimf)
            duthet = + cosaz * u(1) + sinaz * u(2)
            duphi = + sinaz * u(1) - cosaz * u(2)
            IF (PRESENT(dur)) THEN
                dur = -u(3) ! changing sign because both Halo and Aura use
                            ! the Mansinha & Smylie positive-downwards convention,
                            ! but this routine uses a positive-upwards convention.
            END IF
       END IF
       END SUBROUTINE DChange
!-----------------------------------------------------------------

       SUBROUTINE Halo (ucap, z, l, dbot, dtop, & ! inputs
     &                  u) ! output

!   Computes displacements in a uniform elastic half-space
!   of Poisson's ratio 0.25 due to uniform slip (and/or opening) on a
!   rectangular dislocation patch in a vertical plane,
!   whose horizontal and vertical edges form a rectangle.

!   Equations for fault-shear components based on Mansinha and Smylie (1967),
!   J. Geophys. Res., v. 72, p. 4731-4743.

!   Equations for opening/closing (mode-1) components based on Okada [1985],
!   Bull. Seismol. Soc. Am., 74(4), 1135-1154.

       USE DSphere ! provided by Peter Bird as file DSphere.f90
       IMPLICIT NONE

!   Arguments:
       DOUBLE PRECISION,               INTENT(IN)  :: l, dbot, dtop
       DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: ucap, z
       DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: u

       INTEGER :: i
       DOUBLE PRECISION :: d_Okada, del_Okada, dTilde_Okada, epsilon_Okada, factor, I1_Okada, I3_Okada, I5_Okada, l_Okada, &
                         & mu_over_lambda_plus_mu, nu_Okada, p_Okada, q, q_Okada, r, R_Okada, R2_Okada, &
                         & term1, term2, term3, term4, term5, term6, u3_Okada, &
                         & ux_Okada, uy_Okada, uz_Okada, &
                         & w_Okada, x_Okada, XCap_Okada, XCap2_Okada, y_Okada, yTilde_Okada, z_Okada, &
                         & zeta1, zeta3

!   The coordinate system is Cartesian and right-handed,
!   with z3 down, z1 parallel to the fault/crack trace, and z2
!   perpendicular.  The dislocation surface is at z2 = 0,
!   DTOP <= z3 <= DBOT (note that DTOP > 0.), and
!   -L <= z1 <= L.
!   UCAP contains the three components of the displacement of the
!   positive-z2 side relative to the other side;
!   U contains the three components of the displacement
!   at the observation point Z.

       u = 0.0D0 ! all 3 components initialized
       DO 100 i = 1, 4
           factor = 1.0D0
           IF ((i == 2).OR.(i == 3)) factor = -1.0D0
           zeta1 = l
           IF ((i == 3).OR.(i == 4)) zeta1 = -l
           zeta3 = dbot
           IF ((i == 2).OR.(i == 4)) zeta3 = dtop
           r = DSQRT((z(1) - zeta1)**2 + z(2)**2 + (z(3) - zeta3)**2)
           q = DSQRT((z(1) - zeta1)**2 + z(2)**2 + (z(3) + zeta3)**2)
           term1 = z(2) * (z(1) - zeta1) * (2.0D0 / (r * (r + z(3) - zeta3)) &
             &   - (5.0D0 * q + 8.0D0 * zeta3) / (2.0D0 * q * (q + z(3) + zeta3)**2) &
             &   + 4.0D0 * z(3) * zeta3 * (2.0D0 * q + z(3) + zeta3) / &
             &   (q**3 * (q + z(3) + zeta3)**2)) &
             &   + 3.0D0 * DAtanPV((z(1) - zeta1) * (z(3) - zeta3), (r * z(2))) &
             &   - 3.0D0 * DAtanPV((z(1) - zeta1) * (z(3) + zeta3), (q * z(2)))

           term2 = -DLOG(r + z(3) - zeta3) &
             &   + 0.5D0 * DLOG(q + z(3) + zeta3) &
             &   - (5.0D0 * z(3) - 3.0D0 * zeta3) / (2.0D0 * (q + z(3) + zeta3)) &
             &   - 4.0D0 * z(3) * zeta3 / (q * (q + z(3) + zeta3)) &
             &   + z(2)**2 * (2.0D0 / (r * (r + z(3) - zeta3)) &
             &   - (5.0D0 * q + 8.0D0 * zeta3) / (2.0D0 * q * (q + z(3) + zeta3)**2) &
             &   + 4.0D0 * z(3) * zeta3 * (2.0D0 * q + z(3) + zeta3) / &
             &   (q**3 * (q + z(3) + zeta3)**2))

            TERM3=Z(2)*(2.0D0/R-2.0D0/Q+4.0D0*Z(3)*ZETA3/Q**3 &
        &   +3.0D0/(Q+Z(3)+ZETA3)+2.0D0*(Z(03)+3.0D0*ZETA3)/ &
        &   (Q*(Q+Z(3)+ZETA3)))

            TERM4=Z(2)*(2.0D0/R+4.0D0/Q-4.0D0*Z(3)*ZETA3/Q**3 &
        &   -6.0D0*Z(3)/(Q*(Q+Z(3)+ZETA3)))

            TERM5= -DLOG(R+Z(1)-ZETA1)+DLOG(Q+Z(1)-ZETA1) &
        &   +(6.0D0*Z(3)**2+10.0D0*Z(3)*ZETA3)/(Q*(Q+Z(1)-ZETA1)) &
        &   +(6.0D0*Z(3)*(Z(1)-ZETA1))/(Q*(Q+Z(3)+ZETA3)) &
        &   +Z(2)**2*(2.0D0/(R*(R+Z(1)-ZETA1)) &
        &   +4.0D0/(Q*(Q+Z(1)-ZETA1)) &
        &   -4.0D0*Z(3)*ZETA3*(2.0D0*Q+Z(1)-ZETA1)/(Q**3*(Q+Z(1)-ZETA1)**2))

            TERM6=Z(2)*(2.0D0*(Z(3)-ZETA3)/(R*(R+Z(1)-ZETA1)) &
        &   -2.0D0*(Z(3)+2.0D0*ZETA3)/(Q*(Q+Z(1)-ZETA1)) &
        &   -4.0D0*Z(3)*ZETA3*(Z(3)+ZETA3)*(2.0D0*Q+Z(1)-ZETA1)/ &
        &   (Q**3*(Q+Z(1)-ZETA1)**2)) &
        &   +3.0D0*DATANPV((Z(1)-ZETA1)*(Z(3)-ZETA3),(R*Z(2))) &
        &   -3.0D0*DATANPV((Z(1)-ZETA1)*(Z(3)+ZETA3),(Q*Z(2)))

           !Begin with strike-slip component (most common application):
           u(1) = u(1) + factor * term1 * ucap(1) / 37.69911184D0
           u(2) = u(2) + factor * term2 * ucap(1) / 37.69911184D0
           U(3)=U(3)+FACTOR*TERM3*UCAP(1)/37.69911184D0

           !Is there any vertical slip on this dislocation?
           !NOTE that the vertical component U(3) is positive-downwards.
           IF (ucap(3) /= 0.0D0) THEN
               U(1)=U(1)+FACTOR*TERM4*UCAP(3)/37.69911184D0
               U(2)=U(2)+FACTOR*TERM5*UCAP(3)/37.69911184D0
               U(3)=U(3)+FACTOR*TERM6*UCAP(3)/37.69911184D0
           END IF
  100  CONTINUE
       !------end of Mansinha & Smylie [1967] section---------------------------

       !------beginning of Okada [1985, BSSA] section --------------------------
           !Is there any opening/closing (mode-1 component) on this dislocation?
           IF (ucap(2) /= 0.0D0) THEN ! tensile / mode-1 / opening component:

              !Note that Okada uses a right-handed Cartesian coordinate system (like M&S)
              !but that his has z increasing UP, and therefore if we make his x axis (along-strike)
              !equal to the z1 axis of M&S, then Okada's y axis is opposite to the z2 axis of M&S.
              !Another distinction is that M&S have their (z1, z2) origin at a central point on
              !the fault trace (projected to the free surface), but Okada's (x, y) origin
              !is located over the deepest part of the opening patch.

              !Dip of crack/fault:
               del_Okada = Pi_over_2 ! defined in DSphere.f90 (== 90 degrees; only here in Halo)

              !Observer (test) location [which must be on the surface, using formulas of Okada, 1985]:
               x_Okada =  z(1)
               y_Okada = -z(2)
              !N.B. Fortunately, for a vertical patch, the formulas above do NOT require any offset terms
              !     between M&S's (z1, z2) origin point and Okada's (x, y) origin point.
              !     This would not be true for other dips of the fault patch!
               z_Okada =  0.0D0 ! (On free surface. To change this, use formulas of Okada, 1992 instead of Okada, 1985.)

              !Dimensions of dislocation patch:
               d_Okada = dbot
               w_Okada = dbot - dtop ! should be positive.  Note formula would be more complex for non-vertical dip.
               L_Okada = l ! half-length, along strike (same as M&S)

              !Okada dimensional parameters dependent on the above:
               p_Okada = d_Okada ! special simple form for vertical dip, when del_Okada = Pi_over_2
               q_Okada = y_Okada ! only for vertical dip

              !Amount of opening {rate?}:
               u3_Okada = ucap(2)

              !Poisson's ratio term:
               mu_over_lambda_plus_mu = 0.50D0 ! (assumed here in Halo; equivalent to Poisson's ratio of 1/4)

              !Initialize before sum of 4 terms:
               ux_Okada = 0.0D0
               uy_Okada = 0.0D0
               uz_Okada = 0.0D0

               DO 200 i = 1, 4 ! sum over 4 terms (from rectangular integration limits in 2-D):

                  !Integration limits and associated sign-changing factor:
                   factor = 1.0D0 ! unless...
                   IF ((i == 2).OR.(i == 3)) factor = -1.0D0 ! replacing the choice made in the line above
                   epsilon_Okada = x_Okada + L_Okada ! unless...
                   IF (i >= 3) epsilon_Okada = x_Okada - L_Okada ! replacing the choice made in the line above
                   nu_Okada = p_Okada ! unless...
                   IF ((i == 2).OR.(i == 4)) nu_Okada = p_Okada - w_Okada ! replacing the choice made in the line above

                  !Geometric factors dependent on the above:
                   yTilde_Okada = q_Okada ! only for vertical dip
                   dTilde_Okada = nu_Okada
                   R2_Okada = epsilon_Okada**2 + yTilde_Okada**2 + dTilde_Okada**2
                   R_Okada = SQRT(R2_Okada)
                   XCap2_Okada = epsilon_Okada**2 + q_Okada**2
                   XCap_Okada = SQRT(XCap2_Okada)

                   I1_Okada = -0.5D0 * mu_over_lambda_plus_mu * epsilon_Okada*q_Okada/(R_Okada+dTilde_Okada)**2
                   I3_Okada = +0.5D0 * mu_over_lambda_plus_mu * ( (nu_Okada/(R_Okada+dTilde_Okada)) + &
                                                              &   (yTilde_Okada*q_Okada/(R_Okada+dTilde_Okada)**2) - &
                                                              &   LOG(R_Okada+nu_Okada) )
                   I5_Okada = -mu_over_lambda_plus_mu * epsilon_Okada/(R_Okada+dTilde_Okada)

                   ux_Okada = ux_Okada + factor * (u3_Okada/Two_Pi) * ( (q_Okada**2/(R_Okada*(R_Okada+nu_Okada))) - &
                                                                    &   I3_Okada )
                   uy_Okada = uy_Okada + factor * (u3_Okada/Two_Pi) * ( (-dTilde_Okada*q_Okada/(R_Okada*(R_Okada+epsilon_Okada))) - &
                                                                    &   ((epsilon_Okada*q_Okada/(R_Okada*(R_Okada+nu_Okada))) - (ATAN((epsilon_Okada*nu_Okada) / (q_Okada*R_Okada)))) - &
                                                                    &   I1_Okada )
                   uz_Okada = uz_Okada + factor * (u3_Okada/Two_Pi) * ( (yTilde_Okada*q_Okada/(R_Okada*(R_Okada+epsilon_Okada))) - &
                                                                    &   I5_Okada )
200            CONTINUE ! sum over 4 terms

              !Final contribution (reversing signs of two terms, per reversal of axis-directions between M&S versus O):
               u(1) = u(1) + ux_Okada
               u(2) = u(2) - uy_Okada
               u(3) = u(3) - uz_Okada
           END IF ! there is a tensile / mode-1 / opening component of source motion
       !------end of Okada [1985, BSSA] section --------------------------------

       END SUBROUTINE Halo
!-----------------------------------------------------------------

       SUBROUTINE Aura (u4, theta4, x4, l4, dbot4, dtop4, & ! inputs
                      & v4) ! output

!   Computes displacements in an infinite uniform elastic halfspace
!   of Poisson's ratio 0.25 caused by uniform slip over the surface
!   of a buried rectangular dislocation with horizontal and plunging
!   sides.  To calculate displacements caused by more general slip,
!   use this routine many times in the kernel of an integral,
!   assigning the local average slip over each very small rectangle.

!   Equations for shear-displactions based on Mansinha and Smylie (1971),
!   Bull. Seism. Soc. Amer., v. 61, p. 1433-1440.
!   This code has been extended [2020.09] to also model tensile/mode-1/opening
!   components of crack dislocations, according to the solution of
!   Okada [1985], Bull. Seismol. Soc. Amer., v. 75, #4, p. 1135-1154.

       USE DSphere ! provided by Peter Bird as file DSphere.f90
       IMPLICIT NONE

!  Arguments:
       DOUBLE PRECISION,               INTENT(IN)  :: dbot4, dtop4, l4, theta4
       DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: u4, x4
       DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: v4

!  Internal variables:
       INTEGER :: i
       DOUBLE PRECISION :: costh, d_Okada, dbot, del_Okada, deld, dTilde_Okada, dtop, &
                         & epsilon_Okada, factor, h, &
                         & I1_Okada, I3_Okada, I4_Okada, I5_Okada, & ! N.B. I2_Okada (for strike-slip) is not needed.
                         & k, l, L_Okada, &
                         & mu_over_lambda_plus_mu, &
                         & nu_Okada, p_Okada, q, q_Okada, q2, q3, qsquar, &
                         & r, R_Okada, r2, R2_Okada, r3, rsquar, &
                         & secth, sinth, tanth, &
                         & term1, term2, term3, term4, term5, term6, theta, &
                         & ucap, ucap1, u3_Okada, ux_Okada, uy_Okada, uz_Okada, w_Okada, &
                         & x_Okada, xCap_Okada, xCap2_Okada, xi, xi1, xi2, xi3, y_Okada, yTilde_Okada, z_Okada
       DOUBLE PRECISION, DIMENSION(3) :: u(3), v(3), x(3)

      !Copy all input values to new variables, to avoid inadvertent changes!
       DO i = 1, 3
            u(i) = u4(i)
            x(i) = x4(i)
       END DO
       theta = theta4
       l = l4
       dbot = dbot4
       dtop = dtop4

!  U(3) is the dislocation vector (i.e. the total offset vector,
!  sometimes called the Burgers vector),
!  in a Cartesian coordinate system where x1 is parallel to the
!  fault strike, x2 is also horizontal and pointing 90 degrees to the
!  right of +x1 (that is, pointing to the hanging wall,
!  assuming that dip is less than Pi/2 radians!),
!  and x3 is down.  (Take care to make these right-handed.)
!  THETA is the fault dip angle (in radians) from +x2 axis.
!  X(3) gives the observation point in these coordinates.
!  L is the half-length of the rectangle, along strike.
!  DBOT and DTOP are oblique distances in the fault plane (and
!  down the dip) from the surface to the bottom and top of the
!  rectangular dislocation, respectively.
!  Note: Fault does not break the surface, so DTOP > 0.,
!  although it may be small.  If 0.0 is supplied for DTOP4,
!  a small positive depth will be substituted in the calculation.
!  V(3) is the answer: the 3-component displacement at the
!  observation point in same units as U and in x1,x2,x3 coordinates.

      !Fault/crack dip section:
       sinth = SIN(theta)
       costh = COS(theta)
       tanth = sinth / costh
       secth = 1.0D0 / costh
!  NOTE: FAULT MAY NOT BE VERTICAL IN THIS ROUTINE.
!       (To handle vertically-dipping faults, CALL Halo instead.)

       !Dislocation (Burgers) vector section:
       ! ucap1 is the strike-slip component, with sinistral (left-lateral) considered positive
        ucap1 = u(1)
       !ucap is the dip-slip (shear) component, with normal-faulting considered positive:
        ucap = costh * u(2) + sinth * u(3)
       !u3_Okada (not used until the 2nd half of this code, in the Okada [1985] section)
       !   is the tensile / mode-1 / crack-opening component.
       u3_Okada = sinth * u(2) - costh * u(3)
       !  (Of course, crack-closing can be represented by negative crack-opening.)
       !NOTE that Aura would typically execute faster if some of these components
       !     were exactly zero, so that code-sections could be skipped.
       !     However, that rarely happens in REAL*8 floating-point computations!
       !     Extra LOGICAL arguments could be added to control this.

       ! Initialize remote displacement {rate} to zero, before adding terms.
       DO i = 1, 3
           v(i) = 0.0D0
       END DO
       ! Note that v(1:3) follows Mansinha & Smylie [1971] sign-conventions,
       ! with downwards displacement considered positive.

      !---------- begin Mansinha & Smylie [1971] component (shearing only) -----------------
       r2 = x(2) * sinth - x(3) * costh
       r3 = x(2) * costh + x(3) * sinth
       q2 = x(2) * sinth + x(3) * costh
       q3 = -x(2) * costh + x(3) * sinth
       DO 100 i = 1, 4
           factor = 1.0D0
           IF ((i == 2).OR.(i == 3)) factor = -1.0D0
           xi1 = l
           IF (i >= 3) xi1 = -l
           xi = dbot
           IF ((i == 2).OR.(i == 4)) xi = dtop
           xi2 = xi * costh
           xi3 = xi * sinth
           rsquar = (x(1) - xi1)**2 + (x(2) - xi2)**2 + (x(3) - xi3)**2
           qsquar = (x(1) - xi1)**2 + (x(2) - xi2)**2 + (x(3) + xi3)**2
           k = DSQRT(qsquar - (q3 + xi)**2)
           h = DSQRT(qsquar - (x(1) - xi1)**2)
           q = DSQRT(qsquar)
           r = DSQRT(rsquar)

           !Impose arbitrary limit on how close test point can be to the (extended) fault plane;
           !without this limit the equations below may give 0./0. indeterminacy and an ABEND.
           !The following version of the limit was historically used in Dislocation.f90 through 2015,
           !q = MAX(q, 1.00001D0 * ABS(x(1) - xi1))
           !r = MAX(r, 1.00001D0 * ABS(x(1) - xi1))
           !but on 2015.11.03 I found that the second version is preferable, for both
           !Dislocation.f90 (which uses DOUBLE PRECISION in this routine) and for DDislocation.f90:
           q = MAX(q, 1.0000000001D0 * ABS(x(1) - xi1))
           r = MAX(r, 1.0000000001D0 * ABS(x(1) - xi1))

           term1 = (x(1) - xi1) * (2.0D0 * r2 / (r * (r + r3 - xi)) &
         &   - (4.0D0 * q2 - 2.0D0 * x(3) * costh) / (q * (q + q3 + xi)) &
         &   - 3.0D0 * tanth / (q + x(3) + xi3) &
         &   + 4.0D0 * q2 * x(3) * sinth / q**3 &
         &   - 4.0D0 * q2 * q3 * x(3) * sinth * (2.0D0 * q + q3 + xi) / &
         &   (q**3 * (q + q3 + xi)**2)) &
         &   - 6.0D0 * tanth**2 * DAtanPV(((k - q2 * costh) * (q - k) + &
         &   (q3 + xi) * k * sinth), ((x(1) - xi1) * (q3 + xi) * costh)) &
         &   + 3.0D0 * DAtanPV((x(1) - xi1) * (r3 - xi), (r2 * r)) &
         &   - 3.0D0 * DAtanPV((x(1) - xi1) * (q3 + xi), (q2 * q))

           term2 = sinth * (3.0D0 * tanth * secth * Dlog(q + x(3) + xi3) - &
         &   Dlog(r + r3 - xi) - (1.0D0 + 3.0D0 * tanth**2) * Dlog &
         &   (q + q3 + xi)) + 2.0D0 * r2**2 * sinth / (r * (r + r3 - xi)) &
         &   + 2.0D0 * r2 * costh / r - 2.0D0 * sinth * &
         &   (2.0D0 * x(3) * (q2 * costh - q3 * sinth) + q2 * (q2 + x(2) * sinth)) &
         &   / (q * (q + q3 + xi)) - 3.0D0 * tanth * (x(2) - xi2) / &
         &   (q + x(3) + xi3) + 2.0D0 * (q2 * costh - q3 * sinth - &
         &   x(3) * sinth**2) / q + 4.0D0 * q2 * x(3) * sinth * &
         &   ((x(2) - xi2) + q3 * costh) / q**3 - 4.0D0 * q2**2 * q3 * x(3) &
         &   * sinth**2 * (2.0D0 * q + q3 + xi) / (q**3 * (q + q3 + xi)**2)

           TERM3=COSTH*(DLOG(R+R3-XI)+(1.0D0+3.0D0*TANTH**2)*DLOG( &
         &   Q+Q3+XI)-3.0D0*TANTH*SECTH*DLOG(Q+X(3)+XI3)) &
         &   +2.0D0*R2*SINTH/R+2.0D0*SINTH*(Q2+X(2)*SINTH)/Q &
         &   -2.0D0*R2**2*COSTH/(R*(R+R3-XI))+ &
         &   (4.0D0*Q2*X(3)*SINTH**2-2.0D0*(Q2+X(2)*SINTH)*(X(3)+Q3*SINTH))/ &
         &   (Q*(Q+Q3+XI))+ &
         &   4.0D0*Q2*X(3)*SINTH*(X(3)+XI3-Q3*SINTH)/Q**3 &
         &   -4.0D0*Q2**2*Q3*X(3)*COSTH*SINTH*(2.0D0*Q+Q3+XI)/ &
         &   (Q**3*(Q+Q3+XI)**2)

           term4 = (x(2) - xi2) * sinth * (2.0D0 / r + 4.0D0 / q - 4.0D0 * xi3 * x(3) / q**3 &
         &   - 3.0D0 / (q + x(3) + xi3)) - costh * (3.0D0 * Dlog(q + x(3) + xi3) &
         &   + 2.0D0 * (x(3) - xi3) / r + 4.0D0 * (x(3) - xi3) / q + 4.0D0 * xi3 * x(3) * (x(3) + xi3) &
         &   / q**3) + 3.0D0 * (Dlog(q + x(3) + xi3) - sinth * Dlog(q + &
         &   q3 + xi)) / costh + 6.0D0 * x(3) * (costh / q - q2 * sinth / &
         &   (q * (q + q3 + xi)))

         term5 = sinth * (-Dlog(r + x(1) - xi1) + Dlog(q + x(1) - xi1) &
         &   + 4.0D0 * xi3 * x(3) / (q * (q + x(1) - xi1)) + 3.0D0 * (x(1) - xi1) / &
         &   (q + x(3) + xi3) + (x(2) - xi2)**2 * (2.0D0 / (r * (r + x(1) - xi1)) &
         &   + 4.0D0 / (q * (q + x(1) - xi1)) - 4.0D0 * xi3 * x(3) * (2.0D0 * q + x(1) - xi1) / &
         &   (q**3 * (q + x(1) - xi1)**2))) &
         &   - costh * ((x(2) - xi2) * (2.0D0 * (x(3) - xi3) / (r * (r + x(1) - xi1)) &
         &   + 4.0D0 * (x(3) - xi3) / (q * (q + x(1) - xi1)) + 4.0D0 * xi3 * x(3) * &
         &   (x(3) + xi3) * (2.0D0 * q + x(1) - xi1) / (q**3 * (q + x(1) - xi1)**2)) &
         &   + 6.0D0 * DAtanPV((x(1) - xi1) * (x(2) - xi2), ((h + x(3) + xi3) * (q + h))) &
         &   - 3.0D0 * DAtanPV((x(1) - xi1) * (r3 - xi), (r2 * r)) + &
         &   6.0D0 * DAtanPV((x(1) - xi1) * (q3 + xi), (q2 * q))) + &
         &   6.0D0 * (secth * DAtanPV(((k - q2 * costh) * (q - k) + (q3 + xi) * k * sinth) &
         &   , ((x(1) - xi1) * (q3 + xi) * costh)) &
         &   + x(3) * (((sinth**2 - costh**2) * (q3 + xi) + 2.0D0 * q2 * &
         &   costh * sinth) / (q * (q + x(1) - xi1)) + (x(1) - xi1) * sinth**2 &
         &   / (q * (q + q3 + xi))))

           TERM6=SINTH*((X(2)-XI2)*(2.0D0*(X(3)-XI3)/(R*(R+X(1)-XI1)) &
         &   +4.0D0*(X(3)-XI3)/(Q*(Q+X(1)-XI1))-4.0D0*XI3*X(3)* &
         &   (X(3)+XI3)*(2.0D0*Q+X(1)-XI1)/(Q**3*(Q+X(1)-XI1)**2)) &
         &   -6.0D0*DATANPV((X(1)-XI1)*(X(2)-XI2),((H+X(3)+XI3)*(Q+H))) &
         &   +3.0D0*DATANPV((X(1)-XI1)*(R3-XI),(R2*R))- &
         &   6.0D0*DATANPV((X(1)-XI1)*(Q3+XI),(Q2*Q))) &
         &   +COSTH*(DLOG(R+X(1)-XI1)-DLOG(Q+X(1)-XI1)- &
         &   2.0D0*(X(3)-XI3)**2/(R*(R+X(1)-XI1)) &
         &   -4.0D0*((X(3)+XI3)**2-XI3*X(3))/(Q*(Q+X(1)-XI1)) &
         &   -4.0D0*XI3*X(3)*(X(3)+XI3)**2*(2.0D0*Q+X(1)-XI1)/ &
         &   (Q**3*(Q+X(1)-XI1)**2)) &
         &   +6.0D0*X(3)*(COSTH*SINTH*(2.0D0*(Q3+XI)/(Q*(Q+X(1)-XI1)) &
         &   +(X(1)-XI1)/(Q*(Q+Q3+XI))) &
         &   -Q2*(SINTH**2-COSTH**2)/(Q*(Q+X(1)-XI1)))

           deld = dbot - dtop
           v(1) = v(1) + factor * (term1 * ucap1 + term4 * ucap) / 37.69911184D0
           v(2) = v(2) + factor * (term2 * ucap1 + term5 * ucap) / 37.69911184D0
           V(3) = V(3) + FACTOR * (TERM3 * UCAP1 + TERM6 * UCAP) / 37.69911184D0

  100  CONTINUE
      !---------- end Mansinha & Smylie [1971] component (shearing only) -------------------

      !---------- begin Okada [1985] component (tensile/mode-1/opening cracks) -------------
       IF (u3_Okada /= 0.0D0) THEN ! Note this this test is ALMOST ALWAYS passed, and therefore this Okada section almost always runs.
                                   ! However, its contribution may be infinitesimal if the dislocation is basically of shear/fault type.
                                   ! Aura would (typically) run faster if the test were for (ABS(u3_Okada) > n.nnDnn),
                                   ! or if extra LOGICAL arguments were added to control which code sections are run!
            !Note that Okada uses a right-handed Cartesian coordinate system (like M&S)
            !but that his has z increasing UP, and therefore if we make his x axis (along-strike)
            !equal to the z1 axis of M&S, then Okada's y axis is opposite to the z2 axis of M&S.
            !Another distinction is that M&S have their (z1, z2) origin at a central point on
            !the fault trace (projected to the free surface), but Okada's (x, y) origin
            !is located over the deepest part of the opening patch.

            !Poisson's ratio term:
            mu_over_lambda_plus_mu = 0.50D0 ! (assumed here in Aura; equivalent to Poisson's ratio of 1/4)

            !Dip of crack/fault:
            del_Okada = theta
            !Dimensions of dislocation patch:
            d_Okada = dbot * sinth ! measured vertically, not down-dip
            w_Okada = dbot - dtop ! measured down-dip; should be positive.
            L_Okada = l ! half-length, along strike (same as M&S)
                        ! NOTE that Figure 1 of Okada [1985] shows half of the fault-patch with dashed-lines;
                        ! we are adopting the convention that this "half" is part of the source;
                        ! Okada explains in his text how to adjust his formulas [see * below] if this is so.

            !Observer (test) location [which must be on the surface, using formulas of Okada, 1985]:
            x_Okada = +x(1)     ! reflecting parallelism of Okada's x-axis with Mansinha & Smylie's x1-axis
            y_Okada = -x(2) + & ! reflecting anti-parallelism of Okada's y-axis with Mansinha & Smylie's x2-axis;
                    &  d_Okada/tanth  ! also note OFFSET between Okada and M&S origins along the dip-direction!
            z_Okada =  0.0D0    ! (On free surface. To change this, use formulas of Okada, 1992 instead of Okada, 1985.)

            !Distances used in Okada formulas to describe position of observer in a vertical cross-section perpendicular to fault/crack strike:
            p_Okada = y_Okada * costh + d_Okada * sinth ! distance from lower edge of dislocation patch; component measured parallel to the fault/crack
            q_Okada = y_Okada * sinth - d_Okada * costh ! component of distance from fault/crack plane (extended beyond dislocation patch, if necessary)

            !Initialize remote displacement {rate} before sum of 4 terms:
            ux_Okada = 0.0D0
            uy_Okada = 0.0D0
            uz_Okada = 0.0D0

            DO 200 i = 1, 4 ! sum over 4 terms (from rectangular integration limits in 2-D):

               !Integration limit sign-changing factors, and integration limits:
                factor = 1.0D0 ! unless...
                IF ((i == 2).OR.(i == 3)) factor = -1.0D0 ! replacing the choice made in the previous line
                epsilon_Okada = x_Okada + L_Okada ! *This is where we choose the double-length, symmetrical source option.
                IF (i >= 3) epsilon_Okada = x_Okada - L_Okada ! replacing the choice made in the previous line
                nu_Okada = p_Okada ! unless...
                IF ((i == 2).OR.(i == 4)) nu_Okada = p_Okada - w_Okada ! replacing the choice made in the previous line

               !Geometric factors dependent on the above:
                yTilde_Okada = nu_Okada * costh + q_Okada * sinth
                dTilde_Okada = nu_Okada * sinth - q_Okada * costh
                R2_Okada = epsilon_Okada**2 + nu_Okada**2 + q_Okada**2
                !R2_Okada = epsilon_Okada**2 + yTilde_Okada**2 + dTilde_Okada**2 ! (alternate formula given by Okada)
                R_Okada = SQRT(R2_Okada)
                XCap2_Okada = epsilon_Okada**2 + q_Okada**2
                XCap_Okada = SQRT(XCap2_Okada)

               !Cryptic Okada factors/terms:  Note that I2_Okada is not used, because it is associated with strike-slip component of source dislocation.
               !                              Also note that factors must be computed in reverse order (5, 4, 3, 1), due to sequential self-referencing.
                I5_Okada = mu_over_lambda_plus_mu * (2.0D0/costh) * ATAN( (nu_Okada*(XCap_Okada+q_Okada*costh) + XCap_Okada*(R_Okada+XCap_Okada)*sinth) / &
                                                                       &  (epsilon_Okada*(R_Okada+XCap_Okada)*costh) )
                I4_Okada = mu_over_lambda_plus_mu * (1.0D0/costh) * ( LOG(R_Okada+dTilde_Okada) - sinth*LOG(R_Okada+nu_Okada) )
                I3_Okada = mu_over_lambda_plus_mu * ( (1.0D0/costh)*(yTilde_Okada/(R_Okada+dTilde_Okada)) - LOG(R_Okada+nu_Okada) ) + tanth*I4_Okada
                I1_Okada = mu_over_lambda_plus_mu * ( (-1.0D0/costh) * (epsilon_Okada/(R_Okada+dTilde_Okada)) ) - tanth*I5_Okada

               !Final Okada formula for remote displacement due to a finite rectangular opening (mode-1) source dislocation:
                ux_Okada = ux_Okada + factor * (u3_Okada/Two_Pi) * ( (q_Okada**2/(R_Okada*(R_Okada+nu_Okada))) - &
                                                                  &   I3_Okada*(sinth**2) )
                uy_Okada = uy_Okada + factor * (u3_Okada/Two_Pi) * ( (-dTilde_Okada*q_Okada/(R_Okada*(R_Okada+epsilon_Okada))) - &
                                                                 &    sinth*( (epsilon_Okada*q_Okada/(R_Okada*(R_Okada+nu_Okada))) - ATAN((epsilon_Okada*nu_Okada) / (q_Okada*R_Okada)) ) - &
                                                                 &    I1_Okada*(sinth**2) )
                uz_Okada = uz_Okada + factor * (u3_Okada/Two_Pi) * ( (yTilde_Okada*q_Okada/(R_Okada*(R_Okada+epsilon_Okada))) + &
                                                                 &    costh*( (epsilon_Okada*q_Okada/(R_Okada*(R_Okada+nu_Okada))) - ATAN((epsilon_Okada*nu_Okada) / (q_Okada*R_Okada)) ) - &
                                                                 &    I5_Okada*(sinth**2) )
200         CONTINUE ! sum over 4 terms

            !Final contribution to remote displacment {rate}, reversing signs
            !of two terms, per reversal of two axis-directions between M&S versus O:
            v(1) = v(1) + ux_Okada
            v(2) = v(2) - uy_Okada
            v(3) = v(3) - uz_Okada

       END IF ! non-zero opening component
      !---------- end Okada [1985] component (tensile/mode-1/opening cracks) ---------------

      !Copy results to output variables:
       v4(1) = v(1)
       v4(2) = v(2)
       v4(3) = v(3) ! NOTE that this vertical component of remote displacment is positive-downwards,
                    ! according to the Mansinha & Smylie [1971] convention.

       END SUBROUTINE Aura
!------------------------------------------------------------------------------

       SUBROUTINE DOnArc (s, theta1, phi1, theta2, phi2, & ! input
     &                    theta, phi) ! output

!   Computes coordinates (THETA, PHI) =
!  (colatitude, longitude) in radians
!   for a point which lies on the great circle arc from
!  (THETA1, PHI1) to (THETA2, PHI2).
!   Input parameter S expresses the fractional distance,
!   so input S = 0.0D0 returns (THETA1,PHI1) and
!      input S = 1.0D0 returns (THETA2,PHI2).

       IMPLICIT NONE

!   Arguments (see text above):
       REAL*8, INTENT(IN)  :: s, theta1, phi1, theta2, phi2
       REAL*8, INTENT(OUT) :: theta, phi

       REAL*8 :: comple, equat, size
       REAL*8, DIMENSION(3) :: tvec, uvec, uvec1, uvec2
!     These are all unit vectors (in the unit sphere)
!     in a Cartesian coordinate system with
!     x = (0E, 0N), y = (90E, 0N), z = North pole.

       uvec1(1) = DSIN(theta1) * DCOS(phi1)
       uvec1(2) = DSIN(theta1) * DSIN(phi1)
       uvec1(3) = DCOS(theta1)

       uvec2(1) = DSIN(theta2) * DCOS(phi2)
       uvec2(2) = DSIN(theta2) * DSIN(phi2)
       uvec2(3) = DCOS(theta2)

       comple = 1.0D0 - s
       tvec(1) = comple * uvec1(1) + s * uvec2(1)
       tvec(2) = comple * uvec1(2) + s * uvec2(2)
       tvec(3) = comple * uvec1(3) + s * uvec2(3)

       size = DSQRT(tvec(1)**2 + tvec(2)**2 + tvec(3)**2)
       uvec(1) = tvec(1) / size
       uvec(2) = tvec(2) / size
       uvec(3) = tvec(3) / size

       equat = DSQRT(uvec(1)**2 + uvec(2)**2)
       theta = DATAN2(equat, uvec(3))
       phi   = DATan2F(uvec(2), uvec(1))
       END SUBROUTINE DOnArc
!------------------------------------------------------------

    REAL*8 FUNCTION DATan2F(y, x)
        ! Works like ATAN2 but corrects for case of (0.0D0, 0.0D0).
        ! Returns inverse tangent in radians.
        IMPLICIT NONE
        REAL*8, INTENT(IN)    :: y, x
        IF (y == 0.0D0) THEN
            IF (x == 0.0D0) THEN
                DATAN2F = 0.0D0
            ELSE
                DATAN2F = DATAN2(y, x)
            END IF
        ELSE
            DATAN2F = DATAN2(y, x)
        END IF
    END FUNCTION DATan2F
!---------------------------------------------------------------

    REAL*8 FUNCTION DFltLen (phi1, phi2, radius, theta1, theta2)

!      Calculates length of great circle segment between
!      point (THETA1,PHI1) and point (THETA2,PHI2),
!      in physical length units (radians*planet RADIUS).

       IMPLICIT NONE
       REAL*8, INTENT(IN) :: phi1, phi2, radius, theta1, theta2
       DOUBLE PRECISION :: ab
       ab    = DSIN(theta1) * DSIN(theta2) * DCOS(phi1) * DCOS(phi2) + &
     &         DSIN(theta1) * DSIN(theta2) * DSIN(phi1) * DSIN(phi2) + &
     &         DCOS(theta1) * DCOS(theta2)
       ab = DACOS(ab)
       DFltLen = ab * radius
   END FUNCTION DFltLen
!-------------------------------------------------------------------

   REAL*8 FUNCTION DChord (angle1, s, angle2)

!   Returns an angle obtained by interpolation between ANGLE1
!   and ANGLE2.  The interpolation method is not sensitive to
!   possible cycle shifts (of 2*n*PI) between ANGLE1 and ANGLE2.

!   Unit vectors are constructed for ANGLE1 and ANGLE2, and a
!   linear chord is drawn between their tips.

!   Double-precision S is the internal coordinate along the chord;
!   it is dimensionless, with value 0.0 at ANGLE1 and 1.0 at
!   ANGLE2.  (The user may input S values outside this range
!   to get results outside the (smaller) angle between ANGLE1 and
!   ANGLE2, if desired.)  The angle returned is that from the
!   origin to this chord point.

!   This algorithm should work equally well for angles measured
!   either clockwise or counterclockwise from any reference, as
!   long as the usage is consistent.

!   Both the input angles and the result "Chord" are in radians.

       DOUBLE PRECISION, INTENT(IN) :: s
       REAL*8, INTENT(IN) :: angle1, angle2

       REAL*8, DIMENSION(2) :: uvec1, uvec2, vecs

       uvec1(1) = DCOS(angle1)
       uvec1(2) = DSIN(angle1)
       uvec2(1) = DCOS(angle2)
       uvec2(2) = DSIN(angle2)
       vecs(1) = (1.0D0 - s) * uvec1(1) + s * uvec2(1)
       vecs(2) = (1.0D0 - s) * uvec1(2) + s * uvec2(2)
       DChord = DATAN2F(vecs(2), vecs(1))
    END FUNCTION DChord

    DOUBLE PRECISION FUNCTION DAtanPV (y, x)

!   Returns principal value of inverse tangent of (y,x)

       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(IN) :: x, y
       DOUBLE PRECISION :: xx, yy
       IF (x >= 0.0D0) THEN
            xx = x
            yy = y
       ELSE
            xx = -x
            yy = -y
       END IF
       IF((yy == 0.0D0).OR.(xx == 0.0D0)) THEN
          IF(yy == 0.0D0) DAtanPV = 0.0D0
          IF(xx == 0.0D0) DAtanPV = 1.57079632679490D0
       ELSE
          DAtanPV = DATAN2(yy, xx)
       END IF
    END FUNCTION DAtanpv

end module DDislocation


module DIcosahedron
    !
    ! Initial letter "D" indicates a DOUBLE PRECISION (REAL*8) version,
    ! created 2015.02 as part of a systematic upgrade of several of my codes
    ! (Shells, NeoKinema, FiniteMap, & NeoKineMap) to 64-bit precision.
    ! It is intended that the "D" FUNCTIONs and SUBROUTINEs here should be
    ! logically equivalent to those in MODULE Icosahedron,
    ! just more precise.  Naturally, they now all take REAL*8 arguments
    ! in place of the old REAL arguments.
    !
    ! Updated versions copyright 2015 by Peter Bird and the
    ! Regents of the University of California.
    !
    !----------------------------------------------------------------------
    !
    !Tools for creation of global finite element grids by subdivision of the
    !facets of the icosahedron into quasi-equilateral spherical triangles.
    !
    !by Peter Bird, UCLA; written in BASIC ~1982; translated to Fortran 90 2001.11
    !----------------------------------------------------------------------------
      PRIVATE
    !
    ! USER ROUTINES:
    !
                    PUBLIC DMake_Global_Grid
                    PUBLIC DWrite_Global_Grid
    !
    ! UTILITY ROUTINES (called by the user routines):
    !
    !     RECURSIVE SUBROUTINE DDivide
    !               SUBROUTINE DElement_Out
    !               SUBROUTINE DUnitize
    !----------------------------------------------------------------------------

CONTAINS

    SUBROUTINE DMake_Global_Grid (n_slice, &           ! only input(!)
                                & numnod, node_uvec, & ! output: number of nodes, unit vectors of nodes,
                                & numel, nodes)        !         number of elements, element definitions
	   use SharedVars
       !---------------------------------------------------------------
        !Generate a finite element grid of spherical triangles by subdivision of the icosahedron:
        !Level 0 (generated here) has 12 vertices, 30 edges, and 20 triangles.
        !Then, subdivide each face n_slice times, and output triangular elements,
        !  by using SUBR Divide, which calls itself recursively!
       !---------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,                 INTENT(IN)  :: n_slice   ! subdivision level; 0 or higher
        INTEGER,                 INTENT(OUT) :: numnod    ! number of nodes created
        REAL*8, DIMENSION(:,:),  INTENT(OUT) :: node_uvec ! (3, numnod)
        INTEGER,                 INTENT(OUT) :: numel     ! number of elements created
        INTEGER, DIMENSION(:,:), INTENT(OUT) :: nodes     ! (3, numel)
       !---------------------------------------------------------------
        INTEGER                            :: m_slice ! copy of n_slice, allowing it to be changed (counted-down)
        ! INTEGER                            :: facets_done, i, j, k
        INTEGER                            :: facets_done, k
        INTEGER, DIMENSION(3)              :: node_number
        REAL*8, DIMENSION(3)               :: rx, ry, rz, uvec1, uvec2, uvec3
        DOUBLE PRECISION, PARAMETER        :: s = 1.107148719D0
        DOUBLE PRECISION                   :: dot1, dot2, dot3, x1, x2, x3, y1, y2, y3, z1, z2, z3
        DOUBLE PRECISION, DIMENSION(12)    :: lat, lon
        DOUBLE PRECISION, DIMENSION(3)     :: v1, v2, v3
        DOUBLE PRECISION, DIMENSION(3, 12) :: abg  !Cartesian (alpha, beta, gamma) coordinates of these vertices.
        LOGICAL counterclockwise
		character(len=30) :: IFilename
       !---------------------------------------------------------------

       !generate basic form with a vertex (a 5-fold axis) up; highest symmetry axis
        lat(1) = 1.570796327D0
        lon(1) = 0.0D0
        DO i = 2, 6
            lat(i) = lat(1) - s
            lon(i) = (i - 2.0D0) * 1.256637061D0
        END DO
        DO i = 7, 11
            lat(i) = -lat(1) + s
            lon(i) = (i - 7.0D0) * 1.256637061D0 + .628318531D0
        END DO
        lat(12) = -lat(1)
        lon(12) = 0.0D0
        DO i = 1, 12
            abg(1, i) = DCOS(lat(i)) * DCOS(lon(i))
            abg(2, i) = DCOS(lat(i)) * DSIN(lon(i))
            abg(3, i) = DSIN(lat(i))
        END DO
       !-------------------------------------------------------
       !create output file for dumping results as they are found:
	   write(IFilename,"(A,I0,A)") 'Icosahedron_',ThID,'.tmp'
       open(unit=729,file=trim(IFilename),form="UNFORMATTED")
	   !-------------------------------------------------------
       !find all 20 faces and subdivide each into four spherical triangles;
       !WRITE (*, "(' Creating global grid by level-',I1,' subdivision of icosahedron facets:')") n_slice
       !WRITE (*, *) ! advance, because next WRITE will not
        m_slice = n_slice
        facets_done = 0
        DO i = 1, 10
            DO j = (i + 1), 11
                DO k = (j + 1), 12
                    dot1 = abg(1, i) * abg(1, j) + abg(2, i) * abg(2, j) + abg(3, i) * abg(3, j)
                    dot2 = abg(1, j) * abg(1, k) + abg(2, j) * abg(2, k) + abg(3, j) * abg(3, k)
                    dot3 = abg(1, k) * abg(1, i) + abg(2, k) * abg(2, i) + abg(3, k) * abg(3, i)
                    IF ((dot1 > 0.3D0) .AND. (dot2 > 0.3D0) .AND. (dot3 > 0.3D0)) THEN
                        x1 = abg(1, i)
                        x2 = abg(1, j)
                        x3 = abg(1, k)
                        y1 = abg(2, i)
                        y2 = abg(2, j)
                        y3 = abg(2, k)
                        z1 = abg(3, i)
                        z2 = abg(3, j)
                        z3 = abg(3, k)
                        !Note: Divide will call itself, recursively.
                        !Therefore, all inputs are simple numbers (not vectors),
                        !to go on stack as values, not addresses.
                        !Also, note that output is sent to a temporary file, to avoid multiple copies on stack!
                        CALL DDivide(x1, y1, z1, x2, y2, z2, x3, y3, z3, m_slice) ! using copy m_slice = n_slice
                        facets_done = facets_done + 1
                       !WRITE (*, "('+',I8,' facets   out of       20 divided into elements')") facets_done
                    END IF
                END DO
            END DO
        END DO
       !-----------------------------------------------------------------------
       !read binary file, extracting groups of 3 uvecs, and assigning to lists:
        CLOSE (UNIT = 729)
        ! OPEN  (UNIT = 729, FILE = "Icosahedron.tmp", STATUS = "OLD", FORM = "UNFORMATTED")
        OPEN  (UNIT = 729, FILE = trim(IFilename), STATUS = "OLD", FORM = "UNFORMATTED")
        numnod = 0  !initialization
        IF (n_slice == 0) THEN
            numel = 20
        ELSE
            numel = 20 * (4**n_slice)
        END IF
        WRITE (*, *) ! advance, because next WRITE will not
        DO i = 1, numel
            READ (729) rx(1), ry(1), rz(1), rx(2), ry(2), rz(2), rx(3), ry(3), rz(3)
            DO j = 1, 3
                node_number(j) = 0 ! initialization
                k_loop: DO k = 1, numnod
                    IF (rx(j) == node_uvec(1, k)) THEN
                        IF (ry(j) == node_uvec(2, k)) THEN
                            IF (rz(j) == node_uvec(3, k)) THEN
                               !this node is already defined
                                node_number(j) = k
                                EXIT k_loop
                            END IF
                        END IF
                    END IF
                END DO k_loop ! k = 1, numnod
                IF (node_number(j) == 0) THEN
                   !no match was found; define a new node
                    numnod = numnod + 1
                    node_number(j) = numnod
                    node_uvec(1, numnod) = rx(j)
                    node_uvec(2, numnod) = ry(j)
                    node_uvec(3, numnod) = rz(j)
                END IF
               !record this element, after checking for counterclockwise ordering!
               !First, determine side from node "#1" to "#2":
                uvec3(1) = rx(2) - rx(1)
                uvec3(2) = ry(2) - ry(1)
                uvec3(3) = rz(2) - rz(1)
               !Next, determine side from node "#2" to "#3":
                uvec1(1) = rx(3) - rx(2)
                uvec1(2) = ry(3) - ry(2)
                uvec1(3) = rz(3) - rz(2)
                CALL DCross(uvec3, uvec1, uvec2)
               !Note: If numbering is counterclockwise, uvec2 should point up.
                counterclockwise = ((uvec2(1)*rx(1) + uvec2(2)*ry(1) + uvec2(3)*rz(1)) > 0.0D0)
                IF (counterclockwise) THEN
                    nodes(1:3, i) = node_number(1:3)
                ELSE
                    nodes(1, i) = node_number(1)
                    nodes(2, i) = node_number(3)
                    nodes(3, i) = node_number(2)
                END IF
            END DO ! j = 1, 3
            IF (MOD(i, 100) == 0) THEN
               !WRITE (*, "('+',I8,' elements out of ',I8,' scanned for new nodes')") i, numel
            END IF
        END DO ! i = 1, numel
       !WRITE (*, "('+',I8,' elements out of ',I8,' scanned for new nodes')") numel, numel
        CLOSE (UNIT = 729, STATUS = "DELETE")
    END SUBROUTINE DMake_Global_Grid

    SUBROUTINE DWrite_Global_Grid (path, &              ! [Drive:][\Path\] for output (or null)
                                 & n_slice, &
                                 & numnod, node_uvec, &
                                 & numel, nodes)        ! all of these are INTENT(IN); get values
                                                        ! by calling Make_Global_Grid first
        IMPLICIT NONE
        CHARACTER*(*),           INTENT(IN) :: path      ! [Drive:][\Path\] for output (or null)
        INTEGER,                 INTENT(IN) :: n_slice   ! subdivision level; 0 or higher
        INTEGER,                 INTENT(IN) :: numnod    ! number of nodes
        REAL*8, DIMENSION(:,:),  INTENT(IN) :: node_uvec ! (3, numnod)
        INTEGER,                 INTENT(IN) :: numel     ! number of elements
        INTEGER, DIMENSION(:,:), INTENT(IN) :: nodes     ! (3, numel)
       !---------------------------------------------------------------
        CHARACTER*1  :: c1
        CHARACTER*11 :: grid_file
        CHARACTER*80 :: grid_pathfile
        INTEGER :: i
        REAL*8, PARAMETER :: Pi = 3.141592653589793D0
        REAL*8, PARAMETER :: degrees_per_radian = 57.2957795130D0
        REAL*8 :: equat, equat2, lat, lon, phi, theta
        REAL*8, DIMENSION(3) :: uvec
       !---------------------------------------------------------------
        WRITE (c1, "(I1)") n_slice
        grid_file = "global" // c1 // ".feg"
        grid_pathfile = TRIM(path) // grid_file
        WRITE (*, "(' Writing ',A,'...')") TRIM(grid_file)
        OPEN (UNIT = 730, FILE = grid_pathfile) ! unconditional; overwrites any older file
        WRITE (730, "(A)") grid_file
        WRITE (730, "(I6,I6,'     0 1000000 F')") numnod, numnod
        DO i = 1, numnod
            uvec(1:3) = node_uvec(1:3, i)
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           !Note: NOT calling DUvec_2_LonLat, because that is in DSphere or DMap_Projections/DAdobe_Illustrator,
           !      and in many projects would cause confusing duplication of code.  However, steps exactly the same.
            equat2 = uvec(1)*uvec(1) + uvec(2)*uvec(2)
            IF (equat2 == 0.0D0) THEN
                phi = 0.0D0 ! actually undefined; default 0.
                IF (uvec(3) > 0.0D0) THEN
                    theta = 0.0D0 ! N pole
                ELSE
                    theta = Pi    ! S pole
                END IF
            ELSE
                equat = DSQRT(equat2)
                theta = DATAN2(equat, uvec(3))
                phi = DATAN2(uvec(2), uvec(1))
            END IF
            lat = 90.D0 - degrees_per_radian * DABS(theta)
            lat = MAX (lat, -90.0D0)
            lon = degrees_per_radian * phi
            IF (lon > 180.0D0) lon = lon - 360.0D0
            IF (lon <= -180.0D0) lon = lon + 360.0D0
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                WRITE (730, "(I6,F9.3,F8.3,'       0.0       0.0       0.0       0.0')") i, lon, lat
        END DO
        WRITE (730, *) numel
        DO i = 1, numel
            WRITE (730, "(4I6)") i, nodes(1, i), nodes(2, i), nodes(3, i)
        END DO
        WRITE (730, "('     0')")
        CLOSE (UNIT = 730)
        WRITE (*, "('+Writing ',A,'...DONE')") grid_file
    END SUBROUTINE DWrite_Global_Grid

    SUBROUTINE DCross (a_vec, b_vec, c_vec)
      ! vector cross product: a x b = c
        IMPLICIT NONE
        REAL*8, DIMENSION(3), INTENT(IN)  :: a_vec, b_vec
        REAL*8, DIMENSION(3), INTENT(OUT) :: c_vec
        c_vec(1) = a_vec(2)*b_vec(3) - a_vec(3)*b_vec(2)
        c_vec(2) = a_vec(3)*b_vec(1) - a_vec(1)*b_vec(3)
        c_vec(3) = a_vec(1)*b_vec(2) - a_vec(2)*b_vec(1)
    END SUBROUTINE DCross

    RECURSIVE SUBROUTINE DDivide (x1, y1, z1, x2, y2, z2, x3, y3, z3, m_slice)
      !CALLed only by SUBR Make_Global_Grid (but called multiple times) and by itself (ditto).
      !Accepts 3 unit vectors of triangle corners as input, subdivides the
      !triangle m_slice times into 4 triangles, and outputs results to a file.
      !NOTE that Divide calls itself recursively, so all the arguments should
      !be simple numbers passed by value, not by address.
      !To avoid complications, output is dumped into an open binary file (UNIT = 729).
      !-------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(IN)    :: x1, y1, z1, x2, y2, z2, x3, y3, z3
       INTEGER,          INTENT(IN)    :: m_slice
      !-------------------------------------------------------
       INTEGER :: l_slice
       DOUBLE PRECISION :: xe1, xe2, xe3, ye1, ye2, ye3, ze1, ze2, ze3
      !-------------------------------------------------------
       xe1 = 0.5D0 * (x2 + x3)
       xe2 = 0.5D0 * (x3 + x1)
       xe3 = 0.5D0 * (x1 + x2)
       ye1 = 0.5D0 * (y2 + y3)
       ye2 = 0.5D0 * (y3 + y1)
       ye3 = 0.5D0 * (y1 + y2)
       ze1 = 0.5D0 * (z2 + z3)
       ze2 = 0.5D0 * (z3 + z1)
       ze3 = 0.5D0 * (z1 + z2)
       CALL DUnitize(xe1, ye1, ze1)
       CALL DUnitize(xe2, ye2, ze2)
       CALL DUnitize(xe3, ye3, ze3)
       IF (m_slice > 1) THEN
           l_slice = m_slice - 1
           CALL DDivide( x1,  y1,  z1, xe2, ye2, ze2, xe3, ye3, ze3, l_slice)
           CALL DDivide( x2,  y2,  z2, xe1, ye1, ze1, xe3, ye3, ze3, l_slice)
           CALL DDivide( x3,  y3,  z3, xe1, ye1, ze1, xe2, ye2, ze2, l_slice)
           CALL DDivide(xe1, ye1, ze1, xe2, ye2, ze2, xe3, ye3, ze3, l_slice)
       ELSE IF (m_slice == 1) THEN
           CALL DElementOut( x1,  y1,  z1, xe2, ye2, ze2, xe3, ye3, ze3)
           CALL DElementOut( x2,  y2,  z2, xe1, ye1, ze1, xe3, ye3, ze3)
           CALL DElementOut( x3,  y3,  z3, xe1, ye1, ze1, xe2, ye2, ze2)
           CALL DElementOut(xe1, ye1, ze1, xe2, ye2, ze2, xe3, ye3, ze3)
       ELSE ! m_slice == 0
           CALL DElementOut(x1, y1, z1, x2, y2, z2, x3, y3, z3)
       END IF
    END SUBROUTINE DDivide



    SUBROUTINE DElementOut (x1, y1, z1, x2, y2, z2, x3, y3, z3)
       !CALLed only by SUBR DDivide (but called multiple times).
       !Outputs an element defined by 3 unit radius vectors to its vertices.
       !The output (binary) file must already be open, as device 729.
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
       !Write to binary file:
        WRITE (729) x1, y1, z1, x2, y2, z2, x3, y3, z3
    END SUBROUTINE DElementOut



    SUBROUTINE DUnitize (x, y, z)
       !Converts any vector in 3-component double-precision form to a unit vector.
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(INOUT) :: x, y, z
        DOUBLE PRECISION :: r
        r = DSQRT(x * x + y * y + z * z)
        IF (r > 0.0D0) THEN
            x = x / r
            y = y / r
            z = z / r
        ELSE
            x = 1.0D0
            y = 0.0D0
            z = 0.0D0
        END IF
    END SUBROUTINE DUnitize

end module DIcosahedron


module DWeighting

    ! Initial letter "D" indicates a DOUBLE PRECISION (REAL*8) version,
    ! created 2015.07 as part of a systematic upgrade of several of my codes
    ! (Shells, NeoKinema, FiniteMap, & NeoKineMap) to 64-bit precision.
    ! It is intended that the "D" FUNCTIONs and SUBROUTINEs here should be
    ! logically equivalent to those in MODULE Weighting,
    ! just more precise.  Naturally, they now all take REAL*8 arguments
    ! in place of the old REAL arguments.

    ! Computes "weights" (with mean value 1.0)
    ! for a list of "number_of_data" datum positions,
    ! each of which is specified by "theta_radians"
    ! (colatitude, measured from North pole) and
    ! "phi_radians" (longitude, measured East from Greenwich).

    ! Exploits pre-existing code (USE DIcosahedron)
    ! to create a global finite-element grid of equal-sized
    ! spherical triangles.

    ! Then, each datum position is assigned to a finite element,
    ! and counts of total data in each finite element are kept.
    ! The surface area (actually, solid angle, in steradians)
    ! of the element is equally divided among these data.

    ! Finally, surface areas of empty elements are added to
    ! the area associated with the nearest datum.

    ! The result is that weights are proportional to surface areas
    !"associated with" each datum, where "association" is procedurally
    ! defined by this algorithm.

    ! by Peter Bird, UCLA, November 2006.
    ! (c) Copyright 2006 by G. Peter Bird and the
    !     Regents of the University of California;
    !     REAL*8 revision copyright 2015.

    USE DSphere ! Peter Bird's Fortran 90 MODULE file DSphere.f90

    PRIVATE ! <=== Very Important!  We only desire access through
            ! CALL DInitialize_Weighting(subdivision) ! and
            ! CALL DPerform_Weighting(...)            ! ;
            ! we don't want to make all internal variables and arrays visible!
    PUBLIC DInitialize_Weighting
    PUBLIC DIs_Initialization_Needed
    PUBLIC DPerform_Weighting


    ! The following values need to have a "lifetime" that extends
    ! beyond the end of Initialize_Weighting,
    ! so they are available for later calls to Perform_Weighting.
    SAVE
    INTEGER :: n_slice, numNod, numEl
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: nodes
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: neighbor   ! neighbors of each element
    LOGICAL :: initialized = .FALSE.
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: center     ! center uvecs of elements
    REAL*8, DIMENSION(:),   ALLOCATABLE :: a_         ! (plane) areas of elements
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: node_uvec

CONTAINS

SUBROUTINE DInitialize_Weighting (subdivision)
    !Create global finite-element grid and store it.
    !Suggested values of subdivision are 2 ... 5.

    USE DIcosahedron ! Peter Bird's Fortran 90 MODULE file DIcosahedron.f90

    IMPLICIT NONE
    SAVE
    LOGICAL :: chatty
    INTEGER, INTENT(IN) :: subdivision
   !Note that its value is captured and kept here, as n_slice:
    n_slice = subdivision
    numNod = 2 + 10 * (4**subdivision)
    IF (ALLOCATED(node_uvec)) DEALLOCATE(node_uvec)
    ALLOCATE ( node_uvec(3, numNod) )
    numEl = 20 * (4**subdivision)
    IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
    ALLOCATE ( nodes(3, numEl) )

    CALL DMake_Global_Grid (n_slice = subdivision, &                  ! only input(!)
                          & numNod = numNod, node_uvec = node_uvec, & ! output: number of nodes, unit vectors of nodes,
                          & numEl = numEl, nodes = nodes)             !         number of elements, element definitions

    IF (ALLOCATED(a_)) DEALLOCATE(a_)
    ALLOCATE ( a_(numEl) )
    IF (ALLOCATED(center)) DEALLOCATE(center)
    ALLOCATE ( center(3, numEl) )
    IF (ALLOCATED(neighbor)) DEALLOCATE(neighbor)
    ALLOCATE ( neighbor(3, numEl) )

    chatty = .FALSE.
    CALL DLearn_Spherical_Triangles (numEl, nodes, node_uvec, chatty, &
                                   & a_, center, neighbor)

    initialized = .TRUE.

END SUBROUTINE DInitialize_Weighting

LOGICAL FUNCTION DIs_Initialization_Needed()
    IMPLICIT NONE
    DIs_Initialization_Needed = (.NOT.initialized) ! from private to public
END FUNCTION DIs_Initialization_Needed

SUBROUTINE DPerform_Weighting (number_of_data, theta_radians, phi_radians, & ! INTENT(IN)
                             & weights) ! INTENT(OUT)
    USE DSphere ! Peter Bird's Fortran 90 MODULE file DSphere.f90
	use SharedVars
	use ShellSetSubs
    IMPLICIT NONE

    INTEGER :: best
    INTEGER, INTENT(IN) :: number_of_data
    REAL*8 :: distance, trial_distance
    REAL*8, DIMENSION(:), INTENT(IN)  :: theta_radians, phi_radians
    REAL*8, DIMENSION(:), INTENT(OUT) :: weights

    INTEGER :: iEle, n_allocated
    INTEGER, DIMENSION(:),   ALLOCATABLE :: data_in_element    ! one entry for each element
    INTEGER, DIMENSION(:),   ALLOCATABLE :: element_assignment ! one entry for each datum
    LOGICAL :: cold_start, success
    REAL*8 :: element_steradians, &
            & s1, s2, s3, sphere_steradians
    REAL*8, DIMENSION(3) :: uvec
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: uvecs

    IF (.NOT.initialized) THEN
	  write(ErrorMsg,'(A)') "You must CALL Initialize_Weighting(subdivision) before you are allowed to CALL Perform_Weighting()!"
	  call FatalError(ErrorMsg,ThID)
    END IF

    sphere_steradians = 4.0D0 * 3.141592653589793D0
    element_steradians = sphere_steradians / numEl

    ALLOCATE ( uvecs(3, number_of_data) )

    DO j = 1, number_of_data
        CALL DThetaPhi_2_Uvec (theta_radians(j), phi_radians(j), uvec)
        uvecs(1:3, j) = uvec(1:3)
    END DO

    ALLOCATE ( data_in_element(numEl) )
    data_in_element = 0 ! whole array (just initializing)

    ALLOCATE ( element_assignment(number_of_data) )
    element_assignment = 0 ! whole array, just for tidiness

    ! Assign data to elements:
    DO j = 1, number_of_data
        uvec(1:3) = uvecs(1:3, j)
        cold_start = .TRUE.
        CALL DWhich_Spherical_Triangle (uvec, cold_start, &
                                      & numEl, nodes, node_uvec, center, a_, neighbor, &
                                      & success, iEle, s1, s2, s3)
        IF (success) THEN
            element_assignment(j) = iEle
            data_in_element(iEle) = data_in_element(iEle) + 1
        ELSE
	      write(ErrorMsg,'(A,I6,A)') "ERROR: Datum #",j," did not fit into any element of global icosahedral grid."
	      call FatalError(ErrorMsg,ThID)
        END IF
    END DO

    weights = 0.0D0 ! whole array; it will be used to accumulate steradians (until the very last line)

   !Assign element_steradians to one or more data "weights":
    DO i = 1, numEl

        IF (data_in_element(i) > 0) THEN

            n_allocated = 0
                allocating: DO j = 1, number_of_data
                IF (element_assignment(j) == i) THEN
                    weights(j) = weights(j) + element_steradians / data_in_element(i)
                    n_allocated = n_allocated + 1
                    IF (n_allocated == data_in_element(i)) EXIT allocating
                END IF
            END DO allocating

        ELSE ! data_in_element(i) == 0

           !Find nearest datum, and give it the whole solid angle.

           !First, find its center:
            uvec(1) = ( node_uvec(1, nodes(1, i)) + node_uvec(1, nodes(2, i)) + node_uvec(1, nodes(3, i)) ) / 3.0
            uvec(2) = ( node_uvec(2, nodes(1, i)) + node_uvec(2, nodes(2, i)) + node_uvec(2, nodes(3, i)) ) / 3.0
            uvec(3) = ( node_uvec(3, nodes(1, i)) + node_uvec(3, nodes(2, i)) + node_uvec(3, nodes(3, i)) ) / 3.0
            CALL DMake_Uvec(uvec, uvec)

           !Second, search for datum nearest to element center:
            distance = 2.01D0 ! just initializing
            best = 0 ! ditto
            DO j = 1, number_of_data
                trial_distance = DSQRT( (uvec(1) - uvecs(1, j))**2 + &
                                      & (uvec(2) - uvecs(2, j))**2 + &
                                      & (uvec(3) - uvecs(3, j))**2 )
                IF (trial_distance < distance) THEN
                    distance = trial_distance
                    best = j
                END IF
            END DO ! j = 1, numEl
           !Note: "best" is the only important information retained from this search.
            weights(best) = weights(best) + element_steradians

        END IF

    END DO ! i = 1, numEl

    weights(1:number_of_data) = weights(1:number_of_data) * (number_of_data / sphere_steradians)
    !Now, weights are renormalized so that their sum is number_of_data, and their mean value is 1.0D0

    DEALLOCATE ( element_assignment )
    DEALLOCATE ( data_in_element )
    DEALLOCATE ( uvecs )

END SUBROUTINE DPerform_Weighting

    SUBROUTINE DDumb_s123 (element, vector, node, xyz_nod, center, a_, &
                         & s1, s2, s3)
        ! Finds s1, s2, s3 coordinates of position vector "in" element
        ! (whether the point is actually in the element or NOT).
        IMPLICIT NONE
        INTEGER,                   INTENT(IN)  :: element  ! element #
        REAL*8,    DIMENSION(3),   INTENT(IN)  :: vector   ! uvec to point
        INTEGER,   DIMENSION(:,:), INTENT(IN)  :: node     ! element definitions
        REAL*8,    DIMENSION(:,:), INTENT(IN)  :: xyz_nod  ! uvecs of nodes
        REAL*8,    DIMENSION(:,:), INTENT(IN)  :: center   ! uvec of each element (uvec)
        REAL*8,    DIMENSION(:),   INTENT(IN)  :: a_       ! element areas (plane; R == 1.0)
        REAL*8,                    INTENT(OUT) :: s1, s2, s3  ! results
        !- - - - - - - - - - - - - - - - - -  - - - - -  -
        INTEGER :: i1, i2, i3
        REAL*8, DIMENSION(3) :: tv, tvi, tvo, tv1, tv2, v1
        REAL*8 :: d1, dc, t
        IF (element == 0) THEN
            WRITE (*,"(' ERROR: element = 0 in Dumb_s123')")
            CALL DTraceback
        END IF
        i1 = node(1, element)
        i2 = node(2, element)
        i3 = node(3, element)
        !shorten(?) vector to just touch plane element -> v1
        tv1 = center(1:3, element)
        dc = DOT_PRODUCT(vector, tv1)
        IF (dc <= 0.0D0) THEN
            WRITE (*,"(' ERROR: Internal vector >= 90 deg. from element in Dumb_s123')")
            CALL DTraceback
        END IF
        tv2 = xyz_nod(1:3, i1)
        d1 = DOT_PRODUCT(tv2, tv1)
        t = d1 / dc
        v1 = t * vector
        tvi = xyz_nod(1:3,i3) - xyz_nod(1:3,i2)
        tvo = v1(1:3) - xyz_nod(1:3,i3)
        CALL DCross(tvi, tvo, tv)
        s1 = 0.5D0 * DOT_PRODUCT(tv1, tv) / a_(element)
        tvi = xyz_nod(1:3,i1) - xyz_nod(1:3,i3)
        tvo = v1(1:3) - xyz_nod(1:3,i1)
        CALL DCross(tvi, tvo, tv)
        s2 = 0.5D0 * DOT_PRODUCT(tv1, tv) / a_(element)
        s3 = 1.0D0 - s1 - s2
    END SUBROUTINE DDumb_s123

    SUBROUTINE DLearn_Spherical_Triangles (numEl, nodes, node_uvec, chatty, &
                                        & a_, center, neighbor)
       !Creates arrays needed by lookup subr. Which_Spherical_Triangle:
       !  a_       = area of plane triangle below element, when radius == 1.0
       !  center   = uvec pointing to center of element
       !  neighbor = neighboring spherical triangular elements on each
       !             side (or zero if none); note that algorithm
       !             depends on node-location match, not on node-number
       !             match, and therefore ignores intevening faults.
       !These arrays are only meaningful for finite element grids used
       !   with SHELLS and/or RESTORE.
        IMPLICIT NONE
        INTEGER,                   INTENT(IN)  :: numEl ! number of spherical triangle elements
        INTEGER,   DIMENSION(:,:), INTENT(IN)  :: nodes  ! element definitions
        REAL*8,    DIMENSION(:,:), INTENT(IN)  :: node_uvec ! uvecs of nodes
        LOGICAL,                   INTENT(IN)  :: chatty
        REAL*8,    DIMENSION(:),   INTENT(OUT) :: a_
        REAL*8,    DIMENSION(:,:), INTENT(OUT) :: center
        INTEGER,   DIMENSION(:,:), INTENT(OUT) :: neighbor
        INTEGER :: furthest, ia, ib, i1, i2, i3, j, j1, j2, j3, k, l_, m, step_aside
        REAL*8, DIMENSION(3) :: a, b, c, t, u
        IF (chatty) WRITE (*,"(' Learning the spherical triangles...')")
        furthest = (numEl + 1) / 2 ! INTEGER division is intentional
        neighbor = 0 ! whole array, initialized to "no neighbor on this side"
        homes: DO l_ = 1, numEl
           !first, a_
            i1 = nodes(1,l_)
            i2 = nodes(2,l_)
            i3 = nodes(3,l_)
            a = node_uvec(1:3,i2) - node_uvec(1:3,i1)
            b = node_uvec(1:3,i3) - node_uvec(1:3,i2)
            CALL DCross (a, b, c)
            a_(l_) = 0.5D0 * DSQRT(c(1)**2 + c(2)**2 + c(3)**2)
           !second, compute center
            t(1:3) = (node_uvec(1:3,i1)+node_uvec(1:3,i2)+node_uvec(1:3,i3)) / 3.0D0
            CALL DMake_Uvec(t, u)
            center(1:3, l_) = u(1:3)
           !third, find neighbor(?) for each side of element
            sides: DO j = 1, 3  ! 3 sides
                k = 1 + MOD (j, 3)
                ia = nodes(k, l_) ! 1st node along side
                ib = nodes(1 + MOD (k, 3), l_) ! 2nd node along side
                strangers: DO step_aside = 1, furthest
                    m = l_ + step_aside ! I also try -step_aside, below
                    m = 1 + MOD(m-1, numEl) ! wraps around
                    j1 = nodes(1, m)
                    j2 = nodes(2, m)
                    j3 = nodes(3, m)
                    IF (DSame(j1, ib) .AND. DSame(j2, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    ELSE IF (DSame(j2, ib) .AND. DSame(j3, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    ELSE IF (DSame(j3, ib) .AND. DSame(j1, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    END IF
                    m = l_ - step_aside ! I also try +step_aside, above
                    m = 1 + MOD(m-1+numEl, numEl) ! wraps around
                    j1 = nodes(1, m)
                    j2 = nodes(2, m)
                    j3 = nodes(3, m)
                    IF (DSame(j1, ib) .AND. DSame(j2, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    ELSE IF (DSame(j2, ib) .AND. DSame(j3, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    ELSE IF (DSame(j3, ib) .AND. DSame(j1, ia)) THEN
                        neighbor(j, l_) = m
                        EXIT strangers
                    END IF
                END DO strangers
            END DO sides
            IF (chatty) WRITE (*,"('+Learning the spherical triangles...',I6)") l_
        END DO homes
        IF (chatty) WRITE (*,"('+Learning the spherical triangles...DONE    ')")
        CONTAINS
            LOGICAL FUNCTION DSame(i,j)
               ! Are node_uvec #i and #j the same vector?
                INTEGER, INTENT(IN) :: i, j
                REAL*8, PARAMETER :: epsilon = 1.0D-5 ! 64 m on Earth
               !the logic is:
               !Same = (node_uvec(1,i) == node_uvec(1,j)).AND. &
               !     & (node_uvec(2,i) == node_uvec(2,j)).AND. &
               !     & (node_uvec(3,i) == node_uvec(3,j))
               !But, it is written this way for speed:
                IF (ABS(node_uvec(1,i) - node_uvec(1,j)) <= epsilon) THEN
                    IF (ABS(node_uvec(2,i) - node_uvec(2,j)) <= epsilon) THEN
                        IF (ABS(node_uvec(3,i) - node_uvec(3,j)) <= epsilon) THEN
                            DSame = .TRUE.
                        ELSE
                            DSame = .FALSE.
                        END IF
                    ELSE
                        DSame = .FALSE.
                    END IF
                ELSE
                    DSame = .FALSE.
                END IF
            END FUNCTION DSame
    END SUBROUTINE DLearn_Spherical_Triangles

    SUBROUTINE DPull_in(s)
        ! If necessary, adjusts internal coordinates s(1..3) so
        ! that none is negative.
        IMPLICIT NONE
        REAL*8, DIMENSION(3), INTENT(INOUT) :: s
        INTEGER, DIMENSION(1) :: array  ! stupid, to satisfy MINLOC
        REAL*8  :: factor, lowest, highest, medium
        INTEGER :: side, sidea, sideb
        lowest = MINVAL(s)
        IF (lowest < 0.0D0) THEN
            highest = MAXVAL(s)
            medium = 1.00D0 - lowest - highest
            IF (medium > 0.0D0) THEN ! correct to nearest edge
                array = MINLOC(s)
                side = array(1)
                s(side) = 0.0D0
                sidea = 1 + MOD(side, 3)
                sideb = 1 + MOD(sidea, 3)
                factor = 1.00D0 / (1.00D0 - lowest)
                s(sidea) = factor * s(sidea)
                ! s(sideb) = factor * s(sideb) would be logical
                s(sideb) = 1.00D0 - s(sidea)  ! is safer
            ELSE    ! correct to nearest vertex
                array = MAXLOC(s)
                side = array(1)
                s(side) = 1.00D0
                sidea = 1 + MOD(side, 3)
                sideb = 1 + MOD(sidea, 3)
                s(sidea) = 0.0D0
                s(sideb) = 0.0D0
            END IF
        END IF
    END SUBROUTINE DPull_in

    SUBROUTINE DWhich_Spherical_Triangle (b_, cold_start, &
                                        & num_ele, node, xyz_node, center, a_, neighbor, &
                                        & success, iEle, s1, s2, s3)
       !Locates a point (b_, a uvec) in element iEle with internal
       !coordinates (s1, s2, s3) in a SHELLS or RESTORE .feg.
       !and reports success.
       !If (cold_start), makes no use of input iEle, s1, s2, s3
       !If not, uses these values to initialize the search.
       !
       !Note that Learn_Spherical_Triangles can be used to initialize
       !necessary arrays a_, center, and neighbor.
       !
       !Beware of variable name changes (numEl = num_ele, nodes = node, etc.).
        IMPLICIT NONE
        REAL*8,  DIMENSION(3),  INTENT(IN)    :: b_         ! uvec of unknown point
        LOGICAL,                INTENT(IN)    :: cold_start ! mode switch
        INTEGER,                INTENT(IN)    :: num_ele    ! count of elements
        INTEGER, DIMENSION(:,:),INTENT(IN)    :: node       ! element definitions
        REAL*8,  DIMENSION(:,:),INTENT(IN)    :: xyz_node   ! uvecs of nodes
        REAL*8,  DIMENSION(:,:),INTENT(IN)    :: center     ! center uvecs of elements
        REAL*8,  DIMENSION(:),  INTENT(IN)    :: a_         ! (plane) areas of elements
        INTEGER, DIMENSION(:,:),INTENT(IN)    :: neighbor   ! neighbors of each element
        LOGICAL,                INTENT(OUT)   :: success    ! OUTPUT
        INTEGER,                INTENT(INOUT) :: iEle       ! OUTPUT
        REAL*8,                 INTENT(INOUT) :: s1, s2, s3 ! OUTPUT
       !- - - - - - - - - - - - - - - - - - - - - - - - - - - -
        INTEGER                :: back1, back2, back3, i, iet, l_
        REAL*8                 :: r2, r2min, s1t, s2t, s3t
        REAL*8,  DIMENSION(3)  :: s_temp, tv
      ! establish defaults (not found) in case of quick exit
        success = .FALSE.
        IF (cold_start) THEN
            iEle = 0  ! default
            s1 = 0.0D0; s2 = 0.0D0; s3 = 0.0D0 ! default
            !find closest element center to initialize search
            r2min = 4.01D0 ! radians
            DO l_ = 1, num_ele
                r2 = (b_(1) - center(1,l_))**2 +(b_(2) - center(2,l_))**2 +(b_(3) - center(3,l_))**2
                IF (r2 < r2min) THEN
                    r2min = r2
                    iet = l_
                END IF
            END DO
          ! If closest element center is more than 1 radian away, give uP.
            tv = center(1:3, iet)
            IF (DOT_PRODUCT(b_, tv) < 0.540D0) RETURN
        END IF ! cold_start
        ! initialize search memory (with impossible numbers)
        back1 = -1
        back2 = -2
        back3 = -3
        is_it_here: DO
            ! first, check for infinite loop between 2 elements!
            IF (iet == back2) THEN
                ! in loop; force location in one or the other!
                CALL DDumb_s123 (iet, b_, node, xyz_node, center, a_, &
                               & s1t, s2t, s3t)
                s_temp(1) = s1t; s_temp(2) = s2t; s_temp(3) = s3t
                CALL DPull_in(s_temp)
                s1t = s_temp(1); s2t = s_temp(2); s3t = s_temp(3)
                EXIT is_it_here
            ! then, check for infinite loop between 3 elements!
            ELSE IF (iet == back3) THEN
                ! in loop; force location in one or the other!
                CALL DDumb_s123 (iet, b_, node, xyz_node, center, a_, &
                               & s1t, s2t, s3t)
                s_temp(1) = s1t; s_temp(2) = s2t; s_temp(3) = s3t
                CALL DPull_in(s_temp)
                s1t = s_temp(1); s2t = s_temp(2); s3t = s_temp(3)
                EXIT is_it_here
            ELSE
                ! normal operation
                CALL DDumb_s123 (iet, b_, node, xyz_node, center, a_, &
                               & s1t, s2t, s3t)
                IF ((s1t < s2t) .AND. (s1t < s3t)) THEN ! s1 is most negative; most critical
                    IF (s1t >= 0.0D0) THEN
                        EXIT is_it_here ! success
                    ELSE
                        i = neighbor(1, iet)
                        IF (i > 0) THEN
                            back3 = back2
                            back2 = back1
                            back1 = iet
                            iet = i
                            CYCLE is_it_here
                        ELSE
                            RETURN  ! fell off edge of grid
                        ENDIF
                    ENDIF
                ELSE IF ((s2t < s1t) .AND. (s2t < s3t)) THEN    ! s2 is most negative; most critical
                    IF (s2t >= 0.0D0) THEN
                        EXIT is_it_here ! success
                    ELSE
                        i = neighbor(2, iet)
                        IF (i > 0) THEN
                            back3 = back2
                            back2 = back1
                            back1 = iet
                            iet = i
                            CYCLE is_it_here
                        ELSE
                            RETURN  ! fell off edge of grid
                        ENDIF
                    ENDIF
                ELSE    ! s3 is most negative; most critical
                    IF (s3t >= 0.0D0) THEN
                        EXIT is_it_here ! success
                    ELSE
                        i = neighbor(3, iet)
                        IF (i > 0) THEN
                            back3 = back2
                            back2 = back1
                            back1 = iet
                            iet = i
                            CYCLE is_it_here
                        ELSE
                            RETURN  ! fell off edge of grid
                        ENDIF
                    ENDIF
                END IF
            END IF ! in/not in a loop
        END DO is_it_here
        ! successful completion
        iEle = iet
        s1 = s1t
        s2 = s2t
        s3 = s3t
        success = .TRUE.
    END SUBROUTINE DWhich_Spherical_Triangle

end module DWeighting
