! File: interpolation.f90
MODULE interpolation
    ! Interpolation functions by Dalton Kei Sasaki 2023/08/18
    ! Linear interpolation module to be used with f2py
    ! f2py -c -m interpolation interpolation.f90
    
        public :: interp_along_axis
        CONTAINS
    
        REAL FUNCTION  interp_linear(x, x0, y0, x1, y1)
            real, intent(in) :: x, x0, y0, x1, y1
            real :: y
    
            if (x < x0 .OR. x > x1) THEN
                y = 9999.0
            else
                y = (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
            end if
            interp_linear = y
        END FUNCTION interp_linear
    
        FUNCTION find_nearest(x, xs)
            real, intent(in) :: x
            real, dimension(:), intent(in) :: xs
            INTEGER :: i
    
            i = MINLOC(ABS(xs - x),1)
            find_nearest = i
        END FUNCTION find_nearest
    
        FUNCTION Find_Second_Nearest(x, xs)
            real, intent(in) :: x
            real, dimension(:), intent(in) :: xs
            real, dimension(SIZE(xs)) :: xaux
            INTEGER :: i, i1
    
            xaux = xs
            i = find_nearest(x, xaux)
            xaux(i) = 1e10
    
            i1 = find_nearest(x, xaux)
            Find_Second_Nearest = i1
        END FUNCTION Find_Second_Nearest
    
        REAL FUNCTION interp1d(x, xs, ys)
            real, intent(in) :: x
            real, dimension(:), intent(in) :: xs, ys
            INTEGER :: i0, i1
            real :: x0, x1, y0, y1
    
            if (SIZE(xs) /= SIZE(ys)) THEN
                STOP 'IOError'
            end if
    
            i0 = find_nearest(x, xs)
            i1 = Find_Second_Nearest(x, xs)
            if (i0 > i1) THEN
                aux = i0
                i0  = i1
                i1  = aux
            end if
    
    
            x0 = xs(i0)
            x1 = xs(i1)
            y0 = ys(i0)
            y1 = ys(i1)
            interp1d = interp_linear(x, x0, y0, x1, y1)
        END FUNCTION interp1d
    
        FUNCTION interp_along_axis(xs, x0, y0) result(ys)
            real, dimension(:), intent(in) :: xs
            real, dimension(:), intent(in) :: x0
            real, dimension(:), intent(in) :: y0
            real, dimension(SIZE(xs)) :: ys
    
            integer :: i
            real :: y
    
            DO i = 1, SIZE(xs)
                y = interp1d(xs(i), x0, y0)
                ys(i) = y
            END DO
    
        END FUNCTION interp_along_axis
    
        FUNCTION interp3d_along_axis0(zs, x0, y0, z0, varb0,nzs, nx,ny,nz) result(varbs)
            integer, intent(in) :: nx,ny,nz,nzs
            real, dimension(nzs,ny,nx), intent(in) :: zs
            real, dimension(nz,ny,nx), intent(in) :: x0, y0, z0, varb0
            real, dimension(nzs,ny,nx) :: varbs
    
            integer :: i, j, k
    
    
    
            DO i = 1, ny
                DO j = 1, nx
                    varbs(:, i, j) = interp_along_axis(zs(:, i, j), z0(:, i, j), varb0(:, i, j))
                END DO
            END DO
        END FUNCTION interp3d_along_axis0
    
    
    
    END MODULE interpolation
    
    