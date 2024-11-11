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


!def interp_linear(x,x0,y0,x1,y1):
!
!    if (x<x0) or (x>x1):
!        y = np.nan
!    else:
!        y = (y0*(x1-x) + y1*(x-x0))/(x1-x0)
!    return  y
!
!def find_nearest(x, xs):
!    return np.argmin(np.abs(xs-x))
!
!
!def find_second_nearest(x,xs):
!    xaux = xs.copy()
!    i    = find_nearest(x,xs)
!    xaux[i] = 1e10
!    i1 = find_nearest(x,xaux)
!    return i1
!
!def interp1d(x, xs, ys):
!    if len(xs) == len(ys):
!        pass
!    else:
!        raise IOError()
!
!    i0 = find_nearest(x,xs)
!    i1 = find_second_nearest(x,xs)
!
!    if i0 > i1:
!        aux=i1
!        i1 = i0
!        i0 = aux
!
!    x0 = xs[i0]
!    x1 = xs[i1]
!    y0 = ys[i0]
!    y1 = ys[i1]
!    return interp_linear(x,x0,y0,x1,y1)
!    
!
!def interp_along_axis(xs, x0,y0):
!    ys = np.zeros(xs.size) * np.nan
!    for i,x in enumerate(xs):
!        ys[i] = interp1d(x,x0,y0)
!
!    return ys
!        
!
!def interp3d_along_axis0(zs, x0, y0, z0, varb0):
!    shape = x0.shape
!    varbs = np.zeros(zs.shape)
!    if (y0.shape!=shape) or (z0.shape!=shape) or (varb0.shape!=shape):
!        raise IOError
!
!    for i in range(shape[1]):
!        for j in range(shape[2]):
!            varbs[:,i,j] =  interp_along_axis(zs[:,i,j], z0[:,i,j], varb0[:,i,j])
!
!    return varbs
!    
!
!x = np.arange(10)
!y = np.arange(20)
!z = np.arange(5)
!zs = np.arange(0.5,3.5)
!
!shape = (5, 20, 10)
!# Create a random array with values between 0 and 1
!random_array = np.random.rand(*shape)
!
!
!zm, ym, xm = np.meshgrid(z,y,x, indexing='ij')
!
!
!zsm, _, _ = np.meshgrid(zs,y,x, indexing='ij')
!
!a = interp3d_along_axis0(zsm,xm,ym,zm,random_array)
!x = np.arange(10)
!y = np.arange(20)
!z = np.arange(5)
!zs = np.arange(0.5,3.5,0.1)
!
!shape = (5, 20, 10)
!# Create a random array with values between 0 and 1
!random_array = np.random.rand(*shape)
!
!
!zm, ym, xm = np.meshgrid(z,y,x, indexing='ij')
!
!
!zsm, _, _ = np.meshgrid(zs,y,x, indexing='ij')
!
!a = interp3d_along_axis0(zsm,xm,ym,zm,random_array)
!
!plt.close('all')
!fig, ax = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True)
!ax[0].contourf(y,zs, a[:,:,0], vmin=-1, vmax=1)
!ax[1].contourf(y,z, random_array[:,:,0], vmin=-1, vmax=1)
!plt.savefig('bla.png')
!
!# fig, ax = plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
!# ax.contourf(y,z, random_array[:,:,0], vmin=-1, vmax=1)
!# ax.contour(y,zs, a[:,:,0], vmin=-1, vmax=1)
!# plt.savefig('bla.png')