subroutine bicubic_interpolation(l,logR_x,logT_x,log_kappa_x,no_tabular_value_x)

  use blmod, only: logR_array, logT_array, opacity_tables, rpoints, tpoints, &
    op_rows, op_cols
  use parameters
  implicit none

!input:
  integer :: l !gridpoint
  real*8 :: logR_x, logT_x

!output:
  real*8 :: log_kappa_x
  integer :: no_tabular_value_x ! =1 if there is no tabulated value of log kappa
                                !for given log T and log R

!local:
  integer :: i,j
  integer :: ibuffer

  !16 values of the log kappa from the table, used for interpolation
  real*8 :: kim1jm1, kim1j, kim1jp1, kim1jp2
  real*8 :: kijm1,   kij,   kijp1,   kijp2
  real*8 :: kip1jm1, kip1j, kip1jp1, kip1jp2
  real*8 :: kip2jm1, kip2j, kip2jp1, kip2jp2

  real*8 :: dt, dr
  real*8 :: inv_dti, inv_dti1, inv_drj, inv_drj1

  real*8 :: kappa(2,2)
  real*8 :: dkdt(2,2)
  real*8 :: dkdr(2,2)
  real*8 :: ddkdtdr(2,2)

  real*8 :: x,y

  real*8 :: yterm1, yterm2, yterm3, yterm4
  real*8 :: xterm1, xterm2, xterm3, xterm4

  real*8 :: ki, ki1, dkdti, dkdti1

!------------------------------------------------------------------------------

  no_tabular_value_x = 0

  if (logT_x.lt.logT_array(1) .or. logT_x.gt.logT_array(tpoints) &
    .or. logR_x.lt.logR_array(1) .or. logR_x.gt.logR_array(rpoints)) then
    !the point (log T, log R) is out of the table
    no_tabular_value_x = 1
    log_kappa_x = 0.0d0
    goto 255
  endif

  ! i and j are the indices of the left upper corner of the square, where the
  !value (logT,logR) is located
  call map_find_index(tpoints,logT_array,logT_x,ibuffer,i)
  call map_find_index(rpoints,logR_array,logR_x,ibuffer,j)

  if (i.lt.op_cols(j,1) .or. i.gt.(op_cols(j+1,2)-1) &
    .or. j.lt.op_rows(i,1) .or. j.gt.(op_rows(i+1,2)-1)) then
    !the point (log T, log R) is out of the meaningful region of the table
    no_tabular_value_x = 1
    log_kappa_x = 0.0d0
    goto 255
  endif

  kij     = opacity_tables(l,i,j)
  kip1j   = opacity_tables(l,i+1,j)
  kijp1   = opacity_tables(l,i,j+1)
  kip1jp1 = opacity_tables(l,i+1,j+1)

  !values of log kappa in the corners of the square
  kappa(1,1) = kij
  kappa(2,1) = kip1j
  kappa(1,2) = kijp1
  kappa(2,2) = kip1jp1

  !values of the derivatives (D log kappa/D log T) in the corners of the square
  !(constructed from the available log kappa values)
  if(i.eq.op_cols(j,1)) then
    inv_dti1  = 1.0d0/(logT_array(i+2)-logT_array(i))
    kip2j     = opacity_tables(l,i+2,j)
    kip2jp1   = opacity_tables(l,i+2,j+1)
    dkdt(2,1) = (kip2j - kij)*inv_dti1
    dkdt(2,2) = (kip2jp1 - kijp1)*inv_dti1
    dkdt(1,1) = dkdt(2,1)
    dkdt(1,2) = dkdt(2,2)
  else if(i.eq.(op_cols(j+1,2)-1)) then
    inv_dti   = 1.0d0/(logT_array(i+1)-logT_array(i-1))
    kim1j     = opacity_tables(l,i-1,j)
    kim1jp1   = opacity_tables(l,i-1,j+1)
    dkdt(1,1) = (kip1j - kim1j)*inv_dti
    dkdt(1,2) = (kip1jp1 - kim1jp1)*inv_dti
    dkdt(2,1) = dkdt(1,1)
    dkdt(2,2) = dkdt(1,2)
  else
    inv_dti  = 1.0d0/(logT_array(i+1)-logT_array(i-1))
    inv_dti1 = 1.0d0/(logT_array(i+2)-logT_array(i))
    kip2j     = opacity_tables(l,i+2,j)
    kip2jp1   = opacity_tables(l,i+2,j+1)
    kim1j     = opacity_tables(l,i-1,j)
    kim1jp1   = opacity_tables(l,i-1,j+1)
    dkdt(2,1) = (kip2j - kij)*inv_dti1
    dkdt(2,2) = (kip2jp1 - kijp1)*inv_dti1
    dkdt(1,1) = (kip1j - kim1j)*inv_dti
    dkdt(1,2) = (kip1jp1 - kim1jp1)*inv_dti
  end if

  !values of the derivatives (D log kappa/D log R) in the corners of the square
  !(constructed from the available log kappa values)
  if(j.eq.op_rows(i,1)) then
    inv_drj1  = 1.0d0/(logR_array(j+2)-logR_array(j))
    kijp2     = opacity_tables(l,i,j+2)
    kip1jp2   = opacity_tables(l,i+1,j+2)
    dkdr(1,2) = (kijp2 - kij)*inv_drj1
    dkdr(2,2) = (kip1jp2 - kip1j)*inv_drj1
    dkdr(1,1) = dkdr(1,2)
    dkdr(2,1) = dkdr(2,2)
  else if(j.eq.(op_rows(i+1,2)-1)) then
    inv_drj   = 1.0d0/(logR_array(j+1)-logR_array(j-1))
    kijm1     = opacity_tables(l,i,j-1)
    kip1jm1   = opacity_tables(l,i+1,j-1)
    dkdr(1,1) = (kijp1 - kijm1)*inv_drj
    dkdr(2,1) = (kip1jp1 - kip1jm1)*inv_drj
    dkdr(1,2) = dkdr(1,1)
    dkdr(2,2) = dkdr(2,1)
  else
    inv_drj  = 1.0d0/(logR_array(j+1)-logR_array(j-1))
    inv_drj1 = 1.0d0/(logR_array(j+2)-logR_array(j))
    kijp2     = opacity_tables(l,i,j+2)
    kip1jp2   = opacity_tables(l,i+1,j+2)
    kijm1     = opacity_tables(l,i,j-1)
    kip1jm1   = opacity_tables(l,i+1,j-1)
    dkdr(1,2) = (kijp2 - kij)*inv_drj1
    dkdr(2,2) = (kip1jp2 - kip1j)*inv_drj1
    dkdr(1,1) = (kijp1 - kijm1)*inv_drj
    dkdr(2,1) = (kip1jp1 - kip1jm1)*inv_drj
  end if

  !values of the mixed derivatives (D2 log kappa/D log T/D log R)
  !(constructed from the available log kappa values)
  if(j.eq.op_rows(i,1) .and. i.eq.op_cols(j,1)) then
    kip2jp2 = opacity_tables(l,i+2,j+2)
    ddkdtdr(2,2) = (kip2jp2-kip2j-kijp2+kij)*inv_dti1*inv_drj1
    ddkdtdr(1,1) = ddkdtdr(2,2)
    ddkdtdr(1,2) = ddkdtdr(2,2)
    ddkdtdr(2,1) = ddkdtdr(2,2)
  else if(j.eq.op_rows(i,1) &
                    .and. i.ne.op_cols(j,1) .and. i.ne.(op_cols(j+1,2)-1)) then
    kim1jp2 = opacity_tables(l,i-1,j+2)
    kip2jp2 = opacity_tables(l,i+2,j+2)
    ddkdtdr(1,2) = (kip1jp2-kip1j-kim1jp2+kim1j)*inv_dti*inv_drj1
    ddkdtdr(2,2) = (kip2jp2-kip2j-kijp2+kij)*inv_dti1*inv_drj1
    ddkdtdr(1,1) = ddkdtdr(1,2)
    ddkdtdr(2,1) = ddkdtdr(2,2)
  else if(j.eq.op_rows(i,1) .and. i.eq.(op_cols(j+1,2)-1)) then
    kim1jp2 = opacity_tables(l,i-1,j+2)
    ddkdtdr(1,2) = (kip1jp2-kip1j-kim1jp2+kim1j)*inv_dti*inv_drj1
    ddkdtdr(1,1) = ddkdtdr(1,2)
    ddkdtdr(2,1) = ddkdtdr(1,2)
    ddkdtdr(2,2) = ddkdtdr(1,2)
  else if(i.eq.(op_cols(j+1,2)-1) &
                    .and. j.ne.op_rows(i,1) .and. j.ne.(op_rows(i+1,2)-1)) then
    kim1jm1 = opacity_tables(l,i-1,j-1)
    kim1jp2 = opacity_tables(l,i-1,j+2)
    ddkdtdr(1,1) = (kip1jp1-kip1jm1-kim1jp1+kim1jm1)*inv_dti*inv_drj
    ddkdtdr(1,2) = (kip1jp2-kip1j-kim1jp2+kim1j)*inv_dti*inv_drj1
    ddkdtdr(2,1) = ddkdtdr(1,1)
    ddkdtdr(2,2) = ddkdtdr(1,2)
  else if(i.eq.(op_cols(j+1,2)-1) .and. j.eq.(op_rows(i+1,2)-1)) then
    kim1jm1 = opacity_tables(l,i-1,j-1)
    ddkdtdr(1,1) = (kip1jp1-kip1jm1-kim1jp1+kim1jm1)*inv_dti*inv_drj
    ddkdtdr(1,2) = ddkdtdr(1,1)
    ddkdtdr(2,1) = ddkdtdr(1,1)
    ddkdtdr(2,2) = ddkdtdr(1,1)
  else if(j.eq.(op_rows(i+1,2)-1) &
                    .and. i.ne.(op_cols(j+1,2)-1) .and. i.ne.op_cols(j,1)) then
    kim1jm1 = opacity_tables(l,i-1,j-1)
    kip2jm1 = opacity_tables(l,i+2,j-1)
    ddkdtdr(1,1) = (kip1jp1-kip1jm1-kim1jp1+kim1jm1)*inv_dti*inv_drj
    ddkdtdr(2,1) = (kip2jp1-kip2jm1-kijp1+kijm1)*inv_dti1*inv_drj
    ddkdtdr(1,2) = ddkdtdr(1,1)
    ddkdtdr(2,2) = ddkdtdr(2,1)
  else if(j.eq.(op_rows(i+1,2)-1) .and. i.eq.op_cols(j,1)) then
    kip2jm1 = opacity_tables(l,i+2,j-1)
    ddkdtdr(2,1) = (kip2jp1-kip2jm1-kijp1+kijm1)*inv_dti1*inv_drj
    ddkdtdr(1,1) = ddkdtdr(2,1)
    ddkdtdr(1,2) = ddkdtdr(2,1)
    ddkdtdr(2,2) = ddkdtdr(2,1)
  else if(i.eq.op_cols(j,1) &
                    .and. j.ne.(op_rows(i+1,2)-1) .and. j.ne.op_rows(i,1)) then
    kip2jm1 = opacity_tables(l,i+2,j-1)
    kip2jp2 = opacity_tables(l,i+2,j+2)
    ddkdtdr(2,1) = (kip2jp1-kip2jm1-kijp1+kijm1)*inv_dti1*inv_drj
    ddkdtdr(2,2) = (kip2jp2-kip2j-kijp2+kij)*inv_dti1*inv_drj1
    ddkdtdr(1,1) = ddkdtdr(2,1)
    ddkdtdr(1,2) = ddkdtdr(2,2)
  else
    kim1jm1 = opacity_tables(l,i-1,j-1)
    kim1jp2 = opacity_tables(l,i-1,j+2)
    kip2jm1 = opacity_tables(l,i+2,j-1)
    kip2jp2 = opacity_tables(l,i+2,j+2)
    ddkdtdr(1,1) = (kip1jp1-kip1jm1-kim1jp1+kim1jm1)*inv_dti*inv_drj
    ddkdtdr(1,2) = (kip1jp2-kip1j-kim1jp2+kim1j)*inv_dti*inv_drj1
    ddkdtdr(2,1) = (kip2jp1-kip2jm1-kijp1+kijm1)*inv_dti1*inv_drj
    ddkdtdr(2,2) = (kip2jp2-kip2j-kijp2+kij)*inv_dti1*inv_drj1
  end if

    dt = logT_array(i+1) - logT_array(i)
    dr = logR_array(j+1) - logR_array(j)

    !dimensionless coordinates of the given point within the square
    x = (logT_x - logT_array(i))/dt
    y = (logR_x - logR_array(j))/dr

    !terms constructed from the cubic Hermite basis functions:
    yterm1 = y * y * (2.0d0*y - 3.0d0) + 1.0d0
    yterm2 = (1.0d0-y) * (1.0d0-y) * (2.0d0*(1.0d0-y) - 3.0d0) + 1.0d0
    yterm3 = dr * y * ( y * (y - 2.0d0) + 1.0d0)
    yterm4 = dr * (1.0d0-y) * ( (1.0d0-y) * ((1.0d0-y) - 2.0d0) + 1.0d0)

    !interpolation in log R direction
    ki = kappa(1,1)*yterm1 + kappa(1,2)*yterm2 &
                                    + dkdr(1,1)*yterm3 - dkdr(1,2)*yterm4

    ki1 =kappa(2,1)*yterm1 + kappa(2,2)*yterm2 &
                                    + dkdr(2,1)*yterm3 - dkdr(2,2)*yterm4

    dkdti = dkdt(1,1)*yterm1 + dkdt(1,2)*yterm2 &
                                    + ddkdtdr(1,1)*yterm3 - ddkdtdr(1,2)*yterm4

    dkdti1 = dkdt(2,1)*yterm1 + dkdt(2,2)*yterm2 &
                                    + ddkdtdr(2,1)*yterm3 - ddkdtdr(2,2)*yterm4

    !terms constructed from the cubic Hermite basis functions:
    xterm1 = x * x * (2.0d0*x - 3.0d0) + 1.0d0
    xterm2 = (1.0d0-x) * (1.0d0-x) * (2.0d0*(1.0d0-x) - 3.0d0) + 1.0d0
    xterm3 = dt * x * ( x * (x - 2.0d0) + 1.0d0)
    xterm4 = dt * (1.0d0-x) * ( (1.0d0-x) * ((1.0d0-x) - 2.0d0) + 1.0d0)

    !interpolation in log T direction
    log_kappa_x = ki*xterm1 + ki1*xterm2 + dkdti*xterm3 - dkdti1*xterm4

255 continue

end subroutine bicubic_interpolation
