subroutine matrix_arrays(temp_x, lambda_x, inv_kappa_x, eps_x, p_x, lum_x, &
    Aarray_x, Barray_x, Carray_x, Darray_x)

  use blmod, only: delta_mass, delta_cmass, dtime, theta, &
        eps_p, p_p, lum, dedt, dpdt, rho, rho_p, r, Qterm, &
        bomb_heating, Ni_heating
  use parameters
  use physical_constants
  implicit none

!input
  real*8 :: temp_x(imax)
  real*8 :: lambda_x(imax)
  real*8 :: inv_kappa_x(imax)
  real*8 :: eps_x(imax), p_x(imax), lum_x(imax)


!output:
  real*8 :: Aarray_x(imax-1)
  real*8 :: Barray_x(imax-1)
  real*8 :: Carray_x(imax-1)
  real*8 :: Darray_x(imax-1)
  real*8 :: const
  real*8 :: inv_delta_mass
  real*8 :: inv_delta_cmass
  real*8 :: inv_delta_cmass_iM1
  real*8 :: r_iP1_4
  real*8 :: r_4
  real*8 :: lambda_inv_kappa
  real*8 :: temp_3
  real*8 :: inv_rho
  real*8 :: inv_rho_p

!local:
  integer :: i

!------------------------------------------------------------------------------
  !** common constant used in all terms
  const = -theta * dtime * 64.0d0 * pi * pi * a_rad * clite / 3.0d0

  !************* A array **************** coefficients of \delta T(i+1)
  do i=1, imax-2
     inv_delta_mass = 1.0d0 / delta_mass(i)
     inv_delta_cmass = 1.0d0 / delta_cmass(i)
     r_iP1_4 = r(i+1)**4
     r_4 = r(i)**4
     lambda_inv_kappa = lambda_x(i) * inv_kappa_x(i)
     temp_3 = temp_x(i)**3
     inv_rho = 1.0d0 / rho(i)
     inv_rho_p = 1.0d0 / rho_p(i)


     Aarray_x(i) = const*inv_delta_mass * &
          r_iP1_4 * &
          (lambda_x(i+1)*inv_kappa_x(i+1) * &
          temp_x(i+1)**3*inv_delta_cmass )

  !************* B array **************** coefficients of \delta T(i)

     if (i.gt.1) then
       inv_delta_cmass_iM1 = 1.0d0 / delta_cmass(i-1)

       Barray_x(i) = dedt(i) + 0.5d0*dpdt(i)* &
            (inv_rho-inv_rho_p) &

            - const*inv_delta_mass * &
            r_iP1_4  * &
            (lambda_x(i+1)*inv_kappa_x(i+1)* &
           temp_3*inv_delta_cmass ) &

            - const*inv_delta_mass * &
            r_4 * &
            (lambda_inv_kappa* &
            temp_3*inv_delta_cmass_iM1 )

  !************* C array **************** coefficients of \delta T(i-1)
       Carray_x(i) = const*inv_delta_mass * r_4 * &
            (lambda_inv_kappa*temp_x(i-1)**3*inv_delta_cmass_iM1 )

  !************* D array **************** free terms
       Darray_x(i) = &

  !terms from the previous time step
            eps_p(i) - 0.5d0*p_p(i)*(inv_rho-inv_rho_p) &

            - (1-theta)*(lum(i+1)-lum(i))*dtime*inv_delta_mass &

            + Qterm(i) &

  !terms from the intermediate point (with eps_temp, p_temp and temp_temp)
            - ( eps_x(i) + 0.5d0*p_x(i)*(inv_rho-inv_rho_p) &

            + theta*dtime*inv_delta_mass * ( lum_x(i+1)-lum_x(i) ) ) &

  !bomb heating, nickel heating
            + bomb_heating(i)*dtime + Ni_heating(i)*dtime
     end if
  end do

  !********** Outermost and innermost zones **********************************
  !************* A array **************** coefficients of \delta T(i+1)
  Aarray_x(imax-1) = 0.0d0



  !************* B array **************** coefficients of \delta T(i)
  Barray_x(1) = dedt(1) + 0.5d0*dpdt(1)*(1.0d0/rho(1)-1.0d0/rho_p(1)) &
         - const/delta_mass(1) * r(2)**4 * &
        (lambda_x(2)*inv_kappa_x(2)*temp_x(1)**3/delta_cmass(1) )


  Barray_x(imax-1) = & !outer boundary condition L(imax) = L(imax-1)
       dedt(imax-1) + 0.5d0*dpdt(imax-1)* &
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1))


  !************* C array **************** coefficients of \delta T(i-1)

  Carray_x(1) = 0.0d0

  Carray_x(imax-1) = 0.0d0 !outer boundary condition L(imax) = L(imax-1)


  !************* D array **************** free terms


  Darray_x(1) = eps_p(1) - 0.5d0*p_p(1)*(1.0d0/rho(1)-1.0d0/rho_p(1)) &

    - (1-theta)*lum(2)*dtime/delta_mass(1) &

    + Qterm(1) &

    - ( eps_x(1) + 0.5d0*p_x(1)*(1.0d0/rho(1)-1.0d0/rho_p(1)) &

    + theta*dtime/delta_mass(1) * lum_x(2) ) &

    + bomb_heating(1)*dtime + Ni_heating(1)*dtime

  Darray_x(imax-1) = & !outer boundary condition L(imax) = L(imax-1)

       eps_p(imax-1) - 0.5d0*p_p(imax-1)*&
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1)) &

       + Qterm(imax-1) &

       - (eps_x(imax-1)+0.5d0*p_x(imax-1)* &
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1)))&

       + bomb_heating(imax-1)*dtime + Ni_heating(imax-1)*dtime


end subroutine matrix_arrays

