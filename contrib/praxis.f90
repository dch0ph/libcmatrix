function flin ( n, j, l, f, x, nf, v, q0, q1 )
!
!*******************************************************************************
!
!! FLIN is the function of one variable to be minimized by MINNY.
!
!
!  Discussion:
!
!    In fact, what is happening is that the scalar function F(X),
!    where X is an N dimensional vector, is being minimized along a 
!    fixed line.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
  integer n
!
  real f
  real flin
  integer i
  integer j
  real l
  integer nf
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real t(n)
  real v(n,n)
  real x(n)
!
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  if ( j /= 0 ) then
!
!  The search is linear.
!
    do i = 1, n
      t(i) = x(i) + l * v(i,j)
    end do

  else
!
!  The search is along a parabolic space curve.
!
    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

    do i = 1, n
      t(i) = ( qa * q0(i) + qb * x(i) ) + qc * q1(i)
    end do

  end if
!
!  The function evaluation counter NF is incremented.
!
  nf = nf + 1
  flin = f(t,n)

  return
end
subroutine minfit ( m, n, machep, tol, ab, q )
!
!*******************************************************************************
!
!! MINFIT computes the singular value decomposition of an array.
!
!
!  Discussion:
!
!   This is an improved version of minfit (see golub and reinsch, 1969)
!   restricted to m = n, p = 0.
!
!   The singular values of the array ab are returned in q and ab is
!   overwritten with the orthogonal matrix V such that u.diag(q) = ab.v,
!   where u is another orthogonal matrix.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
  integer m
  integer n
!
  real ab(m,n)
  real c
  real e(n)
  real eps
  real f
  real g
  real h
  integer i
  integer ii
  integer j
  integer k
  integer kk
  integer kt
  integer l
  integer l2
  integer ll2
  real machep
  real q(n)
  real s
  real temp
  real tol
  real x
  real y
  real z
!
!  Householder's reduction to bidiagonal form.
!
  if ( n == 1 ) then
    q(1) = ab(1,1)
    ab(1,1) = 1.0
    return
  end if

  eps = machep
  g = 0.0
  x = 0.0

  do i = 1, n

    e(i) = g
    s = 0.0
    l = i + 1

    do j = i, n
      s = s + ab(j,i)**2
    end do

    g = 0.0

    if ( s >= tol ) then

      f = ab(i,i)

      g = sqrt ( s )
      if ( f >= 0.0 ) then
        g = -g
      end if

      h = f * g - s
      ab(i,i) = f - g

      do j = l, n

        f = 0.0
        do k = i, n
          f = f + ab(k,i) * ab(k,j)
        end do
        f = f / h
        do k = i, n
          ab(k,j) = ab(k,j) + f * ab(k,i)
        end do

      end do 

    end if

    q(i) = g

    s = 0.0
    do j = l, n
      s = s + ab(i,j)**2
    end do

    g = 0.0

    if ( s >= tol ) then

      if ( i /= n ) then
        f = ab(i,i+1)
      end if

      g = sqrt ( s )
      if ( f >= 0.0 ) then
        g = - g
      end if

      h = f * g - s

      if ( i /= n ) then

        ab(i,i+1) = f - g

        do j = l, n
          e(j) = ab(i,j) / h
        end do

        do j = l, n

          s = 0.0
          do k = l, n
            s = s + ab(j,k) * ab(i,k)
          end do

          do k = l, n
            ab(j,k) = ab(j,k) + s * e(k)
          end do

        end do

      end if

    end if

    y = abs ( q(i) ) + abs ( e(i) )

  end do

  x = max ( x, y )
!
!  Accumulation of right-hand transformations.
!
  ab(n,n) = 1.0
  g = e(n)
  l = n

  do ii = 2, n

    i = n - ii + 1

    if ( g /= 0.0 ) then

      h = ab(i,i+1) * g

      do j = l, n
        ab(j,i) = ab(i,j) / h
      end do

      do j = l, n

        s = 0.0
        do k = l, n
          s = s + ab(i,k) * ab(k,j)
        end do

        do k = l, n
          ab(k,j) = ab(k,j) + s * ab(k,i)
        end do

      end do

    end if

    ab(i,l:n) = 0.0
    ab(l:n,i) = 0.0
    ab(i,i) = 1.0

    g = e(i)

  end do

  l = i
!
!  Diagonalization of the bidiagonal form.
!
100   continue

  eps = eps * x

  do kk = 1, n

     k = n - kk + 1
     kt = 0

101  continue

     kt = kt + 1

     if ( kt > 30 ) then
       e(k) = 0.0
       write ( *, * ) ' '
       write ( *, * ) 'MINFIT - Fatal error!'
       write ( *, * ) '  QR failed to converge.'
       stop
     end if

     do ll2 = 1, k

       l2 = k - ll2 + 1
       l = l2
       if ( abs ( e(l) ) <= eps ) then
         go to 120
       end if

       if ( l /= 1 ) then
         if ( abs ( q(l-1) ) <= eps ) then
           go to 110
         end if
       end if

     end do
!
!  Cancellation of E(L) if L>1.
!
110  continue

     c = 0.0
     s = 1.0

     do i = l, k

       f = s * e(i)
       e(i) = c * e(i)
       if ( abs ( f ) <= eps ) then
         go to 120
       end if
       g = q(i)
!
!  q(i) = h = sqrt(g*g + f*f).
!
       if ( abs ( f ) < abs ( g ) ) then
         h = abs ( g ) * sqrt ( 1.0 + ( f / g )**2 )
       else if ( f == 0.0 ) then
         h = 0.0
       else
         h = abs ( f ) * sqrt ( 1.0 + ( g / f )**2 )
       end if

       q(i) = h

       if ( h == 0.0 ) then
         g = 1.0
         h = 1.0
       end if

       c = g / h
       s = - f / h

     end do
!
!  Test for convergence.
!
120      continue

     z = q(k)
     if ( l == k ) then
       go to 140
     end if
!
!  Shift from bottom 2*2 minor.
!
     x = q(l)
     y = q(k-1)
     g = e(k-1)
     h = e(k)
     f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y )

     g = sqrt ( f * f + 1.0 )

     if ( f < 0.0 ) then
       temp = f - g
     else
       temp = f + g
     end if

     f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
!
!  Next QR transformation.
!
     c = 1.0
     s = 1.0

     do i = l+1, k

        g = e(i)
        y = q(i)
        h = s * g
        g = g * c

        if ( abs ( f ) < abs ( h ) ) then
          z = abs ( h ) * sqrt ( 1.0 + ( f / h )**2 )
        else if ( f == 0.0 ) then
          z = 0.0
        else
          z = abs ( f ) * sqrt ( 1.0 + ( h / f )**2 )
        end if

        e(i-1) = z

        if ( z == 0.0 ) then
          f = 1.0
          z = 1.0
        end if

        c = f / z
        s = h / z
        f =   x * c + g * s
        g = - x * s + g * c
        h = y * s
        y = y * c

        do j = 1, n
          x = ab(j,i-1)
          z = ab(j,i)
          ab(j,i-1) = x * c + z * s
          ab(j,i) = - x * s + z * c
        end do

        if ( abs ( f ) < abs ( h ) ) then
          z = abs ( h ) * sqrt ( 1.0 + ( f / h ) **2 )
        else if ( f == 0.0 ) then
          z = 0.0
        else
          z = abs ( f ) * sqrt ( 1.0 + ( h / f )**2 )
        end if

        q(i-1) = z

        if ( z == 0.0 ) then
          f = 1.0
          z = 1.0
        end if

        c = f / z
        s = h / z
        f = c * g + s * y
        x = - s * g + c * y

     end do

     e(l) = 0.0
     e(k) = f
     q(k) = x
     go to 101
!
!  Convergence:  Q(K) is made non-negative.
!
140      continue

    if ( z < 0.0) then
      q(k) = - z
      ab(1:n,k) = - ab(1:n,k)
    end if

  end do

  return
end
subroutine minny ( n, j, nits, d2, x1, f1, fk, f, x, t, machep, h, v, q0, q1 )
!
!*******************************************************************************
!
!! MINNY minimizes a scalar function of N variables along a line.
!
!
!  Discussion:
!
!    MINNY minimizes F from x in the direction v(*,j) unless
!    j is less than 1, when a quadratic search is made in the plane
!    defined by q0,q1,x.
!
!    d2 is either zero or an approximation to half f".
!    on entry, x1 is an estimate of the distance from x to the minimum
!    along v(*,j) (or, if j = 0, a curve).  on return, x1 is the distance
!    found.
!    if fk = .true., then f1 is flin(x1).  otherwise x1 and f1 are ignored
!    on entry unless final fx is greater than f1.
!    nits controls the number of times an attempt will be made to halve
!    the interval.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
  real d1
  real d2
  real dmin
  logical dz
  real f
  real f0
  real f1
  real f2
  logical fk
  real flin
  real fm
  real fx
  real h
  integer i
  integer j
  integer k
  real ldt
  real m2
  real m4
  real machep
  integer n
  integer nf
  integer nits
  integer nl
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real s
  real sf1
  real small
  real sx1
  real t
  real t2
  real temp
  real v(n,n)
  real x(n)
  real x1
  real x2
  real xm
!
  external f
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  small = machep**2
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
  sf1 = f1
  sx1 = x1
  k = 0
  xm = 0.0
  fm = fx
  f0 = fx
  dz = ( d2 < machep )
!
!  Find the step size.
! 
  s = 0.0
  do i = 1, n
    s = s + x(i)**2
  end do
  s = sqrt ( s )

  if ( dz ) then
    temp = dmin
  else
    temp = d2
  end if

  t2 = m4 * sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
  s = m4 * s + t
  if ( dz .and. t2 > s ) then
    t2 = s
  end if

  t2 = max ( t2, small )
  t2 = min ( t2, 0.01 * h )

  if ( fk .and. f1 <= fm ) then
    xm = x1
    fm = f1
  end if

  if ( .not. fk .or. abs ( x1 ) < t2 ) then

    if ( x1 >= 0.0 ) then
      temp = 1.0
    else
      temp = - 1.0
    end if

    x1 = temp * t2
    f1 = flin ( n, j, x1, f, x, nf, v, q0, q1 )

  end if

  if ( f1 <= fm ) then
    xm = x1
    fm = f1
  end if
!
!  Evaluate FLIN at another point and estimate the second derivative.
!
4 continue

  if ( dz ) then

    if ( f0 >= f1 ) then
      x2 = 2.0*x1
    else
      x2 = - x1
    end if

    f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )

    if ( f2 <= fm ) then
      xm = x2
      fm = f2
    end if

    d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )

  end if
!
!  Estimate the first derivative at 0.
!
  d1 = ( f1 - f0 ) / x1 - x1 * d2
  dz = .true.
!
!  Predict the minimum.
!
  if ( d2 <= small ) then

    if ( d1 >= 0.0 ) then
      x2 = - h
    else
      x2 = h
    end if

  else

     x2 = ( - 0.5 * d1 ) / d2

  end if

  if ( abs ( x2 ) > h ) then

    if ( x2 <= 0.0 ) then
      x2 = - h
    else
      x2 = h
    end if

  end if
!
!  Evaluate F at the predicted minimum.
!
11    continue

  f2 = flin ( n, j, x2, f, x, nf, v, q0, q1 )
!
!  No success, so try again.
!
  if ( k < nits .and. f2 > f0 ) then
     k = k + 1
     if ( f0 < f1 .and. x1 * x2 > 0.0 ) go to 4
     x2 = 0.5 * x2
     go to 11
  end if
!
!  Increment the one-dimensional search counter.
!
  nl = nl + 1

  if ( f2 > fm ) then
    x2 = xm
  else
    fm = f2
  end if
!
!  Get a new estimate of the second derivative.
!
  if ( abs ( x2 * ( x2 - x1 ) ) > small ) then
    d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
  else
    if ( k > 0 ) then
      d2 = 0.0
    end if
  end if

  d2 = max ( d2, small )

  x1 = x2
  fx = fm

  if ( sf1 < fx ) then
    fx = sf1
    x1 = sx1
  end if
!
!  Update X for linear but not parabolic search.
!
  if ( j /= 0 ) then

    do i = 1, n
      x(i) = x(i) + x1 * v(i,j)
    end do

  end if

  return
end
function praxis ( t0, h0, n, prin, x, f, fmin )
!
!*******************************************************************************
!
!! PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
!
!
!  Discussion:
!
!    PRAXIS returns the minimum of the function F(X,N) of N variables
!    using the principal axis method.  The gradient of the function is
!    not required.
!
!    The approximating quadratic form is
!      q(x') = f(x,n) + (1/2) * (x'-x)-transpose * a * (x'-x)
!    where x is the best estimate of the minimum and 
!      A = inverse(v-transpose) * d * inverse(v)
!   (v(*,*) is the matrix of search directions; d(*) is the array
!   of second differences).  if f has continuous second derivatives
!   near x0, a will tend to the hessian of f at x0 as x approaches x0.
!
!  Modified:
!
!    08 June 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Finding Zeros and Extrema of Functions Without
!    Calculating Derivatives,
!    Stanford University Technical Report STAN-CS-71-198.
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
!    Input, real T0, is a tolerance.  PRAXIS attempts to return 
!    praxis = f(x) such that if x0 is the true local minimum near x, then
!    norm(x-x0) < t0 + squareroot(machep)*norm(x).
!
!    Input, real H0, is the maximum step size.  H0 should be set to about the
!    maximum distance from the initial guess to the minimum.
!    If H0 is set too large or too small, the initial rate of
!    convergence may be slow.
!
!    Input, integer N, is the number of variables upon which
!    the function depends.  N must be at least 2.
!
!    Input, integer PRIN, controls the printing of intermediate results.
!    0, nothing is printed.
!    1, F is printed after every n+1 or n+2 linear minimizations.  
!       final x is printed, but intermediate x is printed only 
!       if n is at most 4.
!    2, the scale factors and the principal values of the approximating 
!       quadratic form are also printed.
!    3, X is also printed after every few linear minimizations.
!    4, the principal vectors of the approximating quadratic form are 
!       also printed.
!
!    Input/output, real X(N), is an array containing on entry a guess 
!    of the point of minimum, on return the estimated point of minimum.
!
!    Input, external F, is the name of the function to be minimized.
!    The function should have the form 
!      FUNCTION F(X,N)
!      INTEGER N
!      REAL F
!      REAL X(N)
!    and accepts X and N as input, returning in F the function value.
!
!    Output, real FMIN, is an estimate of the minimum, used only in printing
!    intermediate results.
!
  integer n
!
  real d(n)
  real df
  real dmin
  real dn
  real dni
  real f
  real f1
  real fmin
  real fx
  real h
  real h0
  integer i
  integer ii
  logical illc
  integer iseed
  integer j
  integer k
  integer k2
  integer kl
  integer kt
  integer ktm
  real large
  real ldfac
  real lds
  real ldt
  real m2
  real m4
  real machep
  integer nl
  integer nf
  real praxis
  integer prin
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real r
  real s
  real scbd
  real sf
  real sl
  real small
  real t
  real t0
  real t2
  real v(n,n)
  real value
  real vlarge
  real vsmall
  real x(n)
  real y(n)
  real z(n)
!
  external f

  data iseed / 1234567 /
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1

  save iseed
!
!  Initialization.
!
  machep = epsilon ( machep )

  small = machep * machep
  vsmall = small * small
  large = 1.0 / small
  vlarge = 1.0 / vsmall
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
!
!     heuristic numbers:
!     if the axes may be badly scaled (which is to be avoided if
!     possible), then set scbd = 10.  otherwise set scbd=1.
!     if the problem is known to be ill-conditioned, set illc = true.
!     otherwise set illc = false.
!     ktm is the number of iterations without improvement before the
!     algorithm terminates.  ktm = 4 is very cautious; usually ktm=1
!     is satisfactory.
!
  scbd = 1.0
  illc = .false.
  ktm = 1

  if ( illc ) then
    ldfac = 0.1
  else
    ldfac = 0.01
  end if

  kt = 0
  nl = 0
  nf = 1
  fx = f(x,n)
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0 * t )
  ldt = h
!
!  The first set of search directions v is the identity matrix.
!
  v(1:n,1:n) = 0.0
  do i = 1, n
    v(i,i) = 1.0
  end do

  d(1) = 0.0
  qd0 = 0.0
  q0(1:n) = x(1:n)
  q1(1:n) = x(1:n)

  if ( prin > 0 ) then
    call print2 ( n, x, prin, fmin )
  end if
!
!  The main loop starts here.
!
40    continue

  sf = d(1)
  d(1) = 0.0
  s = 0.0
!
!  Minimize along the first direction v(*,1).
!  FX must be passed to min by value.
!
  value = fx
  call minny ( n, 1, 2, d(1), s, value, .false., f, x, t, machep, h, v, q0, q1 )

  if ( s <= 0.0 ) then
    do i = 1, n
      v(i,1) = - v(i,1)
    end do
  end if

  if ( sf <= 0.9 * d(1) .or. 0.9 * sf >= d(1) ) then
    d(2:n) = 0.0
  end if
!
!  The inner loop starts here.
!
  do k = 2, n

    y(1:n) = x(1:n)

    sf = fx

    if ( kt > 0 ) then
      illc = .true.
    end if

80  continue

    kl = k
    df = 0.0
!
!  A random step follows (to avoid resolution valleys).
!  PRAXIS assumes that random returns a random number uniformly
!  distributed in (0,1).
!
       if ( illc ) then

         do i = 1, n
           call r_random ( 0.0, 1.0, r )
           s = ( 0.1 * ldt + t2 * 10**kt ) * ( r - 0.5 )
           z(i) = s
           x(1:n) = x(1:n) + s * v(1:n,i)
         end do

         fx = f(x,n)
         nf = nf + 1

       end if
!
!  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
!
       do k2 = k, n

         sl = fx
         s = 0.0
         value = fx
         call minny ( n, k2, 2, d(k2), s, value, .false., f, x, t, machep, &
           h, v, q0, q1 )

         if ( illc ) then
           s = d(k2) * ( ( s + z(k2) )**2 )
         else
           s = sl - fx
         end if

         if ( df <= s ) then
           df = s
           kl = k2
         end if

       end do
!
!  If there was not much improvement on the first try, set
!  ILLC = true and start the inner loop again.
!
       if ( .not. illc ) then
         if ( df < abs ( 100.0 * machep * fx ) ) then
           illc = .true.
           go to 80
         end if
       end if

       if ( k == 2 .and. prin > 1 ) then
         call rvec_print ( n, d, '  The second difference array' )
       end if
!
!  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
!
       do k2 = 1, k-1

         s = 0
         value = fx
         call minny ( n, k2, 2, d(k2), s, value, .false., f, x, t, &
           machep, h, v, q0, q1 )

       end do

       f1 = fx
       fx = sf
       lds = 0

       do i = 1, n
         sl = x(i)
         x(i) = y(i)
         sl = sl - y(i)
         y(i) = sl
         lds = lds + sl**2
       end do

       lds = sqrt ( lds )
!
!  Discard direction V(*,kl).
!
!  If no random step was taken, V(*,KL) is the "non-conjugate"
!  direction along which the greatest improvement was made.
!
       if ( lds > small ) then

         do ii = 1, kl-k
           i = kl - ii
           do j = 1, n
             v(j,i+1) = v(j,i)
           end do
           d(i+1) = d(i)
         end do

         d(k) = 0
         v(1:n,k) = y(1:n) / lds
!
!  Minimize along the new "conjugate" direction V(*,k), which is
!  the normalized vector:  (new x) - (old x).
!
         value = f1
         call minny ( n, k, 4, d(k), lds, value, .true., f, x, t, &
           machep, h, v, q0, q1 )

         if ( lds <= 0.0 ) then
           lds = - lds
           v(1:n,k) = - v(1:n,k)
         end if

       end if

       ldt = ldfac * ldt
       ldt = max ( ldt, lds )

       if ( prin > 0 ) then
         call print2 ( n, x, prin, fmin )
       end if

       t2 = 0.0
       do i = 1, n
         t2 = t2 + x(i)**2
       end do
       t2 = m2 * sqrt(t2) + t
!
!  See whether the length of the step taken since starting the
!  inner loop exceeds half the tolerance.
!
       if ( ldt > 0.5 * t2 ) then
         kt = -1
       end if

       kt = kt + 1

       if ( kt > ktm ) then
         go to 400
       end if

  end do
!
!  The inner loop ends here.
!
!  Try quadratic extrapolation in case we are in a curved valley.
!
171   continue

  call quad ( n, f, x, t, machep, h, v, q0, q1 )

  d(1:n) = 1.0 / sqrt ( d(1:n) )

  dn = 0.0
  do i = 1, n
    dn = max ( dn, d(i) )
  end do

  if ( prin > 3 ) then
    call rmat_print ( v, n, n, n, 'The new direction vectors' )
  end if

  do j = 1, n
    v(1:n,j) = ( d(j) / dn ) * v(1:n,j)
  end do
!
!  Scale the axes to try to reduce the condition number.
!
  if ( scbd > 1.0 ) then

    s = vlarge
    do i = 1, n
      sl = 0.0
      do j = 1, n
        sl = sl + v(i,j)**2
      end do
      z(i) = sqrt ( sl )
      z(i) = max ( z(i), m4 )
      s = min ( s, z(i) )
    end do

    do i = 1, n

      sl = s / z(i)
      z(i) = 1.0 / sl

      if ( z(i) > scbd ) then
        sl = 1.0 / scbd
        z(i) = scbd
      end if

      v(i,1:n) = sl * v(i,1:n)

    end do

  end if
!
!  Calculate a new set of orthogonal directions before repeating
!  the main loop.
!
!  First transpose V for MINFIT:
!
  do i = 2, n
    do j = 1, i-1
      call r_swap ( v(i,j), v(j,i) )
    end do
  end do
!
!  Call MINFIT to find the singular value decomposition of V.
!
!  This gives the principal values and principal directions of the
!  approximating quadratic form without squaring the condition number.
!
  call minfit ( n, n, machep, vsmall, v, d )
!
!  Unscale the axes.
!
  if ( scbd > 1.0 ) then

    do i = 1, n
      v(i,1:n) = z(i) * v(i,1:n)
    end do

    do i = 1, n

      s = 0.0
      do j = 1, n
        s = s + v(j,i)**2
      end do
      s = sqrt ( s )

      d(i) = s * d(i)
      v(1:n,i) = s * v(1:n,i) / s

    end do

  end if

  do i = 1, n

    dni = dn * d(i)

    if ( dni > large ) then
      d(i) = vsmall
    else if ( dni < small ) then
      d(i) = vlarge
    else
      d(i) = 1.0 / dni**2
    end if

  end do
!
!  Sort the eigenvalues and eigenvectors.
!
  call sort ( n, n, d, v )

  dmin = d(n)
  dmin = max ( dmin, small )

  if ( m2 * d(1) > dmin ) then
    illc = .true.
  else
    illc = .false.
  end if

  if ( prin > 1 ) then

    if ( scbd > 1.0 ) then
      call rvec_print ( n, z, '  The scale factors' )
    end if 

    call rvec_print ( n, d, '  Principal values of the quadratic form' )

  end if

  if ( prin > 3 ) then
    call rmat_print ( v, n, n, n, 'The principal axes' )
  end if
!
!  The main loop ends here.
!
  go to 40

400   continue

  if ( prin > 0 ) then
    call rvec_print ( n, x, '  X:' )
  end if

  praxis = fx

  return
end
subroutine print2 ( n, x, prin, fmin )
!
!*******************************************************************************
!
!! PRINT2 prints certain data about the progress of the iteration.
!
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
  integer n
!
  real dmin
  real fmin
  real fx
  integer i
  real ldt
  real ln
  integer nf
  integer nl
  integer prin
  real x(n)
!
  common /global/ fx,ldt,dmin,nf,nl
!
  write ( *, * ) ' '
  write ( *, * ) '  The number of linear searches made is  ', nl
  write ( *, * ) '  The number of function evaluations is  ', nf 
  write ( *, * ) '  The smallest value found is F(X) =     ', fx

  if ( fx <= fmin ) then
    write ( *, 103 ) fmin
  else
    ln = log10 ( fx - fmin )
    write ( *, 102 ) fmin, ln
  end if

  if ( n > 4 .and. prin <= 2 ) then
    return
  end if

  write ( *, 104 ) x(1:n)

  return
102   format (' log (f(x)-',e21.14,') = ',e21.14)
103   format (' log (f(x)-',e21.14,') is undefined.')
104   format (' x is:',e26.14/(e32.14))
end
subroutine quad ( n, f, x, t, machep, h, v, q0, q1 )
!
!*******************************************************************************
!
!! QUAD seeks to minimize the scalar function F along a particular curve.
!
!
!  Discussion:
!
!    The minimizer to be sought is required to lie on a curve defined
!    by Q0, Q1 and X.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
!    Input, integer N, the number of variables.
!
  integer n
!
  real dmin
  real f
  real fx
  real h
  integer i
  real l
  real ldt
  real machep
  integer nf
  integer nl
  real q0(n)
  real q1(n)
  real qa
  real qb
  real qc
  real qd0
  real qd1
  real qf1
  real s
  real t
  real v(n,n)
  real value
  real x(n)
!
  external f
!
  common /global/ fx,ldt,dmin,nf,nl
  common /q/ qa,qb,qc,qd0,qd1,qf1
!
  call r_swap ( fx, qf1 )

  qd1 = 0.0

  do i = 1, n
    s = x(i)
    l = q1(i)
    x(i) = l
    q1(i) = s
    qd1 = qd1 + ( s - l )**2
  end do

  qd1 = sqrt ( qd1 )
  l = qd1
  s = 0.0

  if ( qd0 <= 0.0 .or. qd1 <= 0.0 .or. nl < 3*n*n ) then

    fx = qf1
    qa = 0.0
    qb = qa
    qc = 1.0

  else

    value = qf1

    call minny ( n, 0, 2, s, l, value, .true., f, x, t, &
      machep, h, v, q0, q1 )

    qa = ( l * ( l - qd1 ) ) / ( qd0 * ( qd0 + qd1 ) )
    qb = ( ( l + qd0 ) * ( qd1 - l ) ) / ( qd0 * qd1 )
    qc = ( l * ( l + qd0 ) ) / ( qd1 * ( qd0 + qd1 ) )

  end if

  qd0 = qd1

  do i = 1, n
    s = q0(i)
    q0(i) = x(i)
    x(i) = ( qa * s + qb * x(i) ) + qc * q1(i)
  end do

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  logical, save :: seed = .false.
  real r
  real rhi
  real rlo
  real t
!
  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine rmat_print ( v, lda, m, n, label )
!
!*******************************************************************************
!
!! RMAT_PRINT prints out a matrix.
!
!
!  Discussion:
!
!    The matrix is printed out five columns at a time.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
!    Input, real V(LDA,N), an M by N matrix.
!
!    Input, integer LDA, the leading dimension of the array.
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) LABEL, a label for the matrix.
!
  integer lda
  integer n
!
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) label
  integer m
  real v(lda,n)
!
  write ( *, * ) ' '
  write ( *, '(a)' ) label

  do jlo = 1, n, 5

    jhi = min ( jlo + 4, n )

    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(5g14.6)' ) v(i,jlo:jhi)
    end do

  end do

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector, with an optional title.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  integer n
!
  real a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) title
  end if

  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine sort ( m, n, d, v ) 
!
!*******************************************************************************
!
!! SORT sorts a vector D and adjusts the corresponding columns of a matrix V.
!
!
!  Discussion:
!
!    The routine sorts the elements of D(N) into descending order and moves the
!    corresponding columns of V(N,N).
!    M is the row dimension of V as declared in the calling program.
!
!  Modified:
!
!    15 March 2000
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization with Derivatives,
!    Prentice Hall, 1973.
!
!  Parameters:
!
  integer m
  integer n
!
  real d(n)
  integer i
  integer j
  integer k
  real s
  real v(m,n)
!
  do i = 1, n-1

    k = i
    s = d(i)

    do j = i+1, n
      if ( d(j) > s ) then
        k = j
        s = d(j)
      end if
    end do

    if ( k > i ) then

      d(k) = d(i)
      d(i) = s

      do j = 1, n
        call r_swap ( v(j,i), v(j,k) )
      end do

    end if

  end do

  return
end
