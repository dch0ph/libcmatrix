template<typename T> void kronecker_(Matrix<T>& d, const Matrix<T>& a, int n)
{
  if (n<1)
    throw InvalidParameter("kronecker");

  if (!a) {
    d.identity(n);
    return;
  }

  if (n==1) {
    d=a;
    return;
  }

  const size_t ra=a.rows();
  const size_t ca=a.cols();

  d.create(ra*n,ca*n);
  d=T(0);

  for (size_t i=ra;i--;) {
    const size_t br=i*n;
    for (size_t j=ca;j--;) {
      const T& val=a(i,j);
      if (val!=0.0) {
	const size_t bc=j*n;
	for (size_t k=n;k--;)
	  d(br+k,bc+k)=val;
      }
    }
  }
}

template<typename T> inline void kronecker_(Matrix<T>& d, int n, const Matrix<T>& a, bool transpose =false)
{
  if (n<1)
    throw InvalidParameter("kronecker");

  if (!a) {
    d.identity(n);
    return;
  }

  if (n==1) {
    d=a;
    return;
  }

  size_t ca=a.cols();
  size_t ra=a.rows();
  if (transpose)
    std::swap(ca,ra);

  d.create(n*ra,n*ca);
  d=T(0);

  for (size_t k=n;k--;) {
    const size_t br=k*ra;
    const size_t bc=k*ca;

    for (size_t i=ra;i--;) {
      if (transpose) {
	for (size_t j=ca;j--;)
	  d(br+i,bc+j)=a(j,i);
      }
      else {
	for (size_t j=ca;j--;)
	  d(br+i,bc+j)=a(i,j);
      }
    }
  }
}


template<typename T> inline void kronecker_(Matrix<T>& d, const Matrix<T>& a, const Matrix<T>& b, bool transposeB =false)
{
  if (!a) {
    if (transposeB)
      transpose(d,b);
    else
      d=b;
    return;
  }
  if (!b) {
    d=a;
    return;
  }

  size_t rb=b.rows();
  size_t cb=b.cols();
  if (transposeB)
    std::swap(rb,cb);

  d.create(rb*a.rows(),cb*a.cols());
  d=T(0);

  for (size_t i=a.rows();i--;) {
    const size_t bigr=rb*i;
    for (size_t j=a.cols();j--;) {
      const T& val=a(i,j);
      if (val!=0.0) {
	const size_t bigc=cb*j;
	for (size_t k=rb;k--;) {
	  if (transposeB) {
	    for (size_t l=cb;l--;)
	      d(bigr+k,bigc+l)=val*b(l,k); //!< NB rb and cb have been swapped
	  }
	  else {
	    for (size_t l=cb;l--;)
	      d(bigr+k,bigc+l)=val*b(k,l);
	  }
	}
      }
    }
  }
}
