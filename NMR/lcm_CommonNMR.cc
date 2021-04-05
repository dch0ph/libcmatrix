
template<class T> void gammareduce_(T& dest, const BaseList<T>& source, size_t s, size_t n)
{
  const size_t N(source.size());
  if (n>=N)
    throw BadIndex("gammareduce");
  if (n==0)
    dest=source(s-1);
  else
    unitary_isimtrans(dest,source(s-1),source(n-1));
}

template<class T> void gammareduce_(BaseList<T> dest, const BaseList<T>& source, size_t n)
{
  const size_t N(source.size());
  const size_t nobs(dest.size());
  if (N % nobs)
    throw InvalidParameter("gammareduce: observations must be divisor of number of gamma steps");
  const size_t step(N/nobs);
  size_t s=step;
  for (size_t k=0;k<nobs;k++) {
    gammareduce_(dest(k),source,s,n);
    s+=step;
  }
}


