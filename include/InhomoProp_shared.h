namespace libcmatrix {

  void unitary_isimtrans_identity(cmatrix&, const RCmatrix&, const RCmatrix&);

  template<class T> void unitary_isimtrans_In(cmatrix& dest,const RCmatrix& R,const Matrix<T>& a,const RCmatrix& C,Matrix<T>* tmp =NULL);

  template<class T> void unitary_isimtrans_allownull(Matrix<T>& dest,const rmatrix& R,const Matrix<T>& a,const rmatrix& C,Matrix<T>* tmp =NULL);

 void unitary_isimtrans_In(cmatrix&, const BaseList<double>&, const RCmatrix&);

template<class F> void InhomoHelper_::process(F& obj,rmatrix& A,const cmatrix& sigma0T)
{
  enorm(A,sigma0T);
  if (obj.isdiagonal())
    squeeze(A);
  obj.reset();
}

template<class T> void unitary_isimtrans_In(List<complex>& dest,const RCmatrix& R,const BaseList<T>& a,const RCmatrix& C)
{
  dest.create(a.size());
  switch (R.type()) {
  case RCmatrix::NONE:
    switch (C.type()) {
    case RCmatrix::NONE:
      dest=a;
      return;
    case RCmatrix::REAL:
      multiply(dest,a,C.get_real());
      return;
    case RCmatrix::COMPLEX:
      multiply(dest,a,C.get_complex());
      return;
    default: break;
    }
    break;
  case RCmatrix::REAL:
    if (C.type()==RCmatrix::NONE) {
      transpose_multiply(dest,R.get_real(),a);
      return;
    }
    break;
  case RCmatrix::COMPLEX:
    if (C.type()==RCmatrix::NONE) {
      conj_transpose_multiply(dest,R.get_complex(),a);
      return;
    }
    break;
  }
  throw Failed("unitary_isimtrans: transformation undefined");
 }

template<class F> void InhomoHelper_::process(F& obj,cmatrix& A,const cmatrix& sigma0T,const cmatrix& detectT)
{
  transpose(A,detectT);
  emultiply_ip(A,sigma0T);
  if (obj.isdiagonal())
    squeeze(A);
  obj.reset();
}

template<class F> void InhomoHelper_::process(F& obj,cmatrix& A,const rmatrix& sigma0T, const rmatrix& detectT)
{
  rmatrix tmpr;
  transpose(tmpr,detectT);
  emultiply_ip(tmpr,sigma0T);
  if (obj.isdiagonal())
    squeeze(tmpr);
  A=tmpr;
  obj.reset();
}

struct enormcomplex_ : public ::std::unary_function<complex,complex> {
  complex operator() (const complex& a) { return complex(norm(a)); }
};

struct enormreal_ : public ::std::unary_function<complex,double> {
  complex operator() (double a) { return complex(a*a); }
};

struct enormcomplexip_ {
  void operator() (complex& d) { d=norm(d); }
};

template<class F> void InhomoHelper_::process(F& obj,cmatrix& AU,const cmatrix& sigma0T)
{
  apply(AU,enormcomplex_(),sigma0T);
  if (obj.isdiagonal())
    squeeze(AU);
  obj.reset();
}

template<class F> void InhomoHelper_::process(F& obj,cmatrix& AU,const rmatrix& sigma0T)
{
  apply(AU,enormreal_(),sigma0T);
  if (obj.isdiagonal())
    squeeze(AU);
  obj.reset();
}

template<class F> void InhomoHelper_::process(F& obj,rmatrix& A)
{
  emultiply_ip(A,A);
  if (obj.isdiagonal())
    squeeze(A);
  obj.reset();
}

template<class F> void InhomoHelper_::process(F& obj,cmatrix& A)
{
  apply_ip(enormcomplexip_(),A);
  if (obj.isdiagonal())
    squeeze(A);
  obj.reset();
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A,const cmatrix& sigma0,const cmatrix& det)
{
  if (!aretranspose(sigma0,det)) 
    throw Failed(LCM_INCOMPAT);
  cmatrix sigma0T,detectT;
  unitary_isimtrans_In(sigma0T,obj.row(),sigma0,obj.col(),&detectT);
  unitary_isimtrans_In(detectT,obj.col(),det,obj.row(),&A);
  process(obj,A,sigma0T,detectT);
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A,const rmatrix& sigma0,const rmatrix& det)
{
  if (obj.arebothreal()) {
    if (!aretranspose(sigma0,det))
      throw Failed(LCM_INCOMPAT);
    rmatrix sigma0T,detectT,tmp;
    unitary_isimtrans_allownull(sigma0T,obj.row().get_real(),sigma0,obj.col().get_real(),&tmp);
    unitary_isimtrans_allownull(detectT,obj.col().get_real(),det,obj.row().get_real(),&tmp);
    process(obj,A,sigma0T,detectT);
  }
  else
    observe(obj,A,cmatrix(sigma0),cmatrix(det));
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A, const BaseList<double>& sigma0, const BaseList<double>& det)
{
  if (!arematching(sigma0,det)) 
    throw Failed(LCM_INCOMPAT);
  if (!obj.isdiagonal()) 
    throw Failed("Diagonal operators only valid for diagonal blocks");
  switch (obj.row().type()) {
  case RCmatrix::REAL: {
      rmatrix sigma0T,detectT;
      unitary_isimtrans(sigma0T,sigma0,obj.row().get_real());
      unitary_isimtrans(detectT,det,obj.row().get_real());
      process(obj,A,sigma0T,detectT);
    }
    break;
  case RCmatrix::COMPLEX: {
      cmatrix sigma0T;
      cmatrix& detectT=A;
      unitary_isimtrans(sigma0T,sigma0,obj.row().get_complex());
      unitary_isimtrans(detectT,det,obj.row().get_complex());
      process(obj,A,sigma0T,detectT);
    }
    break;
  default:
    throw Failed("Case not handled");
  }
}

template<class T> void InhomoHelper_::observe(T& obj, rmatrix& A,const cmatrix& sigma0det)
{
  cmatrix sigma0T;
  unitary_isimtrans_In(sigma0T,obj.row(),sigma0det,obj.col());
  process(obj,A,sigma0T);
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A,const cmatrix& sigma0det)
{
  unitary_isimtrans_In(A,obj.row(),sigma0det,obj.col());
  process(obj,A,A);
}

template<class T> void InhomoHelper_::observe(T& obj, rmatrix& A,const rmatrix& sigma0det)
{
  if (obj.arebothreal()) {
    unitary_isimtrans_allownull(A,obj.row().get_real(),sigma0det,obj.col().get_real());
    process(obj,A);
  }
  else {
    cmatrix sigma0T;
    unitary_isimtrans_In(sigma0T,obj.row(),sigma0det,obj.col(),&A);
    process(obj,A,sigma0T);
  }
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A,const rmatrix& sigma0det)
{
  if (obj.arebothreal()) {
    rmatrix sigma0T;
    unitary_isimtrans_allownull(sigma0T,obj.row().get_real(),sigma0det,obj.col().get_real());
    process(obj,A,sigma0T);
  }
  else {
    unitary_isimtrans_In(A,obj.row(),sigma0det,obj.col());
    process(obj,A);
  }
}

template<class T> void InhomoHelper_::observe(T& obj, rmatrix& A,const BaseList<double>& sigma0det)
{
  if (!obj.isdiagonal())
    throw Failed("Diagonal operators only valid for diagonal blocks");
  switch (obj.row().type()) {
  case RCmatrix::REAL: 
    unitary_isimtrans(A,sigma0det,obj.row().get_real());
    process(obj,A);
    break;
  case RCmatrix::COMPLEX: {
      cmatrix sigma0T;
      unitary_isimtrans(sigma0T,sigma0det,obj.row().get_complex());
      process(obj,A,sigma0T);
    }
    break;
  default:
    throw Failed("Case not handled");
  }
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A,const BaseList<double>& sigma0det)
{
  if (!obj.isdiagonal()) 
    throw Failed("Diagonal operators only valid for diagonal blocks");
  unitary_isimtrans_In(A,sigma0det,obj.row());
  process(obj,A,A);
}

template<class T> void InhomoHelper_::observe(T& obj, rmatrix& A)
{
  if (obj.isdiagonal()) 
    throw Failed("observe() only valid for off-diagonal blocks");
  if (obj.arebothreal()) {
    unitary_isimtrans_identity(A,obj.row().get_real(),obj.col().get_real());
    process(obj,A);
  }
  else {
    cmatrix sigma0T;
    unitary_isimtrans_identity(sigma0T,obj.row(),obj.col());
    process(obj,A,sigma0T);
  }
}

template<class T> void InhomoHelper_::observe(T& obj, cmatrix& A)
{
  if (obj.isdiagonal()) 
    throw Failed("observe() only valid for off-diagonal blocks");
  if (obj.arebothreal()) {
    rmatrix AR;
    unitary_isimtrans_identity(AR,obj.row().get_real(),obj.col().get_real());
    process(obj,A,AR);
  }
  else {
    unitary_isimtrans_identity(A,obj.row(),obj.col());
    process(obj,A);
  }
}

} //namespace libcmatrix
