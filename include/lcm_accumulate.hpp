#ifndef LCM_ACCUMULATE_
#define LCM_ACCUMULATE_

namespace {

  template<class Mout,class Min> void accumulate_(Mout& dest,const Min& a, Mout& tmp)
  {
    multiply(tmp,a,dest);
    dest.swap(tmp);
  }
  
  template<class Mout,class Min> inline void accumulate_(Mout& dest,const Min& a,bool& haveflag)
  {
    if (haveflag)
      dest&=a;
    else {
      haveflag=true;
      dest=a;
    }
  }

  template<class Mout,class Min> void accumulate_(Mout& dest,const Min& a, Mout& tmp, bool& haveflag)
  {
    if (haveflag) {
      multiply(tmp,a,dest);
      dest.swap(tmp);
    }
    else {
      haveflag=true;
      dest=a;
    }
  }
  
  template<> void accumulate_(Matrix<complex>& dest,const List<complex>& a, Matrix<complex>& tmp, bool& haveflag)
  {
    if (haveflag) {
      multiply(tmp,a,dest);
      dest.swap(tmp);
    }
    else {
      haveflag=true;
      full(dest,a);
    }
  }

  void accumulate_(BlockedMatrix<complex>& Uacc, const Matrix<complex>& U, Matrix<complex>* Utmp =NULL)
{
  if (!Uacc) {
    Uacc.duplicate_structure(U);
    Uacc.front()=U;
  }
  else {
    assert(Uacc.size()==1);
    Matrix<complex>& Uaccd(Uacc.front());
    if (Utmp) {
      multiply(*Utmp,U,Uaccd);
      Uaccd=*Utmp;
    }
    else
      Uaccd&=U;
  }
}

  void accumulate_(BlockedMatrix<complex>& Uacc, const BaseList<complex>& U) 
  {
    if (!Uacc) {
      Uacc.create(ExplicitList<1,size_t>(U.size()));
      full(Uacc.front(),U);
  }
    else
      Uacc.front()&=U;
  }

  void accumulate_(BlockedMatrix<complex>& Uacc, const BlockedMatrix<complex>& U, int which, BlockedMatrix<complex>& Utmp)//, const matrix_partition_set& part)
  {
    multiply(Utmp,U,Uacc);
    Uacc.swap(Utmp);
  }

  void accumulate_(BlockedMatrix<complex>& Uacc, const BlockedMatrix<complex>& U, int which, BlockedMatrix<complex>* Utmp =NULL)//, const matrix_partition_set& part)
{
  if (which>=0) {
    const size_t uwhich(which);
    if (Utmp) {
      //! rather dodgy - only checks overall size matches
      if (Utmp->raw_storage().size()!=U.raw_storage().size())
	Utmp->duplicate_structure(U);
      accumulate_(Uacc,U(uwhich),&((*Utmp)(uwhich)));//,partp);
    }
    else
      accumulate_(Uacc,U(uwhich));
  }
  else {
    if (!Uacc || (Utmp==NULL))
      Uacc&=U;
    else
      accumulate_(Uacc,U,*Utmp);
  }
}

void accumulate_(BlockedMatrix<complex>& Uacc, const ListList<complex>& U, int which)
{
  if (which>=0)
    accumulate_(Uacc,U(size_t(which)));
  else
    Uacc&=U;
}

}

#endif
