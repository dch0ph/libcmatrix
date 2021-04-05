#ifndef LCM_UnionHolder_h_
#define LCM_UnionHolder_h_

namespace libcmatrix {

  struct Missing2_ {};
  struct Missing3_ {};
  struct Missing4_ {};
  struct Missing5_ {};
  struct Missing6_ {};
  struct Missing7_ {};
  struct Missing8_ {};
  struct Missing9_ {};

  template<class T,class T2> struct unionholder_subcount_ {
    enum { value=1 };
  };
  template<class T> struct unionholder_subcount_<T,T> {
    enum { value=0 };
  };
    
  template<class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_count_ {
    enum { types=1
      +unionholder_subcount_<T2,Missing2_>::value
    +unionholder_subcount_<T3,Missing3_>::value
    +unionholder_subcount_<T4,Missing4_>::value
    +unionholder_subcount_<T5,Missing5_>::value
    +unionholder_subcount_<T6,Missing6_>::value
    +unionholder_subcount_<T7,Missing7_>::value
    +unionholder_subcount_<T8,Missing8_>::value
	   +unionholder_subcount_<T9,Missing9_>::value };
  };

  template<class T,class MissingT> struct union_holder_apply_ {
    template<class F> static typename F::result_type
    apply(const F& func, const T& a) { return func(a); }

    template<class F> static void apply_ip(const F& func, T& a, const T& b) { func(a,b); }
    template<class F> static void apply_ip(const F& func, T& a) { func(a); }
  };

  template<class MissingT> struct union_holder_apply_<MissingT,MissingT> {
    template<class F> static typename F::result_type
    apply(const F&, const MissingT&) { throw InternalError("union_holder_apply_"); }

    template<class F> static void apply_ip(const F&, MissingT&, const MissingT&) { throw InternalError("union_holder_apply_"); }
    template<class F> static void apply_ip(const F&, MissingT&) { throw InternalError("union_holder_apply_"); }
  };

  template<class T,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>  struct unionholder_translate_;

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T1,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=1;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T2,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=2;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T3,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=3;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T4,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=4;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T5,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=5;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T6,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=6;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T7,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=7;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T8,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=8;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_translate_<T9,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    static const size_t index=9;
  };

  template<size_t N,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_;
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<1,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T1 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<2,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T2 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<3,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T3 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<4,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T4 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<5,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T5 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<6,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T6 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<7,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T7 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<8,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T8 type;
  };
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> struct unionholder_reverse_<9,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
    typedef T9 type;
  };
  
  template<class T1,class T2 =Missing2_,class T3 =Missing3_,class T4 =Missing4_,class T5 =Missing5_,class T6 =Missing6_,class T7 =Missing7_,class T8 =Missing8_,class T9 =Missing9_> class CommonHolder {
  public:
    CommonHolder() {}
    explicit CommonHolder(const T1& v) : first_(v) {}
    explicit CommonHolder(const T2& v) : second_(v) {}
    explicit CommonHolder(const T3& v) : third_(v) {}
    explicit CommonHolder(const T4& v) : fourth_(v) {}
    explicit CommonHolder(const T5& v) : fifth_(v) {}
    explicit CommonHolder(const T6& v) : sixth_(v) {}
    explicit CommonHolder(const T7& v) : seventh_(v) {}
    explicit CommonHolder(const T8& v) : eighth_(v) {}
    explicit CommonHolder(const T9& v) : nineth_(v) {}

    static size_t items() { return unionholder_count_<T2,T3,T4,T5,T6,T7,T8,T9>::types; }


    void clear() {
      clear_(first_);
      clear_(second_);
      clear_(third_);
      clear_(fourth_);
      clear_(fifth_);
      clear_(sixth_);
      clear_(seventh_);
      clear_(eighth_);
      clear_(nineth_);
    }

    template<class T> CommonHolder& operator= (const T& v) {
      get(Type2Type<T>())=v;
      return *this;
    }
    
    template<class T> T& operator()(Type2Type<T> in) {
      return get(in);
    }

    template<class T> const T& operator()(Type2Type<T> in) const {
      return get(in);
    }

   template<class T> struct translate : public unionholder_translate_<T,T1,T2,T3,T4,T5,T6,T7,T8,T9> {};
    template<size_t M> struct reverse : public unionholder_reverse_<M,T1,T2,T3,T4,T5,T6,T7,T8,T9> {};

    template<int M> typename reverse<M>::type& operator()(Int2Type<M>) {
      return get(Type2Type<typename reverse<M>::type>());
    }

    template<int M> const typename reverse<M>::type& operator()(Int2Type<M>) const {
      return get(Type2Type<typename reverse<M>::type>());
    }

    T1& get(Type2Type<T1>) { return first_; }
    T2& get(Type2Type<T2>) { return second_; }
    T3& get(Type2Type<T3>) { return third_; }
    T4& get(Type2Type<T4>) { return fourth_; }
    T5& get(Type2Type<T5>) { return fifth_; }
    T6& get(Type2Type<T6>) { return sixth_; }
    T7& get(Type2Type<T7>) { return seventh_; }
    T8& get(Type2Type<T8>) { return eighth_; }
    T9& get(Type2Type<T9>) { return nineth_; }

    const T1& get(Type2Type<T1>) const { return first_; }
    const T2& get(Type2Type<T2>) const { return second_; }
    const T3& get(Type2Type<T3>) const { return third_; }
    const T4& get(Type2Type<T4>) const { return fourth_; }
    const T5& get(Type2Type<T5>) const { return fifth_; }
    const T6& get(Type2Type<T6>) const { return sixth_; }
    const T7& get(Type2Type<T7>) const { return seventh_; }
    const T8& get(Type2Type<T8>) const { return eighth_; }
    const T9& get(Type2Type<T9>) const { return nineth_; }

  protected:
    T1 first_;
    T2 second_;
    T3 third_;
    T4 fourth_;
    T5 fifth_;
    T6 sixth_;
    T7 seventh_;
    T8 eighth_;
    T9 nineth_;

    template<class T> static void clear_(T& v) {
      clear__(v,Int2Type< type_traits<T>::dimensionality >());
    }

  private:
    template<class T,int M> static void clear__(T& v,Int2Type<M>) { v.clear(); }
    template<class T> static void clear__(T&, Int2Type<0>) {}
    
  };

  template<size_t N,class T1,class T2 =Missing2_,class T3 =Missing3_,class T4 =Missing4_,class T5 =Missing5_,class T6 =Missing6_,class T7 =Missing7_,class T8 =Missing8_,class T9 =Missing9_> class UnionHolder 
    : public CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>
  {

    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::clear_;

   size_t type_;

    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::first_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::second_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::third_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::fourth_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::fifth_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::sixth_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::seventh_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::eighth_;
    using CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::nineth_;

 public:

   UnionHolder() : type_(0) {}
    explicit UnionHolder(const T1& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(1) {}
    explicit UnionHolder(const T2& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(2) {}
    explicit UnionHolder(const T3& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(3) {}
    explicit UnionHolder(const T4& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(4) {}
    explicit UnionHolder(const T5& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(5) {}
    explicit UnionHolder(const T6& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(6) {}
    explicit UnionHolder(const T7& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(7) {}
    explicit UnionHolder(const T8& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(8) {}
    explicit UnionHolder(const T9& v) : CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>(v), type_(9) {}

    template<class T> struct translate : public unionholder_translate_<T,T1,T2,T3,T4,T5,T6,T7,T8,T9> {};
    template<size_t M> struct reverse : public unionholder_reverse_<M,T1,T2,T3,T4,T5,T6,T7,T8,T9> {};

    template<int M> explicit UnionHolder(Int2Type<M>) : type_(M) {
      LCM_STATIC_CHECK( (M<=unionholder_count_<T2,T3,T4,T5,T6,T7,T8,T9>::types) , UnionHolder_BadIndex );
    }
    template<class T> explicit UnionHolder(Type2Type<T>) : type_(translate<T>::index) {}
      
    bool operator! () const { return (type_==0); }

    template<typename T> bool istype(Type2Type<T>) const {
      return (type_==translate<T>::index);
    }

    template<typename F> typename F::result_type apply(const F& func) const {
      switch (type_) {
      case 1:
	return func(CommonHolder<T1,T2,T3,T4,T5,T6,T7,T8,T9>::first_);
      case 2:
	return union_holder_apply_<T2,Missing2_>::apply(func,second_);
      case 3:
	return union_holder_apply_<T3,Missing3_>::apply(func,third_);
      case 4:
	return union_holder_apply_<T4,Missing4_>::apply(func,fourth_);
      case 5:
	return union_holder_apply_<T5,Missing5_>::apply(func,fifth_);
      case 6:
	return union_holder_apply_<T6,Missing6_>::apply(func,sixth_);
      case 7:
	return union_holder_apply_<T7,Missing7_>::apply(func,seventh_);
      case 8:
	return union_holder_apply_<T8,Missing8_>::apply(func,eighth_);
      case 9:
	return union_holder_apply_<T9,Missing9_>::apply(func,nineth_);
      }
      throw Undefined("UnionHolder::apply");	
    }

    template<typename F> void apply_ip(const F& func, const UnionHolder<N,T1,T2,T3,T4,T5,T6,T7,T8,T9>& b) {
      if (type_!=b.type())
	throw Mismatch("UnionHolder::apply_ip arguments have storage types");
      switch (type_) {
      case 1:
	func(first_,b(Int2Type<1>()));
	return;
      case 2:
	union_holder_apply_<T2,Missing2_>::apply_ip(func,second_,b(Int2Type<2>()));
	return;
      case 3:
	union_holder_apply_<T3,Missing3_>::apply_ip(func,third_,b(Int2Type<3>()));
	return;
      case 4:
	union_holder_apply_<T4,Missing4_>::apply_ip(func,fourth_,b(Int2Type<4>()));
	return;
      case 5:
	union_holder_apply_<T5,Missing5_>::apply_ip(func,fifth_,b(Int2Type<5>()));
	return;
      case 6:
	union_holder_apply_<T6,Missing6_>::apply_ip(func,sixth_,b(Int2Type<6>()));
	return;
      case 7:
	union_holder_apply_<T7,Missing7_>::apply_ip(func,seventh_,b(Int2Type<7>()));
	return;
      case 8:
	union_holder_apply_<T8,Missing8_>::apply_ip(func,eighth_,b(Int2Type<8>()));
	return;
      case 9:
	union_holder_apply_<T9,Missing9_>::apply_ip(func,nineth_,b(Int2Type<9>()));
	return;
      }
      throw Undefined("UnionHolder::apply_ip");	
    }

    template<typename F> void apply_ip(const F& func) {
      switch (type_) {
      case 1:
	func(first_);
	return;
      case 2:
	union_holder_apply_<T2,Missing2_>::apply_ip(func,second_);
	return;
      case 3:
	union_holder_apply_<T3,Missing3_>::apply_ip(func,third_);
	return;
      case 4:
	union_holder_apply_<T4,Missing4_>::apply_ip(func,fourth_);
	return;
      case 5:
	union_holder_apply_<T5,Missing5_>::apply_ip(func,fifth_);
	return;
      case 6:
	union_holder_apply_<T6,Missing6_>::apply_ip(func,sixth_);
	return;
      case 7:
	union_holder_apply_<T7,Missing7_>::apply_ip(func,seventh_);
	return;
      case 8:
	union_holder_apply_<T8,Missing8_>::apply_ip(func,eighth_);
	return;
      case 9:
	union_holder_apply_<T9,Missing9_>::apply_ip(func,nineth_);
	return;
      }
      throw Undefined("UnionHolder::apply_ip");	
    }

   size_t type() const { return type_; }
   void type(size_t typev) {
     if (type_!=typev) {
       if (typev>N)
	 throw InvalidParameter("UnionHolder::type");
       clear();
       type_=typev;
     }
   }
    template<int M> void type(Int2Type<M>) {
      if (type_!=M) {
	LCM_STATIC_CHECK( M<=N, UnionHolder_BadIndex );
	clear();
	type_=M;
      }
    }
    template<class T> void type(Type2Type<T>) {
      if (type_!=translate<T>::index) {
	clear();
	type_=translate<T>::index;
      }
    }

    //! overrides ::clear method from ::CommonHolder
    void clear() {
      switch (type_) {
      case 1: clear_(first_); break;
      case 2: clear_(second_); break;
      case 3: clear_(third_); break;
      case 4: clear_(fourth_); break;
      case 5: clear_(fifth_); break;
      case 6: clear_(sixth_); break;
      case 7: clear_(seventh_); break;
      case 8: clear_(eighth_); break;
      case 9: clear_(nineth_); break;
      }
      type_=0;
    }

    template<class T> UnionHolder& operator= (const T& v) {
      this->get(Type2Type<T>())=v;
      type(Type2Type<T>()); //only change type if assignment was successful
      return *this;
    }

    UnionHolder& operator= (const UnionHolder& v) {
      if (type_!=v.type_) {
	clear();
	type_=v.type_;
      }
      switch (type_) {
      case 1: first_=v.first_; break;
      case 2: second_=v.second_; break;
      case 3: third_=v.third_; break;
      case 4: fourth_=v.fourth_; break;
      case 5: fifth_=v.fifth_; break;
      case 6: sixth_=v.sixth_; break;
      case 7: seventh_=v.seventh_; break;
      case 8: eighth_=v.eighth_; break;
      case 9: nineth_=v.nineth_; break;
      }
      return *this;
    }

    template<class T> T& set(Type2Type<T> in) {
      type(in);
      return this->get(in);
    }

    template<int M> typename reverse<M>::type& set(Int2Type<M> in) {
      type(in);
      return this->get(Type2Type<typename reverse<M>::type>());
    }

    template<class T> T& operator()(Type2Type<T> in) {
      verify(in);
      return this->get(in);
    }

    template<class T> const T& operator()(Type2Type<T> in) const {
      verify(in);
      return this->get(in);
    }

    template<int M> void verify(Int2Type<M>) const {
      if (type_!=M)
	throw Failed("UnionHolder: wrong type requested");
    }

    template<class T> void verify(Type2Type<T>) const {
      if (type_!=translate<T>::index)
	throw Failed("UnionHolder: wrong type requested");
    }

    template<int M> typename reverse<M>::type& operator()(Int2Type<M> in) {
      verify(in);
      return this->get(Type2Type<typename reverse<M>::type>());
    }

    template<int M> const typename reverse<M>::type& operator()(Int2Type<M> in) const {
      verify(in);
      return this->get(Type2Type<typename reverse<M>::type>());
    }
  };
  
  namespace {
    struct union_holder_print_ {
      union_holder_print_(std::ostream& ostrv) : ostr_(ostrv) {}
      typedef std::ostream& result_type;
      template<class T> std::ostream& operator()(const T& a) const { return ostr_ << a; }
      std::ostream& ostr_;
    };
  }
    
  template<size_t N,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> std::ostream& operator<< (std::ostream& ostr,const UnionHolder<N,T1,T2,T3,T4,T5,T6,T7,T8,T9>& a)
   {
     if (!a)
       return ostr << "<unset>\n";
     else
       return a.apply(union_holder_print_(ostr));
   }

//      switch (a.type()) {
//      case 1: a.print(ostr,a.get(Type2Type<T1>())); break;
//      case 2: a.print(ostr,a.get(Type2Type<T2>())); break;
//      case 3: a.print(ostr,a.get(Type2Type<T3>())); break;
//      case 4: a.print(ostr,a.get(Type2Type<T4>())); break;
//      case 5: a.print(ostr,a.get(Type2Type<T5>())); break;
//      case 6: a.print(ostr,a.get(Type2Type<T6>())); break;
//      case 7: a.print(ostr,a.get(Type2Type<T7>())); break;
//      case 8: a.print(ostr,a.get(Type2Type<T8>())); break;
//      case 9: a.print(ostr,a.get(Type2Type<T9>())); break;
//      default:
//        return ostr << "<unset>\n";
//      }
//      return ostr;
//    }

  template<class R,class C> struct RealComplexHolder : public UnionHolder<2,R,C> {
    //   using UnionHolder<2,R,C>::get;
    //using CommonHolder<R,C>::get;
    //using UnionHolder<2,R,C>::set;

    enum { NONE=0, REAL=1, COMPLEX=2 };
    
    RealComplexHolder() : UnionHolder<2,R,C>() {}
    RealComplexHolder(const C& Cv) : UnionHolder<2,R,C>(Cv) {}
    RealComplexHolder(const R& Rv) : UnionHolder<2,R,C>(Rv) {}

    bool operator!() const { return (this->type()==NONE); }
    bool iscomplex() const { return (this->type()==COMPLEX); }
    bool isreal() const { return (this->type()==REAL); }
    
    C& set_complex() { return this->set(Int2Type<COMPLEX>()); }
    R& set_real() { return this->set(Int2Type<REAL>()); }
    const C& get_complex() const { return this->get(Type2Type<C>()); }
    const R& get_real() const { return this->get(Type2Type<R>()); }
    C& get_complex() { return this->get(Type2Type<C>()); }
    R& get_real() { return this->get(Type2Type<R>()); }
  };

} //namespace libcmatrix

#endif
