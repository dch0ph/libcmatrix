
template<class Op1,class Op2> class unary_composition : public STLSPACE::unary_function<typename Op2::argument_type,typename Op1::result_type>
{
  Op1 op1_;
  Op2 op2_;
 public:
  unary_composition(const Op1& x_, const Op2& y_)
    : op1_(x_), op2_(y_) {}
  typename Op1::result_type operator()(const typename Op2::argument_type& b) {
    return op1_(op2_(b)); }
};

template<class Op1,class Op2> inline unary_composition<Op1,Op2> compose(const Op1& op1, const Op2& op2) {
  return unary_composition<Op1,Op2>(op1,op2); }

template<class Op,class Iter> class unary_iterator_composition : 
public STLSPACE::iterator< 
  typename STLSPACE::iterator_traits<Iter>::iterator_category,typename Op::result_type, typename STLSPACE::iterator_traits<Iter>::difference_type,
   const typename Op::result_type*, const typename Op::result_type&  > 
{
  Op op_;
  Iter iter_;
  typedef typename Op::result_type result_type;
 public:
  unary_iterator_composition(const Op& opv, const Iter& iterv)
    : op_(opv), iter_(iterv) {}
  /*comparison only on iterators to saves time and allows iterators based on different operations to be combined
    (but violates rule: iter1==iter2 => *iter1 == *iter2 ...) */
  template<class Op2> bool operator== (const unary_iterator_composition<Op2,Iter>& v) const {
    return (iter_==v.iter_); }

  template<class Op2> bool operator!= (const unary_iterator_composition<Op2,Iter>& v) const {
    return (iter_!=v.iter_); }

  template<class Op2> bool operator> (const unary_iterator_composition<Op2,Iter>& x) const {
    return (iter_>x.iter_); }

  template<class Op2> bool operator< (const unary_iterator_composition<Op2,Iter>& x) const {
    return (iter_<x.iter_); }

  template<class Op2> bool operator>= (const unary_iterator_composition<Op2,Iter>& x) const {
    return (iter_>=x.iter_); }

  template<class Op2> bool operator<= (const unary_iterator_composition<Op2,Iter>& x) const {
    return (iter_<=x.iter_); }

  unary_iterator_composition& operator++() {
    ++iter_;
    return *this; }

  unary_iterator_composition& operator--() {
    --iter_;
    return *this; }

  unary_iterator_composition operator++(int) {
    const Iter itertmp(iter_);
    ++iter_;
    return unary_iterator_composition<Op,Iter>(op_,itertmp); }

  unary_iterator_composition operator--(int) {
    const Iter itertmp(iter_);
    --iter_;
    return unary_iterator_composition<Op,Iter>(op_,itertmp); }

  unary_iterator_composition& operator+=(typename unary_iterator_composition::difference_type i) {
    iter_+=i;
    return *this; }

  unary_iterator_composition& operator-=(typename unary_iterator_composition::difference_type i) {
    iter_-=i;
    return *this; }

  unary_iterator_composition operator+ (typename unary_iterator_composition::difference_type i) const {
    return unary_iterator_composition<Op,Iter>(op_,iter_+i); }

  unary_iterator_composition operator- (typename unary_iterator_composition::difference_type i) const {
    return unary_iterator_composition<Op,Iter>(op_,iter_-i); }

  template<class Op2> typename unary_iterator_composition::difference_type operator- (const unary_iterator_composition<Op2,Iter>& x) const {
    return iter_-x.iter_; }

  result_type operator[] (typename unary_iterator_composition::difference_type i) const {
    return op_(iter_[i]); }

  result_type operator*() const {
    return op_(*iter_); }
};

template<class Op,class Iter> inline unary_iterator_composition<Op,Iter> compose_unary_iterator(const Op& opv,const Iter& iterv) {
  return unary_iterator_composition<Op,Iter>(opv,iterv); }

/* Creates a cyclic iterator from an arbitrary container
   Resulting iterator is bidirectional (unless container iterator is only forward) */

template<class T> class cyclic_iterator : public STLSPACE::iterator< 
  typename LCM_iterator_cat_< STLSPACE::iterator_traits<typename T::iterator> >::iterator_category,
  LCM_VAL(T) > {
 public:
  cyclic_iterator(T& obj)
    : begin_(obj.begin()),
    cur_(obj.begin()),
    end_(obj.end()) {}

  cyclic_iterator(T& obj,const typename T::iterator & curv) :
    begin_(obj.begin()),
    cur_(curv),
    end_(obj.end()) {}

  typename T::reference operator*() const {
    return *cur_; }

  cyclic_iterator& operator++() {
    if (++cur_==end_)
      cur_=begin_;
    return *this; }

  cyclic_iterator operator++(int) {
    cyclic_iterator tmp=*this;
    return ++tmp; }

  cyclic_iterator& operator--() {
    if (cur_==begin_)
      cur_=end_;
    --cur;
    return *this; }

  cyclic_iterator operator--(int) {
    cyclic_iterator tmp=*this;
    return --tmp; }

 private:
  const typename T::iterator begin_,end_;
  typename T::iterator cur_;
};
     
template<class T> class cyclic_const_iterator : public STLSPACE::iterator< 
  typename LCM_iterator_cat_< STLSPACE::iterator_traits<typename T::const_iterator> >::iterator_category,
  LCM_VAL(T), ptrdiff_t, const LCM_VAL(T)*,const LCM_VAL(T)& > {
 public:
  cyclic_const_iterator(const T& obj) :
    begin_(obj.begin()),
    cur_(obj.begin()),
    end_(obj.end()) {}

  cyclic_const_iterator(const T& obj,const typename T::const_iterator & curv) :
    begin_(obj.begin()),
    cur_(curv),
    end_(obj.end()) {}

  typename T::const_reference operator*() const {
    return *cur_; }

  cyclic_const_iterator& operator++() {
    if (++cur_==end_)
      cur_=begin_;
    return *this; }

  cyclic_const_iterator operator++(int) {
    cyclic_const_iterator tmp=*this;
    return ++tmp; }

  cyclic_const_iterator& operator--() {
    if (cur_==begin_)
      cur_=end_;
    --cur;
    return *this; }

  cyclic_const_iterator operator--(int) {
    cyclic_const_iterator tmp=*this;
    return --tmp; }

 private:
  const typename T::const_iterator begin_,end_;
  typename T::const_iterator cur_;
};
     
template<typename T> struct doesnothing : public STLSPACE::unary_function<T,T> {
  T operator()(const T& a) const {
    return a; }
};

