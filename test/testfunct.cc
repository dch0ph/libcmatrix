#include <cstddef>
#include <algorithm>
#include <functional>
#include "List.h"
#include "ttyio.h"
#include "timer.h"

#define TYPE size_t

using namespace libcmatrix;
using namespace std;

template<typename T,int incr> void pluseq(T& x) { x+=incr; }

template<typename T> class doespluseq // : public unary_function<T&,T&>
{
  T incr;
public:
  doespluseq(const T& incr_) : incr(incr_) {}
  //  T& operator()(T& in) { return in+=incr; }  
  void operator()(T& in) { in+=incr; }  
};

template<typename T> void dofor(BaseList<T> &a)
{
  for_each(a.begin(),a.end(),&pluseq<T,1>);
}

template<typename T> void dofor_funct(BaseList<T> &a,int incr)
{
  for_each(a.begin(),a.end(),doespluseq<T>(incr) );
}

const size_t N=20;

int main(int argc,const char *argv[])
{
  int count=1;
  const size_t N=getint(argc,argv,count,"N? ",20);
  const int reptimes=getint(argc,argv,count,"Repeat times? ",5000000);


  const List<size_t> omylist(slice(1,N));
  List<size_t> mylist=omylist;
  const List<size_t> blist(slice(10,9+N));
  List<size_t> reslist(N);

  timer<> stopwatch;
  int loop;

  cout << "Original list: " << omylist << endl;

//   mylist=omylist;
//   stopwatch.reset();
//   for (loop=reptimes;loop--;) dofor(mylist);
//   if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
//   cout << "for_each(function): " << mylist << endl;

  mylist=omylist;
  stopwatch.reset();
  for (loop=reptimes;loop--;) dofor_funct(mylist,1);
  if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
  cout << "for_each(functional): " << mylist << endl;

  mylist=omylist;
  stopwatch.reset();
  for (loop=reptimes;loop--;) mylist+=1;
  if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
  cout << "+=: " << mylist << endl;

//   mylist=omylist;
//   stopwatch.reset();
//   for (loop=reptimes;loop--;) mylist-=1;
//   if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
//   cout << "-=: " << mylist << endl;

  reslist= TYPE(0);
  stopwatch.reset();
  for (loop=reptimes;loop--;) add(reslist,mylist,blist);
  if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
  cout << "Functional add(A,B,C): " << reslist << endl;

  const TYPE conval=3;

  reslist= TYPE(0);
  stopwatch.reset();
  for (loop=reptimes;loop--;) apply(reslist,std::plus<TYPE>(),mylist,3);
  if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
  cout << "apply(A,add<>,B,3): " << reslist << endl;

  reslist= TYPE(0);
  stopwatch.reset();
  for (loop=reptimes;loop--;) apply(reslist,bind2nd(std::plus<TYPE>(),conval),mylist);
  if (reptimes!=1) cout << "Time taken: " << stopwatch()*1e6/reptimes << " us\n";
  cout << "apply(A,add<,3>,B): " << reslist << endl;  
  
  return 0;
}

