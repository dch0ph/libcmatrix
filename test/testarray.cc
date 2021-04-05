#include "spin_system.h"

int main()
{
  const int nsys=2;
  
  spin_system fred(2);
  spin_system *axs=new spin_system[nsys] (2);
  
  for (int i=0;i<nsys;i++)
    std::cout << axs[i] << std::endl;
  
  delete[] axs;
  return 0;
}
