#include "cmatrix_utils.h"
#include <iostream>

int main()
{
  const int NSTEPS=2;

  const int NS[NSTEPS]={7,8};

  set_seed();

  for (int i=0;i<NSTEPS;i++) {
    const int N=NS[i];
    std::cout << N << " point FTs\n";
    
    List<complex> start(N);

    for (int j=N;j--;) start(j)=complex(random(1.0),random(1.0));
    
    std::cout << "Initial data: " << start << '\n' << std::endl;
    
    const bool usefft=ispowerof2(N);

    std::cout << "Use FFT: " << usefft << '\n' << std::endl;

    List<complex> fftfacs(N),bftfacs(N);
    List<complex> fted,ffted,ftback,fftback;

    ft_create_table(fftfacs,FT_FORWARD);
    ft_create_table(bftfacs,FT_BACKWARD);

    fted=ft(start,FT_FORWARD);
    std::cout << "FTed: " << fted << '\n' << std::endl;

    ft(fted,start,fftfacs,1.0);
    std::cout << "Table FTed: " << fted << '\n' << std::endl;

    if (usefft) {
      ffted=start; fft_ip(ffted,FT_FORWARD);
      std::cout << "FFTed: " << ffted << '\n' << std::endl;
    }

    ftback=ft(fted,FT_BACKWARD); ftback/=double(N);
    std::cout << "orig - FT-1(FT(orig)): " << (start-ftback) << '\n' << std::endl;

    ft(ftback,fted,bftfacs,1.0); ftback/=double(N);
    std::cout << "orig - TFT-1(TFT(orig)): " << (start-ftback) << '\n' << std::endl;

    if (usefft) {
      fftback=fft(ffted,FT_BACKWARD); fftback/=double(N);
      std::cout << "orig - FFT-1(FFT(orig)): " << (start-fftback) << '\n' << std::endl;
    }
  }
  return 0;
}



    
