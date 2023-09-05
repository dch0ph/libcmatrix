FROM ubuntu:22.04 AS build

RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update \
  && apt-get install -y libopenblas-dev mpich bash nano g++ make autoconf less

WORKDIR /libcmatrix

COPY src/ ./src/
COPY include/ ./include/


#RUN CXX=g++ CXXFLAGS="-O3" ./configure --with-openblas --with-sse --with-MPI
#RUN make lib/libcmatrix.a

# docker build . -t libcmatrix
