FROM ubuntu:22.04 AS build

# Leave out mpich for now
RUN DEBIAN_FRONTEND=noninteractive \
  apt-get update \
  && apt-get install -y libopenblas-dev bash nano g++ make autoconf less

WORKDIR /libcmatrix

COPY . .

#RUN CXX=g++ CXXFLAGS="-O3" ./configure --with-openblas --with-sse
#RUN make lib/libcmatrix.a

# docker build . -t libcmatrix
# docker run -it --name libcmatrix libcmatrix