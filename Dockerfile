# Use the latest Ubuntu LTS as base
FROM ubuntu:24.04

LABEL maintainer="paul.hodgkinson@durham.ac.uk"
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    cmake \
    git \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Optional: install Minuit2
WORKDIR /opt
RUN git clone https://github.com/GooFit/Minuit2.git Minuit2-master && \
    mkdir Minuit2-master/build && \
    cd Minuit2-master/build && \
    cmake .. && \
    make && make install
#    sed -i '/link_libraries/s/$/\nlink_libraries("-lstdc++")/' ../CMakeLists.txt && \

# Create working directory
WORKDIR /opt/libcmatrix

# Copy libcmatrix source code into container
# (Assumes you will mount or COPY the source code into this directory)
COPY . .

# Configure with OpenBLAS and Minuit
RUN CPPFLAGS=-I/opt/Minuit2-master/inc LDFLAGS=-I/opt/Minuit2-master/build/lib ./configure --with-openblas  -with-minuit && make

# Optional: build test programs
# WORKDIR /opt/libcmatrix/test
# RUN make all || true  # Allow some test programs to fail

# Set default working directory
WORKDIR /opt/libcmatrix

CMD ["/bin/bash"]

# docker build . -t libcmatrixnew
# docker run -it --name libcmatrixnew libcmatrixnew