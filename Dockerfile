FROM ubuntu:xenial

# Set the working directory to /app
WORKDIR /app

# Install build system and Boost library
RUN apt-get update \
    && apt-get install --yes wget build-essential gcc-multilib libboost-all-dev

# Install GSL
RUN wget -O gsl.tgz ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz \
    && tar -zxf gsl.tgz \
    && mkdir gsl \
    && cd gsl-1.16 \
    && ./configure --prefix=/app/gsl \
    && make \
    && make install
ENV LIBRARY_PATH /app/gsl/lib/
ENV CPLUS_INCLUDE_PATH /app/gsl/include/

# Install SimKernel
RUN apt-get install --yes unzip \
    && wget -O simkernel.zip http://github.com/ChristophKirst/SimKernel/archive/master.zip \
    && unzip simkernel.zip \
    && cd SimKernel-master \
    && make \
    && make install
ENV LIBRARY_PATH /app/gsl/lib/:/app/SimKernel-master/lib/

# Install yaml-cpp
RUN apt-get install --yes cmake \
    && wget -O yaml-cpp.zip https://github.com/jbeder/yaml-cpp/archive/release-0.5.3.zip \
    && unzip yaml-cpp.zip \
    && cd yaml-cpp-release-0.5.3 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make
ENV LIBRARY_PATH /app/gsl/lib/:/app/SimKernel-master/lib/:/app/yaml-cpp-release-0.5.3/build/
ENV LD_LIBRARY_PATH /app/gsl/lib/:/app/SimKernel-master/lib/:/app/yaml-cpp-release-0.5.3/build/
ENV CPLUS_INCLUDE_PATH /app/gsl/include/:/app/gsl/lib/:/app/yaml-cpp-release-0.5.3/include/

# Install Ruby and Rake
RUN apt-get install --yes ruby \
    && gem install rake

# Add the whole repo to the image
COPY . /app

# Build te-causality binaries
RUN cd transferentropy-sim \
    && rake te-extended test

# By default, only execute tests and exit
CMD /app/transferentropy-sim/test
