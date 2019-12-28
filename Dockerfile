# MdMV SW_VAR_DENS_MODEL

# start by building the basic container

FROM alpine:edge

MAINTAINER Mattia de' Michieli Vitturi <demichie@gmail.com>

# install packages from alpine edge repositories

ADD repositories /etc/apk/repositories
RUN apk add --no-cache bash bash-doc bash-completion nano musl-dev m4 zlib-dev git gfortran gdb make curl py3-numpy@community hdf5-dev@testing openmpi-dev@testing lapack-dev@community


# install netcdf libreries

RUN curl -L https://github.com/Unidata/netcdf-c/archive/v4.7.3.tar.gz > v4.7.3.tar.gz \
    && tar -xzvf v4.7.3.tar.gz \
    && cd netcdf-c-4.7.3 \
    && ./configure --enable-remote-fortran-bootstrap --disable-dap --prefix=/usr/local \
    && make install \
    && make build-netcdf-fortran \
    && make install-netcdf-fortran \
    && cd .. \
    && rm v4.7.3.tar.gz \
    && rm -rf netcdf-c-4.7.3


# add user and create group

RUN addgroup -S SWgroup && adduser -S userSW -G SWgroup

# install the code

USER userSW
WORKDIR /home/userSW

RUN curl -LOk https://github.com/demichie/SW_VAR_DENS_MODEL/archive/master.zip \
    && unzip *.zip \
    && cp /home/userSW/SW_VAR_DENS_MODEL-master/TESTS/run_tests.sh . \
    && echo 'cd /home/userSW/SW_VAR_DENS_MODEL-master/TESTS/' | cat - run_tests.sh > temp \
    && mv temp run_tests.sh \
    && chmod +x run_tests.sh \
    && cd /home/userSW/SW_VAR_DENS_MODEL-master \
    && ./configure 
RUN cd /home/userSW/SW_VAR_DENS_MODEL-master \
    && make \
    && make install \
    && cd /home/userSW/SW_VAR_DENS_MODEL-master/UTILS \
    && gfortran -I/usr/local/include p2d_to_netCDF4.f90 -L/usr/local/lib -lnetcdf -lnetcdff -o p2d_to_netCDF4.x \
    && cd /home/userSW/ \
    && rm *.zip

