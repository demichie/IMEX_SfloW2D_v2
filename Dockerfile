# MdMV SW_VAR_DENS_MODEL

# start by building the basic container

FROM alpine:edge

MAINTAINER Mattia de' Michieli Vitturi <demichie@gmail.com>

# install packages from alpine edge repositories

ADD repositories /etc/apk/repositories
RUN apk update \
    && apk add --no-cache bash bash-doc bash-completion nano musl-dev m4 zlib-dev git autoconf automake gfortran gdb make curl py3-numpy@community hdf5-dev@community lapack-dev@community openmpi-dev@community

RUN curl -L http://dl-cdn.alpinelinux.org/alpine/edge/testing/x86_64/netcdf-fortran-dev-4.5.3-r0.apk > netcdf-fortran-dev-4.5.3-r0.apk \
    && curl -L http://dl-cdn.alpinelinux.org/alpine/edge/testing/x86_64/netcdf-fortran-4.5.3-r0.apk > netcdf-fortran-4.5.3-r0.apk \
    && curl -L http://dl-cdn.alpinelinux.org/alpine/edge/community/x86_64/netcdf-4.7.4-r1.apk > netcdf-4.7.4-r1.apk \
    && curl -L http://dl-cdn.alpinelinux.org/alpine/edge/community/x86_64/netcdf-dev-4.7.4-r1.apk > netcdf-dev-4.7.4-r1.apk

RUN apk add --allow-untrusted netcdf-4.7.4-r1.apk
RUN apk add --allow-untrusted netcdf-dev-4.7.4-r1.apk
RUN apk add --allow-untrusted netcdf-fortran-4.5.3-r0.apk
RUN apk add --allow-untrusted netcdf-fortran-dev-4.5.3-r0.apk 

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
    && touch README \
    && autoreconf \
    && ./configure 
RUN cd /home/userSW/SW_VAR_DENS_MODEL-master/src \
    && ls -ald \
    && make 
#    && make install -C /home/userSW/SW_VAR_DENS_MODEL-master 
#    && make 
#    && make install \
#    && cd /home/userSW/SW_VAR_DENS_MODEL-master/UTILS \
#    && gfortran -I/usr/include p2d_to_netCDF4.f90 -L/usr/lib -lnetcdf -lnetcdff -o p2d_to_netCDF4.x \
#    && cd /home/userSW/ \
#    && rm *.zip


