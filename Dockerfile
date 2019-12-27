# MdMV SW_VAR_DENS test

# start by building the basic container

FROM alpine:edge

MAINTAINER Mattia de' Michieli Vitturi <demichie@gmail.com>

RUN apk add --no-cache bash bash-doc bash-completion
RUN apk add --no-cache musl-dev m4 zlib-dev git
RUN apk add --no-cache gfortran gdb make curl

CMD gfortran --version

RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/testing \
  hdf5-dev

RUN curl -L https://github.com/Unidata/netcdf-c/archive/v4.7.3.tar.gz > v4.7.3.tar.gz
RUN tar -xzvf v4.7.3.tar.gz

WORKDIR netcdf-c-4.7.3
RUN ./configure --enable-remote-fortran-bootstrap --disable-dap --prefix=/usr/local
RUN make install
RUN make build-netcdf-fortran
RUN make install-netcdf-fortran
WORKDIR /
RUN rm v4.7.3.tar.gz
RUN rm -rf netcdf-c-4.7.3

RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/testing \
  openmpi-dev

RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/community \
  lapack-dev

RUN addgroup -S SWgroup && adduser -S userSW -G SWgroup
USER userSW
WORKDIR /home/userSW


# install the code
USER userSW
WORKDIR /home/userSW
USER userSW
RUN curl -LOk https://github.com/demichie/SW_VAR_DENS_MODEL/archive/master.zip
USER userSW
RUN unzip *.zip

USER userSW
RUN cp /home/userSW/SW_VAR_DENS_MODEL-master/TESTS/run_tests.sh .
RUN echo 'cd /home/userSW/SW_VAR_DENS_MODEL-master/TESTS/' | cat - run_tests.sh > temp && mv temp run_tests.sh
RUN chmod +x run_tests.sh


WORKDIR /home/userSW/SW_VAR_DENS_MODEL-master
USER userSW
RUN ./configure
RUN make
RUN make install

WORKDIR /home/userSW/SW_VAR_DENS_MODEL-master/UTILS
RUN gfortran -I/usr/local/include p2d_to_netCDF4.f90 -L/usr/local/lib -lnetcdf -lnetcdff -o p2d_to_netCDF4.x

WORKDIR /home/userSW/
USER userSW

RUN rm *.zip

