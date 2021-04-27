# Use phusion/baseimage as base image. To make your builds reproducible, make
# sure you lock down to a specific version, not to `latest`!
# See https://github.com/phusion/baseimage-docker/blob/master/Changelog.md for
# a list of version numbers.
FROM phusion/baseimage:bionic-1.0.0

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# ...put your own build instructions here...
RUN apt-get update
RUN apt-get install -y autoconf automake gfortran nano gdb make curl python3-numpy liblapack-dev libopenmpi-dev libhdf5-dev libnetcdf-dev libnetcdff-dev unzip

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# add user and create group

RUN adduser --disabled-password user_sw
#RUN addgroup -S SWgroup && adduser -S userSW -G SWgroup

# install the code

USER user_sw
WORKDIR /home/user_sw

RUN curl -LOk https://github.com/demichie/SW_VAR_DENS_MODEL/archive/master.zip \
    && unzip *.zip \
    && cp /home/user_sw/SW_VAR_DENS_MODEL-master/TESTS/run_tests.sh . \
    && echo 'cd /home/user_sw/SW_VAR_DENS_MODEL-master/TESTS/' | cat - run_tests.sh > temp \
    && mv temp run_tests.sh \
    && chmod +x run_tests.sh \
    && cd /home/user_sw/SW_VAR_DENS_MODEL-master \
    && touch README \
    && autoreconf \
    && ./configure \
    && make \
    && make install \
    && cd /home/user_sw/SW_VAR_DENS_MODEL-master/UTILS \
    && gfortran -I/usr/include p2d_to_netCDF4.f90 -L/usr/lib -lnetcdf -lnetcdff -o p2d_to_netCDF4.x \
    && cd /home/user_sw/ \
    && rm *.zip
