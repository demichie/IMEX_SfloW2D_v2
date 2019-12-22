# MdMV SW_VAR_DENS test

# start by building the basic container

FROM alpine:edge

MAINTAINER Mattia de' Michieli Vitturi <demichie@gmail.com>

RUN apk add --no-cache bash bash-doc bash-completion
RUN apk add --no-cache musl-dev 
RUN apk add --no-cache gfortran gdb make curl

CMD gfortran --version

RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/testing \
  openmpi-dev

RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/community \
  lapack-dev

RUN addgroup -S SWgroup && adduser -S userSW -G SWgroup
USER userSW

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

WORKDIR /home/userSW/
USER userSW

RUN rm *.zip

