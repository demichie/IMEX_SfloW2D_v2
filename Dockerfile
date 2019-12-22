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
WORKDIR /home/userSW
RUN curl -LOk https://github.com/demichie/SW_VAR_DENS_MODEL/archive/master.zip
WORKDIR /home/userSW/
RUN ls -la

RUN unzip *.zip

WORKDIR /home/userSW/SW_VAR_DENS_MODEL-master
RUN ./configure
RUN make
RUN make install


WORKDIR /home/userSW/
RUN cp /home/userSW/SW_VAR_DENS_MODEL-master/run_tests.sh .
RUN chmod +x run_tests.sh

RUN rm *.zip

