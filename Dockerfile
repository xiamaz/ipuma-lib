FROM graphcore/poplar:2.4.0

COPY ./build  /build

WORKDIR "/"
CMD ["/build/bin/ipusw", "--help"]