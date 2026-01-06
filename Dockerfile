FROM gcc:15 AS build
COPY ./CMakeLists.txt /oink/
COPY ./src/ /oink/src/
COPY ./include/ /oink/include/
COPY ./test/ /oink/test/
COPY ./cmake/ /oink/cmake/

WORKDIR /oink/build/
RUN apt-get update; apt-get install -y cmake
RUN apt-get install -y libboost-all-dev --no-install-recommends
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make solve oink

FROM ubuntu:noble AS runtime
RUN apt-get update && \
    apt-get install --no-install-recommends -y libboost-iostreams1.83.0 libboost-random1.83.0 libboost-filesystem1.83.0 libboost-regex1.83.0 && \
    rm -rf /var/lib/apt/lists/*
COPY --from=build /oink/build/oink /usr/local/bin
COPY ./runfiles.sh /usr/local/bin
RUN chmod +x /usr/local/bin/runfiles.sh
ENTRYPOINT ["runfiles.sh"]