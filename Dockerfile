FROM ubuntu:18.04

RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        clang \
        libgmp3-dev \
        libssl-dev \
        libntl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY . /usr/src/secure-gwas

WORKDIR /usr/src/secure-gwas/code

ENTRYPOINT ["ls"]
