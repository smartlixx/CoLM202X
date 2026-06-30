FROM ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran libopenmpi-dev openmpi-bin \
    libnetcdf-dev libnetcdff-dev liblapack-dev libblas-dev \
    make perl csh git vim \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/*
WORKDIR /workspace
COPY . .
RUN [ -f .github/workflows/create_defineh.bash ] && \
    bash .github/workflows/create_defineh.bash SinglePoint LULC_IGBP URBANOFF Campbell CaMaOFF BGCOFF CROPOFF \
    || true
RUN cd /workspace && make clean 2>/dev/null || true && \
    make ARCH=gnu 2>&1 || (echo "Fallback..." && make 2>&1) || true
RUN useradd -m -s /bin/bash colm && chown -R colm:colm /workspace
USER colm
WORKDIR /workspace
ENTRYPOINT ["/bin/bash"]
