# RTC-Tools with the Base Modelica and FMU routes: see README.md in this folder.
#
#   docker build --pull -t rtc-basemodelica "https://github.com/jgillis/rtc-tools.git#rtc-basemodelica:docker/rtc-basemodelica"
#
# Takes a while on a first build: OpenModelica and CasADi are compiled from source.

# ---------------------------------------------------------------------------------------------
# Toolchain for OpenModelica master without its GUI clients.
# ---------------------------------------------------------------------------------------------
FROM ubuntu:24.04 AS toolchain

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates curl git python3 cmake pkg-config build-essential clang gfortran default-jdk \
      omniidl automake autoconf libtool flex bison zip unzip lsb-release ccache \
      libhwloc-dev uuid-dev libncurses-dev libsqlite3-dev libreadline-dev liblpsolve55-dev \
      libomniorb4-dev libhdf5-dev liblapack-dev libblas-dev libcurl4-gnutls-dev libexpat1-dev \
      libomp-dev libboost-all-dev libffi-dev libminizip-dev zlib1g-dev

# Without a UTF-8 locale Java (ANTLR 3) reads the sources as ASCII and fails on license headers
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# SimulationRuntime/rust needs rustc >= 1.85, more than Ubuntu's cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/rust/cargo PATH=/opt/rust/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --profile minimal --default-toolchain stable --no-modify-path \
 && cargo --version

# ---------------------------------------------------------------------------------------------
# OpenModelica: lsdae-input-derivatives (fmi-ls-dae export, initial residuals, analytic
# derivatives, --fmiInputDerivatives) with the C+Rust --daeMode fix.
# ---------------------------------------------------------------------------------------------
FROM toolchain AS openmodelica

ARG OPENMODELICA_COMMIT=2620460e5b19fea69d176408a48ab4c2685c4a6a

# origin stays upstream: the submodule URLs in .gitmodules are relative to it, and git describe
# (the version of omc) needs its tags.
RUN git init -q /src/OpenModelica && cd /src/OpenModelica \
 && git remote add origin https://github.com/OpenModelica/OpenModelica.git \
 && git remote add fork https://github.com/jgillis/OpenModelica.git \
 && git fetch -q --filter=blob:none --tags origin master \
 && git fetch -q --filter=blob:none fork rtc-basemodelica \
 && git -c advice.detachedHead=false checkout -q ${OPENMODELICA_COMMIT} \
 && git submodule update --init --recursive --filter=blob:none \
      OMCompiler/3rdParty OMCompiler/Compiler/boot/bomc OMCompiler/Examples/ModelicaTutorial

RUN --mount=type=cache,target=/ccache,id=rtc-basemodelica-ccache \
    export CCACHE_DIR=/ccache \
 && cd /src/OpenModelica \
 && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/openmodelica \
      -DOM_ENABLE_GUI_CLIENTS=OFF -DOM_ENABLE_OMSIMULATOR=OFF -DOM_ENABLE_TESTSUITE=OFF \
 && cmake --build build --parallel "$(nproc)" --target install \
 && /opt/openmodelica/bin/omc --version

# ---------------------------------------------------------------------------------------------
# CasADi: fmu3-configuration-mode (Configuration Mode switch, flat residuals, initial residuals)
# with the solvers that the RTC-Tools tests and examples use.
# ---------------------------------------------------------------------------------------------
FROM toolchain AS casadi

ARG CASADI_COMMIT=e5f5df13018085a0fadebc4e532ee88ddcdb71ac

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      python3-dev python3-numpy libpcre2-dev libpcre3-dev yodl byacc \
      coinor-libipopt-dev coinor-libbonmin-dev coinor-libcbc-dev coinor-libclp-dev

# The SWIG fork that CasADi is built with
RUN git clone -q https://github.com/jaeandersson/swig.git /src/swig && cd /src/swig \
 && git -c advice.detachedHead=false checkout -q 2cd61e3435448c525db0331b0a46c193d14c871a \
 && ./autogen.sh && ./configure && make -j"$(nproc)" && make install

RUN git init -q /src/casadi && cd /src/casadi \
 && git remote add origin https://github.com/jgillis/casadi.git \
 && git fetch -q --depth 1 origin ${CASADI_COMMIT} \
 && git -c advice.detachedHead=false checkout -q FETCH_HEAD \
 && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/opt/casadi_lsdae -DPYTHON_PREFIX=/opt/casadi_lsdae \
      -DWITH_PYTHON=ON -DWITH_PYTHON3=ON -DSWIG_EXECUTABLE=/usr/local/bin/swig \
      -DWITH_SELFCONTAINED=ON -DWITH_DEEPBIND=ON \
      -DWITH_GHC_FILESYSTEM=ON -DWITH_BUILD_GHC_FILESYSTEM=ON \
      -DWITH_LIBZIP=ON -DWITH_BUILD_LIBZIP=ON -DWITH_ZLIB=ON -DWITH_BUILD_ZLIB=ON \
      -DWITH_IPOPT=ON -DWITH_BONMIN=ON -DWITH_CBC=ON -DWITH_CLP=ON \
      -DWITH_HIGHS=ON -DWITH_BUILD_HIGHS=ON \
 && cmake --build build --parallel "$(nproc)" --target install

# ---------------------------------------------------------------------------------------------
# The image: omc, CasADi, and editable checkouts of the Python packages in /src
# ---------------------------------------------------------------------------------------------
FROM ubuntu:24.04 AS image

ARG PYMOCA_COMMIT=48d085b3959bc7dd53940709e9f4390c0d8d0e67
ARG CHANNEL_FLOW_COMMIT=656fdd3f64128db0a8dfa1d8f46990a6e94bf3ff
ARG FMI_LS_DAE_COMMIT=f2b7ab7da0afd860a6d2bf9528ad96165da4a66b

# gcc, gfortran and cargo: buildModelFMU compiles the FMU
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates curl git python3 python3-venv python3-dev build-essential gfortran zip unzip \
      liblapack-dev libblas-dev libcurl4-gnutls-dev libexpat1-dev libffi-dev libminizip-dev \
      zlib1g-dev uuid-dev libhdf5-dev \
      coinor-libipopt-dev coinor-libbonmin-dev coinor-libcbc-dev coinor-libclp-dev

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/rust/cargo
ENV PATH=/opt/venv/bin:/opt/openmodelica/bin:/opt/rust/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PYTHONPATH=/opt/casadi_lsdae

COPY --from=openmodelica /opt/openmodelica /opt/openmodelica
COPY --from=casadi /opt/casadi_lsdae /opt/casadi_lsdae
COPY --from=toolchain /opt/rust /opt/rust

RUN python3 -m venv /opt/venv

# The tip of the branch this Dockerfile comes from: a new commit there invalidates the clone
ADD https://api.github.com/repos/jgillis/rtc-tools/git/refs/heads/rtc-basemodelica /tmp/rtc-tools.ref

RUN mkdir -p /src \
 && git clone -q -b rtc-basemodelica https://github.com/jgillis/rtc-tools.git /src/rtc-tools \
 && git clone -q -b rtc-basemodelica https://github.com/jgillis/rtc-tools-channel-flow.git /src/rtc-tools-channel-flow \
 && git -C /src/rtc-tools-channel-flow -c advice.detachedHead=false checkout -q ${CHANNEL_FLOW_COMMIT} \
 && git clone -q -b rtc-basemodelica https://github.com/jgillis/pymoca.git /src/pymoca \
 && git -C /src/pymoca -c advice.detachedHead=false checkout -q ${PYMOCA_COMMIT} \
 && git clone -q -b rtc-basemodelica https://github.com/jgillis/fmi-ls-dae.git /src/fmi-ls-dae \
 && git -C /src/fmi-ls-dae -c advice.detachedHead=false checkout -q ${FMI_LS_DAE_COMMIT}

# --no-deps for our three: rtc-tools declares pymoca == 0.9.* and casadi from PyPI, and these
# checkouts deliberately replace both. The rest of the dependencies are as declared.
# SETUPTOOLS_SCM_PRETEND_VERSION: the forks do not carry the release tags.
RUN pip install -q --no-cache-dir -U pip \
 && pip install -q --no-cache-dir "numpy<2.6" scipy setuptools "antlr4-python3-runtime==4.13.*" \
      "defusedxml>=0.7.0" "netCDF4>=1.6.1" "rtc-tools-standard-library>=0.1.2" \
      "highspy>=1.13.1,<1.15" "pytest>=8.2" \
 && pip install -q --no-cache-dir --no-deps -e /src/pymoca \
 && SETUPTOOLS_SCM_PRETEND_VERSION=1.4.0 pip install -q --no-cache-dir --no-deps -e /src/rtc-tools-channel-flow \
 && SETUPTOOLS_SCM_PRETEND_VERSION=2.8.1 pip install -q --no-cache-dir --no-deps -e /src/rtc-tools

COPY export_fmus.sh export_basemodelica.sh run_tests.sh sitecustomize.py README.md /opt/rtc-basemodelica/

# The exports are build artifacts next to the models, as the drivers and tests expect them.
RUN omc /src/rtc-tools/tests/modelica/install_msl.mos \
 && bash /opt/rtc-basemodelica/export_basemodelica.sh \
 && bash /opt/rtc-basemodelica/export_fmus.sh

WORKDIR /src/rtc-tools
CMD ["bash"]
