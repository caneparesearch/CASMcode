# syntax=docker/dockerfile:1
FROM dengzeyu/casm.0.3.x:latest

RUN apt-get update && \
    apt-get install -y vim make gcc autoconf m4 bash-completion libtool pkg-config automake  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Install CASM dependencies
RUN  conda activate casm_dev && conda remove casm casm-python casm-cpp --force
#RUN git clone -b hengning_Fvib https://github.com/caneparesearch/CASMcode.git /app/CASMcode

#WORKDIR /app/CASMcode
ENV CASM_DIRTY_IS_OK=1
#RUN conda activate casm_dev && bash build.sh
#RUN make install
#RUN conda activate casm_dev && pip install python/casm && bash clean.sh && conda clean -a -y
WORKDIR /
