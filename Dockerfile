FROM conda/miniconda3

MAINTAINER Ignacio Ferres <iferres@pasteur.edu.uy>

RUN apt update -y && apt upgrade -y && apt install -y \
      wget \
      bzip2 \
      zip \
      git \
      libz-dev \
      bc \
      xvfb \
      software-properties-common && \
    apt clean -y && \
    apt autoremove

RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public |  apt-key add - && \
        add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ && \
        apt install -y apt-transport-https && \
        apt update -y && \
        apt install -y \
          adoptopenjdk-12-hotspot \
          openjfx && \
        apt autoremove -y

RUN wget https://github.com/BenjaminAlbrecht84/DAA_Converter/releases/download/v0.9.0/DAA_Converter_v0.9.0.jar

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/long16S/bin:$PATH
