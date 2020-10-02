FROM conda/miniconda3

LABEL authors="Cecilia Salazar <csalazar@pasteur.edu.uy> & Ignacio Ferres <iferres@pasteur.edu.uy>"

RUN apt update -y && apt upgrade -y && apt install -y \
      wget \
      bzip2 \
      zip \
      git \
      libz-dev \
      bc \
      xvfb \
      gnupg \
      software-properties-common && \
    apt clean -y && \
    apt autoremove -y

RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public |  apt-key add - && \
        add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ && \
        apt install -y apt-transport-https && \
        apt update -y && \
        apt install -y \
          adoptopenjdk-12-hotspot \
          openjfx && \
        apt autoremove -y

RUN wget -P /opt https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_18_0.sh && \
	bash /opt/MEGAN_Community_unix_6_18_0.sh -q  && \
	MEGAN="/opt/megan"  && \
	find "$MEGAN"/tools -type f | while read -r file; do b=$(basename "$file"); ln -fs "$file" /usr/local/bin/"$b"; sed -i -e "s@\"\$0\"@\"$file\"@" "$file"; done  && \
	ln -fs "$MEGAN"/MEGAN /usr/local/bin/MEGAN && \
	sed -i -e "s@\"\$0\"@\"$MEGAN/MEGAN\"@" "$MEGAN"/MEGAN  && \
	rm /opt/MEGAN_Community_unix_6_18_0.sh

# RUN wget -P /opt https://github.com/BenjaminAlbrecht84/DAA_Converter/releases/download/v0.9.0/DAA_Converter_v0.9.0.jar

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH $PATH:/usr/local/envs/porefile/bin
