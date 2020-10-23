FROM conda/miniconda3

LABEL authors="Cecilia Salazar <csalazar@pasteur.edu.uy> & Ignacio Ferres <iferres@pasteur.edu.uy>"

# Update base image
RUN apt update -y && apt install -y --no-install-recommends wget procps && apt clean -y && apt autoremove -y

# Copy conda recipe 
COPY environment.yml /

# Update conda and install env
RUN conda env create -f /environment.yml && \
    conda clean -a -y

# Export temporary PATH to allow MEGAN installation (needs openjdk)
ARG PATH=/usr/local/envs/porefile/bin:$PATH

# Install MEGAN6
ARG MEGAN=MEGAN_Community_unix_6_20_5.sh
RUN wget -P /opt https://software-ab.informatik.uni-tuebingen.de/download/megan6/"$MEGAN" && \
    bash /opt/"$MEGAN" -q && rm /opt/"$MEGAN" && echo "-Xmx8G" > /opt/megan/MEGAN.vmoptions

# Export PATH
ENV PATH=/usr/local/envs/porefile/bin:/opt/megan/tools:$PATH

