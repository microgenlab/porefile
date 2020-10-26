FROM conda/miniconda3

LABEL authors="Cecilia Salazar <csalazar@pasteur.edu.uy> & Ignacio Ferres <iferres@pasteur.edu.uy>"

# Copy conda recipe 
COPY environment.yml /

# Update conda and install env
RUN conda env create -f /environment.yml && \
    conda clean -afy

# Export  PATH to allow MEGAN installation (needs openjdk and wget)
ENV PATH=/usr/local/envs/porefile/bin:/opt/megan/tools:$PATH

# Install MEGAN6
ARG MEGAN=MEGAN_Community_unix_6_20_5.sh
RUN wget -P /opt https://software-ab.informatik.uni-tuebingen.de/download/megan6/$MEGAN && \
    bash /opt/$MEGAN -q && rm /opt/$MEGAN && echo "-Xmx8G" > /opt/megan/MEGAN.vmoptions

