FROM docker.io/mambaorg/micromamba:1.3.1

LABEL authors="Cecilia Salazar <csalazar@pasteur.edu.uy> & Ignacio Ferres <iferres@pasteur.edu.uy>"
LABEL org.opencontainers.image.source https://github.com/microgenlab/porefile

# Copy conda recipe 
COPY environment.yml /

# Update conda, create env, and clean 
RUN micromamba create -n porefile -y -f /environment.yml && \
		micromamba clean --all --yes

# Export  PATH to allow MEGAN installation (needs openjdk and wget)
ENV PATH=/opt/conda/envs/porefile/bin:$PATH


