#!/usr/bin/env sh

podman build -t porefile .

# Login to GHCR
# echo $CR_PAT | podman login ghcr.io -u USER --password-stdin

# Push to GHCR
# podman push porefile ghcr.io/microgenlab/porefile:latest
# podman push porefile ghcr.io/microgenlab/porefile:VERSION