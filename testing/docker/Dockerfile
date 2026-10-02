# Image for running the Linux CI gate locally (testing/ci_docker.sh).
# Mirrors what ubuntu-latest provides to .github/workflows/ci.yml: GCC,
# CMake and git (the build embeds the git revision in git_version.h).
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential cmake git \
 && rm -rf /var/lib/apt/lists/*
# The repository is bind-mounted here by ci_docker.sh; it is owned by a
# different uid than root, which git refuses to touch unless told otherwise.
RUN git config --global --add safe.directory /work
WORKDIR /work
