# hash:sha256:dc535bca3157e9d3b0000d952a8445152fd3747e9c1e325a19108b9f6aa58a81
FROM registry.codeocean.com/codeocean/matlab:2020b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        meshlab=2020.03+dfsg1-1 \
    && rm -rf /var/lib/apt/lists/*
