FROM debian:stretch-slim
LABEL maintainer="dirk.winkelhardt@comkom.de"

ARG UID=1000
ARG GID=1000
ARG CRUX_ZIP_URL=https://noble.gs.washington.edu/crux-downloads/crux-3.2/crux-3.2.Linux.x86_64.zip
ARG CRUX_BINARY_IN_ZIP=crux-3.2.Linux.x86_64/bin/crux
ENV UID=$UID
ENV GID=$GID
ENV CRUX_ZIP_URL=$CRUX_ZIP_URL
ENV CRUX_BINARY_IN_ZIP=$CRUX_BINARY_IN_ZIP

# mzML-splitup settings
ENV SPECTRUM_SUFFIX=""

# identification settings
ENV MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTID=2
ENV NUMBER_OF_DECOYS=1000
ENV LOWER_MASS_TOLERANCE_IN_PPM=5
ENV UPPER_MASS_TOLERANCE_IN_PPM=5
ENV THREAD_COUNT=2

WORKDIR /home/max-decoy
COPY . ./src/

# prepare system
RUN apt-get update -y && apt-get install -y gcc curl unzip \
    && groupadd -r -g $GID max-decoy \
	&& useradd -r -M -g max-decoy -u $UID max-decoy \
    && chown -R max-decoy:max-decoy .

USER max-decoy
WORKDIR /home/max-decoy

# prepare user
RUN mkdir /home/max-decoy/results \
    && curl https://sh.rustup.rs -sSf | sh -s -- -y --default-toolchain 1.32.0 --no-modify-path \
    && cd src/ \
    && /home/max-decoy/.cargo/bin/cargo build --release \
    && cp target/release/max_decoy /home/max-decoy \
    && cp .env ~ \
    && cp run_splitup_and_identification.sh ~ \
    && cd /home/max-decoy \
    && rm -r /home/max-decoy/src \
    && curl -o crux.zip $CRUX_ZIP_URL \
    && unzip -j crux.zip $CRUX_BINARY_IN_ZIP \
    && rm crux.zip

CMD [./run_splitup_and_identification.sh]
