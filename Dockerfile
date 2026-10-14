FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN test "$(dpkg --print-architecture)" = amd64 \
    || { echo "Build with --platform linux/amd64: dipcall's release binaries require amd64." >&2; exit 1; }

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential bzip2 ca-certificates curl git zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN git clone --depth 1 --branch 0.25.0 https://github.com/chhylp123/hifiasm.git \
    && make -C hifiasm -j2 \
    && install hifiasm/hifiasm /usr/local/bin/hifiasm \
    && git clone --depth 1 --branch v0.1 https://github.com/lh3/yak.git \
    && make -C yak -j2 \
    && install yak/yak /usr/local/bin/yak

RUN curl -fsSL --retry 3 \
        https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2 \
        -o dipcall.tar.bz2 \
    && mkdir dipcall \
    && tar -xjf dipcall.tar.bz2 -C dipcall --strip-components=1

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    LC_ALL=C \
    LANG=C \
    PYTHONIOENCODING=utf8

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash bcftools ca-certificates gawk make minimap2 parallel perl pigz \
        python3 python3-pip samtools tabix \
    && python3 -m pip install --no-cache-dir svim-asm==1.0.3 vcfpy==0.13.8 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/hifiasm /usr/local/bin/hifiasm
COPY --from=builder /usr/local/bin/yak /usr/local/bin/yak
COPY --from=builder /build/dipcall/ /root/dipcall/
COPY bin/parse_snvs.py /usr/local/bin/parse_snvs.py

RUN chmod 755 /root /usr/local/bin/parse_snvs.py \
    && ln -s /root/dipcall/run-dipcall /usr/local/bin/run-dipcall

WORKDIR /data
CMD ["/bin/bash"]
