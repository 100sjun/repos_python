FROM pytorch/pytorch:2.6.0-cuda12.6-cudnn9-runtime

WORKDIR /app

# 필수 패키지 설치
RUN apt-get clean && rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y --fix-missing curl gcc g++ build-essential tar

# uv 설치
RUN curl -L --insecure https://github.com/astral-sh/uv/releases/download/0.6.9/uv-x86_64-unknown-linux-musl.tar.gz -o uv.tar.gz && \
    tar -xzf uv.tar.gz && \
    mv uv-x86_64-unknown-linux-musl/uv /usr/local/bin/uv && \
    chmod +x /usr/local/bin/uv && \
    rm -rf uv.tar.gz uv-x86_64-unknown-linux-musl

RUN uv --version && uv venv /opt/uv-venv
RUN rm -rf /opt/conda

ENV PATH="/opt/uv-venv/bin:$PATH"
ENV PYTHONPATH="/opt/uv-venv/lib/python3.11/site-packages"
