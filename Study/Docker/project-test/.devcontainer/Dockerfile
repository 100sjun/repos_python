# CUDA 12.6.3 + cuDNN + Ubuntu 22.04 (Runtime only)
FROM nvidia/cuda:12.6.3-cudnn-runtime-ubuntu22.04

# 기본 패키지 설치
RUN apt-get update && apt-get install -y \
    git \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    curl \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Python 심볼릭 링크 생성 (python3 -> python)
RUN ln -sf /usr/bin/python3 /usr/bin/python && \
    ln -sf /usr/bin/pip3 /usr/bin/pip

# 작업 디렉토리 생성
WORKDIR /app

# Python 가상환경 생성 (/opt/venv에 생성)
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# pip 업그레이드
RUN . /opt/venv/bin/activate && \
    pip install --upgrade pip setuptools wheel

# PyTorch GPU 설치 (CUDA 12.6)
RUN . /opt/venv/bin/activate && \
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

# 추가 Python 패키지 설치
RUN . /opt/venv/bin/activate && \
    pip install numpy pandas matplotlib scikit-learn jupyter black flake8 pylint isort

# 환경변수 설정을 .bashrc에 추가하여 터미널에서도 사용 가능하게 함
RUN echo 'export PATH="/opt/venv/bin:$PATH"' >> /root/.bashrc && \
    echo 'source /opt/venv/bin/activate' >> /root/.bashrc

# app 디렉토리에 심볼릭 링크 생성 (선택사항)
RUN ln -sf /opt/venv /app/.venv

# 기본 실행 셸
CMD ["bash"]