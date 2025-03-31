FROM python:3.10-slim

# 기본 패키지 설치
RUN apt-get update && apt-get install -y \
    git curl wget build-essential \ 
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 파이썬 패키지 설치
RUN pip install --upgrade pip \
    && pip install ipython \
    && pip install jupyter

# 작업 디렉토리 설정
WORKDIR /app

# requirements 미리 복사하고 설치
COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt