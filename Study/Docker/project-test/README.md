# PyTorch GPU 개발 환경

CUDA 12.6 기반 PyTorch 개발 환경입니다. DevContainer를 사용하여 쉽게 환경을 구성할 수 있습니다.

## 시작하기

### 필수 요구 사항

- [Docker](https://www.docker.com/products/docker-desktop/)
- [VS Code](https://code.visualstudio.com/)
- [VS Code 용 Dev Containers 확장](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
- NVIDIA GPU와 [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)

### 환경 시작하기

1. VS Code에서 이 프로젝트 폴더를 엽니다
2. 명령 팔레트(F1 또는 Ctrl+Shift+P)에서 "Dev Containers: Reopen in Container"를 선택합니다
3. 컨테이너가 빌드되고 시작될 때까지 기다립니다

## 주요 기능

- CUDA 12.6 + cuDNN 기반 PyTorch
- Python 3 개발 환경
- 데이터 과학 도구 (NumPy, Pandas, Matplotlib, scikit-learn)
- VS Code에 최적화된 Python 개발 환경 (코드 포맷팅, 린팅, 자동완성)
- Jupyter 노트북 지원

## 예제 실행하기

### PyTorch 테스트

기본 PyTorch 설치 및 CUDA 가용성을 확인합니다:

```bash
python test.py
```

### 회귀 모델 예제

간단한 PyTorch 신경망을 사용한 회귀 예제:

```bash
python example.py
```

## 폴더 구조

- `app/`: Python 소스 코드 (컨테이너 내부의 `/app`에 마운트됨)
- `.devcontainer/`: DevContainer 설정 파일
  - `devcontainer.json`: DevContainer 구성
  - `Dockerfile`: 개발 환경 도커 이미지 정의
- `docker-compose.yml`: Docker Compose 설정

## 라이선스

MIT
