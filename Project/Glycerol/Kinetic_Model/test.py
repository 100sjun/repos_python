import torch, torchdiffeq, pathlib
print("torch :", torch.__version__, torch.version.cuda, torch.cuda.is_available())
print("torch path:", pathlib.Path(torch.__file__).parent)
print("torchdiffeq :", torchdiffeq.__version__)
from torchdiffeq import odeint_adjoint  # ← 에러 없으면 성공
