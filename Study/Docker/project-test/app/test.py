import torch

print("🔥 PyTorch version:", torch.__version__)
print("🧠 CUDA available:", torch.cuda.is_available())

if torch.cuda.is_available():
    print("🚀 CUDA version (detected by torch):", torch.version.cuda)
    print("🎯 Current device index:", torch.cuda.current_device())
    print("📛 GPU name:", torch.cuda.get_device_name(torch.cuda.current_device()))
    print("📈 Total GPU memory:", torch.cuda.get_device_properties(0).total_memory / (1024 ** 3), "GB")
else:
    print("⚠️ CUDA is NOT available. CPU-only mode.")