"""
1D 정상상태 열전도 문제 - Newton-Raphson with PyTorch Autograd

물리적 문제:
- 1D 저항 발열체 (길이 L = 1m)
- 열전도도 k = 1 W/m/K
- 경계조건: T(x=0) = 25°C, T(x=L) = 100°C (양단 Dirichlet)
- 선형 발열: q_gen = 0.1 W/m (길이방향 발열, 저항체)

지배방정식 (1D, 정상상태):
    d/dx(k * dT/dx) + q''' = 0
    
여기서 q''' = q_gen / A [W/m³] (체적 발열밀도)
단면적 A를 가정하면: q''' = 0.1 / A

FDM 이산화 (2차 중심차분, 내부 노드 i):
    k * (T_{i+1} - 2*T_i + T_{i-1}) / dx² + q''' = 0

Newton-Raphson:
    T_{n+1} = T_n - J^{-1} * R(T_n)
    
여기서 Jacobian J = ∂R/∂T 는 torch.autograd.functional.jacobian으로 계산
"""

import torch
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, List


# =============================================================================
# 물리적 파라미터
# =============================================================================
N_NODES = 1000                # 노드 수
L = 1.0                      # 도메인 길이 [m]
K = 1.0                      # 열전도도 [W/m/K]
T_LEFT = 25.0                # 왼쪽 경계 온도 [°C]
T_RIGHT = 25.0             # 오른쪽 경계 온도 [°C] (Dirichlet BC)
Q_LINEAR = 1               # 선형 발열량 [W/m]
A_CROSS = 0.001              # 단면적 [m²] (예: 1cm × 1cm)

# 체적 발열밀도 계산 [W/m³]
Q_VOL = Q_LINEAR / A_CROSS   # = 100 W/m³

# 격자 설정
DX = L / (N_NODES - 1)

# Device 설정
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# =============================================================================
# 좌표 생성
# =============================================================================
def create_grid(n_nodes: int, length: float, device: torch.device) -> torch.Tensor:
    """1D 격자 좌표 생성"""
    return torch.linspace(0, length, n_nodes, dtype=torch.float64, device=device)


# =============================================================================
# 잔차 함수 정의 (Residual Function)
# =============================================================================
def residual(T: torch.Tensor, 
                        k: float = K, 
                        q_vol: float = Q_VOL,
                        dx: float = DX,
                        h_left: float = 10.0,    # 입구 열전달계수 [W/m²/K]
                        h_right: float = 10.0,   # 출구 열전달계수 [W/m²/K]
                        T_inf: float = 25.0      # 외부 온도 [°C]
                        ) -> torch.Tensor:
    """
    대류 (Robin) 경계조건
    
    입구 (x=0):  열이 왼쪽으로 빠져나감
                 -k * dT/dx|_{x=0} = h_left * (T[0] - T_inf)
                 (외부가 왼쪽에 있으므로 부호 주의)
                 
    출구 (x=L):  열이 오른쪽으로 빠져나감  
                 -k * dT/dx|_{x=L} = h_right * (T[-1] - T_inf)
    """
    n = T.shape[0]
    R = torch.zeros_like(T)
    
    # =========================================================================
    # 입구 (x=0): 대류 경계조건
    # =========================================================================
    # 왼쪽 외부로 열전달: q = h*(T_surface - T_inf)
    # 에너지 밸런스: k*dT/dx = h*(T[0] - T_inf)  (열이 왼쪽으로 나감)
    # 2차 정확도 전방차분
    dTdx_left = (-3.0*T[0] + 4.0*T[1] - T[2]) / (2.0*dx)
    R[0] = k * dTdx_left - h_left * (T[0] - T_inf)
    
    # =========================================================================
    # 출구 (x=L): 대류 경계조건
    # =========================================================================
    # 오른쪽 외부로 열전달: -k*dT/dx = h*(T_surface - T_inf)
    # 2차 정확도 후방차분
    dTdx_right = (3.0*T[n-1] - 4.0*T[n-2] + T[n-3]) / (2.0*dx)
    R[n-1] = -k * dTdx_right - h_right * (T[n-1] - T_inf)
    
    # =========================================================================
    # 내부 노드: 열전도 + 발열
    # =========================================================================
    T_ip1 = T[2:]
    T_i = T[1:-1]
    T_im1 = T[:-2]
    
    d2T_dx2 = (T_ip1 - 2.0*T_i + T_im1) / (dx**2)
    R[1:-1] = k * d2T_dx2 + q_vol
    
    return R


# =============================================================================
# Newton-Raphson Solver (with Autograd Jacobian)
# =============================================================================
def newton_raphson_solve(
    T_init: torch.Tensor,
    residual_func,
    tol: float = 1e-8,
    max_iter: int = 100,
    verbose: bool = True
) -> Tuple[torch.Tensor, List[float]]:
    """
    Newton-Raphson 방법으로 비선형 시스템 R(T) = 0 해결
    
    알고리즘:
        1. Jacobian J = ∂R/∂T 계산 (autograd 사용)
        2. 선형 시스템 J * ΔT = -R 풀기
        3. 온도 업데이트: T_{n+1} = T_n + ΔT
        4. 수렴 확인: ||R|| < tol
    
    Args:
        T_init: 초기 온도 추정값
        residual_func: 잔차 함수 R(T)
        tol: 수렴 허용 오차
        max_iter: 최대 반복 횟수
        verbose: 출력 여부
    
    Returns:
        T_solution: 수렴된 온도 분포
        residual_history: 잔차 norm 이력
    """
    T = T_init.clone().detach()
    residual_history = []
    
    if verbose:
        print(f"\n{'='*60}")
        print("Newton-Raphson Iteration")
        print(f"{'='*60}")
        print(f"{'Iter':>5} | {'||R||':>15} | {'||ΔT||':>15}")
        print(f"{'-'*5}-+-{'-'*15}-+-{'-'*15}")
    
    for iteration in range(max_iter):
        # 잔차 계산
        R = residual_func(T)
        residual_norm = torch.norm(R).item()
        residual_history.append(residual_norm)
        print(f"residual_norm: {residual_norm}, tol: {tol}")
        # 수렴 확인
        if residual_norm < tol:
            if verbose:
                print(f"{iteration:5d} | {residual_norm:15.6e} | {'CONVERGED':>15}")
                print(f"\n✓ 수렴 완료! (iteration {iteration})")
            break
        
        # =====================================================================
        # Jacobian 계산 (PyTorch autograd 사용)
        # =====================================================================
        # torch.autograd.functional.jacobian(func, inputs)
        # J[i,j] = ∂R[i]/∂T[j]
        J = torch.autograd.functional.jacobian(residual_func, T)
        
        # =====================================================================
        # Newton-Raphson 업데이트
        # J * ΔT = -R  →  ΔT = J^{-1} * (-R)
        # =====================================================================
        delta_T = torch.linalg.solve(J, -R)
        
        # 온도 업데이트
        T = T + delta_T
        
        delta_norm = torch.norm(delta_T).item()
        
        if verbose and (iteration < 10 or iteration % 10 == 0):
            print(f"{iteration:5d} | {residual_norm:15.6e} | {delta_norm:15.6e}")
    
    else:
        if verbose:
            print(f"\n⚠ 최대 반복 횟수 도달 (max_iter = {max_iter})")
    
    return T.detach(), residual_history


# =============================================================================
# 해석해 (Analytical Solution)
# =============================================================================
def analytical_solution(
    x: torch.Tensor,
    k: float = K,
    q_vol: float = Q_VOL,
    T_left: float = T_LEFT,
    T_right: float = T_RIGHT,
    L: float = L
) -> torch.Tensor:
    """
    1D 열전도 방정식의 해석해 (발열 포함)
    
    지배방정식: k * d²T/dx² + q_vol = 0
    
    일반해:
        d²T/dx² = -q_vol/k
        dT/dx = -q_vol/k * x + C1
        T(x) = -q_vol/(2k) * x² + C1*x + C2
    
    경계조건 적용:
        T(0) = T_left  →  C2 = T_left
        T(L) = T_right →  C1 = (T_right - T_left)/L + q_vol*L/(2k)
    
    Args:
        x: 위치 좌표 [m]
        
    Returns:
        T: 해석해 온도 분포 [°C]
    """
    C2 = T_left
    C1 = (T_right - T_left) / L + q_vol * L / (2.0 * k)
    
    T_analytical = -q_vol / (2.0 * k) * x**2 + C1 * x + C2
    return T_analytical


# =============================================================================
# 시각화 함수
# =============================================================================
def plot_results(
    x: torch.Tensor,
    T_numerical: torch.Tensor,
    T_analytical: torch.Tensor,
    residual_history: List[float],
    save_path: str = None
):
    """결과 시각화"""
    
    x_np = x.cpu().numpy()
    T_num_np = T_numerical.cpu().numpy()
    T_ana_np = T_analytical.cpu().numpy()
    error = np.abs(T_num_np - T_ana_np)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 온도 프로파일 비교
    ax1 = axes[0, 0]
    ax1.plot(x_np, T_num_np, 'b-', lw=2, label='Numerical (FDM + NR)')
    ax1.plot(x_np, T_ana_np, 'r--', lw=2, label='Analytical')
    ax1.set_xlabel('Position x [m]')
    ax1.set_ylabel('Temperature T [°C]')
    ax1.set_title('Temperature Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 오차 분포
    ax2 = axes[0, 1]
    ax2.semilogy(x_np, error + 1e-16, 'g-', lw=2)
    ax2.set_xlabel('Position x [m]')
    ax2.set_ylabel('Absolute Error [°C]')
    ax2.set_title(f'Error (Max: {np.max(error):.2e} °C)')
    ax2.grid(True, alpha=0.3)
    
    # 3. Newton-Raphson 수렴
    ax3 = axes[1, 0]
    ax3.semilogy(residual_history, 'ko-', ms=6, lw=2)
    ax3.axhline(y=1e-10, color='r', ls='--', label='Tolerance')
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('Residual Norm ||R||')
    ax3.set_title('Newton-Raphson Convergence')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. 열 플럭스
    ax4 = axes[1, 1]
    dT_dx = np.gradient(T_num_np, x_np)
    q_flux = -K * dT_dx
    ax4.plot(x_np, q_flux, 'm-', lw=2)
    ax4.set_xlabel('Position x [m]')
    ax4.set_ylabel('Heat Flux q [W/m²]')
    ax4.set_title('Heat Flux Profile')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"그래프 저장: {save_path}")
    
    plt.show()


# =============================================================================
# 메인 실행
# =============================================================================
def main():
    print("="*60)
    print("1D Heat Conduction with Resistive Heating")
    print("Newton-Raphson Solver with PyTorch Autograd Jacobian")
    print("="*60)
    
    # 파라미터 출력
    print(f"\n[물리적 파라미터]")
    print(f"  - 도메인 길이 L = {L} m")
    print(f"  - 열전도도 k = {K} W/m/K")
    print(f"  - 선형 발열 q_lin = {Q_LINEAR} W/m")
    print(f"  - 단면적 A = {A_CROSS} m²")
    print(f"  - 체적 발열밀도 q''' = {Q_VOL} W/m³")
    print(f"\n[경계조건]")
    print(f"  - T(x=0) = {T_LEFT} °C (Dirichlet)")
    print(f"  - T(x=L) = {T_RIGHT} °C (Dirichlet)")
    print(f"\n[수치해석 설정]")
    print(f"  - 노드 수 N = {N_NODES}")
    print(f"  - 격자 간격 dx = {DX:.6f} m")
    print(f"  - Device: {DEVICE}")
    
    # 격자 생성
    x = create_grid(N_NODES, L, DEVICE)
    
    # 초기 온도 추정 (선형 분포)
    T_init = torch.linspace(T_LEFT, T_RIGHT, N_NODES, 
                            dtype=torch.float64, device=DEVICE)
    
    # Newton-Raphson 풀이
    T_numerical, res_history = newton_raphson_solve(
        T_init, 
        residual,
        tol=1e-8,
        max_iter=100,
        verbose=True
    )
    
    # 해석해 계산
    T_analytical = analytical_solution(x)
    
    # 결과 검증
    error = torch.abs(T_numerical - T_analytical)
    max_error = torch.max(error).item()
    rms_error = torch.sqrt(torch.mean(error**2)).item()
    
    print(f"\n{'='*60}")
    print("결과 검증")
    print(f"{'='*60}")
    print(f"  최대 오차: {max_error:.6e} °C")
    print(f"  RMS 오차:  {rms_error:.6e} °C")
    print(f"\n[경계조건 확인]")
    print(f"  T(0) = {T_numerical[0].item():.6f} °C (목표: {T_LEFT})")
    print(f"  T(L) = {T_numerical[-1].item():.6f} °C (목표: {T_RIGHT})")
    
    # 최고 온도 위치 (포물선 형태이므로 내부에 최대값 존재 가능)
    T_max = torch.max(T_numerical).item()
    idx_max = torch.argmax(T_numerical).item()
    x_max = x[idx_max].item()
    print(f"\n[최고 온도]")
    print(f"  T_max = {T_max:.4f} °C at x = {x_max:.4f} m")
    
    # 일부 노드 출력
    print(f"\n[온도 프로파일 (일부)]")
    print(f"{'x [m]':>10} | {'T_num [°C]':>12} | {'T_ana [°C]':>12} | {'Error':>12}")
    print("-"*52)
    for i in [0, 24, 49, 74, 99]:
        print(f"{x[i].item():10.4f} | {T_numerical[i].item():12.6f} | "
              f"{T_analytical[i].item():12.6f} | {error[i].item():12.2e}")
    
    # 시각화
    plot_results(x, T_numerical, T_analytical, res_history,
                 save_path='heat_conduction_results.png')
    
    return T_numerical, x


if __name__ == "__main__":
    T_solution, x_grid = main()