"""
2D 이상기체 유동 시뮬레이션 - findiff 사용
튜브 내에서 압력 구배에 의한 기체 유동 모사

지배방정식:
- 연속방정식 (질량 보존): ∂ρ/∂t + ∇·(ρv) = 0
- 모멘텀 방정식: ∂(ρv)/∂t + ∇·(ρvv) = -∇p + μ∇²v
- 에너지 방정식: ∂(ρE)/∂t + ∇·((ρE+p)v) = ∇·(k∇T) + μΦ
- 이상기체 상태방정식: p = ρRT
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # GUI 없는 백엔드 사용
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from findiff import FinDiff
import time

class GasFlowSimulator2D:
    def __init__(self, Lx=2.0, Ly=0.5, nx=100, ny=25, dt=1e-5):
        """
        2D 이상기체 유동 시뮬레이터 초기화
        
        Parameters:
        - Lx, Ly: 계산 영역 크기 [m]
        - nx, ny: 격자 점 개수
        - dt: 시간 간격 [s]
        """
        self.Lx, self.Ly = Lx, Ly
        self.nx, self.ny = nx, ny
        self.dx, self.dy = Lx/(nx-1), Ly/(ny-1)
        self.dt = dt
        
        # 좌표 생성
        self.x = np.linspace(0, Lx, nx)
        self.y = np.linspace(0, Ly, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        # 유체 물성치 (공기, 상온)
        self.R = 287.0      # 기체상수 [J/kg·K]
        self.gamma = 1.4    # 비열비
        self.mu = 1.8e-5    # 동점성 계수 [Pa·s]
        self.k = 0.025      # 열전도도 [W/m·K]
        self.cp = 1005.0    # 정압비열 [J/kg·K]
        
        # findiff 미분 연산자 정의
        self.d_dx = FinDiff(1, self.dx, 1)  # ∂/∂x
        self.d_dy = FinDiff(0, self.dy, 1)  # ∂/∂y
        self.d2_dx2 = FinDiff(1, self.dx, 2)  # ∂²/∂x²
        self.d2_dy2 = FinDiff(0, self.dy, 2)  # ∂²/∂y²
        
        self.initialize_fields()
        
    def initialize_fields(self):
        """초기 조건 및 경계 조건 설정"""
        # 초기 조건 (정지 상태)
        self.rho = np.ones((self.ny, self.nx)) * 1.2  # 밀도 [kg/m³]
        self.u = np.zeros((self.ny, self.nx))         # x 방향 속도 [m/s]
        self.v = np.zeros((self.ny, self.nx))         # y 방향 속도 [m/s]
        self.T = np.ones((self.ny, self.nx)) * 300    # 온도 [K]
        self.p = self.rho * self.R * self.T           # 압력 [Pa]
        
        # 입구 조건 (x=0에서 압력 높음)
        self.p[:, 0] = 1.5 * 101325  # 1.5 atm
        self.rho[:, 0] = self.p[:, 0] / (self.R * self.T[:, 0])
        
        # 출구 조건 (x=Lx에서 압력 낮음)
        self.p[:, -1] = 0.8 * 101325  # 0.8 atm
        self.rho[:, -1] = self.p[:, -1] / (self.R * self.T[:, -1])
        
    def apply_boundary_conditions(self):
        """경계 조건 적용"""
        # 벽면 조건 (no-slip, y=0, y=Ly)
        self.u[0, :] = 0   # 하부 벽면
        self.u[-1, :] = 0  # 상부 벽면
        self.v[0, :] = 0   # 하부 벽면
        self.v[-1, :] = 0  # 상부 벽면
        
        # 온도 경계 조건 (등온 벽면)
        self.T[0, :] = 300   # 하부 벽면
        self.T[-1, :] = 300  # 상부 벽면
        
        # 입구/출구 압력 고정
        self.p[:, 0] = 1.5 * 101325   # 입구
        self.p[:, -1] = 0.8 * 101325  # 출구
        
        # 상태방정식으로 밀도 업데이트
        self.rho[:, 0] = self.p[:, 0] / (self.R * self.T[:, 0])
        self.rho[:, -1] = self.p[:, -1] / (self.R * self.T[:, -1])
    
    def compute_derivatives(self):
        """findiff를 이용한 미분 계산"""
        try:
            # 압력 구배
            dp_dx = self.d_dx(self.p)
            dp_dy = self.d_dy(self.p)
            
            # 속도 구배 (점성 항)
            d2u_dx2 = self.d2_dx2(self.u)
            d2u_dy2 = self.d2_dy2(self.u)
            d2v_dx2 = self.d2_dx2(self.v)
            d2v_dy2 = self.d2_dy2(self.v)
            
            # 온도 구배 (열전도 항)
            d2T_dx2 = self.d2_dx2(self.T)
            d2T_dy2 = self.d2_dy2(self.T)
            
            return dp_dx, dp_dy, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2, d2T_dx2, d2T_dy2
            
        except Exception as e:
            print(f"미분 계산 오류: {e}")
            return [np.zeros_like(self.p) for _ in range(8)]
    
    def solve_momentum(self, dp_dx, dp_dy, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2):
        """모멘텀 방정식 해석 (간단화된 형태)"""
        # 점성력 계산
        visc_u = self.mu * (d2u_dx2 + d2u_dy2)
        visc_v = self.mu * (d2v_dx2 + d2v_dy2)
        
        # 관성항 (간단화)
        inertia_u = self.rho * (self.u * self.d_dx(self.u) + self.v * self.d_dy(self.u))
        inertia_v = self.rho * (self.u * self.d_dx(self.v) + self.v * self.d_dy(self.v))
        
        # 모멘텀 방정식: ρ(Du/Dt) = -∇p + μ∇²u
        du_dt = (-dp_dx + visc_u - inertia_u) / (self.rho + 1e-10)
        dv_dt = (-dp_dy + visc_v - inertia_v) / (self.rho + 1e-10)
        
        # 속도 업데이트
        self.u += self.dt * du_dt
        self.v += self.dt * dv_dt
        
        # 안정성을 위한 속도 제한
        max_vel = 100  # m/s
        self.u = np.clip(self.u, -max_vel, max_vel)
        self.v = np.clip(self.v, -max_vel, max_vel)
    
    def solve_continuity(self):
        """연속방정식 해석 (간단화된 형태)"""
        # ∂ρ/∂t + ∇·(ρv) = 0
        try:
            div_rho_u = self.d_dx(self.rho * self.u) + self.d_dy(self.rho * self.v)
            drho_dt = -div_rho_u
            self.rho += self.dt * drho_dt
            
            # 밀도 양수 유지
            self.rho = np.maximum(self.rho, 0.1)
            
        except Exception as e:
            print(f"연속방정식 해석 오류: {e}")
    
    def solve_energy(self, d2T_dx2, d2T_dy2):
        """에너지 방정식 해석 (간단화된 형태)"""
        # 열전도 항
        heat_conduction = self.k * (d2T_dx2 + d2T_dy2)
        
        # 대류 항 (간단화)
        convection = self.rho * self.cp * (self.u * self.d_dx(self.T) + self.v * self.d_dy(self.T))
        
        # 온도 업데이트
        dT_dt = (heat_conduction - convection) / (self.rho * self.cp + 1e-10)
        self.T += self.dt * dT_dt
        
        # 온도 범위 제한
        self.T = np.clip(self.T, 250, 400)
    
    def update_pressure(self):
        """이상기체 상태방정식으로 압력 업데이트"""
        p_new = self.rho * self.R * self.T
        
        # 압력 변화 제한 (안정성)
        dp_max = 0.1 * np.abs(self.p)
        dp = np.clip(p_new - self.p, -dp_max, dp_max)
        self.p += dp
    
    def time_step(self):
        """한 시간 스텝 계산"""
        # 미분 계산
        dp_dx, dp_dy, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2, d2T_dx2, d2T_dy2 = self.compute_derivatives()
        
        # 방정식 해석
        self.solve_momentum(dp_dx, dp_dy, d2u_dx2, d2u_dy2, d2v_dx2, d2v_dy2)
        self.solve_continuity()
        self.solve_energy(d2T_dx2, d2T_dy2)
        self.update_pressure()
        
        # 경계 조건 적용
        self.apply_boundary_conditions()
    
    def run_simulation(self, max_steps=5000, save_interval=50):
        """시뮬레이션 실행"""
        print("2D 이상기체 유동 시뮬레이션 시작...")
        
        # 결과 저장용 배열
        self.time_history = []
        self.velocity_history = []
        self.pressure_history = []
        self.density_history = []
        self.temperature_history = []
        
        start_time = time.time()
        
        for step in range(max_steps):
            self.time_step()
            
            # 결과 저장
            if step % save_interval == 0:
                current_time = step * self.dt
                self.time_history.append(current_time)
                
                # 속도 크기 계산
                vel_magnitude = np.sqrt(self.u**2 + self.v**2)
                self.velocity_history.append(vel_magnitude.copy())
                self.pressure_history.append(self.p.copy())
                self.density_history.append(self.rho.copy())
                self.temperature_history.append(self.T.copy())
                
                # 진행 상황 출력
                if step % (save_interval * 10) == 0:
                    max_vel = np.max(vel_magnitude)
                    avg_p = np.mean(self.p)
                    print(f"Step {step:5d}: t = {current_time:.3f}s, 최대속도 = {max_vel:.2f} m/s, 평균압력 = {avg_p:.0f} Pa")
        
        elapsed_time = time.time() - start_time
        print(f"시뮬레이션 완료! 소요시간: {elapsed_time:.2f}초")
        
        return np.array(self.velocity_history), np.array(self.pressure_history)

def create_visualization(simulator):
    """시뮬레이션 결과 시각화"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('2D 이상기체 유동 시뮬레이션 (findiff 사용)', fontsize=16)
    
    # 최종 상태 시각화
    vel_magnitude = np.sqrt(simulator.u**2 + simulator.v**2)
    
    # 속도 분포
    im1 = axes[0,0].contourf(simulator.X, simulator.Y, vel_magnitude, levels=20, cmap='viridis')
    axes[0,0].set_title('속도 크기 분포 [m/s]')
    axes[0,0].set_xlabel('x [m]')
    axes[0,0].set_ylabel('y [m]')
    plt.colorbar(im1, ax=axes[0,0])
    
    # 유선도
    axes[0,0].streamplot(simulator.X, simulator.Y, simulator.u, simulator.v, 
                        density=2, color='white')
    
    # 압력 분포
    im2 = axes[0,1].contourf(simulator.X, simulator.Y, simulator.p/1000, levels=20, cmap='plasma')
    axes[0,1].set_title('압력 분포 [kPa]')
    axes[0,1].set_xlabel('x [m]')
    axes[0,1].set_ylabel('y [m]')
    plt.colorbar(im2, ax=axes[0,1])
    
    # 밀도 분포
    im3 = axes[1,0].contourf(simulator.X, simulator.Y, simulator.rho, levels=20, cmap='coolwarm')
    axes[1,0].set_title('밀도 분포 [kg/m³]')
    axes[1,0].set_xlabel('x [m]')
    axes[1,0].set_ylabel('y [m]')
    plt.colorbar(im3, ax=axes[1,0])
    
    # 온도 분포
    im4 = axes[1,1].contourf(simulator.X, simulator.Y, simulator.T, levels=20, cmap='hot')
    axes[1,1].set_title('온도 분포 [K]')
    axes[1,1].set_xlabel('x [m]')
    axes[1,1].set_ylabel('y [m]')
    plt.colorbar(im4, ax=axes[1,1])
    
    plt.tight_layout()
    plt.savefig('2d_gas_flow_result.png', dpi=150, bbox_inches='tight')
    print("결과 저장: 2d_gas_flow_result.png")
    return fig

def create_animation(simulator):
    """애니메이션 생성"""
    if not hasattr(simulator, 'velocity_history'):
        print("시뮬레이션을 먼저 실행해주세요.")
        return None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    def animate(frame):
        ax1.clear()
        ax2.clear()
        
        # 속도 애니메이션
        vel_data = simulator.velocity_history[frame]
        im1 = ax1.contourf(simulator.X, simulator.Y, vel_data, levels=20, cmap='viridis')
        ax1.set_title(f'속도 크기 [m/s] - t = {simulator.time_history[frame]:.3f}s')
        ax1.set_xlabel('x [m]')
        ax1.set_ylabel('y [m]')
        
        # 압력 애니메이션
        pressure_data = simulator.pressure_history[frame]
        im2 = ax2.contourf(simulator.X, simulator.Y, pressure_data/1000, levels=20, cmap='plasma')
        ax2.set_title(f'압력 [kPa] - t = {simulator.time_history[frame]:.3f}s')
        ax2.set_xlabel('x [m]')
        ax2.set_ylabel('y [m]')
        
        return [im1, im2]
    
    ani = FuncAnimation(fig, animate, frames=len(simulator.velocity_history), 
                       interval=100, blit=False, repeat=True)
    
    return ani

if __name__ == "__main__":
    # 시뮬레이션 실행
    print("findiff를 이용한 2D 이상기체 유동 시뮬레이션")
    print("=" * 50)
    
    # 시뮬레이터 생성
    simulator = GasFlowSimulator2D(Lx=2.0, Ly=0.5, nx=80, ny=20, dt=5e-6)
    
    # 시뮬레이션 실행
    vel_history, pressure_history = simulator.run_simulation(max_steps=3000, save_interval=30)
    
    # 결과 시각화
    fig = create_visualization(simulator)
    
    # 애니메이션 생성 (선택사항)
    # ani = create_animation(simulator)
    # if ani:
    #     ani.save('2d_gas_flow.gif', writer='pillow', fps=10)
    #     print("애니메이션 저장 완료: 2d_gas_flow.gif")