"""
2D 이상기체 유동 시뮬레이션 사용 예시
findiff를 이용한 다양한 케이스 실행
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # GUI 없는 백엔드 사용
import matplotlib.pyplot as plt

# 모듈 import 수정
import sys
import os
sys.path.append(os.path.dirname(__file__))

# 파일명에 하이픈이 있어서 직접 import 불가, importlib 사용
import importlib.util
spec = importlib.util.spec_from_file_location("gas_sim", "2D_gas_flow_simulation.py")
gas_sim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gas_sim)

GasFlowSimulator2D = gas_sim.GasFlowSimulator2D
create_visualization = gas_sim.create_visualization
create_animation = gas_sim.create_animation

def example1_basic_flow():
    """예시 1: 기본 압력 구동 유동"""
    print("예시 1: 기본 압력 구동 유동")
    print("-" * 30)
    
    # 기본 설정으로 시뮬레이터 생성
    sim = GasFlowSimulator2D(
        Lx=1.0,      # 길이 1m
        Ly=0.2,      # 높이 0.2m
        nx=60,       # x 방향 격자수
        ny=15,       # y 방향 격자수
        dt=1e-5      # 시간 간격
    )
    
    # 시뮬레이션 실행
    sim.run_simulation(max_steps=2000, save_interval=40)
    
    # 결과 시각화
    fig = create_visualization(sim)
    plt.savefig('basic_flow_result.png', dpi=150, bbox_inches='tight')
    print("결과 저장: basic_flow_result.png")
    
    return sim

def example2_high_pressure_ratio():
    """예시 2: 높은 압력비 유동 (초킹 현상)"""
    print("\n예시 2: 높은 압력비 유동")
    print("-" * 30)
    
    sim = GasFlowSimulator2D(Lx=2.0, Ly=0.3, nx=100, ny=20, dt=5e-6)
    
    # 높은 압력 차이 설정
    sim.p[:, 0] = 5.0 * 101325   # 입구: 5 atm
    sim.p[:, -1] = 0.5 * 101325  # 출구: 0.5 atm
    
    # 밀도 재계산
    sim.rho[:, 0] = sim.p[:, 0] / (sim.R * sim.T[:, 0])
    sim.rho[:, -1] = sim.p[:, -1] / (sim.R * sim.T[:, -1])
    
    sim.run_simulation(max_steps=3000, save_interval=50)
    
    fig = create_visualization(sim)
    plt.savefig('high_pressure_ratio_result.png', dpi=150, bbox_inches='tight')
    print("결과 저장: high_pressure_ratio_result.png")
    
    return sim

def example3_heated_tube():
    """예시 3: 가열된 튜브 유동"""
    print("\n예시 3: 가열된 튜브 유동")
    print("-" * 30)
    
    sim = GasFlowSimulator2D(Lx=1.5, Ly=0.25, nx=80, ny=20, dt=8e-6)
    
    # 벽면 온도를 높게 설정 (가열)
    def apply_heated_boundary(self):
        # 기존 경계 조건
        self.u[0, :] = 0
        self.u[-1, :] = 0
        self.v[0, :] = 0
        self.v[-1, :] = 0
        
        # 가열된 벽면 온도
        self.T[0, :] = 400   # 하부 벽면 가열
        self.T[-1, :] = 350  # 상부 벽면 가열
        
        # 압력 경계 조건
        self.p[:, 0] = 2.0 * 101325
        self.p[:, -1] = 1.0 * 101325
        
        # 밀도 업데이트
        self.rho[:, 0] = self.p[:, 0] / (self.R * self.T[:, 0])
        self.rho[:, -1] = self.p[:, -1] / (self.R * self.T[:, -1])
    
    # 경계조건 함수 교체
    sim.apply_boundary_conditions = lambda: apply_heated_boundary(sim)
    
    sim.run_simulation(max_steps=2500, save_interval=40)
    
    fig = create_visualization(sim)
    plt.savefig('heated_tube_result.png', dpi=150, bbox_inches='tight')
    print("결과 저장: heated_tube_result.png")
    
    return sim

def example4_convergent_nozzle():
    """예시 4: 수렴 노즐 형태"""
    print("\n예시 4: 수렴 노즐 형태")
    print("-" * 30)
    
    sim = GasFlowSimulator2D(Lx=1.0, Ly=0.4, nx=60, ny=25, dt=1e-5)
    
    # 노즐 형태의 경계 생성 (벽면을 막아서 구현)
    for i in range(sim.nx):
        x_pos = sim.x[i]
        # 선형적으로 줄어드는 노즐
        throat_height = 0.4 * (1 - 0.7 * x_pos / sim.Lx)  # 입구 40cm → 출구 12cm
        
        # 상하부 벽면 설정
        for j in range(sim.ny):
            y_pos = sim.y[j]
            if y_pos < (0.4 - throat_height) / 2 or y_pos > (0.4 + throat_height) / 2:
                # 벽면 영역: 속도 = 0
                sim.u[j, i] = 0
                sim.v[j, i] = 0
    
    sim.run_simulation(max_steps=2000, save_interval=35)
    
    fig = create_visualization(sim)
    plt.savefig('convergent_nozzle_result.png', dpi=150, bbox_inches='tight')
    print("결과 저장: convergent_nozzle_result.png")
    
    return sim

def analyze_results(sim, title="시뮬레이션 결과 분석"):
    """결과 분석 및 주요 특성값 계산"""
    print(f"\n{title}")
    print("=" * len(title))
    
    # 최대/평균 속도
    vel_magnitude = np.sqrt(sim.u**2 + sim.v**2)
    max_velocity = np.max(vel_magnitude)
    avg_velocity = np.mean(vel_magnitude)
    
    # 압력 강하
    inlet_pressure = np.mean(sim.p[:, 0])
    outlet_pressure = np.mean(sim.p[:, -1])
    pressure_drop = inlet_pressure - outlet_pressure
    
    # 질량 유량 (입구에서)
    mass_flow_rate = np.sum(sim.rho[:, 0] * sim.u[:, 0] * sim.dy)
    
    # 레이놀즈 수 (근사)
    rho_avg = np.mean(sim.rho)
    u_avg = np.mean(sim.u[sim.u > 0])
    hydraulic_diameter = 4 * sim.Ly  # 평판 채널의 수력직경
    reynolds = rho_avg * u_avg * hydraulic_diameter / sim.mu if u_avg > 0 else 0
    
    # 마하수 (근사)
    T_avg = np.mean(sim.T)
    speed_of_sound = np.sqrt(sim.gamma * sim.R * T_avg)
    mach_number = max_velocity / speed_of_sound if speed_of_sound > 0 else 0
    
    print(f"최대 속도: {max_velocity:.2f} m/s")
    print(f"평균 속도: {avg_velocity:.2f} m/s")
    print(f"압력 강하: {pressure_drop/1000:.2f} kPa")
    print(f"질량 유량: {mass_flow_rate:.6f} kg/s·m")
    print(f"레이놀즈 수: {reynolds:.0f}")
    print(f"마하수: {mach_number:.3f}")
    print(f"음속: {speed_of_sound:.1f} m/s")

def plot_centerline_profiles(sim, title="중심선 프로파일"):
    """중심선에서의 물리량 분포 그래프"""
    center_idx = sim.ny // 2
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(title, fontsize=14)
    
    # 속도 분포
    vel_magnitude = np.sqrt(sim.u**2 + sim.v**2)
    axes[0,0].plot(sim.x, vel_magnitude[center_idx, :], 'b-', linewidth=2)
    axes[0,0].set_xlabel('x [m]')
    axes[0,0].set_ylabel('속도 [m/s]')
    axes[0,0].set_title('중심선 속도 분포')
    axes[0,0].grid(True)
    
    # 압력 분포
    axes[0,1].plot(sim.x, sim.p[center_idx, :]/1000, 'r-', linewidth=2)
    axes[0,1].set_xlabel('x [m]')
    axes[0,1].set_ylabel('압력 [kPa]')
    axes[0,1].set_title('중심선 압력 분포')
    axes[0,1].grid(True)
    
    # 밀도 분포
    axes[1,0].plot(sim.x, sim.rho[center_idx, :], 'g-', linewidth=2)
    axes[1,0].set_xlabel('x [m]')
    axes[1,0].set_ylabel('밀도 [kg/m³]')
    axes[1,0].set_title('중심선 밀도 분포')
    axes[1,0].grid(True)
    
    # 온도 분포
    axes[1,1].plot(sim.x, sim.T[center_idx, :], 'm-', linewidth=2)
    axes[1,1].set_xlabel('x [m]')
    axes[1,1].set_ylabel('온도 [K]')
    axes[1,1].set_title('중심선 온도 분포')
    axes[1,1].grid(True)
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    print("findiff를 이용한 2D 이상기체 유동 시뮬레이션 예시")
    print("=" * 60)
    
    # 예시 실행
    examples = {
        1: example1_basic_flow,
        2: example2_high_pressure_ratio,
        3: example3_heated_tube,
        4: example4_convergent_nozzle
    }
    
    print("실행할 예시를 선택하세요:")
    print("1: 기본 압력 구동 유동")
    print("2: 높은 압력비 유동")
    print("3: 가열된 튜브 유동")
    print("4: 수렴 노즐 형태")
    print("0: 모든 예시 실행")
    
    try:
        choice = int(input("\n선택 (0-4): "))
        
        if choice == 0:
            # 모든 예시 실행
            for i, (num, func) in enumerate(examples.items(), 1):
                print(f"\n{'='*20} 예시 {num} 실행 {'='*20}")
                sim = func()
                analyze_results(sim, f"예시 {num} 결과 분석")
                
                # 중심선 프로파일 그래프
                fig = plot_centerline_profiles(sim, f"예시 {num} 중심선 프로파일")
                plt.savefig(f'example_{num}_profiles.png', dpi=150, bbox_inches='tight')
                print(f"프로파일 저장: example_{num}_profiles.png")
                
        elif choice in examples:
            # 선택한 예시만 실행
            sim = examples[choice]()
            analyze_results(sim)
            
            # 중심선 프로파일 그래프
            fig = plot_centerline_profiles(sim)
            plt.savefig(f'example_{choice}_profiles.png', dpi=150, bbox_inches='tight')
            print(f"프로파일 저장: example_{choice}_profiles.png")
            
            # 애니메이션 생성 옵션
            create_ani = input("\n애니메이션을 생성하시겠습니까? (y/n): ")
            if create_ani.lower() == 'y':
                print("애니메이션 생성 중...")
                ani = create_animation(sim)
                if ani:
                    ani.save(f'example_{choice}_animation.gif', writer='pillow', fps=8)
                    print(f"애니메이션 저장: example_{choice}_animation.gif")
        else:
            print("올바른 번호를 선택해주세요.")
            
    except (ValueError, KeyboardInterrupt):
        print("\n프로그램을 종료합니다.")
    
    print("\n시뮬레이션 완료!")