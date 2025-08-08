"""
간단한 2D 이상기체 유동 테스트
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # GUI 없는 백엔드
import matplotlib.pyplot as plt

# 모듈 import
import importlib.util
spec = importlib.util.spec_from_file_location("gas_sim", "2D_gas_flow_simulation.py")
gas_sim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gas_sim)

def quick_test():
    """빠른 테스트 실행"""
    print("findiff 2D 이상기체 유동 시뮬레이션 테스트")
    print("=" * 50)
    
    # 작은 격자로 빠른 테스트
    simulator = gas_sim.GasFlowSimulator2D(
        Lx=1.0,      # 길이 1m
        Ly=0.2,      # 높이 0.2m  
        nx=40,       # x 방향 격자수 (작게)
        ny=10,       # y 방향 격자수 (작게)
        dt=2e-5      # 시간 간격
    )
    
    print("시뮬레이터 초기화 완료")
    print(f"격자: {simulator.nx} x {simulator.ny}")
    print(f"시간 간격: {simulator.dt}")
    
    # 짧은 시뮬레이션 실행
    print("\n시뮬레이션 실행 중...")
    simulator.run_simulation(max_steps=500, save_interval=50)
    
    # 기본 결과 출력
    vel_magnitude = np.sqrt(simulator.u**2 + simulator.v**2)
    max_vel = np.max(vel_magnitude)
    avg_vel = np.mean(vel_magnitude)
    
    print(f"\n결과:")
    print(f"최대 속도: {max_vel:.2f} m/s")
    print(f"평균 속도: {avg_vel:.2f} m/s")
    print(f"압력 범위: {np.min(simulator.p)/1000:.1f} - {np.max(simulator.p)/1000:.1f} kPa")
    
    # 간단한 그래프 저장
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    
    # 속도 분포
    im1 = ax1.contourf(simulator.X, simulator.Y, vel_magnitude, levels=15, cmap='viridis')
    ax1.set_title('Velocity Magnitude [m/s]')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    plt.colorbar(im1, ax=ax1)
    
    # 압력 분포
    im2 = ax2.contourf(simulator.X, simulator.Y, simulator.p/1000, levels=15, cmap='plasma')
    ax2.set_title('Pressure [kPa]')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    plt.colorbar(im2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('test_result.png', dpi=150, bbox_inches='tight')
    print("\n테스트 결과 저장: test_result.png")
    
    return simulator

if __name__ == "__main__":
    try:
        sim = quick_test()
        print("\n✅ 테스트 성공! 시뮬레이션이 정상 작동합니다.")
    except Exception as e:
        print(f"\n❌ 테스트 실패: {e}")
        import traceback
        traceback.print_exc()