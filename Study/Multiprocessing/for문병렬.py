import multiprocessing as mp
from tqdm import tqdm

def process_function(x):
    # 여기에 병렬로 실행할 작업을 정의합니다
    return x * x

if __name__ == '__main__':
    # 데이터 준비
    data = list(range(1000))

    # 프로세스 풀 생성
    pool = mp.Pool(processes=mp.cpu_count())

    # tqdm을 사용하여 진행 상황을 표시하면서 병렬 처리 실행
    results = list(tqdm(pool.imap(process_function, data), total=len(data)))

    # 풀 종료
    pool.close()
    pool.join()

    print("처리 완료:", results[:10])  # 결과의 일부만 출력