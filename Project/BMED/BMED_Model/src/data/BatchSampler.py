from torch.utils.data import Sampler
import numpy as np

'''
실험 단위로 배치를 구성하는 샘플러
각 실험은 연속적으로 처리되며 시계열 순서가 유지됨
'''
class BatchSampler(Sampler):
    def __init__(self, dataset):
        self.dataset = dataset
        self.exp_indices = dataset.exp_indices
        self.cur_batch = []

    def __iter__(self):
        # 실험 ID 리스트 생성 및 섞기 (실험 순서는 랜덤)
        exp_ids = [exp_id for exp_id, _ in self.exp_indices]
        np.random.shuffle(exp_ids)

        for exp_id in exp_ids:
            # 각 실험마다 해당 실험의 인덱스 가져오기
            cur_indices = None # 변수 초기화
            for e_id, indices in self.exp_indices:
                if e_id == exp_id:
                    cur_indices = indices # 현재 실험셋에 맞는 인덱스 세트 부여
                    break
            
            if cur_indices is None or len(cur_indices) == 0:
                continue
            
            # 하나의 실험 전체를 하나의 배치로 처리
            self.cur_batch = cur_indices
            yield self.cur_batch

    def __len__(self):
        # 대략적인 배치 수 계산
        return len(self.exp_indices)
    
    def get_cur_batch_indices(self):
        '''현재 batch의 index 목록 반환'''
        return self.cur_batch


