#!/usr/bin/env python3
"""
BMED 모델 컴포넌트 테스트 스크립트
이 스크립트는 하이퍼파라미터 최적화 전에 모든 컴포넌트가 올바르게 작동하는지 확인합니다.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

class LSTMEncoder(nn.Module):
    """
    개선된 LSTM 인코더 - device 호환성 및 안정성 강화
    """
    def __init__(self, input_size, hidden_size, num_layers, dropout=0.2):
        super().__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        
        # LSTM layer with improved error handling
        self.lstm = nn.LSTM(
            input_size, hidden_size, num_layers, 
            batch_first=True, dropout=dropout if num_layers > 1 else 0
        )
        
        self.layer_norm = nn.LayerNorm(hidden_size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, seq_len):
        """
        Forward pass with robust device handling
        """
        # 입력 검증
        if x.size(0) != seq_len.size(0):
            raise ValueError(f"Batch size mismatch: input {x.size(0)} vs seq_len {seq_len.size(0)}")
        
        # seq_len을 CPU로 이동하고 정렬
        seq_len_cpu = seq_len.detach().cpu().long()
        
        # 시퀀스 길이가 0인 경우 처리
        if (seq_len_cpu <= 0).any():
            seq_len_cpu = torch.clamp(seq_len_cpu, min=1)
        
        try:
            # 패딩된 시퀀스를 pack하여 효율적 처리
            packed_input = pack_padded_sequence(
                x, seq_len_cpu, batch_first=True, enforce_sorted=False
            )
            packed_output, (hidden, cell) = self.lstm(packed_input)
            
            # 다시 패딩된 형태로 복원 (원래 max_len 길이로)
            lstm_out, output_lengths = pad_packed_sequence(
                packed_output, batch_first=True, total_length=x.size(1)
            )
            
        except Exception as e:
            # Pack/unpack 실패시 일반 LSTM forward pass 사용
            print(f"Warning: Pack/unpack failed ({e}), using standard LSTM forward")
            lstm_out, (hidden, cell) = self.lstm(x)
        
        # Normalization and dropout
        normalized = self.layer_norm(lstm_out)
        return self.dropout(normalized)

class MLPDecoder(nn.Module):
    def __init__(self, hidden_size, output_size, num_layers=2, num_nodes=None, dropout = 0.3):
        super().__init__()

        if num_nodes is None:
            num_nodes = hidden_size
        
        self.layers = nn.ModuleList()

        # 첫 번째 레이어: hidden_size → num_nodes
        self.layers.append(nn.Linear(hidden_size, num_nodes))
        self.layers.append(nn.LayerNorm(num_nodes))
        self.layers.append(nn.ReLU())
        self.layers.append(nn.Dropout(dropout))

        # 중간 은닉층들: num_nodes → num_nodes
        for i in range(num_layers - 1):
            self.layers.append(nn.Linear(num_nodes,num_nodes))
            self.layers.append(nn.LayerNorm(num_nodes))
            self.layers.append(nn.ReLU())
            self.layers.append(nn.Dropout(dropout))

        # 마지막 출력층: num_nodes → output_size
        self.layers.append(nn.Linear(num_nodes, output_size))

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        return x

class StateUpdateLayer(nn.Module):
    """
    개선된 상태 업데이트 레이어 - 물리적 제약 조건 적용
    Bipolar membrane electrodialysis 물리 법칙 기반
    """
    def __init__(self, eps=1e-8):
        super().__init__()
        self.eps = eps  # division by zero 방지
        
    def forward(self, mlp_output, cur_state):
        """
        물리적 제약 조건을 적용하여 다음 상태 계산
        
        Args:
            mlp_output: [batch, seq, 7] - [dVA, dVB, dNALA, dNBLA, dNAK, dNBK, nI]
            cur_state: [batch, seq, 12] - 현재 상태
        """
        # 입력 차원 검증
        if mlp_output.dim() != cur_state.dim():
            raise ValueError(f"Dimension mismatch: mlp_output {mlp_output.shape} vs cur_state {cur_state.shape}")
        
        if cur_state.size(-1) != 12:
            raise ValueError(f"Expected 12 state features, got {cur_state.size(-1)}")
            
        if mlp_output.size(-1) != 7:
            raise ValueError(f"Expected 7 MLP outputs, got {mlp_output.size(-1)}")
        
        # 현재 상태 변수 추출 (차원 유지)
        V = cur_state[..., 0:1]     # 전압 (고정값)
        E = cur_state[..., 1:2]     # 외부 전해질 농도 (고정값)
        VF = cur_state[..., 2:3]    # Feed 부피
        VA = cur_state[..., 3:4]    # Acid 부피
        VB = cur_state[..., 4:5]    # Base 부피
        CFLA = cur_state[..., 5:6]  # Feed LA 농도
        CALA = cur_state[..., 6:7]  # Acid LA 농도
        CBLA = cur_state[..., 7:8]  # Base LA 농도
        CFK = cur_state[..., 8:9]   # Feed K 농도
        CAK = cur_state[..., 9:10]  # Acid K 농도
        CBK = cur_state[..., 10:11] # Base K 농도
        I = cur_state[..., 11:12]   # 전류

        # 물질량 계산 (농도 × 부피)
        NFLA = CFLA * VF; NALA = CALA * VA; NBLA = CBLA * VB
        NFK = CFK * VF; NAK = CAK * VA; NBK = CBK * VB

        # MLP 출력 (변화량)
        dVA = mlp_output[..., 0:1]    # Acid 부피 변화량
        dVB = mlp_output[..., 1:2]    # Base 부피 변화량
        dNALA = mlp_output[..., 2:3]  # Acid LA 물질량 변화량
        dNBLA = mlp_output[..., 3:4]  # Base LA 물질량 변화량
        dNAK = mlp_output[..., 4:5]   # Acid K 물질량 변화량
        dNBK = mlp_output[..., 5:6]   # Base K 물질량 변화량
        nI = mlp_output[..., 6:7]     # 새로운 전류값

        # 새로운 부피 계산 (물질 보존 법칙)
        nVF = VF - dVA - dVB  # Feed에서 빠져나간 부피
        nVA = VA + dVA       # Acid로 들어온 부피
        nVB = VB + dVB       # Base로 들어온 부피
        
        # 새로운 물질량 계산 (물질 보존 법칙)
        nNFLA = NFLA - dNALA - dNBLA  # Feed에서 빠져나간 LA
        nNALA = NALA + dNALA          # Acid로 들어온 LA
        nNBLA = NBLA + dNBLA          # Base로 들어온 LA
        nNFK = NFK - dNAK - dNBK      # Feed에서 빠져나간 K
        nNAK = NAK + dNAK             # Acid로 들어온 K
        nNBK = NBK + dNBK             # Base로 들어온 K
        
        # 부피 제약 조건 적용 (물리적으로 양수 유지)
        nVF = torch.clamp(nVF, min=self.eps)
        nVA = torch.clamp(nVA, min=self.eps)
        nVB = torch.clamp(nVB, min=self.eps)
        
        # 새로운 농도 계산 (농도 = 물질량 / 부피)
        nCFLA = torch.clamp(nNFLA / nVF, min=0)
        nCALA = torch.clamp(nNALA / nVA, min=0)
        nCBLA = torch.clamp(nNBLA / nVB, min=0)
        nCFK = torch.clamp(nNFK / nVF, min=0)
        nCAK = torch.clamp(nNAK / nVA, min=0)
        nCBK = torch.clamp(nNBK / nVB, min=0)
        
        # 전류는 양수 제약
        nI = torch.clamp(nI, min=0)

        # 새로운 상태 조립 (V, E는 고정값이므로 그대로 유지)
        new_state = torch.cat([
            V, E,  # 고정값: 전압, 외부 전해질 농도
            nVF, nVA, nVB,  # 새로운 부피
            nCFLA, nCALA, nCBLA,  # 새로운 LA 농도
            nCFK, nCAK, nCBK,     # 새로운 K 농도
            nI  # 새로운 전류
        ], dim=-1)
        
        return new_state

class BMEDSeq2SeqModel(nn.Module):
    def __init__(self, lstm_params, mlp_params):
        super().__init__()
        self.lstm_encoder = LSTMEncoder(**lstm_params)
        self.mlp_decoder = MLPDecoder(**mlp_params)
        self.mass_balance_layer = StateUpdateLayer()

    def forward(self, x, seq_len):
        # Teacher Forcing: 전체 시퀀스로 다음 상태들 예측
        
        # LSTM으로 시계열 패턴 학습
        lstm_out = self.lstm_encoder(x, seq_len)
        
        # MLP로 상태 변화량 예측
        mlp_out = self.mlp_decoder(lstm_out)
        
        # 물리적 제약 조건 적용하여 다음 상태 계산
        next_states = self.mass_balance_layer(mlp_out, x)
        
        return next_states

def masked_mse_loss(predictions, targets, seq_lengths):
    """
    개선된 마스킹된 MSE 손실 함수 - device 호환성 및 안정성 강화
    """
    # 입력 검증
    if predictions.shape != targets.shape:
        raise ValueError(f"Shape mismatch: predictions {predictions.shape} vs targets {targets.shape}")
    
    if predictions.size(0) != seq_lengths.size(0):
        raise ValueError(f"Batch size mismatch: predictions {predictions.size(0)} vs seq_lengths {seq_lengths.size(0)}")
    
    batch_size, max_len, features = predictions.shape
    
    # seq_lengths를 CPU로 이동하여 arange와 호환되도록 처리
    seq_lengths_cpu = seq_lengths.detach().cpu().long()
    
    # 시퀀스 길이가 0 이하인 경우 처리
    seq_lengths_cpu = torch.clamp(seq_lengths_cpu, min=1, max=max_len)
    
    # 마스크 생성: 실제 시퀀스 길이만큼만 True
    mask = torch.arange(max_len, device='cpu')[None, :] < seq_lengths_cpu[:, None]
    mask = mask.float().to(predictions.device)
    
    # 각 요소별 MSE 계산 (reduction='none')
    loss = F.mse_loss(predictions, targets, reduction='none')
    
    # 마스크 적용하여 패딩 부분 제거
    masked_loss_sum = (loss * mask.unsqueeze(-1)).sum()
    valid_elements = mask.sum() * features
    
    # 0으로 나누기 방지
    if valid_elements == 0:
        return torch.tensor(0.0, device=predictions.device, requires_grad=True)
    
    masked_loss = masked_loss_sum / valid_elements
    
    return masked_loss

def prepare_teacher_forcing_data(input_sequences, seq_lengths):
    """
    Teacher Forcing을 위한 입력-타겟 데이터 준비
    """
    # 입력: 마지막 시점 제외 [:-1]
    inputs = input_sequences[:, :-1, :]
    
    # 타겟: 첫 번째 시점 제외 [1:]  
    targets = input_sequences[:, 1:, :]
    
    # 타겟 시퀀스 길이는 1씩 감소 (마지막 시점 예측 불가)
    target_seq_lengths = torch.clamp(seq_lengths - 1, min=1)
    
    return inputs, targets, target_seq_lengths

def test_model_components():
    """모델의 각 컴포넌트를 개별적으로 테스트"""
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Testing on device: {device}")
    
    # 테스트 데이터 생성
    batch_size = 2
    seq_len = 5
    input_size = 12
    
    # 가짜 입력 데이터 생성
    test_input = torch.randn(batch_size, seq_len, input_size).to(device)
    test_seq_len = torch.tensor([4, 3]).to(device)  # 실제 시퀀스 길이
    
    print(f"Test input shape: {test_input.shape}")
    print(f"Test seq_len: {test_seq_len}")
    
    try:
        # 1. LSTM Encoder 테스트
        print("\n=== Testing LSTM Encoder ===")
        lstm_encoder = LSTMEncoder(input_size=12, hidden_size=64, num_layers=2).to(device)
        lstm_out = lstm_encoder(test_input, test_seq_len)
        print(f"LSTM output shape: {lstm_out.shape}")
        
        # 2. MLP Decoder 테스트
        print("\n=== Testing MLP Decoder ===")
        mlp_decoder = MLPDecoder(hidden_size=64, output_size=7, num_layers=2).to(device)
        mlp_out = mlp_decoder(lstm_out)
        print(f"MLP output shape: {mlp_out.shape}")
        
        # 3. StateUpdateLayer 테스트
        print("\n=== Testing StateUpdateLayer ===")
        state_update = StateUpdateLayer().to(device)
        # StateUpdateLayer는 MLP 출력과 같은 차원의 현재 상태가 필요
        new_state = state_update(mlp_out, test_input)
        print(f"Updated state shape: {new_state.shape}")
        
        # 4. 전체 모델 테스트
        print("\n=== Testing Full Model ===")
        model_params = {
            'lstm': {'input_size': 12, 'hidden_size': 64, 'num_layers': 2, 'dropout': 0.2},
            'mlp': {'hidden_size': 64, 'output_size': 7, 'num_layers': 2, 'num_nodes': 128, 'dropout': 0.3}
        }
        full_model = BMEDSeq2SeqModel(model_params['lstm'], model_params['mlp']).to(device)
        
        # Teacher forcing 데이터 준비
        inputs, targets, target_seq_lengths = prepare_teacher_forcing_data(test_input, test_seq_len)
        print(f"Teacher forcing - inputs: {inputs.shape}, targets: {targets.shape}")
        
        # Forward pass
        predictions = full_model(inputs, target_seq_lengths)
        print(f"Model predictions shape: {predictions.shape}")
        
        # Loss 계산 테스트
        loss = masked_mse_loss(predictions, targets, target_seq_lengths)
        print(f"Loss value: {loss.item()}")
        
        print("\n✅ All components work correctly!")
        return True
        
    except Exception as e:
        print(f"\n❌ Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🔧 BMED 모델 컴포넌트 테스트를 시작합니다...")
    print("="*60)
    
    test_result = test_model_components()
    
    print("="*60)
    if test_result:
        print("🎉 모든 컴포넌트가 정상적으로 작동합니다!")
        print("\n📋 다음 단계 권장사항:")
        print("1. Jupyter 노트북에서 전체 Optuna 하이퍼파라미터 최적화 실행")
        print("2. 최적의 파라미터로 최종 모델 학습")
        print("3. 결과 시각화 및 분석")
        print("\n💡 노트북의 마지막 셀들을 차례대로 실행하시면 됩니다.")
    else:
        print("⚠️ 컴포넌트 테스트에서 오류가 발견되었습니다.")
        print("문제를 해결한 후 다시 테스트해보세요.")