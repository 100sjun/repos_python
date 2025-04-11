import torch
import torch.nn as nn
'''
BMED 사이의 ion migration을 기반으로 하여 각 time step별로 변하는 농도/부피를 계산
'''

class PhysicsLayer(nn.Module):
    def __init__(self):
        super(PhysicsLayer, self).__init__()
        self.time_step = 0.1
    
    def forward(self, migrations, features):
        '''현재 시간 단계 상태 추출'''
        T = features[..., 0:1] # 현재 상태 온도 [C]
        V = features[..., 1:2] # 현재 상태 전압 [V]
        E = features[..., 2:3] # 현재 상태 전해질 농도 [mol/L]
        CFLA = features[..., 3:4] # 현재 상태 Feed LA 농도 [mol/L]
        CFK = features[..., 4:5] # 현재 상태 Feed K+ 농도 [mol/L]
        CALA = features[..., 5:6] # 현재 상태 Acid LA 농도 [mol/L]
        CBK = features[..., 6:7] # 현재 상태 Base K+ 농도 [mol/L]
        VF = features[..., 7:8] # 현재 상태 Feed 부피 [L]
        VA = features[..., 8:9] # 현재 상태 Acid 부피 [L]
        VB = features[..., 9:10] # 현재 상태 Base 부피 [L]

        '''현재 시간 단계 몰 계산'''
        NFLA = CFLA * VF # 현재 상태 Feed LA 몰수 [mol]
        NFK = CFK * VF # 현재 상태 Feed K+ 몰수 [mol]
        NALA = CALA * VA # 현재 상태 Acid LA 몰수 [mol]
        NBK = CBK * VB # 현재 상태 Base K+ 몰수 [mol]

        '''변화량 상태 추출'''
        JLA = migrations[..., 0:1] # 시간당 LA 몰 변화량 [mol/h]
        JK = migrations[..., 1:2] # 시간당 K+ 몰 변화량 [mol/h]
        JVA = migrations[..., 2:3] # 시간당 물 부피 변화량 (Acid) [L/h]
        JVB = migrations[..., 3:4] # 시간당 물 부피 변화량 (Base) [L/h]

        '''time step에서의 변화량 계산'''
        dLA = JLA * self.time_step # LA 몰 변화량 [mol]
        dK = JK * self.time_step # K+ 몰 변화량[mol]
        dVA = JVA * self.time_step # 물 부피 변화량 (Acid) [L]
        dVB = JVB * self.time_step # 물 부피 변화량 (Base) [L]

        '''부피 업데이트'''
        nVF = VF - dVA - dVB # 다음 상태 Feed 부피 [L]
        nVA = VA + dVA # 다음 상태 Acid 부피 [L]
        nVB = VB + dVB # 다음 상태 Base 부피 [L]

        '''몰수 업데이트'''
        nNFLA = NFLA - dLA # 다음 상태 Feed LA 몰수 [mol]
        nNFK = NFK - dK # 다음 상태 Feed K+ 몰수 [mol]
        nNALA = NALA + dLA # 다음 상태 Acid LA 몰수 [mol]
        nNBK = NBK + dK # 다음 상태 Base K+ 몰수 [mol]

        '''농도 업데이트'''
        eps = 1e-6 # 0으로 나누기 방지
        nCFLA = nNFLA / (nVF + eps) # 다음 상태 Feed LA 농도 [mol/L]
        nCFK = nNFK / (nVF + eps) # 다음 상태 Feed K+ 농도 [mol/L] 
        nCALA = nNALA / (nVA + eps) # 다음 상태 Acid LA 농도 [mol/L]
        nCBK = nNBK / (nVB + eps) # 다음 상태 Base K+ 농도 [mol/L]

        '''결과 출력'''
        new_states = torch.cat([
            T, V, E, nCFLA, nCFK, nCALA, nCBK, nVF, nVA, nVB
        ], dim=-1)

        return new_states
        
        
        
        
        
        
        
        
        
        


