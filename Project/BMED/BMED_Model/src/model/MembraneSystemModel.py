from torch import nn
from .PhysicsLayer import PhysicsLayer

'''멤브레인 시스템 모델링을 위한 Physics-Informed feed forward model'''
class MembraneSystemModel(nn.Module):
    def __init__(self, hidden_nodes=64, hidden_layers=3):
        super(MembraneSystemModel, self).__init__()

        # Feed forward network
        layers = []
        
        # input layer
        layers.append(nn.Linear(10, hidden_nodes)) # [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
        layers.append(nn.ReLU())

        # hidden layer
        for _ in range(hidden_layers - 1):
            layers.append(nn.Linear(hidden_nodes, hidden_nodes))
            layers.append(nn.ReLU())

        # output layer
        layers.append(nn.Linear(hidden_nodes, 4)) # [dNLA, dNK, dVA, dVB]
        
        self.migration_predictor = nn.Sequential(*layers)

        # Physics layer
        self.physics_layer = PhysicsLayer()

    def forward(self, x): # hidden 상태는 사용하지 않음
        '''이동량 예측'''
        migrations = self.migration_predictor(x) #[batch_size, 4]

        '''상태 업데이트'''
        new_states = self.physics_layer(migrations, x) #[batch_size, 10]

        return migrations, new_states