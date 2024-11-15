
class CustomModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden = nn.Sequential(
            nn.Linear(5,40),
            nn.ReLU(),
            nn.Linear(40,50),
            nn.ReLU(),
            nn.Linear(50,60),
            nn.ReLU(),
            nn.Linear(60,80),
            nn.ReLU(),
            nn.Linear(80,80),
            nn.ReLU(),
            nn.Linear(80,2)
        )

    def forward(self, x):
        return self.hidden(x)
