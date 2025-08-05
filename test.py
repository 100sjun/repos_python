import numpy as np
import matplotlib.pyplot as plt

def fibonacci(n):
    if n <= 1:
        return n
    return fibonacci(n-1) + fibonacci(n-2)

def plot_fibonacci_sequence(n_terms):
    fib_numbers = [fibonacci(i) for i in range(n_terms)]
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(n_terms), fib_numbers, 'bo-', linewidth=2, markersize=6)
    plt.title(f'First {n_terms} Fibonacci Numbers')
    plt.xlabel('Index')
    plt.ylabel('Fibonacci Value')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return fib_numbers

if __name__ == "__main__":
    n = 15
    print(f"First {n} Fibonacci numbers:")
    fib_sequence = plot_fibonacci_sequence(n)
    print(fib_sequence)
    
    # Calculate golden ratio approximation
    golden_ratio = fib_sequence[-1] / fib_sequence[-2]
    print(f"Golden ratio approximation: {golden_ratio:.6f}")