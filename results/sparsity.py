import numpy as np
import matplotlib.pyplot as plt

def read_matrix(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        matrix = []
        for line in lines:
            row = [float(num) for num in line.split()]
            matrix.append(row)
    return np.array(matrix)

def plot_sparsity_pattern(matrix):
    row, col = np.nonzero(matrix)
    plt.figure(figsize=(10, 10))
    plt.scatter(col, row, s=1) 
    plt.gca().invert_yaxis()  
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Sparsity Pattern of the Kv')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.show()


matrix = read_matrix('result_Kv.txt')
plot_sparsity_pattern(matrix)
