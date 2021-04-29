import matplotlib.pyplot as plt
import numpy as np
import copy
d_x = 0.125
d_time = 0.001
D = 4
L = 8
X = np.arange(-L, L, d_x)

def func():
    N = [10, 100, 1000, 10000]
    for i in range(len(N)):
        prob = [(1 if i == 0 else 0) for i in X]
        for _ in range(N[i]):
            prob_temp = [0 for i in range(len(prob))]
            for j in range(1, len(prob) - 1):
                prob_temp[j] = prob[j] + (D*d_time/(d_x**2)) * (prob[j+ 1] - 2*prob[j] + prob[j-1])
            prob = prob_temp

        plt.plot(X, prob , label = f"N = {N[i]}")
    plt.xlabel("X")
    plt.ylabel("Probability")
    plt.legend()
    plt.show()

def main():
    func()
    d_x = 0.125
    d_y = 0.125
    d_time = 0.001
    D_x = 3
    D_y = 2
    L = 3
    X = np.arange(-L, L, d_x)
    Y = np.arange(-L, L, d_y)
    prob = np.zeros((len(X), len(Y)))
    N = [10, 100, 1000, 10000]
    # print("In progress: ")
    for k in range(len(N)):
        print("N = ", N[k])
        
        for x in range(len(X)):
            for y in range(len(Y)):
                prob[x][y] = 1 if (X[x] == 0 and Y[y] == 0) else 0
                
        for z in range(N[k]):
            prob_temp = np.zeros(prob.shape)
            for i in range(1, len(prob)-1):
                for j in range(1, len(prob[i])-1):
                    xx = D_x / (d_x**2) * (prob[i+1][j] - 2*prob[i][j] + prob[i-1][j])
                    yy = D_y / (d_y**2) * (prob[i][j+1] - 2*prob[i][j] + prob[i][j-1])
                    prob_temp[i][j] = prob[i][j] + d_time* (xx + yy)
            prob = prob_temp
        plt.xlabel("x co-ordinate")
        plt.ylabel("y co-ordinate")
        plt.contour(X, Y, prob, antialiased=False)
        plt.show()
    print()
        
main()