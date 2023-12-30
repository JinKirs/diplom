import matplotlib.pyplot as plt
import numpy as np
sigma = 0.5
sigma1 = 0.3
sigma2 = 0.7
x_a = 0
x_b = 2 * np.pi

def P(t, sigma):
    return np.sinh(2 * sigma * np.pi) / (sigma * sigma * np.pi * np.pi * np.cosh(2 * sigma * (t - np.pi)))

def p_nn():
    t = 0
    p_11 = []
    p_12 = []
    p_22 = []
    while t < 2*np.pi:
        t+= 0.01
        P_11 = P(t, sigma1)
        P_22 = P(t, sigma2)
        P_12 = P(t, (sigma1 + sigma2) / 2)
        p_11.append(P_11)
        p_22.append(P_22)
        p_12.append(P_12)
    return p_11, p_22, p_12

p_11, p_22, p_12 = p_nn()
p_11 = np.array(p_11)
p_22 = np.array(p_22)
p_12 = np.array(p_12)


def solve_quadratic_equation(a, b, c):
    """
    Решение квадратичного уравнения ax^2 + bx + c = 0 для массивов a, b, c.
    Возвращает массив решений.
    """
    discriminant = b**2 - 4*a*c
    sqrt_discriminant = np.sqrt(np.maximum(0, discriminant))  # Избегаем комплексных чисел при отрицательном дискриминанте

    l1 = (-b + sqrt_discriminant) / (2*a)
    l2 = (-b - sqrt_discriminant) / (2*a)

    return l1, l2

# Задаем коэффициенты уравнения
a = 1
b = -p_11-p_22
c = p_11*p_22-p_12**2

l1, l2 = solve_quadratic_equation(a, b, c)
print(l1)
print(l2)

#вывод графика l1
x1 = np.linspace(x_a, x_b, 629)
y1 = l1
fig, ax = plt.subplots(figsize = (12 / 10 * 8, 10 / 12 * 8))
plt.grid()
plt.plot(x1,y1)
plt.title("График λ1")
plt.xlabel("0,2π")
plt.ylabel("λ1")
plt.show()

#вывод графика l2
x2 = np.linspace(x_a, x_b, 629)
y2 = l2
fig, ax = plt.subplots(figsize = (12 / 10 * 8, 10 / 12 * 8))
plt.grid()
plt.plot(x2,y2)
plt.title("График λ2")
plt.xlabel("0,2π")
plt.ylabel("λ2")
plt.show()
