import matplotlib.pyplot as plt
from math import sqrt
import numpy as np


def fourier_coefficients(func, T, N):
    """
    Calculate Fourier coefficients for a given periodic function.

    Parameters:
        func: Function representing the periodic signal.
        T: Period of the signal.
        N: Number of coefficients to calculate.

    Returns:
        Coefficients of the Fourier series.
    """
    coefficients = []

    a0 = (1 / T) * integrate(func, 0, T)
    coefficients.append(a0)

    for n in range(1, N+1):
        an = (2 / T) * integrate(lambda t: func(t) * np.cos(2 * np.pi * n * t / T), 0, T)
        bn = (2 / T) * integrate(lambda t: func(t) * np.sin(2 * np.pi * n * t / T), 0, T)

        print(f"A{n}: {round(an, 3)} \tB{n}: {round(bn, 3)}")

        coefficient = (sqrt(an**2 + bn**2))

        if coefficient > 10e-6:
            coefficients.append((sqrt(an**2 + bn**2)))
        else:
            coefficients.append(0)

    return coefficients


def integrate(func, a, b, N=1_000_000):
    """
    Numerical integration of a function using trapezoidal rule.

    Parameters:
        func: Function to integrate.
        a: Lower limit of integration.
        b: Upper limit of integration.
        N: Number of intervals for integration.

    Returns:
        Approximation of the integral.
    """
    dx = (b - a) / N
    x = np.linspace(a, b, N+1)
    y = func(x)

    return dx * (np.sum(y) - 0.5 * (y[0] + y[-1]))


# Function to create the wave
def rectangle_wave(t, t_u_period, A):
    return np.where(np.mod(t, T) < t_u_period, A, 0)


# Plots the waveform of the given function over given period.
def plot_waveform(func, T, number_of_impulses, title, num_points=1000):

    t = np.linspace(0, T * number_of_impulses, num_points)
    y = func(t)
    plt.plot(t, y)
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()


# Defining the constants needed for the forming of the rectangle wave
A = 1               # Amplitude of the wave
T = 90e-6           # Period of the wave
N = 7               # Number of coefficients to be calculated
tu_period = T/3     # Impulse period

coefficients = fourier_coefficients(lambda t: rectangle_wave(t, tu_period, A), T, N)

print("\nFourier coefficients:")

for n, (An) in enumerate(coefficients[::], start=1):
    print(f"C{n - 1}: {round(An, 3)}")

# Plot the rectangle wave
plot_waveform(lambda t: rectangle_wave(t, tu_period, A), T, 3, 'Rectangle Wave')

# Plot the coefficients
plt.bar(range(len(coefficients)), coefficients, width=0.2)
plt.xlabel('n')
plt.ylabel('C_n$')
plt.title('Fourier Series Coefficients $C_n$')
plt.grid(True)
plt.show()

