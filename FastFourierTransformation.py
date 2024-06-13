from cmath import pi, exp
from math import ceil, log2
from typing import List, Union


def fft(x: List[Union[int, complex]]) -> List[complex]:
    """
    Perform the Fast Fourier Transform (FFT) on a list of numbers.

    Parameters:
    x (List[Union[int, complex]]): A list of integers or complex numbers representing the input signal.

    Returns:
    List[complex]: A list of complex numbers representing the FFT of the input signal.
    """
    N = len(x)
    if N <= 1:
        return x

    # Divide the input into even and odd indices and recursively compute FFT for each
    even = fft(x[0::2])
    odd = fft(x[1::2])

    # Calculate the Twiddle factors for odd indices
    T = [exp(-2j * pi * k / N) * odd[k] for k in range(N // 2)]

    # Combine the results from even and odd indices using the FFT formula
    result = [even[k] + T[k] for k in range(N // 2)] + [even[k] - T[k] for k in range(N // 2)]

    # Return the result after rounding the complex numbers to the specified digits
    return round_complex(result)


# Function to round the real and imaginary parts of complex numbers in the result list
def round_complex(result_list, digits=4):
    return [complex(round(num.real, digits), round(num.imag, digits)) for num in result_list]


# Function to compute the next power of 2 greater than or equal to a given number
def next_power_of_2(n):
    return 1 if n == 0 else 2**ceil(log2(n))


# Function to pad the input list with zeros to make its length a power of 2
def if_pad_input(input_list):
    next_pow2 = next_power_of_2(len(input_list))
    return input_list + [0] * (next_pow2 - len(input_list))


# Function to get input from the user and pad it if needed
def get_input():
    input_str = input("Enter the values (comma separated): ")
    input_list = [int(x) for x in input_str.split(',')]

    return if_pad_input(input_list)


def main():
    input_list = get_input()
    result = fft(input_list)

    print(result)


if __name__ == "__main__":
    main()
