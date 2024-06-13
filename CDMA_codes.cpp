#include <iostream>
#include <vector>
#include <cmath>

using int_2d_matrix = std::vector<std::vector<int>>;

// Constants
constexpr int length_m = 4;
constexpr int voltage_1 = -1;
constexpr int voltage_0 = 1;
const int_2d_matrix h_2 = {{1, 1}, {1, -1}};

#define VAR_TO_STR(x) #x

void print2dVector(const int_2d_matrix& matrix)
{
    for(size_t row = 0; row < matrix.size(); row++)
    {
        for(size_t col = 0; col < matrix[0].size(); col++)
        {
            std::cout << matrix[row][col] << " ";
        }
        std::cout << std::endl;
    } 
    std::cout << std::endl;
}

// Function to calculate tensor product of two matrices
int_2d_matrix tensor_product(const int_2d_matrix& A, const int_2d_matrix& B) 
{
    const size_t m = A.size();
    const size_t n = A[0].size();
    const size_t p = B.size();
    const size_t q = B[0].size();

    int_2d_matrix kronecker(m * p, std::vector<int>(n * q));

    for(size_t i = 0; i < m; i++) 
    {
        for(size_t j = 0; j < n; j++) 
        {
            for(size_t ii = 0; ii < p; ii++) 
            {
                for(size_t jj = 0; jj < q; jj++) 
                {
                    kronecker[i * p + ii][j * q + jj] = A[i][j] * B[ii][jj];
                }
            }
        }
    }

    return kronecker;
}

// Function to calculate Walsh codes
int_2d_matrix walsh_of(size_t k) 
{
    if(k <= 2) 
    {
        return h_2;
    } 
    
    auto curr_h_2 = walsh_of(k - 1);
    return tensor_product(curr_h_2, walsh_of(2));
}

// Function to generate spreading codes
int_2d_matrix channel_sequence_of(int number_of_channels) 
{
    return walsh_of(size_t(std::log2(number_of_channels) + 1));
}

const std::vector<int> convertStringToInts(const std::string& input) 
{
    std::vector<int> result;

    for(char ch : input) 
    {
        result.emplace_back(ch == '0' ? voltage_0 : voltage_1);
    }

    return result;
}

const std::vector<int> get_data_input(void) 
{
    std::string data;

    std::cout << "Enter data: ";
    std::cin >> data;

    std::vector<int> result_data = convertStringToInts(data);

    return result_data; 
}

// Function to perform scalar sum of codes
const std::vector<int> scalar_sum_of_codes(const int_2d_matrix& spreading_codes, const std::vector<int>& data) 
{
    int_2d_matrix mult_data_sp_code(length_m, std::vector<int>(length_m, 0));

    std::cout << VAR_TO_STR(spreading_codes) << ": " << std::endl;
    print2dVector(spreading_codes);

    // Multiplying the spreading code with the current data
    for(size_t row = 0; row < length_m; row++)
    {
        for(size_t col = 0; col < length_m; col++) 
        {
            mult_data_sp_code[row][col] = spreading_codes[row][col] * data[row];
        }
    }

    std::cout << VAR_TO_STR(mult_data_sp_code) << ": " << std::endl;
    print2dVector(mult_data_sp_code);

    std::vector<int> summed_voltages(length_m, 0);

    // Summin g the voltages by columns
    for(size_t row = 0; row < length_m; row++)
    {
        for(size_t col = 0; col < length_m; col++)
        {
            summed_voltages[row] += mult_data_sp_code[col][row];
        }
    }

    std::cout << VAR_TO_STR(summed_voltages) << ": " << std::endl;
    for(int vol : summed_voltages)
    {
        std::cout << vol << " ";
    }
    std::cout << std::endl;

    return summed_voltages;
}

int decode_signal(const std::vector<int>& combined_signal, const std::vector<int>& spreading_code) 
{
    int inner_product = 0;

    for(size_t i = 0; i < length_m; i++) 
    {
        inner_product += combined_signal[i] * spreading_code[i];
    }

    return (inner_product / length_m) < 0 ? 1 : 0;
}

int main() 
{
    int_2d_matrix spreading_codes = channel_sequence_of(length_m);
    
    std::vector<int> combined_signal = scalar_sum_of_codes(spreading_codes, get_data_input());
    std::vector<int> decoded_users;

    // Demonstrate the protocol if all messages are being send at the same time
    for(const auto& spreading_code : spreading_codes) 
    {
        int decoded_user = decode_signal(combined_signal, spreading_code);
        decoded_users.emplace_back(decoded_user);
    }

    for(size_t i = 0; i < decoded_users.size(); i++) 
    {
        std::cout << "Decoded User " << i + 1 << ": " << decoded_users[i] << std::endl;
    }

    return 0;
}
