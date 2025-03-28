#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    // Ensure we have the correct number of arguments
    if (argc < 2) {
        std::cerr << "Usage: cpp_program <rank>" << std::endl;
        return 1;
    }

    // Get the rank from the argument
    int rank = std::atoi(argv[1]);

    // Perform some computation (in this case, just an example with rank)
    int result = rank * rank;  // Example computation: square of the rank

    // Print the result (this will be captured by Python)
    std::cout << "Rank " << rank << " result: " << result << std::endl;

    return 0;
}
