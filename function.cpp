#include <string>

#include "function.h"
// Function to read contents of a file line by line
std::string readInputFromFile(const std::string& filename);

{
    std::ifstream input_file(filename);
    std::string file_contents;

    // Check if the file was opened successfully
    if (input_file.is_open()) {
        // Read the file contents line by line
        std::string line;
        while (std::getline(input_file, line)) {
            file_contents +=
                line + '\n';  // Add newline character after each line
        }
        input_file.close();
    } else {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }

    return file_contents;
}
