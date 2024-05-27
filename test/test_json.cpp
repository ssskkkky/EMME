#include <fstream>
#include <iostream>
#include <sstream>

#include "JsonParser.h"

using namespace util::json;

int main() {
    const char* file_name = "./test/test.json";
    {
        std::ifstream ifs(file_name);
        JsonLexer json_lexer(ifs);
        std::cout << "Tokens:\n";
        while (json_lexer) { std::cout << json_lexer.get_token() << '\n'; }
        std::cout << '\n';
    }
    {
        auto obj = parse_file(file_name);
        double a = obj["a"];
        int b0 = obj["bs"][0];
        std::cout << a << ", " << b0 << '\n';

        try {
            double d = obj["obj"];
        } catch (std::exception& e) {
            std::cout << "Try to get number from object\n";
            std::cout << e.what() << '\n';
        }
        try {
            std::string str = obj["bs"];
        } catch (std::exception& e) {
            std::cout << "Try to get string from array\n";
            std::cout << e.what() << '\n';
        }
        std::cout << '\n';

        std::cout << std::boolalpha;
        std::cout << "Is the whole file an object? " << obj.is_object() << '\n';
        std::cout << "Is the 'bs' property an array? " << obj["bs"].is_array()
                  << "\n\n";

        std::cout << "Print type of every property of the outermost object:\n";
        for (auto& [key, val] : obj.as_object()) {
            std::cout << "    " << key << ": ["
                      << get_value_category_name(val.value_category()) << "]\n";
        }
        std::cout << '\n';

        std::cout << "Unformatted output: " << obj.dump() << "\n\n";
        std::cout << "Formatted output:\n" << obj.pretty_print() << "\n\n";
    }
    return 0;
}
