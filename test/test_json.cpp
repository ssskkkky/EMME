#include <fstream>
#include <iostream>
#include <sstream>

#include "JsonParser.h"

using namespace util::json;

int main() {
    std::ifstream ifs("./test/test.json");
    {
        JsonLexer json_lexer(ifs);
        std::cout << "Tokens:\n";
        while (json_lexer) { std::cout << json_lexer.get_token() << '\n'; }
        std::cout << '\n';
    }
    ifs.clear();
    ifs.seekg(0);  // rewind
    {
        auto obj = JsonParser{JsonParser{ifs}}.parse();
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
    }
    ifs.close();
    return 0;
}
