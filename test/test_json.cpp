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
        JsonLexer json_lexer(ifs);
        auto obj = JsonParser{}.parse(json_lexer);
        double a = obj["a"];
        std::cout << a << '\n';
    }
    ifs.close();
    return 0;
}
