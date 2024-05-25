#include <fstream>
#include <iostream>
#include <sstream>

#include "JsonParser.h"

using namespace util::json;

int main() {
    std::ifstream ifs("./test.json");
    JsonLexer json_lexer(ifs);
    while (json_lexer) { std::cout << json_lexer.get_token() << '\n'; }
    ifs.close();
    return 0;
}
