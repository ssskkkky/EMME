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
        bool bb = obj["primitives"][0].as_boolean();
        std::cout << "a = " << a << "\nbs[0] = " << b0
                  << "\nprimitives[0] = " << std::boolalpha << bb << '\n';

        try {
            double d = obj["abc"];
            static_cast<void>(d);
        } catch (std::exception& e) {
            std::cout << "Try to get an undefined property as number\n";
            std::cout << e.what() << '\n';
        }
        try {
            double d = obj["obj"];
            static_cast<void>(d);
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

        std::cout << "Property a < 42 ? " << (obj["a"] < 42) << '\n';
        obj["a"] += 1;
        std::cout << "Property a after adding 1: " << obj["a"].dump() << '\n';
        obj["a"] = 69;
        std::cout << "Property a after being assigned 69: " << obj["a"].dump()
                  << '\n';
        obj["obj"] = 69.69;
        std::cout << "Property obj after being assigned 69.69: "
                  << obj["obj"].dump() << '\n';
        std::cout << '\n';

        std::cout << "The destruct sequence:\n";
    }
    {  // test syntax error detection
        char json_failure_1[] = "{\"a\":1,,\"b\":2}";
        try {
            parse(json_failure_1);
        } catch (std::exception& e) {
            std::cout << "Try to parse a problemetic json:\n"
                      << json_failure_1 << '\n';
            std::cout << e.what() << '\n';
        }
        char json_failure_2[] = "{\"a\":1,\"b\":2,c:3}";
        try {
            parse(json_failure_2);
        } catch (std::exception& e) {
            std::cout << "Try to parse a problemetic json:\n"
                      << json_failure_2 << '\n';
            std::cout << e.what() << '\n';
        }
    }
    return 0;
}
