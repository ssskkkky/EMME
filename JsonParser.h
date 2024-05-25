#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <istream>
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace util {
namespace json {

enum class ValueCategory {
    Number,
    String,
    Array,
    Object,
};

// forward declarations
struct Object;
struct Array;
struct Number;
struct String;

/**
 * @brief Stores either a JSON object, array, number or string
 *
 */
struct Value {
    const ValueCategory type;
    std::unique_ptr<std::variant<Object, Array, Number, String>> ptr;
};

struct Object {
    std::unordered_map<std::string, Value> content;
};
struct Array {
    std::vector<Value> content;
};
struct Number {
    double content;
};
struct String {
    std::string content;
};

struct JsonLexer {
    enum class TokenName {
        STRING,
        NUMBER,
        BRACE_LEFT = '{',
        BRACE_RIGHT = '}',
        BRACKET_LEFT = '[',
        BRACKET_RIGHT = ']',
        COLON = ':',
        COMMA = ',',
        END_OF_FILE,
    };
    struct Token {
        TokenName name;
        std::string content;
    };

    JsonLexer(std::istream& is);

    Token get_token();

    operator bool() const;

    // any char that can be in a float number
    static bool is_digit(char c);
    static bool is_digit_start(char c);
    // tab, lf, cr or space
    static bool is_whitespace(char c);

   private:
    std::istream& is_;
};

#ifdef EMME_DEBUG
std::ostream& operator<<(std::ostream& os, const JsonLexer::Token& token);
#endif

struct JsonParser {};

}  // namespace json
}  // namespace util

#endif  // JSON_PARSER_H
