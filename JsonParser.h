#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <functional>
#include <istream>
#include <memory>  // unique_ptr
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#ifdef EMME_DEBUG
#include <iostream>
#endif

namespace util {
namespace json {

enum class ValueCategory {
    Number,
    String,
    Array,
    Object,
};

static const char* get_value_category_name(ValueCategory);

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
    Value() = default;

    // Type information in deleter is stored during contruction
    template <typename T>
    Value(ValueCategory cat, T* raw_ptr)
        : ptr(raw_ptr,
              [cat](const void* data) {
                  delete static_cast<const T*>(data);
#ifdef EMME_DEBUG
                  std::cout << get_value_category_name(cat) << " deleted.\n";
#endif
              }),
          value_cat(cat) {
    }

    operator double() const;
    operator std::string() const;

    const Value& operator[](const std::string&);
    const Value& operator[](int);

   private:
    std::unique_ptr<void, std::function<void(void*)>> ptr;
    ValueCategory value_cat;

    void expected_cat(ValueCategory) const;
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
        END_OF_FILE,
        STRING,
        NUMBER,
        BRACE_LEFT = '{',
        BRACE_RIGHT = '}',
        BRACKET_LEFT = '[',
        BRACKET_RIGHT = ']',
        COLON = ':',
        COMMA = ',',
    };
    struct Token {
        TokenName name;
        std::string content;
        // TODO: Add position in case of syntax error
    };

    JsonLexer(std::istream&);

    Token get_token();
    Token peek_token();

    operator bool() const;

   private:
    std::istream& is_;
    Token buffer;
    bool is_buffer_full;
    bool is_buffer_output;

    void read_token_to_buffer();
    // any char that can be in a float number
    static bool is_digit(char c);
    static bool is_digit_start(char c);
    // tab, lf, cr or space
    static bool is_whitespace(char c);
};

std::ostream& operator<<(std::ostream& os, const JsonLexer::Token& token);

struct JsonParser {
    JsonParser(JsonLexer&&);
    Value parse();

   private:
    JsonLexer lexer;

    Value parse_value();
    Value parse_string(const JsonLexer::Token&);
    Value parse_number(const JsonLexer::Token&);
    Value parse_object();
    Value parse_array();

    JsonLexer::Token try_get_and_check(JsonLexer::TokenName);
    JsonLexer::Token try_get_from_lexer(bool = false);
    JsonLexer::Token try_peek_from_lexer();
    void report_syntax_error(const JsonLexer::Token& = {});
};

}  // namespace json
}  // namespace util

#endif  // JSON_PARSER_H
