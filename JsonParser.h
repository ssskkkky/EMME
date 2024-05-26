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
    NumberInt,
    NumberFloat,
    String,
    Array,
    Object,
};

const char* get_value_category_name(ValueCategory);

struct NumberInt {
    int content;
};
struct NumberFloat {
    double content;
};

/**
 * @brief Stores either a JSON object, array, number or string
 *
 */
struct Value {
    using object_container_type = std::unordered_map<std::string, Value>;
    using array_container_type = std::vector<Value>;
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
    // operator int() const;
    operator std::string() const;

    template <typename T>
    T as_number() const requires std::is_arithmetic_v<T> {
        if (value_cat == ValueCategory::NumberFloat) {
            return static_cast<NumberFloat*>(ptr.get())->content;
        } else if (value_cat == ValueCategory::NumberInt) {
            return static_cast<NumberInt*>(ptr.get())->content;
        } else {
            throw std::runtime_error(
                std::string("Incorrect JSON type, require: Number, actually") +
                std::string(get_value_category_name(value_cat)));
        }
        return T{};
    }

    const Value& operator[](const std::string&) const;
    const Value& operator[](int) const;

    const object_container_type& as_object() const;
    const array_container_type& as_array() const;

    std::size_t size() const;
    std::size_t empty() const;

    // unformatted output
    std::string dump() const;

    // formatted output
    std::string pretty_print(std::size_t = 0) const;

   private:
    std::unique_ptr<void, std::function<void(void*)>> ptr;
    ValueCategory value_cat;

    template <typename... Ts>
    void expected_cat(Ts... cats) const
        requires(std::same_as<Ts, ValueCategory>&&...) {
        if (((value_cat != cats) && ...)) {
            std::ostringstream oss;
            if constexpr (sizeof...(Ts) == 1) {
                oss << "Incorrect JSON type, requires: ";
            } else {
                oss << "Incorrect JSON type, requires one of: ";
            }
            ((oss << get_value_category_name(cats) << ", "), ...)
                << "actually: " << get_value_category_name(value_cat);
            throw std::runtime_error(oss.str());
        }
    }

    static void print_space(std::ostream&, std::size_t);
};

struct Object {
    Value::object_container_type content;
};
struct Array {
    Value::array_container_type content;
};
struct String {
    std::string content;
};

struct JsonLexer {
    enum class TokenName {
        END_OF_FILE,
        STRING,
        INTEGER,
        FLOAT,
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
    Token buffer{};
    bool is_buffer_full{};
    bool is_buffer_output{};

    void read_token_to_buffer();
    // any char that can be in a float number
    static bool is_digit(char c);
    static bool is_digit_start(char c);
    // tab, lf, cr or space
    static bool is_whitespace(char c);
};

#ifdef EMME_DEBUG
std::ostream& operator<<(std::ostream& os, const JsonLexer::Token& token);
#endif

struct JsonParser {
    JsonParser(JsonLexer&&);
    Value parse();

   private:
    JsonLexer lexer;

    Value parse_value();
    Value parse_string(const JsonLexer::Token&);
    Value parse_int(const JsonLexer::Token&);
    Value parse_float(const JsonLexer::Token&);
    Value parse_object();
    Value parse_array();

    JsonLexer::Token try_get_and_check(JsonLexer::TokenName);
    JsonLexer::Token try_get_from_lexer(bool = false);
    JsonLexer::Token try_peek_from_lexer();
    void report_syntax_error(const JsonLexer::Token& = {});
};

Value parse(std::istream&);
Value parse_file(std::string);

}  // namespace json
}  // namespace util

#endif  // JSON_PARSER_H
