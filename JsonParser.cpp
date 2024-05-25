#include <cstdlib>  // atof
#include <ostream>
#include <sstream>

#include "JsonParser.h"

namespace util {
namespace json {

const char* get_value_category_name(ValueCategory val_cat) {
#define PROCESS_CAT_NAME(p) \
    case (p):               \
        return #p;          \
        break;
    switch (val_cat) {
        PROCESS_CAT_NAME(ValueCategory::Number)
        PROCESS_CAT_NAME(ValueCategory::String)
        PROCESS_CAT_NAME(ValueCategory::Array)
        PROCESS_CAT_NAME(ValueCategory::Object)
    }
#undef PROCESS_CAT_NAME
    return "";  // unreachable
}

Value::operator double() const {
    expected_cat(ValueCategory::Number);
    return static_cast<Number*>(ptr.get())->content;
}
Value::operator std::string() const {
    expected_cat(ValueCategory::String);
    return static_cast<String*>(ptr.get())->content;
}
const Value& Value::operator[](const std::string& key) {
    expected_cat(ValueCategory::Object);
    return static_cast<Object*>(ptr.get())->content[key];
}
const Value& Value::operator[](int idx) {
    expected_cat(ValueCategory::Array);
    return static_cast<Array*>(ptr.get())->content.at(idx);
}

void Value::expected_cat(ValueCategory cat) const {
    if (value_cat != cat) {
        std::ostringstream oss;
        oss << "Incorrect JSON type, require: " << get_value_category_name(cat)
            << ", actually: " << get_value_category_name(value_cat);
        throw std::runtime_error(oss.str());
    }
}

JsonLexer::JsonLexer(std::istream& is) : is_(is) {}

JsonLexer::Token JsonLexer::get_token() {
    if (!is_buffer_full || is_buffer_output) { read_token_to_buffer(); }
    is_buffer_output = true;
    return buffer;
}

JsonLexer::Token JsonLexer::peek_token() {
    if (!is_buffer_full || is_buffer_output) { read_token_to_buffer(); }
    is_buffer_output = false;
    return buffer;
}

void JsonLexer::read_token_to_buffer() {
    char c;
    // skip all whitespaces
    while (is_.get(c) && is_whitespace(c)) {}
    buffer.content.clear();
    if (!is_) {
        buffer.name = TokenName::END_OF_FILE;
    } else if (c == static_cast<char>(TokenName::BRACE_LEFT) ||
               c == static_cast<char>(TokenName::BRACE_RIGHT) ||
               c == static_cast<char>(TokenName::BRACKET_LEFT) ||
               c == static_cast<char>(TokenName::BRACKET_RIGHT) ||
               c == static_cast<char>(TokenName::COLON) ||
               c == static_cast<char>(TokenName::COMMA)) {
        buffer.name = static_cast<TokenName>(c);
    } else if (c == '"') {
        // a string
        buffer.name = TokenName::STRING;
        while (is_.get(c) && c != '"') { buffer.content.push_back(c); }
    } else if (is_digit_start(c)) {
        // a number
        buffer.name = TokenName::NUMBER;
        do { buffer.content.push_back(c); } while (is_.get(c) && is_digit(c));
        is_.unget();
    }
    is_buffer_full = true;
}
JsonLexer::operator bool() const {
    return !is_.fail();
}
bool JsonLexer::is_digit(char c) {
    return c == '.' || c == 'E' || c == 'e' || is_digit_start(c);
}
bool JsonLexer::is_digit_start(char c) {
    return (c >= '0' && c <= '9') || c == '-' || c == '+';
}
bool JsonLexer::is_whitespace(char c) {
    return c == '\t' || c == '\n' || c == '\r' || c == ' ';
}

std::ostream& operator<<(std::ostream& os, const JsonLexer::Token& token) {
    return os << "{ Name: " << ([&token] {
#define PROCESS_TOKEN_NAME(p) \
    case (p):                 \
        return #p;            \
        break;
               switch (token.name) {
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::STRING)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::NUMBER)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACE_LEFT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACE_RIGHT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACKET_LEFT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACKET_RIGHT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::COLON)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::COMMA)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::END_OF_FILE)
               }
#undef PROCESS_TOKEN_NAME
               return "(No such name)";
           })()
              << ", Content: " << token.content << " }";
}

Value JsonParser::parse(JsonLexer& lexer) {
    if (!lexer) { report_syntax_error(); }  // token list is empty
    auto value = parse_value(lexer);
    try_get_from_lexer(lexer, true);  // Check if lexer ends
    return value;
}

Value JsonParser::parse_value(JsonLexer& lexer) {
    auto token = try_get_from_lexer(lexer);
    switch (token.name) {
        case JsonLexer::TokenName::STRING:
            return parse_string(token);
        case JsonLexer::TokenName::NUMBER:
            return parse_number(token);
        case JsonLexer::TokenName::BRACE_LEFT:
            return parse_object(lexer);
        case JsonLexer::TokenName::BRACKET_LEFT:
            return parse_array(lexer);
        default:
            report_syntax_error(token);
            return Value{};  // unreachable
    }
}

Value JsonParser::parse_string(const JsonLexer::Token& token) {
    if (token.name != JsonLexer::TokenName::STRING) {
        report_syntax_error(token);
    }
    return {ValueCategory::String, new String{token.content}};
}

Value JsonParser::parse_number(const JsonLexer::Token& token) {
    if (token.name != JsonLexer::TokenName::NUMBER) {
        report_syntax_error(token);
    }
    return {ValueCategory::Number,
            new Number{std::atof(token.content.c_str())}};
}

Value JsonParser::parse_object(JsonLexer& lexer) {
    auto obj = new Object;
    auto token = try_peek_from_lexer(lexer);
    // empty object
    if (token.name == JsonLexer::TokenName::BRACE_RIGHT) {
        return {ValueCategory::Object, obj};
    }
    while (true) {
        auto key = try_get_and_check(lexer, JsonLexer::TokenName::STRING);
        try_get_and_check(lexer, JsonLexer::TokenName::COLON);
        obj->content.emplace(key.content, parse_value(lexer));  // value
        token = try_get_from_lexer(lexer);
        if (token.name == JsonLexer::TokenName::BRACE_RIGHT) { break; }
        if (token.name != JsonLexer::TokenName::COMMA) {
            report_syntax_error(token);
        }
    }
    return {ValueCategory::Object, obj};
}

Value JsonParser::parse_array(JsonLexer& lexer) {
    auto arr = new Array;
    auto token = try_peek_from_lexer(lexer);
    // empty array
    if (token.name == JsonLexer::TokenName::BRACKET_RIGHT) {
        return {ValueCategory::Array, arr};
    }
    while (true) {
        arr->content.emplace_back(parse_value(lexer));
        token = try_get_from_lexer(lexer);
        if (token.name == JsonLexer::TokenName::BRACKET_RIGHT) { break; }
        if (token.name != JsonLexer::TokenName::COMMA) {
            report_syntax_error(token);
        }
    }
    return {ValueCategory::Array, arr};
}

JsonLexer::Token JsonParser::try_get_and_check(
    JsonLexer& lexer,
    JsonLexer::TokenName expected_token_name) {
    auto token = lexer.get_token();
    if (token.name != expected_token_name) { report_syntax_error(token); }
    return token;
}

JsonLexer::Token JsonParser::try_get_from_lexer(JsonLexer& lexer,
                                                bool end_expected) {
    auto token = lexer.get_token();
    if (end_expected ^ (token.name == JsonLexer::TokenName::END_OF_FILE)) {
        report_syntax_error(token);
    }
    return token;
}

JsonLexer::Token JsonParser::try_peek_from_lexer(JsonLexer& lexer) {
    auto token = lexer.peek_token();
    if (token.name == JsonLexer::TokenName::END_OF_FILE) {
        report_syntax_error(token);
    }
    return token;
}

void JsonParser::report_syntax_error(const JsonLexer::Token& token) {
    std::ostringstream oss;
    oss << "Syntax error in JSON file at token ";
    oss << token;
    throw std::runtime_error(oss.str());
}

}  // namespace json
}  // namespace util
