#include "JsonParser.h"

#include <cstdlib>  // atof, atoi
#include <cstring>  // strcmp, strlen
#include <fstream>  // ifstream
#include <sstream>  // ostringstream

namespace util {
namespace json {

const char* get_value_category_name(ValueCategory val_cat) {
#define PROCESS_CAT_NAME(p) \
    case (p):               \
        return #p;          \
        break;
    switch (val_cat) {
        PROCESS_CAT_NAME(ValueCategory::Null)
        PROCESS_CAT_NAME(ValueCategory::NumberInt)
        PROCESS_CAT_NAME(ValueCategory::NumberFloat)
        PROCESS_CAT_NAME(ValueCategory::Boolean)
        PROCESS_CAT_NAME(ValueCategory::String)
        PROCESS_CAT_NAME(ValueCategory::Array)
        PROCESS_CAT_NAME(ValueCategory::Object)
        PROCESS_CAT_NAME(ValueCategory::TypedArrayComplexFloat)
        PROCESS_CAT_NAME(ValueCategory::TypedArrayComplexDouble)
    }
#undef PROCESS_CAT_NAME
    return "";  // unreachable
}

Value::Value() : value_cat{} {}

Value::Value(const Value& other) : value_cat{} {
    // value_cat is initialized as Null to prevent accessing non allocated
    // address
    *this = other;
}

Value& Value::operator=(const Value& other) {
    other.expected_cat(ValueCategory::Boolean, ValueCategory::NumberFloat,
                       ValueCategory::NumberInt, ValueCategory::String);
    // modify content directly with "plain value"
    if (value_cat == other.value_cat) {
#define PROCESS_TYPE(t)                                      \
    case (ValueCategory::t):                                 \
        static_cast<t*>(ptr.get())->content =                \
            static_cast<const t*>(other.ptr.get())->content; \
        break;
        switch (value_cat) {
            PROCESS_TYPE(Boolean)
            PROCESS_TYPE(NumberFloat)
            PROCESS_TYPE(NumberInt)
            PROCESS_TYPE(String)
            default:
                break;
        }
#undef PROCESS_TYPE
    } else {
        *this = other.clone();
    }
    return *this;
}

Value::operator double() const {
    expected_cat(ValueCategory::NumberFloat, ValueCategory::NumberInt);
    if (value_cat == ValueCategory::NumberFloat) {
        return static_cast<NumberFloat*>(ptr.get())->content;
    } else {
        return static_cast<NumberInt*>(ptr.get())->content;
    }
}

bool Value::as_boolean() const {
    return static_cast<Boolean*>(ptr.get())->content;
};

std::string& Value::as_string() {
    expected_cat(ValueCategory::String);
    return static_cast<String*>(ptr.get())->content;
}
const std::string& Value::as_string() const {
    expected_cat(ValueCategory::String);
    return static_cast<String*>(ptr.get())->content;
}

Value::operator std::string() const {
    expected_cat(ValueCategory::String);
    return static_cast<String*>(ptr.get())->content;
}

Value& Value::operator[](const std::string& key) {
    return as_object()[key];
}
Value& Value::operator[](const char* key) {
    return as_object()[key];
}

const Value& Value::at(const std::string& key) const {
    return as_object().at(key);
}
const Value& Value::at(std::size_t idx) const {
    return as_array().at(idx);
}

Value& Value::at(const std::string& key) {
    return as_object().at(key);
}
Value& Value::at(std::size_t idx) {
    return as_array().at(idx);
}

const Value::object_container_type& Value::as_object() const {
    expected_cat(ValueCategory::Object);
    return static_cast<Object*>(ptr.get())->content;
}
const Value::array_container_type& Value::as_array() const {
    expected_cat(ValueCategory::Array);
    return static_cast<Array*>(ptr.get())->content;
}

Value::object_container_type& Value::as_object() {
    expected_cat(ValueCategory::Object);
    return static_cast<Object*>(ptr.get())->content;
}
Value::array_container_type& Value::as_array() {
    expected_cat(ValueCategory::Array);
    return static_cast<Array*>(ptr.get())->content;
}

bool Value::is_object() const {
    return value_cat == ValueCategory::Object;
}
bool Value::is_array() const {
    return value_cat == ValueCategory::Array;
}
bool Value::is_number() const {
    return value_cat == ValueCategory::NumberInt ||
           value_cat == ValueCategory::NumberFloat;
};
bool Value::is_string() const {
    return value_cat == ValueCategory::String;
}
bool Value::is_boolean() const {
    return value_cat == ValueCategory::Boolean;
}

std::size_t Value::size() const {
    expected_cat(ValueCategory::Object, ValueCategory::Array,
                 ValueCategory::TypedArrayComplexFloat,
                 ValueCategory::TypedArrayComplexDouble);
    switch (value_cat) {
        case ValueCategory::Object:
            return static_cast<Object*>(ptr.get())->content.size();
        case ValueCategory::Array:
            return static_cast<Array*>(ptr.get())->content.size();
        case ValueCategory::TypedArrayComplexFloat:
            return static_cast<TypedArray<std::complex<float>>*>(ptr.get())
                ->content.size();
        case ValueCategory::TypedArrayComplexDouble:
            return static_cast<TypedArray<std::complex<double>>*>(ptr.get())
                ->content.size();
        default:
            return 0;  // unreachable
    }
}
std::size_t Value::empty() const {
    expected_cat(ValueCategory::Object, ValueCategory::Array);
    if (value_cat == ValueCategory::Object) {
        return static_cast<Object*>(ptr.get())->content.empty();
    } else {
        return static_cast<Array*>(ptr.get())->content.empty();
    }
}

ValueCategory Value::value_category() const {
    return value_cat;
}

std::string Value::dump() const {
    std::ostringstream oss;
    auto dump_typed_array =
        [&oss]<typename T>(const std::vector<T>& typed_array) {
            oss << '[';
            for (std::size_t i = 0; i < typed_array.size(); ++i) {
                const auto& val = typed_array[i];
                if (i) { oss << ','; }
                oss << '[' << val.real() << ',' << val.imag() << ']';
            }
            oss << ']';
        };

    switch (value_cat) {
        case ValueCategory::Null:
            oss << "null";
            break;
        case ValueCategory::Boolean:
            oss << (as_boolean() ? "true" : "false");
            break;
        case ValueCategory::NumberInt:
            oss << as_number<int>();
            break;
        case ValueCategory::NumberFloat:
            oss << as_number<double>();
            break;
        case ValueCategory::String:
            oss << '"' << static_cast<std::string>(*this) << '"';
            break;
        case ValueCategory::Object:
            oss << '{';
            for (const auto& [key, val] : as_object()) {
                oss << '"' << key << '"' << ':' << val.dump() << ',';
            }
            if (!empty()) { oss.seekp(-1, std::ios_base::cur); }
            oss << '}';
            break;
        case ValueCategory::Array:
            oss << '[';
            for (std::size_t i = 0; i < size(); ++i) {
                if (i) { oss << ','; }
                oss << operator[](i).dump();
            }
            oss << ']';
            break;
        case ValueCategory::TypedArrayComplexFloat:
            dump_typed_array(as_typed_array<std::complex<float>>());
            break;
        case ValueCategory::TypedArrayComplexDouble:
            dump_typed_array(as_typed_array<std::complex<double>>());
            break;
    }
    return oss.str();
}

std::string Value::pretty_print(std::size_t indent) const {
    std::ostringstream oss;
    auto print_typed_array =
        [&oss, indent]<typename T>(const std::vector<T>& typed_array) {
            oss << "[\n";
            for (std::size_t i = 0; i < typed_array.size(); ++i) {
                const auto& val = typed_array[i];
                if (i) { oss << ",\n"; }
                print_space(oss, indent + 4);
                oss << '[' << val.real() << ", " << val.imag() << ']';
            }
            if (typed_array.empty()) {
                oss.seekp(-1, std::ios_base::cur);
                oss << ' ';
            } else {
                oss << '\n';
                print_space(oss, indent);
            }
            oss << ']';
        };

    switch (value_cat) {
        case ValueCategory::Null:
        case ValueCategory::NumberInt:
        case ValueCategory::NumberFloat:
        case ValueCategory::Boolean:
        case ValueCategory::String:
            oss << dump();
            break;
        case ValueCategory::Object:
            oss << "{\n";
            for (const auto& [key, val] : as_object()) {
                print_space(oss, indent + 4);
                oss << '"' << key << '"' << ": " << val.pretty_print(indent + 4)
                    << ",\n";
            }
            if (empty()) {
                oss.seekp(-1, std::ios_base::cur);
                oss << ' ';
            } else {
                oss.seekp(-2, std::ios_base::cur);
                oss << '\n';
                print_space(oss, indent);
            }
            oss << '}';
            break;
        case ValueCategory::Array:
            oss << "[\n";
            for (std::size_t i = 0; i < size(); ++i) {
                if (i) { oss << ",\n"; }
                print_space(oss, indent + 4);
                oss << operator[](i).pretty_print(indent + 4);
            }
            if (empty()) {
                oss.seekp(-1, std::ios_base::cur);
                oss << ' ';
            } else {
                oss << '\n';
                print_space(oss, indent);
            }
            oss << ']';
            break;
        case ValueCategory::TypedArrayComplexFloat:
            print_typed_array(as_typed_array<std::complex<float>>());
            break;
        case ValueCategory::TypedArrayComplexDouble:
            print_typed_array(as_typed_array<std::complex<double>>());
            break;
    }
    return oss.str();
}

Value Value::clone() const {
    switch (value_cat) {
        case ValueCategory::NumberInt:
            return {value_cat, new NumberInt{as_number<int>()}};
            break;
        case ValueCategory::NumberFloat:
            return {value_cat, new NumberFloat{as_number<double>()}};
            break;
        case ValueCategory::Boolean:
            return {value_cat, new Boolean{as_boolean()}};
            break;
        case ValueCategory::String:
            return {value_cat, new String{as_string()}};
            break;
        case ValueCategory::Array: {
            auto arr = new Array;
            for (const auto& value : as_array()) {
                arr->content.push_back(value.clone());
            }
            return {value_cat, arr};
        }
        case ValueCategory::Object: {
            auto obj = new Object;
            for (const auto& [key, val] : as_object()) {
                obj->content.emplace(key, val.clone());
            }
            return {value_cat, obj};
        }
        case ValueCategory::TypedArrayComplexFloat:
            return {value_cat, new TypedArray<std::complex<float>>{
                                   as_typed_array<std::complex<float>>()}};
        case ValueCategory::TypedArrayComplexDouble:
            return {value_cat, new TypedArray<std::complex<double>>{
                                   as_typed_array<std::complex<double>>()}};
        case ValueCategory::Null:
            break;
    }
    return {value_cat, new Null{}};
}

void Value::print_space(std::ostream& os, std::size_t indent) {
    for (std::size_t i = 0; i < indent; ++i) { os << ' '; }
}

Value Value::create_object() {
    return {ValueCategory::Object, new Object{}};
}
Value Value::create_array(std::size_t n) {
    return {ValueCategory::Array, new Array{array_container_type{n}}};
}

JsonLexer::JsonLexer(std::istream& is, std::string file_name)
    : is_(is), filename(file_name), row(1), col(1) {}

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

const std::string& JsonLexer::get_filename() const {
    return filename;
}

void JsonLexer::read_token_to_buffer() {
    char c;
    // skip all whitespaces
    while (is_.get(c) && is_whitespace(c)) {
        ++col;
        if (c == '\n') {
            ++row;
            col = 1;
        }
    }
    buffer.content.clear();
    if (!is_) {
        buffer.name = TokenName::END_OF_FILE;
        buffer.col = 0;
    } else if (c == static_cast<char>(TokenName::BRACE_LEFT) ||
               c == static_cast<char>(TokenName::BRACE_RIGHT) ||
               c == static_cast<char>(TokenName::BRACKET_LEFT) ||
               c == static_cast<char>(TokenName::BRACKET_RIGHT) ||
               c == static_cast<char>(TokenName::COLON) ||
               c == static_cast<char>(TokenName::COMMA)) {
        buffer.name = static_cast<TokenName>(c);
        buffer.content.push_back(c);
        buffer.col = col;
        ++col;
    } else if (c == '"') {
        // a string
        buffer.name = TokenName::STRING;
        buffer.col = col;
        ++col;
        while (is_.get(c) && c != '"') {
            buffer.content.push_back(c);
            ++col;
        }
        ++col;
    } else if (is_digit_start(c)) {
        // a number
        bool is_float{};
        buffer.col = col;
        do {
            buffer.content.push_back(c);
            is_float |= c == '.';
            ++col;
        } while (is_.get(c) && is_digit(c));
        buffer.name = is_float ? TokenName::FLOAT : TokenName::INTEGER;
        is_.unget();
    } else if (c == 't' || c == 'f') {
        // a boolean
        buffer.name = TokenName::PRIMITIVE;
        char expected[6] = "true";
        if (c == 'f') { std::strcpy(expected, "false"); }
        char tmp[6] = {c};
        is_.readsome(tmp + 1, std::strlen(expected) - 1);
        if (std::strcmp(tmp, expected) == 0) {
            buffer.content = expected;
        } else {
            report_lexical_error();
        }
    } else if (c == 'n') {
        // null
        buffer.name = TokenName::PRIMITIVE;
        char tmp[5] = {c};
        is_.readsome(tmp + 1, 3);
        if (std::strcmp(tmp, "null") == 0) {
            buffer.content = "null";
        } else {
            report_lexical_error();
        }
    } else {
        report_lexical_error();
    }
    buffer.row = row;
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

void JsonLexer::report_lexical_error() const {
    std::ostringstream oss;
    oss << filename << ':' << row << ':' << col
        << ": error: unrecognized token";
    throw std::runtime_error(oss.str());
}

#ifdef EMME_DEBUG
std::ostream& operator<<(std::ostream& os, const JsonLexer::Token& token) {
    return os << "{ Name: " << ([&token] {
#define PROCESS_TOKEN_NAME(p) \
    case (p):                 \
        return #p;            \
        break;
               switch (token.name) {
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::PRIMITIVE)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::STRING)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::INTEGER)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::FLOAT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACE_LEFT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACE_RIGHT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACKET_LEFT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::BRACKET_RIGHT)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::COLON)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::COMMA)
                   PROCESS_TOKEN_NAME(JsonLexer::TokenName::END_OF_FILE)
               }
#undef PROCESS_TOKEN_NAME
               return "(No such name)";  // unreachable
           })()
              << ", Content: '" << token.content << "', position: ("
              << token.row << ", " << token.col << ')' << " }";
}
#endif  // EMME_DEBUG

JsonParser::JsonParser(JsonLexer&& json_lexer) : lexer(std::move(json_lexer)) {}

Value JsonParser::parse() {
    if (!lexer) { report_syntax_error(); }  // token list is empty
    auto value = parse_value();
    try_get_from_lexer(true);  // Check if lexer ends
    return value;
}

Value JsonParser::parse_value() {
    auto token = try_get_from_lexer();
    switch (token.name) {
        case JsonLexer::TokenName::STRING:
            return parse_string(token);
        case JsonLexer::TokenName::INTEGER:
            return parse_int(token);
        case JsonLexer::TokenName::FLOAT:
            return parse_float(token);
        case JsonLexer::TokenName::BRACE_LEFT:
            return parse_object();
        case JsonLexer::TokenName::BRACKET_LEFT:
            return parse_array();
        case JsonLexer::TokenName::PRIMITIVE:
            return parse_primitive(token);
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

Value JsonParser::parse_int(const JsonLexer::Token& token) {
    if (token.name != JsonLexer::TokenName::INTEGER) {
        report_syntax_error(token);
    }
    return {ValueCategory::NumberInt,
            new NumberInt{std::atoi(token.content.c_str())}};
}
Value JsonParser::parse_float(const JsonLexer::Token& token) {
    if (token.name != JsonLexer::TokenName::FLOAT) {
        report_syntax_error(token);
    }
    return {ValueCategory::NumberFloat,
            new NumberFloat{std::atof(token.content.c_str())}};
}
Value JsonParser::parse_primitive(const JsonLexer::Token& token) {
    if (token.name != JsonLexer::TokenName::PRIMITIVE) {
        report_syntax_error(token);
    }
    if (token.content.front() == 'n') {
        return {ValueCategory::Null, new Null{}};
    }
    return {ValueCategory::Boolean, new Boolean{token.content.front() == 't'}};
}

Value JsonParser::parse_object() {
    auto obj = new Object;
    auto token = try_peek_from_lexer();
    // empty object
    if (token.name == JsonLexer::TokenName::BRACE_RIGHT) {
        try_get_from_lexer();
        return {ValueCategory::Object, obj};
    }
    while (true) {
        auto key = try_get_and_check(JsonLexer::TokenName::STRING);
        try_get_and_check(JsonLexer::TokenName::COLON);
        obj->content.emplace(key.content, parse_value());  // value
        token = try_get_from_lexer();
        if (token.name == JsonLexer::TokenName::BRACE_RIGHT) { break; }
        if (token.name != JsonLexer::TokenName::COMMA) {
            report_syntax_error(token);
        }
    }
    return {ValueCategory::Object, obj};
}

Value JsonParser::parse_array() {
    auto arr = new Array;
    auto token = try_peek_from_lexer();
    // empty array
    if (token.name == JsonLexer::TokenName::BRACKET_RIGHT) {
        try_get_from_lexer();
        return {ValueCategory::Array, arr};
    }
    while (true) {
        arr->content.emplace_back(parse_value());
        token = try_get_from_lexer();
        if (token.name == JsonLexer::TokenName::BRACKET_RIGHT) { break; }
        if (token.name != JsonLexer::TokenName::COMMA) {
            report_syntax_error(token);
        }
    }
    return {ValueCategory::Array, arr};
}

JsonLexer::Token JsonParser::try_get_and_check(
    JsonLexer::TokenName expected_token_name) {
    auto token = lexer.get_token();
    if (token.name != expected_token_name) { report_syntax_error(token); }
    return token;
}

JsonLexer::Token JsonParser::try_get_from_lexer(bool end_expected) {
    auto token = lexer.get_token();
    if (end_expected ^ (token.name == JsonLexer::TokenName::END_OF_FILE)) {
        report_syntax_error(token);
    }
    return token;
}

JsonLexer::Token JsonParser::try_peek_from_lexer() {
    auto token = lexer.peek_token();
    if (token.name == JsonLexer::TokenName::END_OF_FILE) {
        report_syntax_error(token);
    }
    return token;
}

void JsonParser::report_syntax_error(const JsonLexer::Token& token) {
    std::ostringstream oss;
    oss << lexer.get_filename() << ':' << token.row << ':' << token.col
        << ": error: unexpected content '" << token.content << '\'';
#ifdef EMME_DEBUG
    oss << ". The token is " << token;
#endif
    throw std::runtime_error(oss.str());
}

Value parse(std::istream& is) {
    return JsonParser{is}.parse();
}

Value parse(std::string str) {
    std::stringstream ss;
    ss.str(str);
    return parse(ss);
}

Value parse_file(std::string filename) {
    std::ifstream ifs(filename);  // Its destructor will close the file
    if (!ifs) { throw std::runtime_error("File " + filename + " not found"); }
    return JsonParser{JsonLexer{ifs, filename}}.parse();
}

}  // namespace json
}  // namespace util
