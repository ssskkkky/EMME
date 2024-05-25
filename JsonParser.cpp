#include "JsonParser.h"

#ifdef EMME_DEBUG
#include <ostream>
#endif

namespace util {
namespace json {

JsonLexer::JsonLexer(std::istream& is) : is_(is) {}

JsonLexer::Token JsonLexer::get_token() {
    char c;
    Token token;
    // skip all whitespaces
    while (is_.get(c) && is_whitespace(c)) {}
    if (!is_) {
        token.name = TokenName::END_OF_FILE;
    } else if (c == static_cast<char>(TokenName::BRACE_LEFT) ||
               c == static_cast<char>(TokenName::BRACE_RIGHT) ||
               c == static_cast<char>(TokenName::BRACKET_LEFT) ||
               c == static_cast<char>(TokenName::BRACKET_RIGHT) ||
               c == static_cast<char>(TokenName::COLON) ||
               c == static_cast<char>(TokenName::COMMA)) {
        token.name = static_cast<TokenName>(c);
    } else if (c == '"') {
        // a string
        token.name = TokenName::STRING;
        while (is_.get(c) && c != '"') { token.content.push_back(c); }
    } else if (is_digit_start(c)) {
        // a number
        token.name = TokenName::NUMBER;
        do { token.content.push_back(c); } while (is_.get(c) && is_digit(c));
        is_.unget();
    }
    return token;
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

#ifdef EMME_DEBUG
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
#endif
}  // namespace json
}  // namespace util
