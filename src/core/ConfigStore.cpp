// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ConfigStore.cpp
// =============================================================================
#include "core/ConfigStore.hpp"

#include <array>
#include <cctype>
#include <charconv>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace simall::core
{

// ----------------------------- ConfigValue ----------------------------------

namespace
{
const ConfigArray kEmptyArr;
const ConfigObject kEmptyObj;
}

bool ConfigValue::asBool(bool def) const
{
    if (auto* p = std::get_if<bool>(&v_))
        return *p;
    if (auto* p = std::get_if<std::int64_t>(&v_))
        return *p != 0;
    return def;
}
std::int64_t ConfigValue::asInt(std::int64_t def) const
{
    if (auto* p = std::get_if<std::int64_t>(&v_))
        return *p;
    if (auto* p = std::get_if<double>(&v_))
        return static_cast<std::int64_t>(*p);
    if (auto* p = std::get_if<bool>(&v_))
        return *p ? 1 : 0;
    return def;
}
double ConfigValue::asDouble(double def) const
{
    if (auto* p = std::get_if<double>(&v_))
        return *p;
    if (auto* p = std::get_if<std::int64_t>(&v_))
        return static_cast<double>(*p);
    return def;
}
std::string ConfigValue::asString(std::string_view def) const
{
    if (auto* p = std::get_if<std::string>(&v_))
        return *p;
    return std::string(def);
}
const ConfigArray& ConfigValue::asArray() const
{
    if (auto* p = std::get_if<ConfigArray>(&v_))
        return *p;
    return kEmptyArr;
}
const ConfigObject& ConfigValue::asObject() const
{
    if (auto* p = std::get_if<ConfigObject>(&v_))
        return *p;
    return kEmptyObj;
}

// ----------------------------- JSON parser ----------------------------------
//
// Tight, self-contained JSON-1 parser (RFC 8259 subset). Sufficient for
// configuration files: supports objects, arrays, strings, numbers (int or
// double), bool, null. Throws std::runtime_error on syntax errors.

namespace
{

class JsonParser
{
public:
    explicit JsonParser(std::string_view text) : t_(text) {}

    ConfigValue parse()
    {
        skipWs();
        ConfigValue v = parseValue();
        skipWs();
        if (i_ != t_.size())
            throw std::runtime_error("JSON: trailing data at " + pos());
        return v;
    }

private:
    void skipWs()
    {
        while (i_ < t_.size()) {
            char c = t_[i_];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
                ++i_;
            else if (c == '/' && i_ + 1 < t_.size() && t_[i_ + 1] == '/') { // line comment ext
                while (i_ < t_.size() && t_[i_] != '\n')
                    ++i_;
            } else
                break;
        }
    }
    char peek()
    {
        if (i_ >= t_.size())
            throw std::runtime_error("JSON: unexpected EOF");
        return t_[i_];
    }
    char get()
    {
        if (i_ >= t_.size())
            throw std::runtime_error("JSON: unexpected EOF");
        return t_[i_++];
    }
    std::string pos() { return "offset " + std::to_string(i_); }

    ConfigValue parseValue()
    {
        skipWs();
        char c = peek();
        if (c == '{')
            return parseObject();
        if (c == '[')
            return parseArray();
        if (c == '"')
            return ConfigValue{parseString()};
        if (c == 't' || c == 'f')
            return ConfigValue{parseBool()};
        if (c == 'n') {
            expectLit("null");
            return ConfigValue{};
        }
        return parseNumber();
    }

    void expectLit(std::string_view lit)
    {
        if (t_.size() - i_ < lit.size() || t_.compare(i_, lit.size(), lit) != 0)
            throw std::runtime_error("JSON: expected '" + std::string(lit) + "' at " + pos());
        i_ += lit.size();
    }

    bool parseBool()
    {
        if (peek() == 't') {
            expectLit("true");
            return true;
        }
        expectLit("false");
        return false;
    }

    std::string parseString()
    {
        if (get() != '"')
            throw std::runtime_error("JSON: expected string at " + pos());
        std::string out;
        while (true) {
            char c = get();
            if (c == '"')
                return out;
            if (c == '\\') {
                char e = get();
                switch (e) {
                case '"':
                    out += '"';
                    break;
                case '\\':
                    out += '\\';
                    break;
                case '/':
                    out += '/';
                    break;
                case 'b':
                    out += '\b';
                    break;
                case 'f':
                    out += '\f';
                    break;
                case 'n':
                    out += '\n';
                    break;
                case 'r':
                    out += '\r';
                    break;
                case 't':
                    out += '\t';
                    break;
                case 'u': {
                    if (t_.size() - i_ < 4)
                        throw std::runtime_error("JSON: bad \\u");
                    unsigned cp = 0;
                    for (int k = 0; k < 4; ++k) {
                        char h = get();
                        cp <<= 4;
                        if (h >= '0' && h <= '9')
                            cp |= unsigned(h - '0');
                        else if (h >= 'a' && h <= 'f')
                            cp |= unsigned(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F')
                            cp |= unsigned(h - 'A' + 10);
                        else
                            throw std::runtime_error("JSON: bad hex");
                    }
                    // Encode as UTF-8.
                    if (cp < 0x80)
                        out += char(cp);
                    else if (cp < 0x800) {
                        out += char(0xC0 | (cp >> 6));
                        out += char(0x80 | (cp & 0x3F));
                    } else {
                        out += char(0xE0 | (cp >> 12));
                        out += char(0x80 | ((cp >> 6) & 0x3F));
                        out += char(0x80 | (cp & 0x3F));
                    }
                    break;
                }
                default:
                    throw std::runtime_error("JSON: bad escape");
                }
            } else {
                out += c;
            }
        }
    }

    ConfigValue parseNumber()
    {
        std::size_t start = i_;
        if (peek() == '-')
            ++i_;
        bool isDouble = false;
        while (i_ < t_.size()) {
            char c = t_[i_];
            if (c == '.' || c == 'e' || c == 'E') {
                isDouble = true;
                ++i_;
            } else if ((c >= '0' && c <= '9') || c == '+' || c == '-') {
                ++i_;
            } else
                break;
        }
        std::string_view tok = t_.substr(start, i_ - start);
        if (isDouble) {
            try {
                return ConfigValue{std::stod(std::string(tok))};
            } catch (...) {
                throw std::runtime_error("JSON: bad number " + std::string(tok));
            }
        }
        std::int64_t n = 0;
        auto res = std::from_chars(tok.data(), tok.data() + tok.size(), n);
        if (res.ec != std::errc{})
            throw std::runtime_error("JSON: bad integer " + std::string(tok));
        return ConfigValue{n};
    }

    ConfigValue parseArray()
    {
        if (get() != '[')
            throw std::runtime_error("JSON: expected [");
        ConfigArray out;
        skipWs();
        if (peek() == ']') {
            ++i_;
            return ConfigValue{std::move(out)};
        }
        while (true) {
            out.push_back(parseValue());
            skipWs();
            char c = get();
            if (c == ']')
                return ConfigValue{std::move(out)};
            if (c != ',')
                throw std::runtime_error("JSON: expected ',' or ']'");
            skipWs();
        }
    }

    ConfigValue parseObject()
    {
        if (get() != '{')
            throw std::runtime_error("JSON: expected {");
        ConfigObject out;
        skipWs();
        if (peek() == '}') {
            ++i_;
            return ConfigValue{std::move(out)};
        }
        while (true) {
            skipWs();
            std::string key = parseString();
            skipWs();
            if (get() != ':')
                throw std::runtime_error("JSON: expected ':'");
            ConfigValue val = parseValue();
            out.emplace(std::move(key), std::move(val));
            skipWs();
            char c = get();
            if (c == '}')
                return ConfigValue{std::move(out)};
            if (c != ',')
                throw std::runtime_error("JSON: expected ',' or '}'");
        }
    }

    std::string_view t_;
    std::size_t i_ = 0;
};

void deepMerge(ConfigValue& dst, const ConfigValue& src)
{
    if (src.isObject() && dst.isObject()) {
        auto& d = std::get<ConfigObject>(dst.raw());
        const auto& s = std::get<ConfigObject>(src.raw());
        for (auto const& [k, v] : s) {
            auto it = d.find(k);
            if (it == d.end())
                d.emplace(k, v);
            else
                deepMerge(it->second, v);
        }
    } else {
        dst = src;
    }
}

ConfigValue coerceCli(std::string_view raw)
{
    if (raw == "true")
        return ConfigValue{true};
    if (raw == "false")
        return ConfigValue{false};
    if (raw == "null")
        return ConfigValue{};
    // try int
    {
        std::int64_t n = 0;
        auto res = std::from_chars(raw.data(), raw.data() + raw.size(), n);
        if (res.ec == std::errc{} && res.ptr == raw.data() + raw.size())
            return ConfigValue{n};
    }
    // try double
    try {
        std::size_t consumed = 0;
        double d = std::stod(std::string(raw), &consumed);
        if (consumed == raw.size())
            return ConfigValue{d};
    } catch (...) {
    }
    return ConfigValue{std::string(raw)};
}

void writeAtPath(ConfigObject& root, std::string_view dotted, ConfigValue value)
{
    std::size_t pos = 0;
    ConfigObject* cur = &root;
    while (true) {
        std::size_t dot = dotted.find('.', pos);
        std::string key{dotted.substr(pos, dot - pos)};
        if (dot == std::string_view::npos) {
            (*cur)[std::move(key)] = std::move(value);
            return;
        }
        auto [it, inserted] = cur->try_emplace(key, ConfigValue{ConfigObject{}});
        if (!it->second.isObject())
            it->second = ConfigValue{ConfigObject{}};
        cur = &std::get<ConfigObject>(it->second.raw());
        pos = dot + 1;
    }
}

} // namespace

// ----------------------------- ConfigStore Impl -----------------------------

struct ConfigStore::Impl
{
    std::array<ConfigValue, static_cast<std::size_t>(ConfigLayer::_Count)> layers;

    Impl()
    {
        for (auto& l : layers)
            l = ConfigValue{ConfigObject{}};
    }

    ConfigValue effective() const
    {
        ConfigValue out = ConfigValue{ConfigObject{}};
        for (auto const& l : layers)
            deepMerge(out, l);
        return out;
    }
};

ConfigStore::ConfigStore() : p_(std::make_unique<Impl>()) {}
ConfigStore::~ConfigStore() = default;

void ConfigStore::setLayer(ConfigLayer layer, ConfigValue root)
{
    if (!root.isObject())
        root = ConfigValue{ConfigObject{}};
    p_->layers[static_cast<std::size_t>(layer)] = std::move(root);
}

void ConfigStore::mergeLayer(ConfigLayer layer, const ConfigValue& overlay)
{
    auto& dst = p_->layers[static_cast<std::size_t>(layer)];
    if (!dst.isObject())
        dst = ConfigValue{ConfigObject{}};
    deepMerge(dst, overlay);
}

bool ConfigStore::loadJsonFile(ConfigLayer layer, const std::filesystem::path& file)
{
    std::ifstream in(file, std::ios::binary);
    if (!in)
        return false;
    std::ostringstream ss;
    ss << in.rdbuf();
    try {
        JsonParser parser{ss.str()};
        ConfigValue v = parser.parse();
        mergeLayer(layer, v);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

void ConfigStore::applyCli(int argc, const char* const* argv)
{
    auto& dst = p_->layers[static_cast<std::size_t>(ConfigLayer::Cli)];
    if (!dst.isObject())
        dst = ConfigValue{ConfigObject{}};
    auto& obj = std::get<ConfigObject>(dst.raw());
    for (int i = 0; i < argc; ++i) {
        std::string_view a{argv[i]};
        if (a.size() < 3 || a[0] != '-' || a[1] != '-')
            continue;
        a.remove_prefix(2);
        std::size_t eq = a.find('=');
        if (eq == std::string_view::npos) {
            writeAtPath(obj, a, ConfigValue{true});
        } else {
            std::string_view k = a.substr(0, eq);
            std::string_view v = a.substr(eq + 1);
            writeAtPath(obj, k, coerceCli(v));
        }
    }
}

std::optional<ConfigValue> ConfigStore::get(std::string_view dottedPath) const
{
    ConfigValue eff = p_->effective();
    if (!eff.isObject())
        return std::nullopt;
    const ConfigObject* obj = &eff.asObject();
    std::size_t pos = 0;
    ConfigValue last;
    while (true) {
        std::size_t dot = dottedPath.find('.', pos);
        std::string key{dottedPath.substr(pos, dot - pos)};
        auto it = obj->find(key);
        if (it == obj->end())
            return std::nullopt;
        if (dot == std::string_view::npos)
            return it->second;
        if (!it->second.isObject())
            return std::nullopt;
        obj = &it->second.asObject();
        pos = dot + 1;
    }
}

bool ConfigStore::getBool(std::string_view path, bool def) const
{
    auto v = get(path);
    return v ? v->asBool(def) : def;
}
std::int64_t ConfigStore::getInt(std::string_view path, std::int64_t def) const
{
    auto v = get(path);
    return v ? v->asInt(def) : def;
}
double ConfigStore::getDouble(std::string_view path, double def) const
{
    auto v = get(path);
    return v ? v->asDouble(def) : def;
}
std::string ConfigStore::getString(std::string_view path, std::string_view def) const
{
    auto v = get(path);
    return v ? v->asString(def) : std::string(def);
}

ConfigValue ConfigStore::effective() const
{
    return p_->effective();
}

} // namespace simall::core
