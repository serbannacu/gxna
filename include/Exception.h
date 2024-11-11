#pragma once

#include <exception>
#include <sstream>
#include <string>

namespace gxna {

// Lightweight exception class.
// Store and print an error message with some minimal formatting.

class Exception : public std::exception {
 public:
    explicit Exception(const std::string& msg)
        : m_msg(msg)
    {}

    virtual const char *what() const noexcept {
        return m_msg.c_str();
    }

    // Append string representation of value to exception message.
    template<class T>
    Exception& operator<<(const T& val) {
        std::ostringstream oss;
        oss << val;
        m_msg += oss.str();
        return *this;
    }

 private:
    std::string m_msg;
};

}  // namespace gxna
