#pragma once

#include <exception>
#include <string>

namespace gxna {

class Exception : public std::exception {
public:
    Exception(const std::string& msg) : m_msg(msg) {}

    virtual const char *what() const noexcept { return m_msg.c_str(); }

private:
    std::string m_msg;
};

} // namespace gxna
