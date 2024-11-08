#pragma once

#include <fstream>
#include <iostream>

namespace gxna {

// Basic mechanism to log to a file.
// As long as the LogGuard object exists, output to std::clog will be
// redirected to the file.

// The ctor creates an ofstream and redirects to it, the dtor destroys it
// and restores clog to its original buffer. This ensure that clog will
// not try to write to an invalid buffer after the ofstream no longer exists.

class LogGuard {
 public:
    LogGuard(const char *filename) {
        m_os = new std::ofstream(filename);
        if (*m_os) {
            m_buf = std::clog.rdbuf();
            std::clog.rdbuf(m_os->rdbuf());  // redirect clog
        }
        else {  // stream opening failed
            delete m_os;
            m_os = nullptr;
        }
    }

    ~LogGuard() {
        if (m_os) {
            std::clog.rdbuf(m_buf);  // restore to original buffer
            delete m_os;
        }
    }

    LogGuard(const LogGuard&) = delete;
    LogGuard& operator=(const LogGuard&) = delete;  // non copyable

 private:
    std::streambuf *m_buf = nullptr;
    std::ofstream *m_os;
};

}  // gxna
