#ifndef ARGUMENT_EXCEPTION
#define ARGUMENT_EXCEPTION

#include <stdlib.h>
#include <string>

class ArgumentException: public std::exception {
private:
public:
    std::string message;
    ArgumentException(std::string msg) : message(msg) {
        this->message = msg;
    }
};
#endif