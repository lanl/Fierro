#ifndef ARGUMENT_EXCEPTION
#define ARGUMENT_EXCEPTION

#include <stdlib.h>
#include <string>

class ArgumentException: public std::exception {
private:
    std::string message;
public:
    ArgumentException(std::string msg) : message(msg) {
        this->message = msg;
    }
    char* what() {
        return (char*)("Argument Error: " + this->message).c_str();
    }
};
#endif