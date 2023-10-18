#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#ifdef DEBUG
    #define Address(str) std::cout << "Address:" << hex << long(this) << str << std::endl;
#else
    #define Address(str)
#endif

#endif // __UTILS_HPP__