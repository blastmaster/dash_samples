#include <iostream>
#include <libdash.h>

/**
 * utility function to return all elements of an array as a std::string
 * separated by " ".
 */
template<class T>
std::string dump(dash::Array<T>& arr)
{
    std::stringstream ss;
    for (auto it : arr) {
        ss << (T)it << " ";
    }
    ss << std::endl;
    return ss.str();
}

/**
 * get the border elements from a given dash::Array and write them in the given
 * border array. So each unit has his borders in separate local storage.
 */
template<class T>
void borders_from_values(dash::Array<T>& b, dash::Array<T>& v, int offset=1)
{
    if (b.lbegin() != b.begin())
        *(b.lbegin()) = *(v.lbegin() + (offset - 1));
    if (b.lend() != b.end())
        *(--b.lend()) = *(v.lend() - offset);
}

/**
 * swap border value with neighboor unit
 */
template<class T>
void borders_iexchange(dash::Array<T>& b)
{
    b.barrier();
    if (b.lend() != b.end()) {
        T tmp = *(--b.lend());
        *(--b.lend()) = *(b.lend());
        *(b.lend()) = tmp;
    }
    b.barrier();
}
