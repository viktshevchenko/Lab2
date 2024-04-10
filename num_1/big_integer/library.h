#ifndef BIG_INTEGER_LIBRARY_H
#define BIG_INTEGER_LIBRARY_H

#include <corecrt.h>
#include <vector>
#include <string>

class big_integer final {

private:
    int _oldest_digit;
    unsigned int *_other_digits;

public:
    big_integer(int *digits, size_t digits_count);
    big_integer(std::vector<int> const &digits);
    big_integer(std::string const &value, size_t base);

    big_integer(big_integer const &other);
    big_integer operator=(big_integer const &other);
    ~big_integer();

    big_integer &operator+=(big_integer const &other);
    big_integer operator+(big_integer const &other) const;

    big_integer &operator-=(big_integer const &other);
    big_integer operator-(big_integer const &other) const;

    big_integer &operator*=(big_integer const &other);
    big_integer operator*(big_integer const &other) const;

    big_integer &operator/=(big_integer const &other);
    big_integer operator/(big_integer const &other) const;

    big_integer &operator%=(big_integer const &other);
    big_integer operator%(big_integer const &other) const;

    bool operator==(big_integer const &other) const;
    bool operator!=(big_integer const &other) const;

    bool operator<(big_integer const &other) const;
    bool operator<=(big_integer const &other) const;

    bool operator>(big_integer const &other) const;
    bool operator>=(big_integer const &other) const;

    big_integer operator~();
    big_integer operator&(big_integer const &other) const;
    big_integer operator|(big_integer const &other) const;
    big_integer operator^(big_integer const &other) const;

    big_integer operator<<(size_t sift_value) const;
    big_integer operator>>(size_t sift_value) const;

    //friend << and friend >>

};


#endif //BIG_INTEGER_LIBRARY_H
