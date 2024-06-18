#include "big_integer.h"
#include <iostream>
#include <cstring>
#include <cmath>

void big_integer::clear() {
    _oldest_digit = 0;
    delete[] _other_digits;
    _other_digits = nullptr;
}

void big_integer::copy(big_integer const &other) {
    _oldest_digit = other._oldest_digit;
    _other_digits = nullptr;
    if (other._other_digits == nullptr)
        return;
    _other_digits = new unsigned int [*other._other_digits];
    std::memcpy(_other_digits, other._other_digits, sizeof(unsigned int) * (*other._other_digits));
}

void big_integer::initialize(int const *digits, size_t digits_count) {
    if (digits == nullptr)
        throw std::logic_error("pointer to digits array must not be nullptr");
    if (digits_count == 0)
        throw std::logic_error("digits array length must  be GT 0");
    unsigned int count_of_additional_zeros = 0;
    while (*digits == 0)
        count_of_additional_zeros++;
    _oldest_digit = digits[digits_count - 1];
    _other_digits = (digits_count == 1) ? nullptr : new unsigned int[digits_count - count_of_additional_zeros];
    if (_other_digits == nullptr)
        return;
    *_other_digits = (unsigned int)digits_count - count_of_additional_zeros;
    std::memcpy(_other_digits + 1, digits + count_of_additional_zeros, sizeof(unsigned int) * (digits_count - 1 - count_of_additional_zeros));
}

void big_integer::initialize(std::vector<int> const &digits, size_t digits_count) {
    _other_digits = nullptr;
    if (digits.empty() || digits_count == 0)
        throw std::logic_error("std::vector<int> of digits should not be empty");
    _oldest_digit = digits[digits_count - 1];
    if (digits_count == 1)
        return;
    _other_digits = new unsigned int[digits_count];
    *_other_digits = (unsigned int)digits_count;
    for (auto i = 0; i < digits_count - 1; ++i)
        _other_digits[1 + i] = *reinterpret_cast<unsigned int const *>(&digits[i]);
}

void big_integer::initialize(std::string const &value, size_t base) {
    if (value.empty())
        throw std::logic_error("string with value must not be nullptr");
    if (base < 2 || base > 36)
        throw std::invalid_argument("Invalid base");

    int const zero[] = {0};
    big_integer number(zero, 1);
    unsigned int digit = 0, pow = 0, num = 0, last_digit = 0;

    for (auto i = size(value); i > 0; i--) {
        digit = 0;
        char c = value[i - 1];
        if (c >= '0' && c <= '9')
            digit = c - '0';
        else if (c >= 'A' && c <= 'Z')
            digit = c - 'A' + 10;
        else if (c == '-')
            continue;
        else
            throw std::invalid_argument("Invalid character in input string");

        last_digit = digit * std::pow(base, pow++);
        num += last_digit;
        if (num < num - last_digit) {
            int const digit_arr[] = {static_cast<int>(num - last_digit)};
            big_integer big_int_digit(digit_arr, 1);
            number += big_int_digit;
            num = last_digit;
        }
    }

    if (value[0] == '-')
        number.change_sign();

    *this = number;
}

big_integer::big_integer(int const *digits, size_t digits_count) {
    initialize(digits, digits_count);
}

big_integer::big_integer(std::vector<int> const &digits) {
    initialize(digits, digits.size());
}

big_integer::big_integer(std::string const &value, size_t base) {
    initialize(value, base);
}

big_integer::big_integer(big_integer const &other) {
    copy(other);
}

big_integer::~big_integer() {
    clear();
}

big_integer &big_integer::operator=(big_integer const &other){
    if (this != &other){
        clear();
        copy(other);
    }
    return *this;
}

inline bool big_integer::is_equal_to_zero() const noexcept {
    return _oldest_digit == 0 && _other_digits == nullptr;
}

inline int big_integer::sign() const noexcept {
    if (is_equal_to_zero())
        return 0;

    return 1 - (static_cast<int>((*reinterpret_cast<unsigned int const *>(&_oldest_digit) >> ((sizeof(int) << 3) - 1))) << 1);
}

inline int big_integer::get_digits_count() const noexcept {
    return static_cast<int>(_other_digits == nullptr ? 1 : *_other_digits);
}

inline unsigned int big_integer::get_digit(int position) const noexcept {
    if (_other_digits == nullptr)
        return position == 0 ? _oldest_digit : 0;
    int const digits_count = get_digits_count();
    if (position < digits_count - 1)
        return _other_digits[position + 1];
    if (position == digits_count - 1)
        return _oldest_digit;
    return 0;
}

big_integer &big_integer::change_sign() {
    _oldest_digit ^= (1 << ((sizeof(int) << 3) - 1));
    return *this;
}

void big_integer::print_byte(std::ostream &stream, unsigned char byte_value) {
    for(int i = 0; i < 8; i++)
        stream << ((byte_value >> (7 - i)) & 1);
}

void big_integer::dump_int_value(std::ostream &stream, int value) {
    auto *p = (unsigned char *)&value;
    for (int i = 0; i < sizeof(int); i++) {
        print_byte(stream, *p++);
        stream << ' ';
    }
}

big_integer &big_integer::operator+=(big_integer const &other) {
    if (other.is_equal_to_zero())
        return *this;

    if (is_equal_to_zero())
        return *this = other;

    if (sign() == -1) {
        change_sign();
        *this += -other;
        return change_sign();
    }

    if (other.sign() == -1)
        return *this -= -other;

    auto const first_value_digits_count = get_digits_count();
    auto const second_value_digits_count = other.get_digits_count();
    auto const digits_count = std::max(first_value_digits_count, second_value_digits_count);

    unsigned int operation_result = 0;

    constexpr int shift = sizeof(unsigned int) << 2;
    constexpr int mask = (1 << shift) - 1;

    std::vector<int> result_digits(digits_count + 1);

    for (int i = 0; i < digits_count; ++i) {
        unsigned int first_value_digit = get_digit(i);
        unsigned int second_value_digit = other.get_digit(i);
        result_digits[i] = 0;

        for (int j = 0; j < 2; ++j) {
            operation_result += (first_value_digit & mask) + (second_value_digit & mask);
            first_value_digit >>= shift;
            second_value_digit >>= shift;
            *reinterpret_cast<unsigned int *>(&result_digits[i]) |= ((operation_result & mask) << shift * j);
            operation_result >>= shift;
        }
    }

    auto result_digits_count = result_digits.size();

    if (operation_result == 1) {
        if ((*reinterpret_cast<unsigned int *>(&result_digits[digits_count - 1]) >> ((sizeof(unsigned int) << 3) - 1)) == 0)
            *reinterpret_cast<unsigned int *>(&result_digits[digits_count - 1]) |= (1u << ((sizeof(unsigned int) << 3) - 1));
        else
            result_digits.back() = 1;
    }
    else if ((*reinterpret_cast<unsigned int *>(&result_digits[digits_count - 1]) >> ((sizeof(unsigned int) << 3) - 1)) == 0)
        --result_digits_count;

    clear();
    initialize(result_digits, result_digits_count);

    return *this;
}

big_integer big_integer::operator+(big_integer const &other) const {
    return big_integer(*this) += other;
}

big_integer &big_integer::operator-=(big_integer const &other) {
    if (other.is_equal_to_zero())
        return *this;

    if (is_equal_to_zero())
        return *this = -other;

    if (*this == other) {
        int const zero = 0;
        return *this = big_integer(&zero, 1);
    }

    if (sign() == -1 && other.sign() == -1)
        return (this->change_sign() -= -other).change_sign();

    if (sign() == -1 && other.sign() == 1)
        return (this->change_sign() += other).change_sign();

    if (sign() == 1 && other.sign() == -1)
        return *this += -other;

    if (*this < other)
        return (*this = other - *this).change_sign();

    bool flag = false;
    std::vector<int> digits(get_digits_count());
    for(unsigned i = 0; i < this->get_digits_count(); i++) {
        auto temp1 = get_digit(i);
        auto temp2 = other.get_digit(i);
        auto result = 0;

        if (flag) {
            temp1--;
            if (temp1 == UINT_MAX) {
                temp1 = 0;
                flag = true;
            }
            else
                flag = false;
        }
        else if (temp1 < temp2) {
            result = temp1 + UINT_MAX - temp2;
            flag = true;
        }
        else
            result = temp1 - temp2;
            flag = false;

        digits.push_back(result);
    }

    clear();
    initialize(digits, digits.size());

    return *this;
}

big_integer big_integer::operator-(big_integer const &other) const {
    return big_integer(*this) -= other;
}

big_integer big_integer::operator-() const {
    return big_integer(*this).change_sign();
}


big_integer &big_integer::operator*=(big_integer const &other) {
    if (other.is_equal_to_zero())
        return *this = other;

    if (is_equal_to_zero())
        return *this;

    auto need_to_change_sign = false;
    if (sign() == -1) {
        change_sign();
        need_to_change_sign = !need_to_change_sign;
    }

    if (other.sign() == -1) {
        need_to_change_sign = !need_to_change_sign;
        *this *= -other;
    }

    auto const first_value_digits_count = get_digits_count();
    auto const second_value_digits_count = other.get_digits_count();

    constexpr int shift = sizeof(unsigned int) << 2;
    constexpr int mask = (1 << shift) - 1;

    big_integer multiplication_result(std::vector<int> { 0 });
    unsigned int single_digits_multiplication_result[2];
    single_digits_multiplication_result[1] = 0;
    int const *single_digits_multiplication_result_int_pointer = reinterpret_cast<int const *>(single_digits_multiplication_result);

    for (int i = 0; i < first_value_digits_count; ++i) {
        unsigned int const first_value_digit_loword = get_digit(i) & mask;
        unsigned int const first_value_digit_hiword = get_digit(i) >> shift;

        for (int j = 0; j < second_value_digits_count; ++j) {
            unsigned int const second_value_digit_loword = other.get_digit(j) & mask;
            unsigned int const second_value_digit_hiword = other.get_digit(j) >> shift;

            single_digits_multiplication_result[0] = first_value_digit_loword * second_value_digit_loword;
            multiplication_result += big_integer(single_digits_multiplication_result_int_pointer, 1 + (single_digits_multiplication_result[0] >> ((sizeof(unsigned int) << 3) - 1))) <<= ((sizeof(unsigned int) << 2) * ((i + j) << 1));

            single_digits_multiplication_result[0] = first_value_digit_loword * second_value_digit_hiword;
            multiplication_result += big_integer(single_digits_multiplication_result_int_pointer, 1 + (single_digits_multiplication_result[0] >> ((sizeof(unsigned int) << 3) - 1))) <<= ((sizeof(unsigned int) << 2) * (((i + j) << 1) + 1));

            single_digits_multiplication_result[0] = first_value_digit_hiword * second_value_digit_loword;
            multiplication_result += big_integer(single_digits_multiplication_result_int_pointer, 1 + (single_digits_multiplication_result[0] >> ((sizeof(unsigned int) << 3) - 1))) <<= ((sizeof(unsigned int) << 2) * (((i + j) << 1) + 1));

            single_digits_multiplication_result[0] = first_value_digit_hiword * second_value_digit_hiword;
            multiplication_result += big_integer(single_digits_multiplication_result_int_pointer, 1 + (single_digits_multiplication_result[0] >> ((sizeof(unsigned int) << 3) - 1))) <<= ((sizeof(unsigned int) << 2) * (((i + j) << 1) + 2));
        }
    }

    *this = multiplication_result;

    if (need_to_change_sign)
        change_sign();

    return *this;
}

big_integer big_integer::operator*(big_integer const &other) const {
    return big_integer(*this) *= other;
}


big_integer &big_integer::operator<<=(size_t shift_value) {
    if (is_equal_to_zero() || shift_value == 0)
        return *this;

    auto value_sign = sign();
    if (value_sign == -1)
        change_sign();

    auto const added_by_shift_at_other_digits_digits_count = shift_value / (sizeof(unsigned int) << 3);
    shift_value %= (sizeof(unsigned int) << 3);

    auto added_by_shift_at_oldest_digit_digits_count = 0;
    if (_oldest_digit != 0) {
        unsigned int oldest_digit = *reinterpret_cast<unsigned int *>(&_oldest_digit);
        int oldest_value_bit_index = 0;
        while (oldest_digit != 1) {
            oldest_digit >>= 1;
            ++oldest_value_bit_index;
        }

        if (oldest_value_bit_index + shift_value > (sizeof(int) << 3) - 2)
            ++added_by_shift_at_oldest_digit_digits_count;
    }

    if (added_by_shift_at_oldest_digit_digits_count != 0 || added_by_shift_at_other_digits_digits_count != 0) {
        auto const added_digits_count = added_by_shift_at_oldest_digit_digits_count + added_by_shift_at_other_digits_digits_count;

        if (_other_digits == nullptr) {
            _other_digits = new unsigned int[added_digits_count + 1];
            *_other_digits = added_digits_count + 1;
            std::memset(_other_digits + 1, 0, sizeof(unsigned int) * (added_digits_count - 1));
            if (added_by_shift_at_oldest_digit_digits_count != 0) {
                _other_digits[*_other_digits - 1] = _oldest_digit;
                _oldest_digit = 0;
            }
            else
                _other_digits[*_other_digits - 1] = 0;
        }
        else {
            auto *new_digits = new unsigned int[added_digits_count + *_other_digits];
            std::memset(new_digits + 1, 0, sizeof(unsigned int) * added_digits_count);
            if (added_by_shift_at_oldest_digit_digits_count != 0) {
                new_digits[added_digits_count + *_other_digits - 1] = _oldest_digit;
                _oldest_digit = 0;
            }
            std::memcpy(new_digits + 1 + added_by_shift_at_other_digits_digits_count, _other_digits + 1, sizeof(unsigned int) * (*_other_digits - 1));
            *new_digits = *_other_digits + added_digits_count;

            clear();
            _other_digits = new_digits;
        }
    }

    if (shift_value != 0) {
        auto const digits_count = get_digits_count();
        unsigned int part_to_move_to_next_digit = 0;
        for (auto i = 0; i < digits_count; ++i) {
            auto digit_value = get_digit(i);
            auto *digit_address = i == digits_count - 1 ? reinterpret_cast<unsigned int *>(&_oldest_digit) : _other_digits + 1 + i;
            *digit_address <<= shift_value;
            *digit_address |= part_to_move_to_next_digit;
            part_to_move_to_next_digit = digit_value >> ((sizeof(unsigned int) << 3) - shift_value);
        }
    }

    if (value_sign == -1)
        change_sign();

    return *this;
}

big_integer big_integer::operator<<(size_t shift_value) const {
    return big_integer(*this) <<= shift_value;
}


bool big_integer::operator<(big_integer const &other) const {
    if (is_equal_to_zero() && other.is_equal_to_zero())
        return false;

    if (is_equal_to_zero())
        return other.sign() == 1;

    if (other.is_equal_to_zero())
        return sign() == -1;

    if (sign() == -1 && other.sign() == -1)
        return !(-*this < -other);

    if (sign() == -1 && other.sign() == 1)
        return true;

    if (sign() == 1 && other.sign() == -1)
        return false;

    if (_other_digits == nullptr && other._other_digits == nullptr)
        return _oldest_digit < other._oldest_digit;

    if (_other_digits != nullptr && other._other_digits == nullptr)
        return false;

    if (_other_digits == nullptr && other._other_digits != nullptr)
        return true;

    if (*_other_digits != *other._other_digits)
        return *_other_digits < *other._other_digits;

    for (auto i = *_other_digits - 1; i >= 0; --i) {
        auto this_digit = get_digit(i);
        auto other_digit = other.get_digit(i);

        if (this_digit != other_digit)
            return this_digit < other_digit;
    }

    return false;
}

bool big_integer::operator<=(big_integer const &other) const {
    return *this < other || *this == other;
}


bool big_integer::operator==(big_integer const &other) const {
    if (this == &other)
        return true;

    if (_oldest_digit != other._oldest_digit)
        return false;

    if (_other_digits == nullptr && other._other_digits == nullptr)
        return true;

    if (_other_digits == nullptr ^ other._other_digits == nullptr)
        return false;

    if (*_other_digits != *other._other_digits)
        return false;

    for (auto i = 0; i < *_other_digits; ++i) {
        if (_other_digits[1 + i] != other._other_digits[1 + i])
            return false;
    }

    return true;
}

bool big_integer::operator!=(big_integer const &other) const {
    return !(*this == other);
}

int main() {
    std::cout << 135 + UINT_MAX - 200;
}