#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <cassert>


constexpr unsigned int K = 10; // bitness of fragment. Value = 0..2^k
constexpr unsigned int M = 10; // amount of fragments (Y)
constexpr unsigned int NUMBERS_AMOUNT = 1 << K;
constexpr unsigned int NUMBER_MAX = NUMBERS_AMOUNT - 1;
constexpr unsigned int EMPTY_NUMBER = NUMBERS_AMOUNT;
constexpr unsigned int YS_AMOUNT = M;
constexpr unsigned int XS_AMOUNT = 2 * M - 1;
constexpr unsigned int FUNCS_AMOUNT = XS_AMOUNT;


// ----------------------------------------------------------------------------
struct Number {
    unsigned int value;

    Number() noexcept : value(EMPTY_NUMBER) {}
    Number(unsigned int value) noexcept : value(value) {}

    inline bool empty() const noexcept { return value == EMPTY_NUMBER; };
    inline unsigned int operator()() const noexcept { return value; }
    inline bool operator==(const Number& other) const noexcept { return value == other(); };
    inline bool operator!=(const Number& other) const noexcept { return value != other(); };
};

inline Number XOR(const Number& a, const Number& b) {
    return a() ^ b();
}

std::ostream& operator<<(std::ostream& stream, const Number& number) {
    if (number.empty()) stream << "*"; else stream << number();
    return stream;
}
// ----------------------------------------------------------------------------
using VectorX = std::array<Number, XS_AMOUNT>;
using VectorY = std::array<Number, YS_AMOUNT>;

bool operator==(const VectorX& v1, const VectorX& v2) {
    for (unsigned int i = 0; i < XS_AMOUNT; i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
};

bool operator!=(const VectorX& v1, const VectorX& v2) {
    return !(operator==(v1, v2));
}

std::ostream& operator<<(std::ostream& stream, const VectorX& vector) {
    if (vector.size() == 0) {
        return stream;
    }

    stream << "{";
    for (unsigned int i = 0; i < vector.size() - 1; i++) {
        stream << vector[i] << ",";
    }
    stream << vector[vector.size() - 1];
    stream << "}";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const VectorY& vector) {
    if (vector.size() == 0) {
        return stream;
    }

    stream << "{";
    for (unsigned int i = 0; i < vector.size() - 1; i++) {
        stream << vector[i] << ",";
    }
    stream << vector[vector.size() - 1];
    stream << "}";

    return stream;
}
// ----------------------------------------------------------------------------
unsigned int generateRandomUniformInt(unsigned int low, unsigned int high) {
    const auto range = high - low + 1;
    return rand() % range + low;
}

Number generateRandomNumber() {
    return generateRandomUniformInt(0U, NUMBER_MAX);
}

VectorX generateRandomVectorX() {
    VectorX vector;
    for (Number& number : vector) {
        number = generateRandomNumber();
    }

    return vector;
}
VectorY generateRandomVectorY() {
    VectorY vector;
    for (Number& number : vector) {
        number = generateRandomNumber();
    }

    return vector;
}
// ----------------------------------------------------------------------------
struct Func {
    public:
        static constexpr unsigned int SIZE = NUMBERS_AMOUNT;
    private:
        std::array<Number, NUMBERS_AMOUNT> map;
        unsigned int emptyAmount;

    public:
        Func() : emptyAmount(SIZE) {
            for (unsigned int i = 0; i < map.size(); i++) {
                map[i] = Number(EMPTY_NUMBER);
            }
        }

        void set(const Number& arg, const Number& value) {
            emptyAmount += value.empty() ? 1 : -1;
            map[arg()] = value;;
        }
        const Number& at(const Number& arg) const {
            return map[arg()];
        }

        bool emptyAt(const Number& arg) const {
            return (map[arg()].empty());
        }

        unsigned int getEmptyAmount() const {
            return emptyAmount;
        }
};

using Funcs = std::array<Func*, FUNCS_AMOUNT>;

std::ostream& operator<<(std::ostream& stream, const Funcs& funcs) {
    stream << "Funcs:" << std::endl;
    for (unsigned int i = 0; i < Func::SIZE; i++) {
        for (const Func* func : funcs) {
            stream << i << " -> " << func->at(i) << "\t\t";
        }
        stream << std::endl;
    }

    return stream;
}
// ----------------------------------------------------------------------------


int main() {
    srand(410);
    std::cout << "--------------------BEGIN----------------------" << std::endl << std::endl;
    Funcs funcs;
    for (unsigned int i = 0; i < funcs.size(); i++) funcs[i] = new Func();

    const VectorY y = generateRandomVectorY();
    std::cout << ":> Generated Y = " << y << std::endl;
    std::vector<VectorX> xs;

    const auto allFuncsSufficient = [](const std::array<Func*, FUNCS_AMOUNT> funcs){
        for (unsigned int j = 0; j < funcs.size(); j++) {
            if (funcs[j]->getEmptyAmount() < M - 1) {
                return false;
            }
        }
        return true;
    };
    const auto pickRandomEmpty = [](const Func* func, const Number& forbiddenIndex){
        const unsigned int targetIndex = generateRandomUniformInt(0, func->getEmptyAmount() - 1);
        unsigned int iterator = 0;
        for (unsigned int k = 0; ; k++) {
            if (k == Func::SIZE) k = 0; // allow looping to ensure that we find targetIndex-s element
            if (func->at(k).empty() && Number(k) != forbiddenIndex) {
                if (iterator == targetIndex) {
                    return Number(k);
                } else {
                    iterator++;
                }
            }
        }
    };
    const auto pickRandomNotEmpty = [](const Func* func){
        if (func->getEmptyAmount() == Func::SIZE) {
            return Number(EMPTY_NUMBER);
        }

        const unsigned int targetIndex = generateRandomUniformInt(0, Func::SIZE - func->getEmptyAmount() - 1);
        unsigned int iterator = 0;
        for (unsigned int k = 0; ; k++) {
            if (!func->at(k).empty()) {
                if (iterator == targetIndex) {
                    return Number(k);
                } else {
                    iterator++;
                }
            }
        }
    };
    const auto pickRandomTarget = [](const Func* func, const Number& target){
        std::vector<unsigned int> indices;
        for (unsigned int k = 0; k < Func::SIZE; k++) {
            if (func->at(k) == target) {
                indices.push_back(k);
            }
        }
        return (indices.size() == 0) ?
            Number(EMPTY_NUMBER) :
            Number(indices[generateRandomUniformInt(0, indices.size() - 1)]);
    };
    do {
        std::cout << "Funcs sufficient " << xs.size() << std::endl;

        // Step 2
        unsigned int i = M - 1; // 1-based indexing

        // Step 3
        std::array<Number, M+1> v_last; // input into last row
        for (unsigned int j = 0; j < v_last.size(); j++) {
            v_last[j] = pickRandomEmpty(funcs[j], EMPTY_NUMBER);
            assert(!v_last[j].empty() && "No empty elements in the last row");
        }

        std::vector<Number> old_v; // i
        for (unsigned int i = 0; i < v_last.size(); i++) old_v.push_back(v_last[i]); // initial values
        for (; i >= 2;) { // Step 4, 7
            i--; // Step 4
            // Step 5
            Number v0 = pickRandomNotEmpty(funcs[0]);
            if (v0.empty()) {
                v0 = pickRandomEmpty(funcs[0], v_last[0]);
                assert(!v0.empty() && "Failed to find empty element");
                funcs[0]->set(v0, generateRandomNumber()); // w0
            }

            // Step 6
            std::vector<Number> new_v; // i-1
            new_v.push_back(v0);
            for (unsigned int j = 1; j <= 2*M-i-1; j++) {
                const Number prevV = new_v[new_v.size() - 1];
                const Number prevW = funcs[j - 1]->at(prevV);
                const Number w = XOR(prevW, old_v[j - 1]);

                Number v = pickRandomTarget(funcs[j], w);
                if (v.empty()) { // then pick any from empty
                    v = pickRandomEmpty(funcs[j], v_last[j]);
                    assert(!v.empty() && "No empty element found");
                    funcs[j]->set(v, w);
                }
                new_v.push_back(v);
            }

            old_v = new_v;
        } // Step 7

        // Step 8
        VectorX x;
        for (unsigned int j = 0; j < XS_AMOUNT; j++) x[j] = old_v[j]; // old_v === new_v now

        // Step 9
        funcs[0]->set(v_last[0], generateRandomNumber());

        // Step 10
        for (unsigned int j = 1; j < M + 1; j++) {
            const Number prevW = funcs[j - 1]->at(v_last[j - 1]);
            const Number w = XOR(XOR(x[j - 1], y[j - 1]), prevW);
            funcs[j]->set(v_last[j], w);
        }

        // Step 11 is integrated into Func logic
        xs.push_back(x);
    } while (allFuncsSufficient(funcs)); // Step 12

    std::cout << ":> Done. Funcs insufficient!" << std::endl;
    std::cout << "Generated " << xs.size() << " vectors" << std::endl;


    std::cout << std::endl << "---------------------END-----------------------" << std::endl;
    return 0;
}
