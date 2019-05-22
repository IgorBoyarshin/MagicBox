#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <cassert>
#include <algorithm>


constexpr static unsigned int K = 16; // bitness of fragment. Value = 0..2^k
constexpr static unsigned int M = 7; // amount of fragments (Y)
constexpr static unsigned int NUMBERS_AMOUNT = 1 << K;
constexpr static unsigned int NUMBER_MAX = NUMBERS_AMOUNT - 1;
constexpr static unsigned int EMPTY_NUMBER = NUMBERS_AMOUNT;
constexpr static unsigned int YS_AMOUNT = M;
constexpr static unsigned int XS_AMOUNT = 2 * M - 1;
constexpr static unsigned int FUNCS_AMOUNT = XS_AMOUNT;
constexpr static bool         DO_CHECK_UNIQUENESS = false;
constexpr static bool         VERBOSE = false;
// ----------------------------------------------------------------------------
struct Number {
    unsigned int value;

    constexpr explicit Number()                   noexcept : value(EMPTY_NUMBER) {}
    constexpr explicit Number(unsigned int value) noexcept : value(value)        {}

    constexpr bool         empty()                         const noexcept { return value == EMPTY_NUMBER; };
    constexpr unsigned int operator()()                    const noexcept { return value; }
    constexpr bool         operator==(const Number& other) const noexcept { return value == other(); };
    constexpr bool         operator!=(const Number& other) const noexcept { return value != other(); };
};

constexpr Number XOR(const Number& a, const Number& b) noexcept {
    return Number(a() ^ b());
}

std::ostream& operator<<(std::ostream& stream, const Number& number) noexcept {
    if (number.empty()) stream << "*"; else stream << number();
    return stream;
}
// ----------------------------------------------------------------------------
using ArrayX = std::array<Number, XS_AMOUNT>;
using ArrayY = std::array<Number, YS_AMOUNT>;

constexpr bool operator==(const ArrayX& v1, const ArrayX& v2) noexcept {
    for (unsigned int i = 0; i < XS_AMOUNT; i++) {
        if (v1[i] != v2[i]) return false;
    }
    return true;
};

constexpr bool operator!=(const ArrayX& v1, const ArrayX& v2) noexcept {
    return !(operator==(v1, v2));
}

template <typename T>
constexpr bool unique(const T& x, const std::vector<T>& vs) noexcept {
    for (const T& v : vs) {
        if (x == v) return false;
    }
    return true;
}

constexpr ArrayY yFromX(const ArrayX& x) noexcept {
    ArrayY y;
    for (unsigned int i = 0; i < YS_AMOUNT; i++) y[i] = x[i];
    return y;
}

std::ostream& operator<<(std::ostream& stream, const ArrayX& vector) noexcept {
    if (vector.size() == 0) return stream;

    stream << "{";
    for (unsigned int i = 0; i < vector.size() - 1; i++) {
        stream << vector[i] << ",";
    }
    stream << vector[vector.size() - 1];
    stream << "}";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const ArrayY& vector) noexcept {
    if (vector.size() == 0) return stream;

    stream << "{";
    for (unsigned int i = 0; i < vector.size() - 1; i++) {
        stream << vector[i] << ",";
    }
    stream << vector[vector.size() - 1];
    stream << "}";

    return stream;
}
// ----------------------------------------------------------------------------
// Range: [low; high]
unsigned int generateRandomUniformInt(unsigned int low, unsigned int high) noexcept {
    const auto range = high - low + 1;
    return std::rand() % range + low;
}

Number generateRandomNumber() noexcept {
    return Number(generateRandomUniformInt(0U, NUMBER_MAX));
}

ArrayY generateRandomArrayY() noexcept {
    ArrayY vector;
    for (Number& number : vector) number = generateRandomNumber();
    return vector;
}
// ----------------------------------------------------------------------------
struct Func {
    public:
        static constexpr unsigned int SIZE = NUMBERS_AMOUNT;
    private:
        // std::array<> doesn't work on large numbers because initializes on stack
        // std::array<Number, NUMBERS_AMOUNT> map;
        Number* map;
        std::vector<Number>* cache;
        unsigned int emptyAmount;

    public:
        Func() : map  (new Number             [NUMBERS_AMOUNT]()), // default initialization
                 cache(new std::vector<Number>[NUMBERS_AMOUNT]()), // default initialization
                 emptyAmount(SIZE) {}

        ~Func() {
            delete[] map;
            delete[] cache;
        }

        void set(const Number& arg, const Number& value) noexcept {
            if (map[arg()].empty() != value.empty()) {
                emptyAmount += value.empty() ? 1 : -1;
            }

            if (value.empty()) {
                const auto key = map[arg()];
                auto& items = cache[key()];
                const auto index = std::find(items.cbegin(), items.cend(), arg);
                items.erase(index);
            } else {
                const auto key = value;
                cache[key()].push_back(arg);
            }

            map[arg()] = value;
        }

        constexpr Number at(unsigned int arg) const noexcept {
            return map[arg];
        }

        constexpr Number at(const Number& arg) const noexcept {
            return map[arg()];
        }

        constexpr bool emptyAt(const Number& arg) const noexcept {
            return (map[arg()].empty());
        }

        constexpr unsigned int getEmptyAmount() const noexcept {
            return emptyAmount;
        }

        constexpr const std::vector<Number>& cacheEntryFor(const unsigned int arg) const noexcept {
            return cache[arg];
        }
};

using Funcs = std::array<Func*, FUNCS_AMOUNT>;

std::ostream& operator<<(std::ostream& stream, const Funcs& funcs) noexcept {
    stream << "Funcs:" << std::endl;
    for (unsigned int i = 0; i < Func::SIZE; i++) {
        for (const Func* func : funcs) {
            stream << i << " -> " << func->at(i) << "\t\t";
        }
        stream << std::endl;
    }

    return stream;
}

ArrayY forward(const Funcs& funcs, const ArrayX& startX) noexcept {
    ArrayX x = startX;
    for (unsigned int width = XS_AMOUNT; width > YS_AMOUNT; width--) {
        for (unsigned int j = 0; j < width; j++) {
            x[j] = funcs[j]->at(x[j]);
        }
        for (unsigned int j = 0; j < width - 1; j++) {
            x[j] = XOR(x[j], x[j+1]);
        }
    }
    for (unsigned int j = 0; j < M; j++) {
        x[j] = XOR(x[j], startX[j]);
    }

    return yFromX(x);
}

bool hashesMatch(const Funcs& funcs, const std::vector<ArrayX>& xs, const ArrayY& y) noexcept {
    for (const ArrayX& x : xs) {
        if (forward(funcs, x) != y) return false;
    }
    return true;
}
// ----------------------------------------------------------------------------
int main() {
    srand(410);
    std::cout << "--------------------BEGIN----------------------" << std::endl << std::endl;
    std::cout << "Running for K=" << K << " and M=" << M << std::endl;
    Funcs funcs;
    for (auto& func : funcs) func = new Func(); // auto <=> Func*&

    const ArrayY y = generateRandomArrayY();
    if (VERBOSE) std::cout << ":> Generated Y = " << y << std::endl;
    std::vector<ArrayX> xs;

    // Each Func must be able to provide at least additional (M - 1) empty cells
    // for upcoming iteration
    const auto allFuncsSufficient = [](const Funcs& funcs) {
        return std::all_of(funcs.cbegin(), funcs.cend(),
            [](Func* func){ return func->getEmptyAmount() >= M - 1; });
    };
    const auto pickRandomEmpty = [](const Func* func, const Number& forbiddenNumber) {
        if (func->getEmptyAmount() == 0) return Number(EMPTY_NUMBER);
        const unsigned int startingIndex = generateRandomUniformInt(0, Func::SIZE - 1);
        for (unsigned int k = startingIndex;; k++) {
            if (k == Func::SIZE) k = 0;
            if (k == forbiddenNumber()) continue;
            if (func->at(k).empty()) return Number(k);
        }
    };
    const auto pickRandomNotEmpty = [](const Func* func){
        if (func->getEmptyAmount() == Func::SIZE) return Number(EMPTY_NUMBER);
        const unsigned int startingIndex = generateRandomUniformInt(0, Func::SIZE - 1);
        for (unsigned int k = startingIndex;; k++) {
            if (k == Func::SIZE) k = 0;
            if (!func->at(k).empty()) return Number(k);
        }
    };
    const auto pickRandomTarget = [](const Func* func, const Number& target){
        const auto& indices = func->cacheEntryFor(target());
        return indices.empty() ?
            Number(EMPTY_NUMBER) :
            indices[generateRandomUniformInt(0, indices.size() - 1)];
    };

    unsigned int duplicatesCount = 0;
    do {
        if (VERBOSE) std::cout << ":> Funcs sufficient ";

        // Step 2
        unsigned int i = M - 1; // 1-based indexing

        // Step 3
        std::array<Number, M+1> v_last; // input into last row
        for (unsigned int j = 0; j < v_last.size(); j++) {
            v_last[j] = pickRandomEmpty(funcs[j], Number(EMPTY_NUMBER));
            assert(!v_last[j].empty() && "No empty elements in the last row");
        }

        std::vector<Number> old_v; // i
        old_v.reserve(v_last.size());
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
            new_v.reserve(2 * M - i);
            new_v.push_back(v0);
            for (unsigned int j = 1; j <= 2*M-i-1; j++) {
                const Number prevV = new_v[j - 1];
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
        ArrayX x;
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
        if (DO_CHECK_UNIQUENESS && !unique(x, xs)) {
            std::cout << "--------- Duplicate found ----------" << std::endl;
            duplicatesCount++;
        } else {
            xs.push_back(x);
            std::cout << "--- New X generated (" << xs.size() << ") ---" << '\r';
        }
    } while (allFuncsSufficient(funcs)); // Step 12

    std::cout << std::endl;
    if (VERBOSE) std::cout << ":> Done. Funcs insufficient!" << std::endl;
    std::cout << "Generated " << xs.size() << " vectors" << std::endl;
    if (DO_CHECK_UNIQUENESS) std::cout << duplicatesCount << " duplicates discarded" << std::endl;

    // if (!hashesMatch(funcs, xs, y)) {
    //     std::cout << ">>> ERROR: Hash does not match!" << std::endl;
    // }

    std::cout << std::endl << "---------------------END-----------------------" << std::endl;
    return 0;
}
