#include <iostream>
#include <array>
#include <string>
#include <random>


static const unsigned int NUMBER_BITS = 3;
static const unsigned int NUMBER_MAX = (1 << NUMBER_BITS) - 1;
static const unsigned int NUMBERS_IN_VECTOR = 3;


// [low, high]
template<typename T = unsigned int>
T getRandomUniformInt(T low, T high) {
    static std::random_device rd;
    static std::mt19937 e2(rd());
    static std::uniform_int_distribution<T> dist(low, high);

    return dist(e2);
}


typedef unsigned int Number;
typedef std::array<Number, NUMBERS_IN_VECTOR> Vector;
// ----------------------------------------------------------------------------
struct Func {
    public:
        static const unsigned int SIZE = NUMBER_MAX + 1;
        // num_in -> num_out
        Number map[SIZE];
    private:
        static const Number EMPTY_NUMBER = SIZE;

    public:
        Func() {
            for (unsigned int i = 0; i < SIZE; i++) map[i] = EMPTY_NUMBER;
        }

        Number apply(Number arg) const {
            return map[arg];
        }

        bool emptyAt(Number number) const {
            return (map[number] == EMPTY_NUMBER);
        }

        static bool emptyNumber(Number number) {
            return (number == EMPTY_NUMBER);
        }

        Number& operator[](Number index) {
            return map[index];
        }

        const Number& operator[](Number index) const {
            return map[index];
        }
};
// ----------------------------------------------------------------------------
class MagicBox {
    public:
        static const unsigned int FUNCS_AMOUNT = 2 * NUMBERS_IN_VECTOR - 1;

        std::array<Func, FUNCS_AMOUNT> funcs;

        Func& operator[](unsigned int index) {
            return funcs[index];
        }

        const Func& operator[](unsigned int index) const {
            return funcs[index];
        }

    public:
        MagicBox() {}

        Vector apply(const Vector& input) const {
            const unsigned int LAYER_SIZE_START = FUNCS_AMOUNT;
            Number layer[LAYER_SIZE_START];
            for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) layer[i] = input[i];
            for (unsigned int i = 0; i < NUMBERS_IN_VECTOR - 1; i++) layer[NUMBERS_IN_VECTOR + i] = input[i];

            for (unsigned int layerSize = LAYER_SIZE_START; layerSize > NUMBERS_IN_VECTOR; layerSize--) {
                for (unsigned int funcIndex = 0; funcIndex < layerSize; funcIndex++) {
                    Number& number = layer[funcIndex];
                    if (funcs[funcIndex].emptyAt(number)) {
                        std::cout << ":> Usage of empty func entry." << std::endl;
                    }

                    number = funcs[funcIndex].apply(number);
                }
                for (unsigned int funcIndex = 0; funcIndex < layerSize - 1; funcIndex++) {
                    layer[funcIndex] = layer[funcIndex] ^ layer[funcIndex + 1];
                }
            }

            for (unsigned int funcIndex = 0; funcIndex < NUMBERS_IN_VECTOR; funcIndex++) {
                layer[funcIndex] = layer[funcIndex] ^ input[funcIndex];
            }

            Vector output;
            for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) output[i] = layer[i];

            return output;
        }
};
// ----------------------------------------------------------------------------
void fillRemainingRandom(MagicBox& magicBox) {
    for (Func& func : magicBox.funcs) {
        for (Number& number : func.map) {
            if (Func::emptyNumber(number)) {
                number = getRandomUniformInt(static_cast<unsigned int>(0), Func::SIZE - 1);
            }
        }
    }
}


void print(const MagicBox& magicBox) {
    for (unsigned int i = 0; i < Func::SIZE; i++) {
        for (const Func& func : magicBox.funcs) {
            std::cout << i << " -> ";
            const std::string element = (func.emptyAt(i) ?
                    std::string("_") : std::to_string(func.map[i]));
            std::cout << element << "\t\t";
        }
        std::cout << std::endl;
    }
}


Vector generateRandomVector() {
    Vector vector;
    for (Number& number : vector) {
        number = getRandomUniformInt((unsigned int)0, NUMBER_MAX);
    }

    return vector;
}
// ----------------------------------------------------------------------------
/*
 * Input: empty MB
 * Output: filled MB, returns the target Y
 * TODO
 */
/* Vector algorithm1(MagicBox& magicBox) { */
/*     const Vector Y = generateRandomVector(); */
/*  */
/*     while (true) { */
/*         const Vector X = generateRandomVector(); */
/*         unsigned int */
/*     } */
/*  */
/*     return Y; */
/* } */
// ----------------------------------------------------------------------------
int main() {
    std::cout << "--------------------BEGIN----------------------" << std::endl;

    MagicBox mb;
    mb[0][1] = 7;
    mb[0][2] = 4;
    mb[1][1] = 7;
    mb[1][2] = 5;
    mb[2][2] = 6;
    mb[2][3] = 4;
    mb[3][1] = 6;
    mb[3][3] = 0;
    mb[4][2] = 5;

    Vector n {1, 2, 3};
    for (Number nn : mb.apply(n)) {
        std::cout << nn << " ";
    }
    std::cout << std::endl;

    // mb[1][2] = 3;
    // fillRandomRemaining(mb);
    print(mb);


    std::cout << "---------------------END-----------------------" << std::endl;
    return 0;
}
