#include <iostream>
#include <array>
#include <string>
#include <random>


static const unsigned int NUM_BITS = 3;
static const unsigned int NUMS_AMOUNT = 3;


// [low, high]
template<typename T = unsigned int>
T getRandomUniformInt(T low, T high) {
    static std::random_device rd;
    static std::mt19937 e2(rd());
    static std::uniform_int_distribution<T> dist(low, high);

    return dist(e2);
}


typedef std::array<unsigned int, NUMS_AMOUNT> Number;


struct Func {
    public:
        static const unsigned int SIZE = 1 << NUM_BITS;
        // num_in -> num_out
        unsigned int map[SIZE];
    private:
        static const unsigned int EMPTY_ELEMENT = SIZE;

    public:
        Func() {
            for (unsigned int i = 0; i < SIZE; i++) map[i] = EMPTY_ELEMENT;
        }

        bool emptyAt(unsigned int index) const {
            return (map[index] == EMPTY_ELEMENT);
        }

        static bool emptyElement(unsigned int elem) {
            return (elem == EMPTY_ELEMENT);
        }

        unsigned int& operator[](unsigned int index) {
            return map[index];
        }

        const unsigned int& operator[](unsigned int index) const {
            return map[index];
        }
};


class MagicBox {
    public:
        static const unsigned int FUNCS_AMOUNT = 2 * NUMS_AMOUNT - 1;

        std::array<Func, FUNCS_AMOUNT> funcs;

    public:
        MagicBox() {

        }
};


void fillRandomRemaining(MagicBox& magicBox) {
    for (Func& func : magicBox.funcs) {
        for (unsigned int& element : func.map) {
            if (Func::emptyElement(element)) {
                element = getRandomUniformInt(static_cast<unsigned int>(0), Func::SIZE - 1);
            }
        }
    }
}


void print(const MagicBox magicBox) {
    for (unsigned int i = 0; i < Func::SIZE; i++) {
        for (const Func& func : magicBox.funcs) {
            std::cout << i << " -> ";
            const std::string element = (func.emptyAt(i) ? std::string("*") : std::to_string(func.map[i]));
            std::cout << element << "\t\t";
        }
        std::cout << std::endl;
    }
}


int main() {
    std::cout << "--------------------BEGIN----------------------" << std::endl;

    MagicBox mb;
    fillRandomRemaining(mb);
    print(mb);


    std::cout << "---------------------END-----------------------" << std::endl;
    return 0;
}
