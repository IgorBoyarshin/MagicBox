#include <iostream>
#include <array>
#include <utility>
#include <variant>
#include <functional>
#include <vector>
#include <set>
#include <string>
#include <random>
#include <stack>
#include <cassert>


static const unsigned int NUMBERS_IN_VECTOR = 3;
static const unsigned int NUMBER_BITS = 3;
static const unsigned int NUMBER_MAX = (1 << NUMBER_BITS) - 1;
static const unsigned int EMPTY_NUMBER = NUMBER_MAX + 1;
static const unsigned int FUNCS_AMOUNT = 2 * NUMBERS_IN_VECTOR - 1;


typedef unsigned int Number;
typedef std::array<Number, NUMBERS_IN_VECTOR> Vector;


struct Coord1 {
    public:
        unsigned int x;

        Coord1(unsigned int x) : x(x) {}
};
struct Coord2 {
    public:
        unsigned int x;
        unsigned int y;

        Coord2(unsigned int x, unsigned int y) : x(x), y(y) {}
};
struct Coord3 {
    public:
        unsigned int x;
        unsigned int y;
        unsigned int z;

        Coord3(unsigned int x, unsigned int y, unsigned int z) : x(x), y(y), z(z) {}
};


// [low, high]
template<typename T = unsigned int>
T getRandomUniformInt(T low, T high) {
    static std::random_device rd;
    static std::mt19937 e2(rd());
    static std::uniform_int_distribution<T> dist(low, high);

    return dist(e2);
}


Number generateRandomNumber() {
    return getRandomUniformInt((unsigned int)0, NUMBER_MAX);
}

Vector generateRandomVector() {
    Vector vector;
    for (Number& number : vector) {
        number = generateRandomNumber();
    }

    return vector;
}
// ----------------------------------------------------------------------------
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

        Number apply(const Number& arg) const {
            return map[arg];
        }

        bool emptyAt(const Number& number) const {
            return (map[number] == EMPTY_NUMBER);
        }

        Number& operator[](const Number& index) {
            return map[index];
        }

        const Number& operator[](const Number& index) const {
            return map[index];
        }
};

bool emptyNumber(const Number& number) {
    return (number == EMPTY_NUMBER);
}

// ----------------------------------------------------------------------------
// (row, column) of Func
typedef std::array<Func, FUNCS_AMOUNT> Functions;
// ----------------------------------------------------------------------------
Number XOR(const Number& a, const Number& b) {
    return a ^ b;
}
// ----------------------------------------------------------------------------
struct Cell {
    public:
        Number value;

        Cell() : value(EMPTY_NUMBER) {}

        Cell(Number value) : value(value) {}

        bool empty() const {
            return value == EMPTY_NUMBER;
        }

        // TODO: check that works as expected
        // Cell& operator=(Number& newValue) {
        //     value = newValue;
        //     return *this;
        // }
};

struct Block {
    private:
        Func& func;

        Number& operator[](unsigned int index) {
            return at(index);
        }

        Number& at(unsigned int index) {
            if (index == 0) {
                return in;
            } else if (index == 1) {
                return out;
            } else {
                std::cout << ":> [Block.operator[]()]: Unsupported index used: "
                    << index << "." << std::endl;
            }
        }


    public:
        Number in;
        Number out;

        Block(Func& func) : func(func), in(EMPTY_NUMBER), out(EMPTY_NUMBER) {}

        bool empty() const {
            return (in == EMPTY_NUMBER) || (out == EMPTY_NUMBER);
        }

        bool full() const {
            return (in != EMPTY_NUMBER) && (out != EMPTY_NUMBER);
        }

        const Number& operator[](unsigned int index) const {
            if (index == 0) {
                return in;
            } else if (index == 1) {
                return out;
            } else {
                std::cout << ":> [Block.operator[]()]: Unsupported index used: "
                    << index << "." << std::endl;
            }
        }

        // Allow modifications only throught here
        void update(unsigned int index, const Number& number) {
            if (full()) { // the only case when func is already set
                func[in] = EMPTY_NUMBER; // unset
            }
            at(index) = number;
            if (full()) {
                func[in] = number;
            }
        }
};

class Snapshot {
    public:
        // Real physical Funcs
        std::array<Func, FUNCS_AMOUNT>& funcs;

        // Abstract snapshot Funcs
        std::vector<std::vector<Block>> rows;
        std::array<Cell, NUMBERS_IN_VECTOR> buffers;
        std::array<Cell, NUMBERS_IN_VECTOR> xs;
        const std::array<Cell, NUMBERS_IN_VECTOR> ys;

    private:
        void init() {
            const unsigned int AMOUNT_OF_ROWS = NUMBERS_IN_VECTOR - 1;
            rows.reserve(AMOUNT_OF_ROWS);
            for (unsigned int rowIndex = 0; rowIndex < AMOUNT_OF_ROWS; rowIndex++) {
                const unsigned int AMOUNT_OF_BLOCKS_IN_ROW = (2 * NUMBERS_IN_VECTOR - 1) - rowIndex;
                std::vector<Block> row;
                row.reserve(AMOUNT_OF_BLOCKS_IN_ROW);

                for (unsigned int blockIndex = 0; blockIndex < AMOUNT_OF_BLOCKS_IN_ROW; blockIndex++) {
                    row.emplace_back(funcs[blockIndex]);
                }

                rows.push_back(row);
            }
        }

        static std::array<Cell, NUMBERS_IN_VECTOR> convertVectorToCells(const Vector& vector) {
            std::array<Cell, NUMBERS_IN_VECTOR> cells;
            for (unsigned int i = 0; i < vector.size(); i++) {
                const Number& number = vector[i];
                cells[i] = Cell(number);
            }

            return cells;
        }

    public:
        Snapshot(
                std::array<Func, FUNCS_AMOUNT>& funcs,
                const Vector& ys)
                : funcs(funcs), ys(convertVectorToCells(ys)) {
            init();
        }

        bool existsEmpty() const {
            for (const auto& row : rows) {
                for (const Block& block : row) {
                    if (block.empty()) {
                        return true;
                    }
                }
            }
            for (const Cell& cell : buffers) {
                if (cell.empty()) {
                    return true;
                }
            }
            for (const Cell& x : xs) {
                if (x.empty()) {
                    return true;
                }
            }

            return false;
        }
};
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
enum LocationType {
    LocationType_X, LocationType_Buffer, LocationType_Func
};

struct Location {
    public:
        const LocationType type;
        const std::variant<Coord1, Coord3> coord;

        Location(LocationType type, const std::variant<Coord1, Coord3>& coord)
            : type(type), coord(coord) {}

        void place(const Number& number, Snapshot& snapshot) const {
            switch (type) {
                case LocationType_X:
                    assert(std::holds_alternative<Coord1>(coord));
                    snapshot.xs[std::get<Coord1>(coord).x] = number;
                    break;
                case LocationType_Buffer:
                    assert(std::holds_alternative<Coord1>(coord));
                    snapshot.buffers[std::get<Coord1>(coord).x] = number;
                    break;
                case LocationType_Func: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const auto& coordinates = std::get<Coord3>(coord);
                        snapshot.rows[coordinates.x][coordinates.y].update(coordinates.z, number);
                    }
                    break;
                default:
                    /* std::cout << ":> [Location.place()]: Usage of unhandled enum." << std::endl; */
                    // TODO: check that works as expected
                    assert(false && ":> [Location.place()]: Usage of unhandled enum.");
            }
        }
};
// ----------------------------------------------------------------------------
// Unconditional
struct Action {
    private:
        // TODO: as of now, every Action holds the same redundant field. Fix.
        Snapshot& snapshot;

        const Location location;
        const Number value;

        void undo() {
            location.place(EMPTY_NUMBER, snapshot);
        }

    public:
        Action(
                Snapshot& snapshot,
                const Location& location,
                const Number& value)
            : snapshot(snapshot), location(location), value(value) {
                location.place(value, snapshot);
            }

        ~Action() {
            undo();
        }

        const Number& getValue() const {
            return value;
        }
};


// Conditional, random
struct Step {
    private:
        static const unsigned int ALLOWED_RETRIES = 4;

        // TODO: make const Location
        Location baseLocation;
        std::vector<Action> actions;
        std::set<Number> triedNumbers;

    public:
        Step(const Location& location) : baseLocation(location) {}

        /* Step(const Action& baseAction) : baseAction(baseAction) { */
        /*     triedNumbers.insert(baseAction.getValue()); */
        /* } */

        /* Step(const Action& baseAction, const std::vector<Action>& actions) : baseAction(baseAction), actions(actions) { */
        /*     triedNumbers.insert(baseAction.getValue()); */
        /* } */

        // inline unsigned int amountOfTries() const {
        //     return usedOptions.size();
        // }
        
        void apply(Snapshot& snapshot, const Number& number) {
            triedNumbers.insert(number);
            actions.clear();
            actions.emplace_back(snapshot, baseLocation, number);
        }

        void applyActions(const std::vector<Action>& newActions) {
        // HERE
            actions.insert(actions.end(), newActions.begin(), newActions.end());
        }

        void fail() {
            actions.clear();
        }

        bool belowLimit() {
            return (triedNumbers.size() < ALLOWED_RETRIES);
        }

        std::optional<Number> getUntriedNumber() const {
            for (unsigned int i = 0; i <= NUMBER_MAX; i++) {
                const Number number(i);
                if (triedNumbers.find(number) == triedNumbers.end()) {
                    // This number hasn't beed tried
                    return {number};
                }
            }

            return std::nullopt;
        }
};
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
class Stack {
    private:
        static const unsigned int ALLOWED_RETRIES = 4;
        std::stack<Step> steps;
        std::stack<unsigned int> stepsTries;

    public:
        Stack() {}

        void push(Step step) {
            steps.push(step);
            stepsTries.push(0);
        }

        void emplace(const Location& location) {
            steps.emplace(location);
            stepsTries.push(0);
        }

        std::optional<Step> pop() {
            if (steps.size() == 0) {
                return std::nullopt;
            }

            const Step step = steps.top();
            steps.pop();
            stepsTries.pop();

            return {step};
        }

        // HERE
        /* std::optional<Step&> top() { */
        /*     if (steps.size() > 0) { */
        /*         return {steps.top()}; */
        /*     } */
        /*  */
        /*     return std::nullopt; */
        /*  */
        /*     // TODO: why doesn't work??? */
        /*     #<{(| return (steps.size() > 0) ? ({steps.top()}) : std::nullopt; |)}># */
        /* } */

        const Step& top() const {
            return steps.top();
        }

        bool belowLimit() {
            return (stepsTries.top() < ALLOWED_RETRIES);
        }

        // bool retry() {
        //     if (++stepsTries.top() < ALLOWED_RETRIES) {
        //
        //     } else {
        //         stack.pop();
        //         stepsTries.pop();
        //         return retry();
        //     }
        // }

        void applyActions(const std::vector<Action>& actions) {
            steps.top().applyActions(actions);
        }

};
// ----------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& stream, const Vector& vector) {
    if (vector.size() == 0) {
        return stream;
    }

    stream << "(";
    for (unsigned int i = 0; i < vector.size() - 1; i++) {
        stream << vector[i] << ",";
    }
    stream << vector[vector.size() - 1];
    stream << ")";

    return stream;
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
typedef std::function<std::optional<Location>(unsigned int index, const Snapshot& snapshot, const Functions& funcs)> StrategyFunction;

class Strategy {
    private:
        StrategyFunction function;
        unsigned int nextIndex;
    public:
        Strategy (StrategyFunction function) : function(function), nextIndex(0) {}

        std::optional<Location> get(unsigned int index, const Snapshot& snapshot, const Functions& funcs) {
            nextIndex = index + 1;
            return function(index, snapshot, funcs);
        }

        std::optional<Location> getNext(const Snapshot& snapshot, const Functions& funcs) {
            return get(nextIndex, snapshot, funcs);
        }

        bool revert() {
            // this strategy ain't gonna be used anymore if nextIndex was 0
            // upon call => don't care about the value anymore.
            /* return (nextIndex-- > 0); */

            // This solution doesn't spoil nextIndex
            if (nextIndex > 0) {
                nextIndex--;
                return true;
            }

            return false;
        }
};
// ----------------------------------------------------------------------------
std::optional<Coord2> findEmptyFunc(const Snapshot& snapshot) {
    unsigned int rowIndex = 0;
    unsigned int columnIndex = 0;
    for (const auto& row : snapshot.rows) {
        for (const Block& block : row) {
            if (block.empty()) {
                return std::optional<Coord2>{{rowIndex, columnIndex}};
            }
            columnIndex++;
        }
        rowIndex++;
    }

    return std::nullopt;
}
// ----------------------------------------------------------------------------
class MagicBox {
public:
    Functions funcs;

    Func& operator[](unsigned int index) {
        return funcs[index];
    }

    const Func& operator[](unsigned int index) const {
        return funcs[index];
    }

public:
    MagicBox() {}


    void construct() {

    }

    void work() {
        const Vector y = generateRandomVector();
        std::vector<Vector> xs;
        while (generateNewX(xs, y));
        outputResult(xs, y);
    }

    bool generateNewX(std::vector<Vector>& xs, const Vector& y) {
        Snapshot snapshot(funcs, y);
        Stack stepsStack;
        Strategy strategy([](
                    __attribute__((unused)) unsigned int index,
                   __attribute__((unused))  const Snapshot& snapshot,
                   __attribute__((unused)) const Functions& funcs) -> std::optional<Location> {
            if (index < NUMBERS_IN_VECTOR) { // first three insertions
                return {{LocationType_X, index}};
            } else {
                if (const std::optional<Coord2> emptyFunc = findEmptyFunc(snapshot)) {
                    const Coord2& funcCoord = *emptyFunc;
                    const unsigned int inFieldCoord = 0;
                    const Coord3 location = Coord3(funcCoord.x, funcCoord.y, inFieldCoord);
                    return {{LocationType_Func, location}};
                } else {
                    return std::nullopt;
                }
            }
        });


        /* while (const std::optional<Location> locationOpt = strategy.getNext(snapshot, funcs)) { */
        /*     const Location location = *locationOpt; */
        /*  */
        /*     // stepsStack.emplace(location); */
        /*     // std::optional<Step&> stepOpt = stepsStack.top(); */
        /*     std::optional<Step> stepOpt{location}; */
        /*     do { */
        /*         Step& step = *stepOpt; */
        /*         const std::optional<Number> numberOpt = step.getUntriedNumber(); */
        /*         const std::optional<std::vector<Action>> actionsOpt = */
        /*             (numberOpt) ? (step.apply(snapshot, *numberOpt), distributePokeFrom(snapshot, funcs, location)) */
        /*                         : (std::nullopt); */
        /*         if (numberOpt && actionsOpt) { */
        /*             step.applyActions(*actionsOpt); */
        /*             stepsStack.push(step); */
        /*  */
        /*             stepOpt = std::nullopt; */
        /*         } else { */
        /*             step.fail(); // inc tries */
        /*             while (stepOpt && !stepsStack.belowLimit()) { */
        /*                 #<{(| stepOpt = stepsStack.pop(); // std::nullopt if stack became empty |)}># */
        /*                 if (stepOpt) { */
        /*                     stepOpt->fail(); // child failed => he failed */
        /*                 } */
        /*                 strategy.revert(); */
        /*             } */
        /*  */
        /*             // By this line we've either reverted to a valid Step on Stack */
        /*             // or made stepOpt == std::nullopt if the stack delpeted */
        /*             if (!stepOpt) { */
        /*                 return false; // failed to construct Vector<X> */
        /*             } */
        /*         } */
        /*     } while (stepOpt); */
        /* } */





        /* while (const std::optional<Location> locationOpt = strategy.getNext(snapshot, funcs)) { */
        /*     const Location location = *locationOpt; */
        /*  */
        /*     std::optional<Step> stepOpt{location}; */
        /*     do { */
        /*         Step& step = *stepOpt; */
        /*         const std::optional<Number> numberOpt = step.getUntriedNumber(); */
        /*         const std::optional<std::vector<Action>> actionsOpt = */
        /*             (numberOpt) ? (step.apply(*number), distributePokeFrom(snapshot, funcs, location)) */
        /*                         : (std::nullopt); */
        /*         if (numberOpt && actionsOpt) { */
        /*             step.applyActions(*actionsOpt); */
        /*             stepsStack.push(step); */
        /*  */
        /*             stepOpt = std::nullopt; */
        /*         } else { */
        /*             step.fail(); // inc tries */
        /*             while (stepOpt && !stepsStack.belowLimit()) { */
        /*                 stepOpt = stepsStack.pop(); // std::nullopt if stack became empty */
        /*                 if (stepOpt) { */
        /*                     stepOpt->fail(); // child failed => he failed */
        /*                 } */
        /*                 strategy.revert(); */
        /*             } */
        /*  */
        /*             // By this line we've either reverted to a valid Step on Stack */
        /*             // or made stepOpt == std::nullopt if the stack delpeted */
        /*             if (!stepOpt) { */
        /*                 return false; // failed to construct Vector<X> */
        /*             } */
        /*         } */
        /*     } while (stepOpt); */
        /* } */



        /* while (const std::optional<Location> locationOpt = strategy.getNext(snapshot, funcs)) { */
        /*     const Location location = *locationOpt; */
        /*  */
        /*     #<{(| TODO: Apply Number here |)}># */
        /*     // Either random number or deduced from context */
        /*     if (const std::optional<Number> numberOpt = stepsStack.top().getUntriedNumber()) { */
        /*         const Number number = *numberOpt; */
        /*         Step step(Action(snapshot, location, number)); */
        /*  */
        /*  */
        /*         if (const auto actionsOpt = distributePokeFrom(snapshot, funcs, location)) { */
        /*             stepsStack.applyActions(*actionsOpt); */
        /*             #<{(| stepsStack.push(Step(*actionsOpt)); // A successful step |)}># */
        /*         } else { */
        /*             // An unconditional conflict arose down the line as the result of */
        /*             // applying the Number at Location => revert back. */
        /*             std::cout << ":> Conflict here" << std::endl; */
        /*  */
        /*             stepsStack.retry(); */
        /*             strategy.revert(); */
        /*         } */
        /*     } else { */
        /*         // Current step depleted => revert */
        /*     } */
        /*  */
        /* } */




        return false;
    }


    // returns std::nullopt if a conflict arose
    std::optional<std::vector<Action>> distributePokeFrom(Snapshot& snapshot, Functions& funcs, const Location& location) {
        /* switch(type) { */
        /*     case : */
        /*         poke1 poke2 */
                // get from poke()
        // }

        /* return summary poke */
        return {};
    }


    // returns std::nullopt if a conflict arose
    std::optional<std::vector<Action>> poke(Snapshot& snapshot, Functions& funcs, const Location& location) {
        // if (!empty) {
        // return {{}};

        std::vector<Action> actions;

        /* switch (location.type) { */
        /*     case LocationType_X: { */
        /*             assert(std::holds_alternative<Coord1>(location.coord)); */
        /*             const unsigned int xIndex = std::get<Coord1>(location.coord).x; */
        /*  */
        /*             const Cell& cell = snapshot.xs[xIndex]; */
        /*             if (!cell.empty()) { */
        /*                 // assert(false && "[MagicBox.generateNewX()]: X Cell not empty."); // remove me later */
        /*                 std::cout << ":> [MagicBox.generateNewX()]: X Cell not empty at: " << xIndex << " (" << cell.value << ")." << std::endl; */
        /*                 return std::nullopt; */
        /*             } */
        /*             const Number numberToPut = generateRandomNumber(); // TODO: rand here?? */
        /*             snapshot.xs[xIndex].value = numberToPut; */
        /*  */
        /*             const std::optional<std::vector<Action>> pokeBuffer = poke(); */
        /*  */
        /*  */
        /*             const std::optional<std::vector<Action>> pokeFunc = poke(); */
        /*             if (xIndex < NUMBERS_IN_VECTOR - 1) { */
        /*                 const std::optional<std::vector<Action>> pokeFuncSecondary = poke(); */
        /*             } */
        /*         } */
        /*         break; */
        /*     case LocationType_Buffer: */
        /*         // TODO */
        /*         break; */
        /*     case LocationType_Func: */
        /*         // TODO */
        /*         break; */
        /* } */

        return {actions};
    }



    void outputResult(const std::vector<Vector>& xs, const Vector& y) {
        std::cout << ":> For Y " << y << " the following Xs were generated (" << xs.size() << "):" << std::endl;
        for (const Vector& x : xs) {
            std::cout << x << std::endl;
        }
    }


    /* Vector apply(const Vector& input) const { */
    /*     const unsigned int LAYER_SIZE_START = FUNCS_AMOUNT; */
    /*     Number layer[LAYER_SIZE_START]; */
    /*     for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) layer[i] = input[i]; */
    /*     for (unsigned int i = 0; i < NUMBERS_IN_VECTOR - 1; i++) layer[NUMBERS_IN_VECTOR + i] = input[i]; */
    /*  */
    /*     for (unsigned int layerSize = LAYER_SIZE_START; layerSize > NUMBERS_IN_VECTOR; layerSize--) { */
    /*         for (unsigned int funcIndex = 0; funcIndex < layerSize; funcIndex++) { */
    /*             Number& number = layer[funcIndex]; */
    /*             if (funcs[funcIndex].emptyAt(number)) { */
    /*                 std::cout << ":> Usage of empty func entry." << std::endl; */
    /*             } */
    /*  */
    /*             number = funcs[funcIndex].apply(number); */
    /*         } */
    /*         for (unsigned int funcIndex = 0; funcIndex < layerSize - 1; funcIndex++) { */
    /*             layer[funcIndex] = layer[funcIndex] ^ layer[funcIndex + 1]; */
    /*         } */
    /*     } */
    /*  */
    /*     for (unsigned int funcIndex = 0; funcIndex < NUMBERS_IN_VECTOR; funcIndex++) { */
    /*         layer[funcIndex] = layer[funcIndex] ^ input[funcIndex]; */
    /*     } */
    /*  */
    /*     Vector output; */
    /*     for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) output[i] = layer[i]; */
    /*  */
    /*     return output; */
    /* } */
};
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
/* void fillRemainingRandom(MagicBox& magicBox) { */
/*     for (Func& func : magicBox.funcs) { */
/*         for (Number& number : func.map) { */
/*             if (Func::emptyNumber(number)) { */
/*                 number = getRandomUniformInt(static_cast<unsigned int>(0), Func::SIZE - 1); */
/*             } */
/*         } */
/*     } */
/* } */


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
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
int main() {
    std::cout << "--------------------BEGIN----------------------" << std::endl;

    MagicBox mb;
    // mb.work();
    /* mb[0][1] = 7; */
    /* mb[0][2] = 4; */
    /* mb[1][1] = 7; */
    /* mb[1][2] = 5; */
    /* mb[2][2] = 6; */
    /* mb[2][3] = 4; */
    /* mb[3][1] = 6; */
    /* mb[3][3] = 0; */
    /* mb[4][2] = 5; */
    /*  */
    /* Vector n {1, 2, 3}; */
    /* for (Number nn : mb.apply(n)) { */
    /*     std::cout << nn << " "; */
    /* } */
    /* std::cout << std::endl; */

    // mb[1][2] = 3;
    // fillRandomRemaining(mb);
    // print(mb);


    std::cout << "---------------------END-----------------------" << std::endl;
    return 0;
}
