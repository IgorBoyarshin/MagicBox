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
#include <stdlib.h>


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
// template<typename T = unsigned int>
// T getRandomUniformInt(T low, T high) {
//     static std::random_device rd;
//     static std::mt19937 e2(rd());
//     static std::uniform_int_distribution<T> dist(low, high);
//
//     return dist(e2);
// }


unsigned int getRandomUniformInt(unsigned int low, unsigned int high) {
    const auto range = high - low;
    return rand() % range + low;
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

static const unsigned int BLOCK_IN = 0;
static const unsigned int BLOCK_OUT = 1;
struct Block {
    private:
        Func& func;

        Number& operator[](unsigned int index) {
            return at(index);
        }

        Number& at(unsigned int index) {
            if (index == BLOCK_IN) {
                return in;
            } else if (index == BLOCK_OUT) {
                return out;
            } else {
                std::cout << ":> [Block.operator[]()]: Unsupported index used: "
                    << index << "." << std::endl;
            }
        }


        Number in;
        Number out;

    public:
        Block(Func& func) : func(func), in(EMPTY_NUMBER), out(EMPTY_NUMBER) {}

        bool existsEmpty() const {
            return !full();
        }

        bool empty() const {
            return (in == EMPTY_NUMBER) && (out == EMPTY_NUMBER);
        }

        bool emptyIn() const {
            return (in == EMPTY_NUMBER);
        }

        bool emptyOut() const {
            return (out == EMPTY_NUMBER);
        }

        bool full() const {
            return (in != EMPTY_NUMBER) && (out != EMPTY_NUMBER);
        }

        const Number& operator[](unsigned int index) const {
            if (index == BLOCK_IN) {
                return in;
            } else if (index == BLOCK_OUT) {
                return out;
            } else {
                std::cout << ":> [Block.operator[]()]: Unsupported index used: "
                    << index << "." << std::endl;
            }
        }

        const Number& In() const {
            return in;
        }

        const Number& Out() const {
            return out;
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
                    if (block.existsEmpty()) {
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
        LocationType type;
        std::variant<Coord1, Coord3> coord;

        Location(const LocationType& type, const std::variant<Coord1, Coord3>& coord)
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

        bool empty(Snapshot& snapshot) const {
            switch (type) {
                case LocationType_X:
                    assert(std::holds_alternative<Coord1>(coord));
                    return emptyNumber(snapshot.xs[std::get<Coord1>(coord).x].value);
                    break;
                case LocationType_Buffer:
                    assert(std::holds_alternative<Coord1>(coord));
                    return emptyNumber(snapshot.buffers[std::get<Coord1>(coord).x].value);
                    break;
                case LocationType_Func: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const auto& coordinates = std::get<Coord3>(coord);
                        const Block& block = snapshot.rows[coordinates.x][coordinates.y];
                        return emptyNumber(block[coordinates.z]);
                    }
                    break;
                default:
                    assert(false && ":> [Location.place()]: Usage of unhandled enum.");
            }
        }
};
// ----------------------------------------------------------------------------
// Unconditional
struct Action {
    private:
        // TODO: as of now, every Action holds the same redundant field. Fix.
        // Snapshot& snapshot;

        Location location;
        Number value;

    public:
        Action(
                const Location& location,
                const Number& value)
            : location(location), value(value) {}

        void undo(Snapshot& snapshot) {
            location.place(EMPTY_NUMBER, snapshot);
        }

        void apply(Snapshot& snapshot) {
            location.place(value, snapshot);
        }

        const Number& getValue() const {
            return value;
        }
};


// Conditional, random
struct Step {
    private:
        static const unsigned int ALLOWED_RETRIES = 4;

        Location baseLocation;
        std::vector<Action> actions;
        std::set<Number> triedNumbers;

    public:
        Step(const Location& location) : baseLocation(location) {}

        void apply(Snapshot& snapshot, const Number& number) {
            triedNumbers.insert(number);
            actions.clear();
            actions.emplace_back(baseLocation, number);
        }

        void applyActions(const std::vector<Action>& newActions) {
            actions.insert(actions.end(), newActions.begin(), newActions.end());
        }

        void fail() {
            actions.clear();
        }

        bool belowLimit() {
            return (triedNumbers.size() < ALLOWED_RETRIES);
        }

        void cleanup(Snapshot& snapshot) {
            for (Action& action : actions) {
                action.undo(snapshot);
            }
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

        std::optional<std::reference_wrapper<Step>> top() {
            if (steps.size() > 0) {
                return {steps.top()};
            }

            return std::nullopt;

            // TODO: why doesn't work???
            /* return (steps.size() > 0) ? ({steps.top()}) : std::nullopt; */
        }

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
            if (block.existsEmpty()) {
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


        while (const std::optional<Location> locationOpt = strategy.getNext(snapshot, funcs)) {
            const Location location = *locationOpt;

            std::optional<Step> stepOpt{location};
            do {
                Step& step = *stepOpt;
                const std::optional<Number> numberOpt = step.getUntriedNumber();
                const std::optional<std::vector<Action>> actionsOpt =
                    (numberOpt) ? (poke(snapshot, location, *numberOpt))
                                : (std::nullopt);
                /* const std::optional<std::vector<Action>> actionsOpt = */
                /*     (numberOpt) ? (step.apply(snapshot, *numberOpt), poke(snapshot, location)) */
                /*                 : (std::nullopt); */
                if (numberOpt && actionsOpt) {
                    step.applyActions(*actionsOpt);
                    stepsStack.push(step);

                    stepOpt = std::nullopt;
                } else {
                    step.fail(); // inc tries
                    while (stepOpt && !stepsStack.belowLimit()) {
                        stepOpt = stepsStack.pop(); // std::nullopt if stack became empty
                        if (stepOpt) {
                            stepOpt->fail(); // child failed => he failed
                        }
                        strategy.revert();
                    }

                    // By this line we've either reverted to a valid Step on Stack
                    // or made stepOpt == std::nullopt if the stack delpeted
                    if (!stepOpt) {
                        return false; // failed to construct Vector<X>
                    }
                }
            } while (stepOpt);
        }

        return true;
    }


    std::optional<std::vector<Action>> reactOnFuncUpdate(Snapshot& snapshot, const Coord2& coord) {
        Block& block = snapshot.rows[coord.x][coord.y];

        if (block.empty()) {
            return {std::vector<Action>{}}; // nothing new for this place
        }
        const unsigned int functionIndex = coord.x;
        if (block.full()) {
            if (snapshot.funcs[functionIndex].apply(block.In()) == block.Out()) {
                return {std::vector<Action>{}}; // complies, eveyrthing good
            } else {
                return std::nullopt; // doesn't comply
            }
        } else if (block.emptyOut()) {
            if (const Number funcRes = snapshot.funcs[functionIndex].apply(block.In()); !emptyNumber(funcRes)) {
                // can deduce Out from function => try placing Out
                std::vector<Action> actions;
                if (const bool successfullPoke = insertActions(snapshot, actions, LocationType_Func, Coord3(coord.x, coord.y, BLOCK_OUT), funcRes); !successfullPoke) return std::nullopt;
                return actions;
            }
        }

        return {std::vector<Action>{}}; // block.In is empty => nothing is determined yet
    }


    bool insertActions(Snapshot& snapshot, std::vector<Action>& actions, const LocationType& locationType, std::variant<Coord1, Coord3> coord, const Number& value) {
        if (const auto& actionsOpt = poke(snapshot, Location(locationType, coord), value)) {
            actions.insert(actions.end(), actionsOpt->begin(), actionsOpt->end());
            return true;
        }

        return false;
    }


    // returns std::nullopt if a conflict arose
    // returns {{}} (an ampty vector) if there was nothing to do
    // TODO: create lambda for extractiong Actions from pokes
    std::optional<std::vector<Action>> poke(Snapshot& snapshot, const Location& location, const Number& number) {
        std::vector<Action> actions;
        actions.emplace_back(location, number); // self
        actions[0].apply(snapshot);

        switch (location.type) {
            // is to poked set only from the outside
            case LocationType_X:
                if (const auto successfullPoke = pokeX(snapshot, actions, location, number);
                        !successfullPoke) return std::nullopt;
                break;
            case LocationType_Buffer:
                if (const auto successfullPoke = pokeBuffer(snapshot, actions, location, number);
                        !successfullPoke) return std::nullopt;
                break;
            case LocationType_Func: {
                    assert(std::holds_alternative<Coord3>(location.coord));
                    const Coord3 coord = std::get<Coord3>(location.coord);
                    switch (coord.z) {
                        case BLOCK_IN:
                            if (const auto successfullPoke = pokeFuncIn(snapshot, actions, location, number);
                                    !successfullPoke) return std::nullopt;
                            break;
                        case BLOCK_OUT:
                            if (const auto successfullPoke = pokeFuncOut(snapshot, actions, location, number);
                                    !successfullPoke) return std::nullopt;
                            break;
                    }
                }
                break;
        }

        return {actions};
    }


    bool pokeX(Snapshot& snapshot, std::vector<Action> actions, const Location& location, const Number& number) {
        assert(std::holds_alternative<Coord1>(location.coord));
        const unsigned int index = std::get<Coord1>(location.coord).x;

        // Main In
        assert(snapshot.rows[0][index].emptyIn() && ":> IN not empty when X got poked.");
        if (const bool successfullPoke = insertActions(snapshot, actions,
                    LocationType_Func,
                    Coord3(0, index, BLOCK_IN),
                    number);
                !successfullPoke) return false;

        // Secondary In
        if (index < NUMBERS_IN_VECTOR - 1) {
            assert(snapshot.rows[0][index + NUMBERS_IN_VECTOR].emptyIn() && ":> Secondary IN not empty when X got poked.");
            if (const bool successfullPoke = insertActions(snapshot, actions,
                        LocationType_Func,
                        Coord3(0, index + NUMBERS_IN_VECTOR, BLOCK_IN),
                        number);
                    !successfullPoke) return false;
        }

        // Buffer
        assert(snapshot.buffers[index].empty() && ":> Buffer not empty when X got poked.");
        if (const bool successfullPoke = insertActions(snapshot, actions,
                    LocationType_Buffer,
                    Coord1(index),
                    XOR(number, snapshot.ys[index].value));
                !successfullPoke) return false;

        return true;
    }

    bool pokeBuffer(Snapshot& snapshot, std::vector<Action> actions, const Location& location, const Number& number) {
        assert(std::holds_alternative<Coord1>(location.coord));
        const unsigned int index = std::get<Coord1>(location.coord).x;

        if (!snapshot.xs[index].empty() && (snapshot.xs[index].value != number)) {
                // The value we're trying to set to Buffer does not equal the current value of X
                return false;
        }
        snapshot.xs[index].value = number; // just set X, but never poke it

        // Main In
        assert(snapshot.rows[0][index].emptyIn() && ":> IN not empty when Buffer got poked and X was empty.");
        if (const bool successfullPoke = insertActions(snapshot, actions,
                    LocationType_Func,
                    Coord3(0, index, BLOCK_IN),
                    number);
                !successfullPoke) return false;

        // Secondary In
        if (index < NUMBERS_IN_VECTOR - 1) {
            assert(snapshot.rows[0][index + NUMBERS_IN_VECTOR].emptyIn() && ":> Secondary IN not empty when Buffer got poked and X was empty.");
            if (const bool successfullPoke = insertActions(snapshot, actions,
                        LocationType_Func,
                        Coord3(0, index + NUMBERS_IN_VECTOR, BLOCK_IN),
                        number);
                    !successfullPoke) return false;
        }

        // Outs
        const unsigned int lastRowIndex = snapshot.rows.size() - 1;
        const Block& block1Const = snapshot.rows[lastRowIndex][index];
        const Block& block2Const = snapshot.rows[lastRowIndex][index + 1];
        Block& block1 = snapshot.rows[lastRowIndex][index];
        Block& block2 = snapshot.rows[lastRowIndex][index + 1];
        if (block1.emptyOut() && block2.emptyOut()) {
            if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                return false;
            }
        }
        if (block1.emptyOut() && !block2.emptyOut()) {
            if (const bool successfullPoke = insertActions(snapshot, actions,
                        LocationType_Func,
                        Coord3(lastRowIndex, index, BLOCK_OUT),
                        XOR(number, block2.Out()));
                    !successfullPoke) return false;
        } else if (block1.emptyOut() && !block2.emptyOut()) {
            if (const bool successfullPoke = insertActions(snapshot, actions,
                        LocationType_Func,
                        Coord3(lastRowIndex, index + 1, BLOCK_OUT),
                        XOR(number, block1.Out()));
                    !successfullPoke) return false;
        }

        return true;
    }

    bool pokeFuncIn(Snapshot& snapshot, std::vector<Action> actions, const Location& location, const Number& number) {
        assert(std::holds_alternative<Coord3>(location.coord));
        const Coord3 coord = std::get<Coord3>(location.coord);
        const unsigned int rowIndex = coord.x;
        const unsigned int columnIndex = coord.y;

        if (rowIndex == 0) {
            const bool isMainIn = columnIndex < NUMBERS_IN_VECTOR;
            const unsigned int xIndex =
                (isMainIn) ? columnIndex
                           : (columnIndex - NUMBERS_IN_VECTOR);
            if (snapshot.xs[xIndex].empty()) {
                // Set X
                snapshot.xs[xIndex].value = number; // just set X, but never poke it

                // Set INs
                const unsigned int otherColumnIndex =
                    (isMainIn) ? (columnIndex + NUMBERS_IN_VECTOR)
                               : (columnIndex - NUMBERS_IN_VECTOR);
                if (const bool successfullPoke = insertActions(snapshot, actions,
                            LocationType_X,
                            Coord1(otherColumnIndex),
                            number);
                        !successfullPoke) return false;

                // Set Buffer
                if (const bool successfullPoke = insertActions(snapshot, actions,
                            LocationType_Buffer,
                            Coord1(columnIndex),
                            number);
                        !successfullPoke) return false;
            }
        } else {
            // Outs
            const Block& block1Const = snapshot.rows[rowIndex - 1][columnIndex];
            const Block& block2Const = snapshot.rows[rowIndex - 1][columnIndex + 1];
            Block& block1 = snapshot.rows[rowIndex - 1][columnIndex];
            Block& block2 = snapshot.rows[rowIndex - 1][columnIndex + 1];
            if (block1.emptyOut() && block2.emptyOut()) {
                if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                    return false;
                }
            }
            if (block1.emptyOut() && !block2.emptyOut()) {
                if (const bool successfullPoke = insertActions(snapshot, actions,
                            LocationType_Func,
                            Coord3(rowIndex - 1, columnIndex, BLOCK_OUT),
                            XOR(number, block2.Out()));
                        !successfullPoke) return false;
            } else if (block1.emptyOut() && !block2.emptyOut()) {
                if (const bool successfullPoke = insertActions(snapshot, actions,
                            LocationType_Func,
                            Coord3(rowIndex - 1, columnIndex + 1, BLOCK_OUT),
                            XOR(number, block1.Out()));
                        !successfullPoke) return false;
            }
        }

        return true;
    }

    bool pokeFuncOut(Snapshot& snapshot, std::vector<Action> actions, const Location& location, const Number& number) {
        return true;
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
    srand (314);

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
