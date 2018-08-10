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


static const bool DO_LOG = false;
#define LOG(x) if(DO_LOG) x


static const unsigned int NUMBERS_IN_VECTOR = 3;
static const unsigned int NUMBER_BITS = 3;
static const unsigned int NUMBER_MAX = (1 << NUMBER_BITS) - 1;
static const unsigned int EMPTY_NUMBER = NUMBER_MAX + 1;
static const unsigned int FUNCS_AMOUNT = 2 * NUMBERS_IN_VECTOR - 1;


typedef unsigned int Number;
typedef std::array<Number, NUMBERS_IN_VECTOR> Vector;

// std::ostream& operator<<(std::ostream& stream, const Number& number) {
//     if (number == EMPTY_NUMBER) {
//         stream << "_";
//     } else {
//         stream << (0 + static_cast<unsigned int>(number));
//     }
//
//     return stream;
// }


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

std::ostream& operator<<(std::ostream& stream, const Coord1& coord) {
    stream << "C1(" << coord.x << ")";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const Coord2& coord) {
    stream << "C2(" << coord.x << "," << coord.y << ")";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const Coord3& coord) {
    stream << "C3(" << coord.x << "," << coord.y << "," << coord.z << ")";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const std::variant<Coord1, Coord2, Coord3>& coordVar) {
    if (std::holds_alternative<Coord1>(coordVar)) {
        stream << std::get<Coord1>(coordVar);
    } else if (std::holds_alternative<Coord2>(coordVar)) {
        stream << std::get<Coord3>(coordVar);
    } else if (std::holds_alternative<Coord3>(coordVar)) {
        stream << std::get<Coord3>(coordVar);
    }

    return stream;
}


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
    const auto range = high - low + 1;
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

        Func() {
            for (unsigned int i = 0; i < SIZE; i++) map[i] = EMPTY_NUMBER;
        }

        Number apply(const Number& arg) const {
            return map[arg];
        }

        Number at(const Number& arg) const {
            return map[arg];
        }

        bool set(const Number& arg) const {
            return (map[arg] != EMPTY_NUMBER);
        }

        bool emptyAt(const Number& number) const {
            return (map[number] == EMPTY_NUMBER);
        }

        void update(const Number& arg, const Number& value) {
            map[arg] = value;
        }

        // Number& operator[](const Number& index) {
        //     return map[index];
        // }

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

        Cell(const Number& value) : value(value) {}

        bool empty() const {
            return emptyNumber(value);
        }
};


static const unsigned int BLOCK_IN = 0;
static const unsigned int BLOCK_OUT = 1;
struct Block {
    private:
        Func& func;
        Number in;
        Number out;

        Number& operator[](unsigned int index) {
            return at(index);
        }

        Number& at(unsigned int index) {
            assert(index == 0 || index == 1);
            return (index == BLOCK_IN) ? in : out;
        }

        const Number& at(unsigned int index) const {
            assert(index == 0 || index == 1);
            return (index == BLOCK_IN) ? in : out;
        }


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
            return at(index);
        }

        const Number& In() const {
            return in;
        }

        const Number& Out() const {
            return out;
        }

        // Allow modifications only throught here
        void update(unsigned int index, const Number& number) {
            /* if (full()) { // the only case when func is already set */
            /*     func.update(in, EMPTY_NUMBER); // unset */
            /* } */
            at(index) = number;
            /* if (full()) { */
            /*     func.update(in, number); */
            /* } */
        }
};

class Snapshot {
    public:
        // Real physical Funcs
        std::array<Func, FUNCS_AMOUNT>& funcs;

        // Abstract snapshot Funcs
        std::vector<std::vector<Block>> blocks;
        std::array<Cell, NUMBERS_IN_VECTOR> buffers;
        std::array<Cell, NUMBERS_IN_VECTOR> xs;
        const std::array<Cell, NUMBERS_IN_VECTOR> ys;

    private:
        void init() {
            const unsigned int AMOUNT_OF_ROWS = NUMBERS_IN_VECTOR - 1;
            blocks.reserve(AMOUNT_OF_ROWS);
            for (unsigned int rowIndex = 0; rowIndex < AMOUNT_OF_ROWS; rowIndex++) {
                const unsigned int AMOUNT_OF_BLOCKS_IN_ROW = (2 * NUMBERS_IN_VECTOR - 1) - rowIndex;
                std::vector<Block> row;
                row.reserve(AMOUNT_OF_BLOCKS_IN_ROW);

                for (unsigned int blockIndex = 0; blockIndex < AMOUNT_OF_BLOCKS_IN_ROW; blockIndex++) {
                    row.emplace_back(funcs[blockIndex]);
                }

                blocks.push_back(row);
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
            for (const auto& row : blocks) {
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


Vector extractXs(const Snapshot& snapshot) {
    Vector vector;
    for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) {
        const Cell& cell = snapshot.xs[i];
        vector[i] = cell.value;
    }

    return vector;
}


std::ostream& operator<<(std::ostream& stream, const Snapshot& snapshot) {
    stream << "+++++++Snap+++++" << std::endl;
    /* stream << "Xs: "; */
    for (const auto& x : snapshot.xs) {
        stream << x.value << " ";
    }
    stream << std::endl;
    stream << std::endl;
    for (const auto& row : snapshot.blocks) {
        for (const auto& block : row) {
            stream << block.In() << " ";
        }
        stream << std::endl;
        for (const auto& block : row) {
            stream << block.Out() << " ";
        }
        stream << std::endl;
        stream << std::endl;
    }
    for (const auto& b : snapshot.buffers) {
        stream << b.value << " ";
    }
    stream << std::endl;
    for (const auto& y : snapshot.ys) {
        stream << y.value << " ";
    }
    stream << std::endl;
    stream << "++++++++++++++++" << std::endl;

    return stream;
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
enum LocationType {
    LocationType_X, LocationType_Buffer, LocationType_Block, LocationType_Func
};

struct Location {
    public:
        LocationType type;
        std::variant<Coord1, Coord2, Coord3> coord;

        Location(const LocationType& type, const std::variant<Coord1, Coord2, Coord3>& coord)
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
                case LocationType_Block: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const auto& coordinates = std::get<Coord3>(coord);
                        snapshot.blocks[coordinates.x][coordinates.y].update(coordinates.z, number);
                    }
                    break;
                case LocationType_Func: {
                        assert(std::holds_alternative<Coord2>(coord));
                        const auto& coordinates = std::get<Coord2>(coord);
                        snapshot.funcs[coordinates.x].update(coordinates.y, number);
                    }
                    break;
            }
        }

        const Number& value(Snapshot& snapshot) const {
            switch (type) {
                case LocationType_X:
                    assert(std::holds_alternative<Coord1>(coord));
                    return snapshot.xs[std::get<Coord1>(coord).x].value;
                    break;
                case LocationType_Buffer:
                    assert(std::holds_alternative<Coord1>(coord));
                    return snapshot.buffers[std::get<Coord1>(coord).x].value;
                    break;
                case LocationType_Block: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const auto& coordinates = std::get<Coord3>(coord);
                        const Block& block = snapshot.blocks[coordinates.x][coordinates.y];
                        return block[coordinates.z];
                    }
                    break;
                case LocationType_Func: {
                        assert(std::holds_alternative<Coord2>(coord));
                        const auto& coordinates = std::get<Coord2>(coord);
                        return snapshot.funcs[coordinates.x][coordinates.y];
                    }
                    break;
            }

            return EMPTY_NUMBER;
        }

        bool empty(Snapshot& snapshot) const {
            return (emptyNumber(value(snapshot)));
        }
};

std::ostream& operator<<(std::ostream& stream, const LocationType& locationType) {
    switch (locationType) {
        case LocationType_X:
            stream << "LocationType_X";
            break;
        case LocationType_Buffer:
            stream << "LocationType_Buffer";
            break;
        case LocationType_Block:
            stream << "LocationType_Block";
            break;
        case LocationType_Func:
            stream << "LocationType_Func";
            break;
    }

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Location& location) {
    stream << "Loc(" << location.type << ", " << location.coord << ")";

    return stream;
}

// ----------------------------------------------------------------------------
// Unconditional
struct Action {
    private:
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

        friend std::ostream& operator<<(std::ostream& stream, const Action& action);
};

std::ostream& operator<<(std::ostream& stream, const Action& action) {
    stream << "Act(" << action.location << ", " << action.value << ")";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const std::vector<Action>& v) {
    for (const Action& a : v) {
        stream << a << "; ";
    }

    return stream;
}


void undoActions(Snapshot& snapshot, std::vector<Action>& actions) {
    for (Action& action : actions) {
        action.undo(snapshot);
    }
}


// Conditional, random
struct Step {
    private:
        static const unsigned int ALLOWED_RETRIES = 4;

        Location baseLocation;
        std::vector<Action> actions;
        std::set<Number> triedNumbers;

    public:
        Step(const Location& location) : baseLocation(location) {}

        void useNumber(const Number& number) {
            triedNumbers.insert(number);
            // cleanup(snapshot);
            // actions.emplace_back(baseLocation, number);
        }

        void applyActions(const std::vector<Action>& newActions) {
            actions.insert(actions.end(), newActions.begin(), newActions.end());
        }

        void fail(Snapshot& snapshot) {
            cleanup(snapshot);
        }

        bool belowLimit() {
            return (triedNumbers.size() < ALLOWED_RETRIES);
        }

        void cleanup(Snapshot& snapshot) {
            for (Action& action : actions) {
                action.undo(snapshot);
            }
            actions.clear();
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

        friend std::ostream& operator<<(std::ostream& stream, const Step& step);
};

std::ostream& operator<<(std::ostream& stream, const Step& step) {
    stream << "Step(" << step.baseLocation << "," << step.actions << ")";

    return stream;
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
class Stack {
    private:
        std::stack<Step> steps;

    public:
        Stack() {}

        void push(Step step) {
            steps.push(step);
        }

        void emplace(const Location& location) {
            steps.emplace(location);
        }

        std::optional<Step> pop() {
            if (steps.size() == 0) {
                return std::nullopt;
            }

            const Step step = steps.top();
            steps.pop();

            return {step};
        }

        std::optional<std::reference_wrapper<Step>> top() {
            if (steps.size() > 0) {
                return {steps.top()};
            }

            return std::nullopt;
        }

        /* const Step& top() const { */
        /*     return steps.top(); */
        /* } */

        bool belowLimit() {
            return (steps.top().belowLimit());
        }

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

        // bool revert() {
        //     if (nextIndex > 0) {
        //         nextIndex--;
        //         return true;
        //     }
        //
        //     return false;
        // }
};
// ----------------------------------------------------------------------------
std::optional<Coord2> findEmptyBlock(const Snapshot& snapshot) {
    unsigned int rowIndex = 0;
    unsigned int columnIndex = 0;
    for (const auto& row : snapshot.blocks) {
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
void dumpFuncs(const Functions& funcs) {
    std::cout << "Funcs:" << std::endl;
    for (const Func& f : funcs) {
        for (unsigned int i = 0; i <= NUMBER_MAX; i++) {
            std::cout << f.at(i) << " ";
        }
        std::cout << std::endl;
    }
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

private:
    bool numberHasSecondaryIn(unsigned int index) const {
        return (index < NUMBERS_IN_VECTOR - 1);
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
        LOG(std::cout << ">> Genning new." << std::endl;)
        Snapshot snapshot(funcs, y);
        Stack stepsStack;
        Strategy strategy([](
                    unsigned int index,
                    const Snapshot& snapshot,
                    const Functions& funcs)
                -> std::optional<Location> {
            if (index < NUMBERS_IN_VECTOR) { // first three insertions, index is 0-based
                return {{LocationType_X, index}};
            } else {
                if (const std::optional<Coord2> emptyBlock = findEmptyBlock(snapshot)) {
                    const Coord2& blockCoord = *emptyBlock;
                    return {{LocationType_Block,
                             Coord3(
                                    blockCoord.x,
                                    blockCoord.y,
                                    snapshot.blocks[blockCoord.x][blockCoord.y].emptyIn() ? BLOCK_IN : BLOCK_OUT)}};
                } else {
                    return std::nullopt;
                }
            }
        });

        while (const auto locationOpt = strategy.getNext(snapshot, funcs)) {
            const Location location = *locationOpt;
            std::cout << ">> Got next strategy: " << location << "." << std::endl;

            std::optional<Step> stepOpt{location};
            do {
                /* LOG(std::cout << "Iteration: " << iter++ << std::endl;) */
                LOG(std::cout << ">> Inside " << (*stepOpt) << std::endl;)
                Step& step = *stepOpt;
                const std::optional<Number> numberOpt = step.getUntriedNumber();
                if (numberOpt) {
                    std::cout << "With Number " << (*numberOpt) << std::endl;
                } else {
                    std::cout << "With Number " << "*" << std::endl;
                }
                LOG(dumpFuncs(snapshot.funcs);)
                const std::optional<std::vector<Action>> actionsOpt =
                    (numberOpt) ? (stepOpt->useNumber(*numberOpt), poke(snapshot, location, *numberOpt))
                                : (std::nullopt);
                if (numberOpt && actionsOpt) {
                    std::cout << ">> Success " << std::endl;
                    std::cout << snapshot << std::endl;
                    LOG(dumpFuncs(snapshot.funcs);)
                    step.applyActions(*actionsOpt);
                    stepsStack.push(step);

                    stepOpt = std::nullopt;
                } else {
                    std::cout << ">> Fail " << std::endl;
                    /* step.fail(snapshot); // inc tries. Never need to empty Actions here */
                    while (stepOpt && !stepsStack.belowLimit()) {
                        stepOpt = stepsStack.pop(); // std::nullopt if stack became empty
                        if (stepOpt) {
                            stepOpt->fail(snapshot); // child failed => parent failed. Empty parent's Actions
                        }
                    }

                    // By this line we've either reverted to a valid Step on Stack
                    // or made stepOpt == std::nullopt if the stack depleted
                    if (!stepOpt) {
                        return false; // failed to construct Vector<X>
                    }
                }
            } while (stepOpt);
            // (!stepOpt) => a Step succeeded and we need to generate a subsequent
            // Location (or finish, if no more is to be done).
        }

        // Everything is swell, extract the Xs from snapshot
        xs.push_back(extractXs(snapshot));

        return false; // TODO: change back to true
    }


    bool reactOnFuncUpdateAndInsert(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& coord) {
        LOG(std::cout << ">>> Inside reactOnFunc for " << coord << std::endl;)
        Block& block = snapshot.blocks[coord.x][coord.y];

        if (block.empty()) {
            LOG(std::cout << "both empty" << std::endl;)
            return true; // nothing new for this place
        }
        const unsigned int functionIndex = coord.y;
        if (block.full()) {
            LOG(std::cout << "both full, ";)
            LOG(std::cout << (snapshot.funcs[functionIndex].apply(block.In()) == block.Out() ? "matched" : "no match");)
            LOG(std::cout << std::endl;)
            /* std::cout << block.In() << " : " << block.Out() << std::endl; */
            /* dumpFuncs(snapshot.funcs); */
            return snapshot.funcs[functionIndex].apply(block.In()) == block.Out();
        } else if (block.emptyOut()) {
            if (const Number funcRes = snapshot.funcs[functionIndex].apply(block.In());
                    !emptyNumber(funcRes)) {
                // can deduce Out from function => try placing Out
                LOG(std::cout << "Func set for Out. placing Out into Block";)
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Block,
                            Coord3(coord.x, coord.y, BLOCK_OUT),
                            funcRes);
                        !successfulPoke) return false;
            }
            LOG(std::cout << "Func not set for Out";)
        }

        LOG(std::cout << ">>> Succ reactOnFunc " << std::endl;)
        return true; // block.In is empty => nothing is determined yet
    }


    bool pokeAndInsert(
            Snapshot& snapshot,
            std::vector<Action>& actions,
            const LocationType& locationType,
            std::variant<Coord1, Coord2, Coord3> coord,
            const Number& value) {
        LOG(std::cout << "Inside pokeAndInsert" << std::endl;)
        if (const auto& actionsOpt = poke(snapshot, Location(locationType, coord), value)) {
            actions.insert(actions.end(), actionsOpt->begin(), actionsOpt->end());
            LOG(std::cout << "Succ done pokeAndInsert" << std::endl;)
            return true;
        }

        LOG(std::cout << "Failed pokeAndInsert" << std::endl;)
        LOG(dumpFuncs(snapshot.funcs);)
        return false;
    }


    // returns std::nullopt if a conflict arose
    std::optional<std::vector<Action>> poke(Snapshot& snapshot, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside Poke " << std::endl;)
        std::vector<Action> actions;
        actions.emplace_back(location, number); // self
        actions[0].apply(snapshot);
        LOG(std::cout << actions << std::endl;)
        LOG(dumpFuncs(snapshot.funcs);)

        switch (location.type) {
            case LocationType_X: // is to poked set only from the outside
                if (const auto successfulPoke = pokeX(snapshot, actions, location, number);
                        !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType_Buffer:
                if (const auto successfulPoke = pokeBuffer(snapshot, actions, location, number);
                        !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType_Block:
                assert(std::holds_alternative<Coord3>(location.coord));
                switch (const Coord3 coord = std::get<Coord3>(location.coord);
                        coord.z) {
                    case BLOCK_IN:
                        if (const auto successfulPoke = pokeFuncIn(snapshot, actions, location, number);
                                !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                        break;
                    case BLOCK_OUT:
                        if (const auto successfulPoke = pokeFuncOut(snapshot, actions, location, number);
                                !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                        break;
                }
                break;
            case LocationType_Func:
                assert(false && "Why would you poke a Func???");
                break;
        }

        LOG(std::cout << ">>> Succ done Poke " << std::endl;)
        return {actions};
    }


    bool pokeX(Snapshot& snapshot, std::vector<Action>& actions, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside PokeX " << std::endl;)
        assert(std::holds_alternative<Coord1>(location.coord));
        const unsigned int index = std::get<Coord1>(location.coord).x;

        // Main In
        assert(snapshot.blocks[0][index].emptyIn() && ":> IN not empty when X got poked.");
        if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                    LocationType_Block,
                    Coord3(0, index, BLOCK_IN),
                    number);
                !successfulPoke) return false;

        // Secondary In
        if (numberHasSecondaryIn(index)) {
            assert(snapshot.blocks[0][index + NUMBERS_IN_VECTOR].emptyIn() && ":> Secondary IN not empty when X got poked.");
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType_Block,
                        Coord3(0, index + NUMBERS_IN_VECTOR, BLOCK_IN),
                        number);
                    !successfulPoke) return false;
        }

        // Buffer
        assert(snapshot.buffers[index].empty() && ":> Buffer not empty when X got poked.");
        if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                    LocationType_Buffer,
                    Coord1(index),
                    XOR(number, snapshot.ys[index].value));
                !successfulPoke) return false;

        LOG(std::cout << ">>> Succ done PokeX" << std::endl;)
        return true;
    }

    bool pokeBuffer(Snapshot& snapshot, std::vector<Action>& actions, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside PokeBuffer " << std::endl;)
        assert(std::holds_alternative<Coord1>(location.coord));
        const unsigned int index = std::get<Coord1>(location.coord).x;
        const Number numberX = XOR(number, snapshot.ys[index].value);

        if (!snapshot.xs[index].empty() && (numberX != snapshot.xs[index].value)) {
            // The value we're trying to set to Buffer does not comply with the current value of X
            return false;
        }
        if (snapshot.xs[index].empty()) {
            // just set X, but never poke it
            actions.emplace_back(Location(LocationType_X, Coord1(index)), numberX);
            (actions.end() - 1)->apply(snapshot);
        }

        // Main In
        if (snapshot.blocks[0][index].emptyIn()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType_Block,
                        Coord3(0, index, BLOCK_IN),
                        numberX);
                    !successfulPoke) return false;
        }

        // Secondary In
        if (numberHasSecondaryIn(index)) {
            if (snapshot.blocks[0][index + NUMBERS_IN_VECTOR].emptyIn()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Block,
                            Coord3(0, index + NUMBERS_IN_VECTOR, BLOCK_IN),
                            numberX);
                        !successfulPoke) return false;
            }
        }

        // Outs
        const unsigned int lastRowIndex = snapshot.blocks.size() - 1;
        const Block& block1Const = snapshot.blocks[lastRowIndex][index];
        const Block& block2Const = snapshot.blocks[lastRowIndex][index + 1];
        Block& block1 = snapshot.blocks[lastRowIndex][index];
        Block& block2 = snapshot.blocks[lastRowIndex][index + 1];
        if (!block1Const.emptyOut() && !block2Const.emptyOut()) {
            if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                return false;
            }
        }
        if (block1.emptyOut() && !block2.emptyOut()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType_Block,
                        Coord3(lastRowIndex, index, BLOCK_OUT),
                        XOR(number, block2.Out()));
                    !successfulPoke) return false;
        } else if (block1.emptyOut() && !block2.emptyOut()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType_Block,
                        Coord3(lastRowIndex, index + 1, BLOCK_OUT),
                        XOR(number, block1.Out()));
                    !successfulPoke) return false;
        }

        LOG(std::cout << ">>> Succ done PokeBuffer" << std::endl;)
        return true;
    }

    bool pokeFuncIn(Snapshot& snapshot, std::vector<Action>& actions, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside PokeIn " << std::endl;)
        assert(std::holds_alternative<Coord3>(location.coord));
        const Coord3 coord = std::get<Coord3>(location.coord);
        const unsigned int rowIndex = coord.x;
        const unsigned int columnIndex = coord.y;

        // Self Out
        Block& block = snapshot.blocks[rowIndex][columnIndex];
        const auto& funcRes = snapshot.funcs[columnIndex].apply(number);
        if (emptyNumber(block.Out())) {
            if (!emptyNumber(funcRes)) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Block,
                            Coord3(rowIndex, columnIndex, BLOCK_OUT),
                            funcRes);
                        !successfulPoke) return false;
            }
        } else {
            if (emptyNumber(funcRes)) {
                if (const auto successfulSet = setFuncAndPokeColumnAndInsert(snapshot, actions,
                            Coord2(rowIndex, columnIndex));
                        !successfulSet) return false;
            } else {
                if (block.Out() != funcRes) {
                    return false;
                }
            }
        }

        if (rowIndex == 0) {
            const bool isMainIn = columnIndex < NUMBERS_IN_VECTOR;
            const unsigned int xIndex =
                (isMainIn) ? columnIndex
                           : (columnIndex - NUMBERS_IN_VECTOR);
            if (snapshot.xs[xIndex].empty()) {
                // Set X (just set, never poke)
                actions.emplace_back(Location(LocationType_X, Coord1(xIndex)), number);
                (actions.end() - 1)->apply(snapshot);

                // Set INs
                const unsigned int otherColumnIndex =
                    (isMainIn) ? (columnIndex + NUMBERS_IN_VECTOR)
                               : (columnIndex - NUMBERS_IN_VECTOR);
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_X,
                            Coord1(otherColumnIndex),
                            number);
                        !successfulPoke) return false;

                // Set Buffer
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Buffer,
                            Coord1(columnIndex),
                            number);
                        !successfulPoke) return false;
            } else {
                if (snapshot.xs[xIndex].value != number) {
                    return false;
                }
            }
        } else { // rowIndex >= 1
            // foreign Outs above
            const Block& block1Const = snapshot.blocks[rowIndex - 1][columnIndex];
            const Block& block2Const = snapshot.blocks[rowIndex - 1][columnIndex + 1];
            Block& block1 = snapshot.blocks[rowIndex - 1][columnIndex];
            Block& block2 = snapshot.blocks[rowIndex - 1][columnIndex + 1];
            if (!block1Const.emptyOut() && !block2Const.emptyOut()) {
                if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                    return false;
                }
            }
            if (block1.emptyOut() && !block2.emptyOut()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Block,
                            Coord3(rowIndex - 1, columnIndex, BLOCK_OUT),
                            XOR(number, block2.Out()));
                        !successfulPoke) return false;
            } else if (block1.emptyOut() && !block2.emptyOut()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType_Block,
                            Coord3(rowIndex - 1, columnIndex + 1, BLOCK_OUT),
                            XOR(number, block1.Out()));
                        !successfulPoke) return false;
            }
        }

        LOG(std::cout << ">>> Succ done PokeIn" << std::endl;)
        return true;
    }

    bool pokeFuncOut(Snapshot& snapshot, std::vector<Action>& actions, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside PokeOut " << std::endl;)
        assert(std::holds_alternative<Coord3>(location.coord));
        const Coord3 coord = std::get<Coord3>(location.coord);
        const unsigned int rowIndex = coord.x;
        const unsigned int columnIndex = coord.y;

        // Self In
        if (const auto& block = snapshot.blocks[rowIndex][columnIndex];
                !block.emptyIn()) {
            if (const auto& funcRes = snapshot.funcs[columnIndex].at(block.In());
                    !emptyNumber(funcRes)) {
                if (funcRes != number) {
                    return false;
                }
            } else {
                std::cout << "still here 0" << std::endl;
                if (const auto successfulSet = setFuncAndPokeColumnAndInsert(snapshot, actions,
                            Coord2(rowIndex, columnIndex));
                        !successfulSet) return false;
            }
        }
        // too unlikely that func is full and exactly 1 input has us as output
        // (then we could derive In), so don't even bother here.


        const auto processPairAndInsert = [&](std::vector<Action>& actions, Location location1, Location location2, const Number& number) {
            if (location1.empty(snapshot)) {
                if (!location2.empty(snapshot)) {
                    if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                                location1.type,
                                location1.coord,
                                XOR(location2.value(snapshot), number));
                            !successfulPoke) return false;
                }
            } else {
                if (location2.empty(snapshot)) {
                    if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                                location2.type,
                                location2.coord,
                                XOR(location1.value(snapshot), number));
                            !successfulPoke) return false;
                } else {
                    if (XOR(location1.value(snapshot), location2.value(snapshot)) != number) {
                        return false;
                    }
                }
            }

            return true;
        };

        // Connections through XOR
        const auto rowsAmount = snapshot.blocks.size() - 1;
        const auto& row = snapshot.blocks[rowsAmount];
        const bool bottom = rowIndex == rowsAmount;
        if (columnIndex < row.size() - 1) { // left tail
            if (const bool successfulProcess = processPairAndInsert(
                        actions,
                        Location(LocationType_Block, Coord3(rowIndex, columnIndex + 1, BLOCK_OUT)),
                        bottom ? Location(LocationType_Buffer, Coord1(columnIndex))
                               : Location(LocationType_Block, Coord3(rowIndex + 1, columnIndex, BLOCK_IN)),
                        number);
                    !successfulProcess) return false;
        }
        if (columnIndex >= 1) { // right tail
            if (const bool successfulProcess = processPairAndInsert(
                        actions,
                        Location(LocationType_Block, Coord3(rowIndex, columnIndex - 1, BLOCK_OUT)),
                        bottom ? Location(LocationType_Buffer, Coord1(columnIndex - 1))
                               : Location(LocationType_Block, Coord3(rowIndex + 1, columnIndex - 1, BLOCK_IN)),
                        number);
                    !successfulProcess) return false;
        }

        LOG(std::cout << "Succ done PokeOut" << std::endl;)
        return true;
    }


    bool setFuncAndPokeColumnAndInsert(Snapshot& snapshot, std::vector<Action>& actions, Coord2 blockCoord) {
        LOG(std::cout << ">>> Inside setFuncAndPokeAndInsert " << std::endl;)
        const auto& block = snapshot.blocks[blockCoord.x][blockCoord.y];
        assert(snapshot.funcs[blockCoord.y].emptyAt(block.In()));

        actions.emplace_back(Location(LocationType_Func, Coord2(blockCoord.y, block.In())), block.Out()); // set func
        (actions.end() - 1)->apply(snapshot);
        /* snapshot.funcs[blockCoord.y].update(block.In(), block.Out()); // set func */

        // func for this column has beed altered => poke all Blocks in this column
        for (unsigned int row = 0; row < snapshot.blocks.size(); row++) {
            // case for row == rowIndex is handled correctly in reactOnFuncUpdateAndInsert, so no worries
            if (const auto updateSuccessful = reactOnFuncUpdateAndInsert(snapshot, actions, Coord2(row, blockCoord.y));
                    !updateSuccessful) return false;
        }

        LOG(std::cout << ">>> Succ done setFuncAndPokeAndInsert " << std::endl;)
        return true;
    }




    void outputResult(const std::vector<Vector>& xs, const Vector& y) {
        std::cout << std::endl;
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
    mb.work();


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
