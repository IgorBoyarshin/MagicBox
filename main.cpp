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
#include <cstring>


// ----------------------------------------------------------------------------
static const bool DO_LOG = false;
#define LOG(x) if(DO_LOG) {x}
// ----------------------------------------------------------------------------
static const unsigned int NUMBERS_IN_VECTOR = 4;
static const unsigned int NUMBER_BITS = 5;

static const unsigned int NUMBERS_AMOUNT = (1 << NUMBER_BITS);
static const unsigned int NUMBER_MAX = NUMBERS_AMOUNT - 1;
static const unsigned int EMPTY_NUMBER = NUMBER_MAX + 1;
static const unsigned int FUNCS_AMOUNT = 2 * NUMBERS_IN_VECTOR - 1;
// ----------------------------------------------------------------------------
struct Number {
    unsigned int value;

    Number() : value(EMPTY_NUMBER) {}
    Number(unsigned int value) : value(value) {}

    inline bool empty() const { return value == EMPTY_NUMBER; };
    inline unsigned int operator()() const { return value; }
    inline bool operator==(const Number& other) const { return value == other(); };
    inline bool operator!=(const Number& other) const { return value != other(); };
};

std::ostream& operator<<(std::ostream& stream, const Number& number) {
    if (number.empty()) {
        stream << "*";
    } else {
        stream << number();
    }

    return stream;
}

typedef std::array<Number, NUMBERS_IN_VECTOR> Vector;

bool full(const Vector& vector) {
    for (const Number& number : vector) {
        if (number.empty()) {
            return false;
        }
    }

    return true;
}


inline Number XOR(const Number& a, const Number& b) {
    return a() ^ b();
}

bool operator==(const Vector& v1, const Vector& v2) {
    for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }

    return true;
};

bool operator!=(const Vector& v1, const Vector& v2) {
    return !(operator==(v1, v2));
}

bool unique(const Vector& vector, const std::vector<Vector>& xs) {
    for (const Vector& x : xs) {
        if (vector == x) {
            return false;
        }
    }

    return true;
}

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
// [low, high]
// template<typename T = unsigned int>
// T generateRandomUniformInt(T low, T high) {
//     static std::random_device rd;
//     static std::mt19937 e2(rd());
//     static std::uniform_int_distribution<T> dist(low, high);
//
//     return dist(e2);
// }

unsigned int generateRandomUniformInt(unsigned int low, unsigned int high) {
    const auto range = high - low + 1;
    return rand() % range + low;
}

Number generateRandomNumber() {
    return generateRandomUniformInt((unsigned int)0, NUMBER_MAX);
}

Vector generateRandomVector() {
    Vector vector;
    for (Number& number : vector) {
        number = generateRandomNumber();
    }

    return vector;
}
// ----------------------------------------------------------------------------
struct Func {
    public:
        static const unsigned int SIZE = NUMBERS_AMOUNT;

    private:
        // num_in -> num_out
        Number map[SIZE];

    public:
        Func() {
            for (unsigned int i = 0; i < SIZE; i++) map[i] = Number(EMPTY_NUMBER);
        }

        bool setAt(const Number& arg) const {
            return (!map[arg()].empty());
        }

        bool emptyAt(const Number& arg) const {
            return (map[arg()].empty());
        }

        void update(const Number& arg, const Number& value) {
            map[arg()] = value;
        }

        const Number& at(const Number& arg) const {
            return map[arg()];
        }

        const Number& operator[](const Number& arg) const {
            return map[arg()];
        }
};

typedef std::array<Func, FUNCS_AMOUNT> Functions;

std::ostream& operator<<(std::ostream& stream, const Functions& funcs) {
    stream << "Funcs:" << std::endl;
    for (unsigned int i = 0; i < Func::SIZE; i++) {
        for (const Func& func : funcs) {
            stream << i << " -> " << func.at(i) << "\t\t";
        }
        stream << std::endl;
    }

    return stream;
}
// ----------------------------------------------------------------------------
enum class BlockCell {
    IN, OUT
};

std::ostream& operator<<(std::ostream& stream, const BlockCell& blockCell) {
    switch (blockCell) {
        case BlockCell::IN:
            stream << "IN";
            break;
        case BlockCell::OUT:
            stream << "OUT";
            break;
    }

    return stream;
}

struct Block {
    private:
        Number in;
        Number out;

        Number& at(BlockCell blockCell) {
            return (blockCell == BlockCell::IN) ? in : out;
        }

        const Number& at(BlockCell blockCell) const {
            return (blockCell == BlockCell::IN) ? in : out;
        }

    public:
        Block() : in(EMPTY_NUMBER), out(EMPTY_NUMBER) {}

        bool existsEmpty() const {
            return !full();
        }

        bool empty() const {
            return (in.empty() && out.empty());
        }

        bool emptyIn() const {
            return in.empty();
        }

        bool emptyOut() const {
            return out.empty();
        }

        bool full() const {
            return (!in.empty() && !out.empty());
        }

        const Number& operator[](BlockCell blockCell) const {
            return at(blockCell);
        }

        const Number& In() const {
            return in;
        }

        const Number& Out() const {
            return out;
        }

        void update(BlockCell blockCell, const Number& number) {
            at(blockCell) = number;
        }
};
// ----------------------------------------------------------------------------
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
        BlockCell z;

        Coord3(unsigned int x, unsigned int y, BlockCell z) : x(x), y(y), z(z) {}
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
        stream << std::get<Coord2>(coordVar);
    } else if (std::holds_alternative<Coord3>(coordVar)) {
        stream << std::get<Coord3>(coordVar);
    }
    return stream;
}
// ----------------------------------------------------------------------------
class Snapshot {
    public:
        Functions& funcs;

        std::vector<std::vector<Block>> blocks;
        std::array<Number, NUMBERS_IN_VECTOR> buffers;
        std::array<Number, NUMBERS_IN_VECTOR> xs;
        const std::array<Number, NUMBERS_IN_VECTOR> ys;

    private:
        void init() {
            const unsigned int AMOUNT_OF_ROWS = NUMBERS_IN_VECTOR - 1;
            blocks.reserve(AMOUNT_OF_ROWS);
            for (unsigned int rowIndex = 0; rowIndex < AMOUNT_OF_ROWS; rowIndex++) {
                const unsigned int AMOUNT_OF_BLOCKS_IN_ROW = (2 * NUMBERS_IN_VECTOR - 1) - rowIndex;
                std::vector<Block> row;
                row.reserve(AMOUNT_OF_BLOCKS_IN_ROW);

                for (unsigned int blockIndex = 0; blockIndex < AMOUNT_OF_BLOCKS_IN_ROW; blockIndex++) {
                    row.emplace_back();
                }

                blocks.push_back(row);
            }
        }

    public:
        Snapshot(Functions& funcs, const Vector& ys) : funcs(funcs), ys(ys) {
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
            for (const Number& buffer : buffers) {
                if (buffer.empty()) {
                    return true;
                }
            }
            for (const Number& x : xs) {
                if (x.empty()) {
                    return true;
                }
            }

            return false;
        }
};


bool full(const Snapshot& snapshot) {
    if (!full(snapshot.xs     )) return false;
    if (!full(snapshot.ys     )) return false;
    if (!full(snapshot.buffers)) return false;

    for (const auto& row : snapshot.blocks) {
        for (const auto& block : row) {
            if (block.existsEmpty()) {
                return false;
            }
        }
    }

    return true;
}

std::optional<Coord2> findEmptyBlock(const Snapshot& snapshot) {
    unsigned int rowIndex = 0;
    for (const auto& row : snapshot.blocks) {
        unsigned int columnIndex = 0;
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


std::ostream& operator<<(std::ostream& stream, const Snapshot& snapshot) {
    stream << "+++++++Snap+++++" << std::endl;
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
enum class LocationType {
    X, Buffer, Block, Func
};

struct Location {
    public:
        LocationType type;
        std::variant<Coord1, Coord2, Coord3> coord;

        Location(const LocationType& type, const std::variant<Coord1, Coord2, Coord3>& coord)
            : type(type), coord(coord) {}

        void place(const Number& number, Snapshot& snapshot) const {
            switch (type) {
                case LocationType::X:
                    assert(std::holds_alternative<Coord1>(coord));
                    snapshot.xs[std::get<Coord1>(coord).x] = number;
                    break;
                case LocationType::Buffer:
                    assert(std::holds_alternative<Coord1>(coord));
                    snapshot.buffers[std::get<Coord1>(coord).x] = number;
                    break;
                case LocationType::Block: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const Coord3& coords = std::get<Coord3>(coord);
                        snapshot.blocks[coords.x][coords.y].update(coords.z, number);
                    }
                    break;
                case LocationType::Func: {
                        assert(std::holds_alternative<Coord2>(coord));
                        const Coord2& coords = std::get<Coord2>(coord);
                        snapshot.funcs[coords.x].update(coords.y, number);
                    }
                    break;
            }
        }

        const Number& value(Snapshot& snapshot) const {
            switch (type) {
                case LocationType::X:
                    assert(std::holds_alternative<Coord1>(coord));
                    return snapshot.xs[std::get<Coord1>(coord).x];
                    break;
                case LocationType::Buffer:
                    assert(std::holds_alternative<Coord1>(coord));
                    return snapshot.buffers[std::get<Coord1>(coord).x];
                    break;
                case LocationType::Block: {
                        assert(std::holds_alternative<Coord3>(coord));
                        const Coord3& coords = std::get<Coord3>(coord);
                        return snapshot.blocks[coords.x][coords.y][coords.z];
                    }
                    break;
                case LocationType::Func: {
                        assert(std::holds_alternative<Coord2>(coord));
                        const Coord2& coords = std::get<Coord2>(coord);
                        return snapshot.funcs[coords.x][coords.y];
                    }
                    break;
            }

            assert(false && "Switch didn't handle all cases.");
            return snapshot.xs[0]; // invalid number!!! (just to return smth)
        }

        bool empty(Snapshot& snapshot) const {
            return (value(snapshot).empty());
        }
};

std::ostream& operator<<(std::ostream& stream, const LocationType& locationType) {
    switch (locationType) {
        case LocationType::X:
            stream << "LocationType::X";
            break;
        case LocationType::Buffer:
            stream << "LocationType::Buffer";
            break;
        case LocationType::Block:
            stream << "LocationType::Block";
            break;
        case LocationType::Func:
            stream << "LocationType::Func";
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

void undoActions(Snapshot& snapshot, std::vector<Action>& actions) {
    for (Action& action : actions) {
        action.undo(snapshot);
    }
}

std::ostream& operator<<(std::ostream& stream, const Action& action) {
    stream << "Act(" << action.location << ", " << action.value << ")";

    return stream;
}
std::ostream& operator<<(std::ostream& stream, const std::vector<Action>& actions) {
    for (const Action& action : actions) {
        stream << action << "; ";
    }

    return stream;
}
// ----------------------------------------------------------------------------
// Conditional, random
struct Step {
    private:
        static const unsigned int ALLOWED_RETRIES = NUMBERS_AMOUNT; // <= NUMBERS_AMOUNT

        Location baseLocation;
        std::vector<Action> actions;
        bool triedNumbers[NUMBERS_AMOUNT];
        unsigned int triedNumbersAmount;

    public:
        Step(const Location& location) : baseLocation(location), triedNumbersAmount(0) {
            memset(triedNumbers, false, NUMBERS_AMOUNT * sizeof(bool));
        }

        void useNumber(const Number& number) {
            if (!triedNumbers[number()]) triedNumbersAmount++;
            triedNumbers[number()] = true;
        }

        void applyActions(const std::vector<Action>& newActions) {
            actions.insert(actions.end(), newActions.begin(), newActions.end());
        }

        void fail(Snapshot& snapshot) {
            for (Action& action : actions) {
                action.undo(snapshot);
            }
            actions.clear();
        }

        Location getLocation() const {
            return baseLocation;
        }

        bool belowLimit() const {
            return (triedNumbersAmount < ALLOWED_RETRIES);
        }

        std::optional<Number> getUntriedNumber() const {
            if (triedNumbersAmount >= ALLOWED_RETRIES) return std::nullopt;

            std::vector<Number> untriedNumbers;
            for (unsigned int i = 0; i < NUMBERS_AMOUNT; i++) {
                if (!triedNumbers[i]) untriedNumbers.emplace_back(i);
            }
            return { untriedNumbers[generateRandomUniformInt(0, untriedNumbers.size() - 1)] };
        }

        friend std::ostream& operator<<(std::ostream& stream, const Step& step);
};

std::ostream& operator<<(std::ostream& stream, const Step& step) {
    stream << "Step(" << step.baseLocation << "," << step.actions << "," << step.triedNumbersAmount << ")";

    return stream;
}
// ----------------------------------------------------------------------------
class Stack {
    private:
        std::stack<Step> steps;

    public:
        void push(const Step& step) {
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
};
// ----------------------------------------------------------------------------
typedef std::function<std::optional<Location>(unsigned int index, const Snapshot& snapshot, const Functions& funcs)> StrategyFunction;

class Strategy {
    private:
        const StrategyFunction function;
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
            if (nextIndex > 0) {
                nextIndex--;
                return true;
            }

            return false;
        }
};
// ----------------------------------------------------------------------------
class MagicBox {
public:
    Functions funcs;

private:
    bool hasSecondaryIn(unsigned int index) const {
        return (index < NUMBERS_IN_VECTOR - 1);
    }

public:
    std::vector<Vector> work(const Vector& y) {
        std::vector<Vector> xs;
        while (generateNewX(xs, y));
        outputResult(xs, y);

        return xs;
    }

    bool generateNewX(std::vector<Vector>& xs, const Vector& y) {
        static auto counter = 0;
        LOG(std::cout << ">> Generating new X(" << counter++ << ")" << std::endl;)
        Snapshot snapshot(funcs, y);
        Stack stepsStack;
        Strategy strategy([](
                    unsigned int index,
                    const Snapshot& snapshot,
                    const Functions& funcs)
                -> std::optional<Location> {
            if (index < NUMBERS_IN_VECTOR) { // first three insertions
                return { {LocationType::X, {index}} };
            } else {
                if (const std::optional<Coord2> emptyBlock = findEmptyBlock(snapshot)) {
                    const Coord2& blockCoord = *emptyBlock;
                    return {{LocationType::Block,
                             Coord3(
                                    blockCoord.x,
                                    blockCoord.y,
                                    snapshot.blocks[blockCoord.x][blockCoord.y].emptyIn()
                                        ? BlockCell::IN
                                        : BlockCell::OUT)}};
                } else {
                    return std::nullopt;
                }
            }
        });

        while (const auto locationOpt = strategy.getNext(snapshot, funcs)) {
            const Location location = *locationOpt;
            LOG(std::cout << ">> Got next strategy: " << location << "." << std::endl;)

            std::optional<Step> stepOpt{location};
            do {
                Step& step = *stepOpt;
                const std::optional<Number> numberOpt = step.getUntriedNumber();

                {
                    LOG(std::cout << ">> Inside " << (*stepOpt) << std::endl;)
                    LOG(if (numberOpt) { std::cout << "With Number " << (*numberOpt) << std::endl; } else { std::cout << "With Number " << "*" << std::endl; })
                    LOG(std::cout << snapshot << std::endl;)
                }

                const std::optional<std::vector<Action>> actionsOpt =
                    (numberOpt) ? (stepOpt->useNumber(*numberOpt),
                                    poke(snapshot, step.getLocation(), *numberOpt))
                                : (std::nullopt);
                if (numberOpt && actionsOpt) {
                    {
                        LOG(std::cout << ">> Success " << std::endl;)
                        LOG(std::cout << snapshot << std::endl;)
                        LOG(std::cout << snapshot.funcs;)
                    }

                    step.applyActions(*actionsOpt);
                    stepsStack.push(step);
                    stepOpt = std::nullopt;
                } else {
                    LOG(std::cout << ">> Fail " << std::endl;)
                    while (stepOpt && !stepOpt->belowLimit()) {
                        stepOpt = stepsStack.pop(); // returns std::nullopt if stack became empty
                        if (stepOpt) {
                            stepOpt->fail(snapshot); // child failed => parent failed. Empty parent's Actions
                        }
                        strategy.revert();
                    }

                    // By this line we've either reverted to a valid Step on Stack
                    // or made stepOpt == std::nullopt if the stack depleted
                    if (!stepOpt) {
                        return false; // failed to construct Vector<X>
                    }
                }
            } while (stepOpt);
            // (!stepOpt) => a Step succeeded and we need to generate a subsequent
            // Location (or finish, if nothing more is to be done).
        }

        // Everything is swell, extract the Xs from snapshot
        if (!full(snapshot)) {
            std::cout << "  Snapshot not full?!" << std::endl;
            return false;
        }
        if (!unique(snapshot.xs, xs)) {
            LOG(std::cout << "  Not unique!!" << std::endl;)
            return false;
        }
        xs.push_back(snapshot.xs);

        LOG(std::cout << snapshot.funcs;)

        return true;
    }


    bool reactOnFuncUpdateAndInsert(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& coord) {
        LOG(std::cout << ">>> Inside reactOnFunc for " << coord << std::endl;)
        assert(snapshot.blocks.size() > coord.x && snapshot.blocks[coord.x].size() > coord.y);
        Block& block = snapshot.blocks[coord.x][coord.y];

        if (block.empty()) {
            LOG(std::cout << "both empty" << std::endl;)
            return true; // nothing new for this place
        }
        const unsigned int functionIndex = coord.y;
        if (block.full()) {
            LOG(std::cout << "both full, ";)
            LOG(std::cout << (snapshot.funcs[functionIndex].at(block.In()) == block.Out() ? "matched" : "no match");)
            LOG(std::cout << std::endl;)
            return snapshot.funcs[functionIndex].at(block.In()) == block.Out();
        } else if (block.emptyOut()) {
            if (const Number funcRes = snapshot.funcs[functionIndex].at(block.In());
                    !funcRes.empty()) {
                // can deduce Out from function => try placing Out
                LOG(std::cout << "Func set for Out. placing Out into Block";)
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(coord.x, coord.y, BlockCell::OUT),
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
        return false;
    }


    // returns std::nullopt if a conflict arose
    std::optional<std::vector<Action>> poke(Snapshot& snapshot, const Location& location, const Number& number) {
        LOG(std::cout << ">>> Inside Poke " << std::endl;)
        std::vector<Action> actions;
        actions.emplace_back(location, number); // self
        actions[0].apply(snapshot);
        LOG(std::cout << actions << std::endl;)
        LOG(std::cout << snapshot << std::endl;)
        LOG(std::cout << snapshot.funcs;)

        switch (location.type) {
            case LocationType::X: // is to poked set only from the outside
                if (const auto successfulPoke = pokeX(snapshot, actions, location, number);
                        !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType::Buffer:
                if (const auto successfulPoke = pokeBuffer(snapshot, actions, location, number);
                        !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType::Block:
                assert(std::holds_alternative<Coord3>(location.coord));
                switch (const Coord3 coord = std::get<Coord3>(location.coord);
                        coord.z) {
                    case BlockCell::IN:
                        if (const auto successfulPoke = pokeFuncIn(snapshot, actions, location, number);
                                !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                        break;
                    case BlockCell::OUT:
                        if (const auto successfulPoke = pokeFuncOut(snapshot, actions, location, number);
                                !successfulPoke) return (undoActions(snapshot, actions), std::nullopt);
                        break;
                }
                break;
            case LocationType::Func:
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
                    LocationType::Block,
                    Coord3(0, index, BlockCell::IN),
                    number);
                !successfulPoke) return false;

        // Secondary In
        if (hasSecondaryIn(index)) {
            assert(snapshot.blocks[0][index + NUMBERS_IN_VECTOR].emptyIn() && ":> Secondary IN not empty when X got poked.");
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(0, index + NUMBERS_IN_VECTOR, BlockCell::IN),
                        number);
                    !successfulPoke) return false;
        }

        // Buffer
        /* assert(snapshot.buffers[index].empty() && ":> Buffer not empty when X got poked."); */
        if (snapshot.buffers[index].empty()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType::Buffer,
                        Coord1(index),
                        XOR(number, snapshot.ys[index].value));
                    !successfulPoke) return false;
        } else {
            if (XOR(number, snapshot.ys[index].value) != snapshot.buffers[index].value) {
                // doesn't comply
                return false;
            }
        }

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
            actions.emplace_back(Location(LocationType::X, Coord1(index)), numberX);
            (actions.end() - 1)->apply(snapshot);
        }

        // Main In
        if (snapshot.blocks[0][index].emptyIn()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(0, index, BlockCell::IN),
                        numberX);
                    !successfulPoke) return false;
        }

        // Secondary In
        if (hasSecondaryIn(index)) {
            if (snapshot.blocks[0][index + NUMBERS_IN_VECTOR].emptyIn()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(0, index + NUMBERS_IN_VECTOR, BlockCell::IN),
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
                        LocationType::Block,
                        Coord3(lastRowIndex, index, BlockCell::OUT),
                        XOR(number, block2.Out()));
                    !successfulPoke) return false;
        } else if (block1.emptyOut() && !block2.emptyOut()) {
            if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(lastRowIndex, index + 1, BlockCell::OUT),
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
        const auto& funcRes = snapshot.funcs[columnIndex].at(number);
        if (block.Out().empty()) {
            if (!funcRes.empty()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(rowIndex, columnIndex, BlockCell::OUT),
                            funcRes);
                        !successfulPoke) return false;
            }
        } else {
            if (funcRes.empty()) {
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
                actions.emplace_back(Location(LocationType::X, Coord1(xIndex)), number);
                (actions.end() - 1)->apply(snapshot);

                // Set INs
                const unsigned int otherColumnIndex =
                    (isMainIn) ? (columnIndex + NUMBERS_IN_VECTOR)
                               : (columnIndex - NUMBERS_IN_VECTOR);
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::X,
                            Coord1(otherColumnIndex),
                            number);
                        !successfulPoke) return false;

                // Set Buffer
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::Buffer,
                            Coord1(columnIndex),
                            number);
                        !successfulPoke) return false;
            } else {
                if (snapshot.xs[xIndex] != number) {
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
                            LocationType::Block,
                            Coord3(rowIndex - 1, columnIndex, BlockCell::OUT),
                            XOR(number, block2.Out()));
                        !successfulPoke) return false;
            } else if (block1.emptyOut() && !block2.emptyOut()) {
                if (const bool successfulPoke = pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(rowIndex - 1, columnIndex + 1, BlockCell::OUT),
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
                    !funcRes.empty()) {
                if (funcRes != number) {
                    return false;
                }
            } else {
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
            if (!processPairAndInsert(
                    actions,
                    Location(LocationType::Block, Coord3(rowIndex, columnIndex + 1, BlockCell::OUT)),
                    bottom ? Location(LocationType::Buffer, Coord1(columnIndex))
                           : Location(LocationType::Block, Coord3(rowIndex + 1, columnIndex, BlockCell::IN)),
                    number)) return false;
        }
        if (columnIndex >= 1) { // right tail
            if (!processPairAndInsert(
                    actions,
                    Location(LocationType::Block, Coord3(rowIndex, columnIndex - 1, BlockCell::OUT)),
                    bottom ? Location(LocationType::Buffer, Coord1(columnIndex - 1))
                           : Location(LocationType::Block, Coord3(rowIndex + 1, columnIndex - 1, BlockCell::IN)),
                    number)) return false;
        }

        LOG(std::cout << "Succ done PokeOut" << std::endl;)
        return true;
    }


    bool setFuncAndPokeColumnAndInsert(Snapshot& snapshot, std::vector<Action>& actions, Coord2 blockCoord) {
        LOG(std::cout << ">>> Inside setFuncAndPokeAndInsert " << std::endl;)
        const auto& block = snapshot.blocks[blockCoord.x][blockCoord.y];
        assert(snapshot.funcs[blockCoord.y].emptyAt(block.In()));

        // set func
        actions.emplace_back(
                Location(LocationType::Func, Coord2(blockCoord.y, block.In().value)),
                block.Out());
        (actions.end() - 1)->apply(snapshot);

        // func for this column has beed altered => poke all Blocks in this column
        for (unsigned int row = 0; row < ((blockCoord.y <= NUMBERS_IN_VECTOR) ? (snapshot.blocks.size()) : (snapshot.blocks.size() + NUMBERS_IN_VECTOR - blockCoord.y)); row++) {
            // case for row == rowIndex is handled correctly in reactOnFuncUpdateAndInsert, so no worries
            if (const auto updateSuccessful = reactOnFuncUpdateAndInsert(snapshot, actions, Coord2(row, blockCoord.y));
                    !updateSuccessful) return false;
        }

        LOG(std::cout << ">>> Succ done setFuncAndPokeAndInsert " << std::endl;)
        return true;
    }




    void outputResult(const std::vector<Vector>& xs, const Vector& y) {
        std::cout << ":> For Y " << y << " the following Xs were generated (" << xs.size() << "):" << std::endl;
        for (const Vector& x : xs) {
            std::cout << x << std::endl;
        }
    }


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

                number = funcs[funcIndex].at(number);
            }
            for (unsigned int funcIndex = 0; funcIndex < layerSize - 1; funcIndex++) {
                layer[funcIndex] = XOR(layer[funcIndex], layer[funcIndex + 1]);
            }
        }

        for (unsigned int funcIndex = 0; funcIndex < NUMBERS_IN_VECTOR; funcIndex++) {
            layer[funcIndex] = XOR(layer[funcIndex], input[funcIndex]);
        }

        Vector output;
        for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) output[i] = layer[i];

        return output;
    }
};
// ----------------------------------------------------------------------------
void fillRemainingRandom(MagicBox& magicBox) {
    for (Func& func : magicBox.funcs) {
        for (unsigned int i = 0; i < Func::SIZE; i++) {
            if (func.emptyAt(i)) {
                // TODO ???
                func.update(i, generateRandomUniformInt(0, Func::SIZE - 1));
            }
        }
    }
}
// ----------------------------------------------------------------------------
int main() {
    std::cout << "--------------------BEGIN----------------------" << std::endl << std::endl;
    srand (309);

    MagicBox mb;
    const Vector y = generateRandomVector();
    const auto xs = mb.work(y);
    fillRemainingRandom(mb);


    std::cout << std::endl << "Checking all possible Xs..." << std::endl;
    unsigned int counter = 0;
    std::array<unsigned int, NUMBERS_IN_VECTOR> nextToUse;
    for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) nextToUse[i] = 0;
    while (true) {
        // Construct X
        Vector x{};
        for (unsigned int i = 0; i < NUMBERS_IN_VECTOR; i++) {
            x[i] = nextToUse[i];
        }

        // Check
        if (mb.apply(x) == y) {
            std::cout << x << " -> " << y << std::endl;
            counter++;
        }

        // Prepare for next iteration
        nextToUse[NUMBERS_IN_VECTOR - 1]++;
        for (unsigned int i = NUMBERS_IN_VECTOR - 1; i > 0; i--) {
            if (nextToUse[i] == NUMBERS_AMOUNT) {
                nextToUse[i-1]++;
                nextToUse[i] = 0;
            }
        }
        if (nextToUse[0] == NUMBERS_AMOUNT) {
            break;
        }
    }
    std::cout << "The model has (" << counter << ") Xs total" << std::endl;


    std::cout << std::endl << "---------------------END-----------------------" << std::endl;
    return 0;
}
