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
#include <cmath>
#include <iomanip>
#include <algorithm>


// 16-16: Ff=0.5 [4075] in 2:03
// 14-18: Ff=0.5 [18500] in 25:30
// ----------------------------------------------------------------------------
static const bool DO_LOG = false;
#define LOG(x) if(DO_LOG) {x}

static const bool DO_CHECK = false;
// ----------------------------------------------------------------------------
static const unsigned int NUMBERS_IN_VECTOR = 8;
static const unsigned int NUMBER_BITS = 9;

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

/* typedef std::array<Number, NUMBERS_IN_VECTOR> Vector; */
using Vector = std::array<Number, NUMBERS_IN_VECTOR>;

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

inline float interpolate(float a, float b, float f) {
    return a + f * (b - a);
}
inline float avg(long long* arr, unsigned int length) {
    long long sum = 0LL;
    long long max = 0LL;
    unsigned int size = 0;
    for (unsigned int i = 0; i < length; i++) {
        if (arr[i] > max) max = arr[i];
        if (arr[i] != 0LL) size++;
        sum += arr[i];
    }
    if (size == length) {
        sum -= max;
        size--;
    }

    return (1.0f * sum) / size;
}
// ----------------------------------------------------------------------------
// [low, high]
// template<typename T = unsigned int>
// T generateRandomUniformInt(T low, T high) {
//     static std::random_device rd;
//     static std::seed_seq seed{1, 2, 3, 300};
//     static std::mt19937 e2(seed);
//     static std::uniform_int_distribution<T> dist(low, high);
//
//     return dist(e2);
// }

unsigned int generateRandomUniformInt(unsigned int low, unsigned int high) {
    const auto range = high - low + 1;
    return rand() % range + low;
}

Number generateRandomNumber() {
    return generateRandomUniformInt(0U, NUMBER_MAX);
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
        Number* map;

        std::vector<Number> usedCells;

    public:
        Func() {
            map = new Number[SIZE];
            for (unsigned int i = 0; i < SIZE; i++) map[i] = Number(EMPTY_NUMBER);
        }

        virtual ~Func() {
            delete[] map;
        }

        bool setAt(const Number& arg) const {
            return (!map[arg()].empty());
        }

        bool emptyAt(const Number& arg) const {
            return (map[arg()].empty());
        }

        void update(const Number& arg, const Number& value) {
            map[arg()] = value;

            if (value.empty()) {
                if (const auto it = std::find(usedCells.begin(), usedCells.end(), arg());
                        it != usedCells.end()) {
                    usedCells.erase(it);
                }
            } else {
                usedCells.emplace_back(arg());
            }
        }

        std::optional<Number> getRandomUsed() const {
            if (usedCells.empty()) return std::nullopt;
            return { usedCells[generateRandomUniformInt(0, usedCells.size() - 1)] };
        }

        const Number& at(const Number& arg) const {
            return map[arg()];
        }

        const Number& operator[](const Number& arg) const {
            return map[arg()];
        }
};

/* typedef std::array<Func, FUNCS_AMOUNT> Functions; */
using Functions = std::array<Func, FUNCS_AMOUNT>;

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
        stream << x << " ";
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
        stream << b << " ";
    }
    stream << std::endl;
    for (const auto& y : snapshot.ys) {
        stream << y << " ";
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
        std::set<unsigned int> triedNumbers;

    public:
        Step(const Location& location)
            : baseLocation(location) {
        }

        void useNumber(const Number& number) {
            triedNumbers.insert(number());
        }

        void applyActions(const std::vector<Action>& newActions) {
            actions.insert(actions.end(), newActions.begin(), newActions.end());
        }

        void revert(Snapshot& snapshot) {
            for (Action& action : actions) {
                action.undo(snapshot);
            }
            actions.clear();
        }

        Location getLocation() const {
            return baseLocation;
        }

        bool belowLimit(float Ff, float Fs) const {
            const float min = interpolate(NUMBERS_AMOUNT / NUMBER_BITS, NUMBER_BITS, std::sqrt(Ff));
            const float max = interpolate(NUMBERS_AMOUNT, NUMBERS_AMOUNT / NUMBER_BITS, Ff);
            const unsigned int retries = interpolate(min, max, 1.0f - std::sqrt(std::sqrt(1.0f - Fs)));
            return triedNumbers.size() < retries;
            /* return (triedNumbers.size() < ALLOWED_RETRIES); */
        }

        std::optional<Number> getUntriedNumber(float Ff, float Fs, const Functions& funcs) const {
            if (!belowLimit(Ff, Fs)) return std::nullopt;

            /* if (baseLocation.type == LocationType::Block) { */
            /*     assert(std::holds_alternative<Coord3>(baseLocation.coord)); */
            /*     const unsigned int funcIndex = std::get<Coord3>(baseLocation.coord).y; */
            /*     if (generateRandomUniformInt(0,10) < 8) { */
            /*         const auto numberOpt = funcs[funcIndex].getRandomUsed(); */
            /*         if (numberOpt && (triedNumbers.find((*numberOpt)()) == triedNumbers.end())) { */
            /*             #<{(| std::cout << "Reusing" << std::endl; |)}># */
            /*             return numberOpt; */
            /*         } */
            /*     } */
            /* } */

            const unsigned int untriedCount = NUMBERS_AMOUNT - triedNumbers.size();
            const unsigned int index = generateRandomUniformInt(0U, untriedCount - 1);
            for (unsigned int i = 0, count = 0; i < NUMBERS_AMOUNT; i++) {
                if (triedNumbers.find(i) != triedNumbers.end()) continue;
                if (count++ == index) {
                    return { Number(i) };
                }
            }
            return std::nullopt;
        }

        friend std::ostream& operator<<(std::ostream& stream, const Step& step);
};

std::ostream& operator<<(std::ostream& stream, const Step& step) {
    stream << "Step(" << step.baseLocation << "," << step.actions << "," << step.triedNumbers.size() << ")";

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

        void undoAll(Snapshot& snapshot) {
            while (auto stepOpt = pop()) {
                stepOpt->revert(snapshot);
            }
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
/* typedef std::function<std::optional<Location>(unsigned int index, const Snapshot& snapshot, const Functions& funcs)> StrategyFunction; */
using StrategyFunction = std::function<std::optional<Location>(unsigned int index, const Snapshot& snapshot, const Functions& funcs)>;

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

    MagicBox() {
        for (unsigned int i = 0; i < prevStepsDone.size(); i++) prevStepsDone[i] = 0;
        /* for (unsigned int i = 0; i < lastTriesCount; i++) { */
        /*     lastTries[i] = 0; */
        /* } */
    }

private:
    bool hasSecondaryIn(unsigned int index) const {
        return (index < NUMBERS_IN_VECTOR - 1);
    }

    const unsigned int allowedSubsequentDropsCount = 4;
    unsigned int subsequentDrops = 0;
    bool lastWasDropped = false;
    std::array<unsigned int, 3> prevStepsDone;
    /* static const unsigned int lastTriesCount = NUMBER_BITS; */
    /* unsigned int lastTryIndex = 0; */
    /* long long lastTries[lastTriesCount]; */

    /* float getAverageTries() { */
    /*     return avg(lastTries, lastTriesCount); */
    /* } */

    float calculateFf(const Functions& funcs) const {
        unsigned int fill = 0;
        for (const Func& func : funcs) {
            for (unsigned int i = 0; i < Func::SIZE; i++) {
                if (!func.emptyAt(i)) fill++;
            }
        }

        return 1.0f * fill / (FUNCS_AMOUNT * Func::SIZE);
    }

    float calculateFs(const Snapshot& snapshot) const {
        unsigned int fill = 0;
        unsigned int size = 0;
        for (const Number& n : snapshot.xs) {
            if (!n.empty()) fill++;
            size++;
        }
        for (const Number& n : snapshot.buffers) {
            if (!n.empty()) fill++;
            size++;
        }
        for (const auto& row : snapshot.blocks) {
            for (const Block& b : row) {
                if (!b.emptyIn()) fill++;
                if (!b.emptyOut()) fill++;
                size++;
                size++;
            }
        }

        return 1.0f * fill / size;
    }

public:
    std::vector<Vector> work(const Vector& y) {
        std::vector<Vector> xs;
        while (generateNewX(xs, y));
        // outputResult(xs, y);

        return xs;
    }

    bool generateNewX(std::vector<Vector>& xs, const Vector& y) {;
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
                                    blockCoord.x, // row
                                    blockCoord.y, // column
                                    snapshot.blocks[blockCoord.x][blockCoord.y].emptyIn()
                                        ? BlockCell::IN
                                        : BlockCell::OUT)}};
                } else {
                    return std::nullopt;
                }
            }
        });

        const float Ff = calculateFf(funcs);
        unsigned int stepsDone = 0;
        const unsigned int prevStepsAvg = [](const decltype(prevStepsDone)& prevStepsDone){
            unsigned int avg = 0;
            for (unsigned int i = 0; i < prevStepsDone.size(); i++) {
                if (prevStepsDone[i] != 0) {
                    avg += prevStepsDone[i];
                }
            }
            return avg / prevStepsDone.size();
        }(prevStepsDone);
        while (const auto locationOpt = strategy.getNext(snapshot, funcs)) {
            const Location location = *locationOpt;
            LOG(std::cout << ">> Got next strategy: " << location << "." << std::endl;)

            std::optional<Step> stepOpt = {Step(location)};
            do {
                stepsDone++;
                if (stepsDone > 5 * prevStepsAvg && prevStepsAvg > 0) {
                    if (lastWasDropped) {
                        subsequentDrops++;
                    } else {
                        lastWasDropped = true;
                        subsequentDrops = 1;
                    }
                    std::cout << ":> Timed out(" << stepsDone << "). Dropping("
                        << subsequentDrops << ")" << std::endl;

                    if (subsequentDrops > allowedSubsequentDropsCount) {
                        std::cout << "Too subsequent many drops. Terminating" << std::endl;
                        return false;
                    } else {
                        stepsStack.undoAll(snapshot); // there would remain useless leftovers in funcs otherwise
                        return true;
                    }
                }

                const float Fs = calculateFs(snapshot);
                Step& step = *stepOpt;
                const std::optional<Number> numberOpt = step.getUntriedNumber(Ff, Fs, funcs);

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
                    while (!stepOpt->belowLimit(calculateFf(funcs), calculateFs(snapshot))) {
                        stepOpt = stepsStack.pop(); // returns std::nullopt if stack became empty
                        if (!stepOpt) { // stack depleted
                            std::cout << ":> Finished because the stack has depleted." << std::endl;
                            return false; // failed to construct Vector<X>
                        }
                        stepOpt->revert(snapshot); // child failed => parent failed. Empty parent's Actions
                        strategy.revert();
                    }
                }
            } while (stepOpt);
            // (!stepOpt) => a Step succeeded and we need to generate a subsequent
            // Location (or finish, if nothing more is to be done).
        }

        for (unsigned int i = 0; i < prevStepsDone.size() - 1; i++) prevStepsDone[i+1] = prevStepsDone[i];
        prevStepsDone[0] = stepsDone;

        // Everything is swell, extract the Xs from snapshot
        assert(full(snapshot) && "  Snapshot not full?!");
        if (!unique(snapshot.xs, xs)) {
            static unsigned int duplicatesCounter = 0;
            std::cout << ":> Duplicate generated (" << ++duplicatesCounter << "). Skipping." <<  std::endl;
            static const unsigned int duplicatesLimit = 20;
            if (duplicatesCounter > duplicatesLimit) {
                std::cout << ":> Finished because duplicates limit has beed exceeded." << std::endl;
                return false;
            } else {
                stepsStack.undoAll(snapshot); // there would remain useless leftovers in funcs otherwise
                return true;
            }
        }
        xs.push_back(snapshot.xs);

        LOG(std::cout << snapshot.funcs;)

        static auto co = 0;
        std::cout << "New x(" << ++co << ")"
            << " Ff={" << calculateFf(snapshot.funcs) << "}"
            << " Steps=[" << stepsDone << "]"
            << std::endl;

        lastWasDropped = false;
        return true;
    }


    bool reactOnFuncUpdateAndInsert(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& coord) {
        LOG(std::cout << ">>> Inside reactOnFunc for " << coord << std::endl;)
        const unsigned int row = coord.x;
        const unsigned int column = coord.y;
        assert(row < snapshot.blocks.size() && column < snapshot.blocks[coord.x].size());
        Block& block = snapshot.blocks[row][column];

        if (block.empty()) {
            LOG(std::cout << "both empty" << std::endl;)
            return true; // nothing to change for this place
        }

        const unsigned int functionIndex = column;
        if (const Number funcRes = snapshot.funcs[functionIndex].at(block.In());
                block.full()) {
            LOG(std::cout << "both full, ";)
            LOG(std::cout << (snapshot.funcs[functionIndex].at(block.In()) == block.Out() ? "matched" : "no match");)
            LOG(std::cout << std::endl;)
            return (funcRes == block.Out());
        } else if (block.emptyOut()) {
            if (!funcRes.empty()) {
                // can deduce Out from function => try placing Out
                LOG(std::cout << "Func set for Out. placing Out into Block";)
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(row, column, BlockCell::OUT),
                            funcRes)) return false;
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
            const std::variant<Coord1, Coord2, Coord3>& coord,
            const Number& value) {
        LOG(std::cout << "Inside pokeAndInsert " << locationType << std::endl;)
        if (const auto& actionsOpt = poke(snapshot, Location(locationType, coord), value)) {
            LOG(std::cout << "Succ done pokeAndInsert " << locationType << std::endl;)
            actions.insert(actions.end(), actionsOpt->begin(), actionsOpt->end());
            return true;
        }

        LOG(std::cout << "Failed pokeAndInsert " << locationType << std::endl;)
        return false;
    }


    // returns std::nullopt if a conflict arose
    // Undoes all actions that it caused in case of a conflict
    std::optional<std::vector<Action>> poke(Snapshot& snapshot, const Location& location, const Number& number) {
        std::vector<Action> actions;
        actions.emplace_back(location, number); // self
        actions[0].apply(snapshot);
        {
            LOG(std::cout << ">>> Inside Poke " << std::endl;)
            LOG(std::cout << "Just done: " << actions << std::endl;)
            LOG(std::cout << snapshot << std::endl;)
            LOG(std::cout << snapshot.funcs;)
        }

        switch (location.type) {
            case LocationType::X: // is to be poked only from the outside
                assert(std::holds_alternative<Coord1>(location.coord));
                if (!pokeX(snapshot, actions, std::get<Coord1>(location.coord), number))
                    return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType::Buffer:
                assert(std::holds_alternative<Coord1>(location.coord));
                if (!pokeBuffer(snapshot, actions, std::get<Coord1>(location.coord), number))
                    return (undoActions(snapshot, actions), std::nullopt);
                break;
            case LocationType::Block:
                assert(std::holds_alternative<Coord3>(location.coord));
                switch (const Coord3 coord = std::get<Coord3>(location.coord);
                        coord.z) {
                    case BlockCell::IN:
                        if (!pokeFuncIn(snapshot, actions, Coord2(coord.x, coord.y), number))
                            return (undoActions(snapshot, actions), std::nullopt);
                        break;
                    case BlockCell::OUT:
                        if (!pokeFuncOut(snapshot, actions, Coord2(coord.x, coord.y), number))
                            return (undoActions(snapshot, actions), std::nullopt);
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


    bool pokeX(Snapshot& snapshot, std::vector<Action>& actions, const Coord1& coord, const Number& number) {
        LOG(std::cout << ">>> Inside PokeX " << std::endl;)
        const unsigned int index = coord.x;

        // Main In
        if (const Block& block = snapshot.blocks[0][index];
                block.emptyIn()) {
            if (!pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(0, index, BlockCell::IN),
                        number)) return false;
        } else if (block.In() != number) return false;

        // Secondary In
        if (hasSecondaryIn(index)) {
            if (const Block& block = snapshot.blocks[0][index + NUMBERS_IN_VECTOR];
                    block.emptyIn()) {
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(0, index + NUMBERS_IN_VECTOR, BlockCell::IN),
                            number)) return false;
            } else if (block.In() != number) return false;
        }

        // Buffer
        const Number expectedBufferValue = XOR(number, snapshot.ys[index].value);
        if (const Number& buffer = snapshot.buffers[index];
                buffer.empty()) {
            if (!pokeAndInsert(snapshot, actions,
                        LocationType::Buffer,
                        Coord1(index),
                        expectedBufferValue)) return false;
        } else if (expectedBufferValue != buffer.value) return false;

        LOG(std::cout << ">>> Succ done PokeX" << std::endl;)
        return true;
    }

    bool pokeBuffer(Snapshot& snapshot, std::vector<Action>& actions, const Coord1& coord, const Number& number) {
        LOG(std::cout << ">>> Inside PokeBuffer " << std::endl;)
        const unsigned int index = coord.x;
        const Number expectedValueX = XOR(number, snapshot.ys[index]);

        // X
        if (const Number& numberX = snapshot.xs[index];
                numberX.empty()) {
            // just set X, but never poke it
            actions.emplace_back(Location(LocationType::X, Coord1(index)), expectedValueX);
            (actions.end() - 1)->apply(snapshot);
            /* if (!canContinue(snapshot)) return false; */
        } else if (expectedValueX != numberX) return false;

        // Main In
        if (const Number& in = snapshot.blocks[0][index].In();
                in.empty()) {
            if (!pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(0, index, BlockCell::IN),
                        expectedValueX)) return false;
        } else if (in != expectedValueX) return false;

        // Secondary In
        if (hasSecondaryIn(index)) {
            if (const Number& in = snapshot.blocks[0][index + NUMBERS_IN_VECTOR].In();
                    in.empty()) {
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(0, index + NUMBERS_IN_VECTOR, BlockCell::IN),
                            expectedValueX)) return false;
            } else if (in != expectedValueX) return false;
        }

        // Outs
        const unsigned int lastRowIndex = snapshot.blocks.size() - 1;
        const Block& block1Const = snapshot.blocks[lastRowIndex][index];
        const Block& block2Const = snapshot.blocks[lastRowIndex][index + 1];
        if (!block1Const.emptyOut() && !block2Const.emptyOut()) {
            if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                return false;
            }
        } else if (block1Const.emptyOut() && !block2Const.emptyOut()) {
            if (!pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(lastRowIndex, index, BlockCell::OUT),
                        XOR(number, block2Const.Out()))) return false;
        } else if (!block1Const.emptyOut() && block2Const.emptyOut()) {
            if (!pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(lastRowIndex, index + 1, BlockCell::OUT),
                        XOR(number, block1Const.Out()))) return false;
        }

        LOG(std::cout << ">>> Succ done PokeBuffer" << std::endl;)
        return true;
    }

    bool pokeFuncIn(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& coord, const Number& number) {
        LOG(std::cout << ">>> Inside PokeIn " << std::endl;)
        const unsigned int rowIndex = coord.x;
        const unsigned int columnIndex = coord.y;

        // Self Out
        const auto& funcRes = snapshot.funcs[columnIndex].at(number);
        if (const Number& out = snapshot.blocks[rowIndex][columnIndex].Out();
                out.empty()) {
            if (!funcRes.empty() && !pokeAndInsert(snapshot, actions,
                        LocationType::Block,
                        Coord3(rowIndex, columnIndex, BlockCell::OUT),
                        funcRes)) return false;
        } else {
            if (funcRes.empty()) {
                if (!setFuncAndPokeColumnAndInsert(snapshot, actions, coord)) return false;
            } else if (funcRes != out) return false;
        }

        if (rowIndex == 0) {
            const bool mainIn = columnIndex < NUMBERS_IN_VECTOR;
            const unsigned int xIndex = (mainIn) ? columnIndex
                                                 : (columnIndex - NUMBERS_IN_VECTOR);
            // X
            if (snapshot.xs[xIndex].empty()) {
                // Set X (just set, never poke)
                actions.emplace_back(Location(LocationType::X, Coord1(xIndex)), number);
                (actions.end() - 1)->apply(snapshot);
            } else if (snapshot.xs[xIndex] != number) return false;

            // Secondary In
            if (hasSecondaryIn(xIndex)) {
                const unsigned int otherColumnIndex =
                    (mainIn) ? (columnIndex + NUMBERS_IN_VECTOR)
                             : (columnIndex - NUMBERS_IN_VECTOR);
                if (const Number& in = snapshot.blocks[0][otherColumnIndex].In();
                        in.empty()) {
                    if (!pokeAndInsert(snapshot, actions,
                                LocationType::Block,
                                Coord3(0, otherColumnIndex, BlockCell::IN),
                                number)) return false;
                } else if (in != number) return false;
            }

            // Buffer
            const Number expectedValueBuffer = XOR(number, snapshot.ys[xIndex]);
            if (const Number& buffer = snapshot.buffers[xIndex];
                    buffer.empty()){
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Buffer,
                            Coord1(xIndex),
                            expectedValueBuffer)) return false;
            } else if (buffer != expectedValueBuffer) return false;
        } else { // rowIndex >= 1
            // foreign Outs above self
            const Block& block1Const = snapshot.blocks[rowIndex - 1][columnIndex];
            const Block& block2Const = snapshot.blocks[rowIndex - 1][columnIndex + 1];
            if (!block1Const.emptyOut() && !block2Const.emptyOut()) {
                if (XOR(block1Const.Out(), block2Const.Out()) != number) {
                    return false;
                }
            } else if (block1Const.emptyOut() && !block2Const.emptyOut()) {
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(rowIndex - 1, columnIndex, BlockCell::OUT),
                            XOR(number, block2Const.Out()))) return false;
            } else if (!block1Const.emptyOut() && block2Const.emptyOut()) {
                if (!pokeAndInsert(snapshot, actions,
                            LocationType::Block,
                            Coord3(rowIndex - 1, columnIndex + 1, BlockCell::OUT),
                            XOR(number, block1Const.Out()))) return false;
            }
        }

        LOG(std::cout << ">>> Succ done PokeIn" << std::endl;)
        return true;
    }

    bool pokeFuncOut(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& coord, const Number& number) {
        LOG(std::cout << ">>> Inside PokeOut " << std::endl;)
        const unsigned int rowIndex = coord.x;
        const unsigned int columnIndex = coord.y;

        // Self In
        if (const Number& in = snapshot.blocks[rowIndex][columnIndex].In();
                !in.empty()) {
            if (const auto& funcRes = snapshot.funcs[columnIndex].at(in);
                    !funcRes.empty()) {
                if (funcRes != number) return false;
            } else if (!setFuncAndPokeColumnAndInsert(snapshot, actions, coord)) return false;
        }
        // too unlikely that func is full and exactly 1 input has us as output
        // (then we could derive In), so don't even bother here.

        const static auto processPairAndInsert = [&](std::vector<Action>& actions,
                const Location& location1, const Location& location2, const Number& number) {
            if (location1.empty(snapshot)) {
                if (!location2.empty(snapshot)) {
                    if (!pokeAndInsert(snapshot, actions,
                                location1.type,
                                location1.coord,
                                XOR(location2.value(snapshot), number))) return false;
                }
            } else {
                if (location2.empty(snapshot)) {
                    if (!pokeAndInsert(snapshot, actions,
                                location2.type,
                                location2.coord,
                                XOR(location1.value(snapshot), number))) return false;
                } else if (XOR(location1.value(snapshot), location2.value(snapshot)) != number) return false;
            }

            return true;
        };

        // Connections through XOR
        const unsigned int rowSize = snapshot.blocks[rowIndex].size();
        const unsigned int lastRowIndex = snapshot.blocks.size() - 1;
        const bool bottom = rowIndex == lastRowIndex;
        if (columnIndex < rowSize - 1) { // left tail
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


    bool setFuncAndPokeColumnAndInsert(Snapshot& snapshot, std::vector<Action>& actions, const Coord2& blockCoord) {
        LOG(std::cout << ">>> Inside setFuncAndPokeAndInsert " << std::endl;)
        const auto& block = snapshot.blocks[blockCoord.x][blockCoord.y];
        assert(snapshot.funcs[blockCoord.y].emptyAt(block.In()));

        // set func
        actions.emplace_back(
                Location(LocationType::Func, Coord2(blockCoord.y, block.In().value)),
                block.Out());
        (actions.end() - 1)->apply(snapshot);

        // func for this column has beed altered => poke all Blocks in this column
        const unsigned int rowsAmount = (blockCoord.y <= NUMBERS_IN_VECTOR)
            ? (snapshot.blocks.size())
            : (snapshot.blocks.size() + NUMBERS_IN_VECTOR - blockCoord.y);
        for (unsigned int row = 0; row < rowsAmount; row++) {
            // case for row == rowIndex is handled correctly in reactOnFuncUpdateAndInsert, so no worries
            if (!reactOnFuncUpdateAndInsert(snapshot, actions, Coord2(row, blockCoord.y))) return false;
        }

        LOG(std::cout << ">>> Succ done setFuncAndPokeAndInsert " << std::endl;)
        return true;
    }


    void outputResult(const std::vector<Vector>& xs, const Vector& y) {
        std::cout << ":> For Y " << y << " the following Xs were generated (" << xs.size() << "):" << std::endl;
        unsigned int counter = 0;
        for (const Vector& x : xs) {
            std::cout << "[" << ++counter << "]: " << x << std::endl;
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
                assert(!funcs[funcIndex].emptyAt(number) && ":> Usage of empty func entry.");
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
void fillRemainingRandom(Functions& funcs) {
    for (Func& func : funcs) {
        for (unsigned int i = 0; i < Func::SIZE; i++) {
            if (func.emptyAt(i)) {
                func.update(i, generateRandomUniformInt(0U, NUMBER_MAX));
            }
        }
    }
}
// ----------------------------------------------------------------------------
int main() {
    std::cout << "--------------------BEGIN----------------------" << std::endl << std::endl;
    srand(410);


    MagicBox mb;
    const Vector y = generateRandomVector();
    const auto xs = mb.work(y);
    fillRemainingRandom(mb.funcs);

    if (DO_CHECK) {
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
                std::cout << "[" << ++counter<< "]: " << x << " -> " << y << std::endl;
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
    }


    std::cout << std::endl << "---------------------END-----------------------" << std::endl;
    return 0;
}
