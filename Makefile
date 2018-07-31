# The name of the main file and executable
mainFileName = main
# Files that have .h and .cpp versions
classFiles =
# Files that only have the .h version
justHeaderFiles =
# Compilation flags
OPTIMIZATION_FLAG = -O0
COMPILER_GLAGS = -Wall
LINKER_FLAGS =


# Auxiliary
filesObj = $(addsuffix .o, $(mainFileName) $(classFiles))
filesH = $(addsuffix .h, $(classFiles) $(justHeaderFiles))


all: cleanExe $(mainFileName)


# Compiler
%.o: %.cpp $(filesH)
	g++ $(COMPILER_GLAGS) $(OPTIMIZATION_FLAG) -c $<

# Linker
$(mainFileName): $(filesObj)
	g++ $(COMPILER_GLAGS) $(OPTIMIZATION_FLAG) $(LINKER_FLAGS) $^ -o $@


# Utils
clean:
	rm -f *.o $(mainFileName)

cleanExe:
	rm -f $(mainFileName)
