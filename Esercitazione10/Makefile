CXX       = mpic++
CXXFLAGS  = -std=c++11 -O2
SRCS      = random.cpp utils.cpp ga.cpp main.cpp
OBJS      = $(SRCS:.cpp=.o)
EXEC      = tsp_mpi

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS)

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp

utils.o: utils.cpp utils.h random.h
	$(CXX) $(CXXFLAGS) -c utils.cpp

ga.o: ga.cpp ga.h utils.h random.h
	$(CXX) $(CXXFLAGS) -c ga.cpp

main.o: main.cpp ga.h utils.h random.h
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(OBJS) $(EXEC) results_* cities_* best_path_*
