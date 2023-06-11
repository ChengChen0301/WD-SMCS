CXX = g++  
CXXFLAGS='-std=c++11'
LIBS= -I/opt/software/boost-1-58-0/include/ -I/usr/include/gsl -L/opt/software/boost-1-58-0/lib/ -L/usr/lib64/ -lm -lgsl -lgslcblas
#OMPFLAGS = -fopenmp


prom = main
deps = simulation.h vec.h variables.h 
obj = main.o simulation.o vec.o variables.o


$(prom): $(obj)
	$(CXX) ${CPP_FLAGS} -o $(prom) $(obj)

debug: DEBUG = -DDEBUG
debug: smcs

%.o: %.c $(deps)
	$(CXX) -c $< -o $@



