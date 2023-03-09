CFLAGS = ${shell pkg-config opencv --cflags}
LDFLAGS = ${shell pkg-config opencv --libs}

SOURCES = smcs.cpp simulation.cpp simulation.h variables.cpp variables.h vec.cpp vec.h
OBJECTS = smcs.o simulation.o variables.o vec.o
TARGET = smcs

$(TARGET):$(OBJECTS)
	g++ -std=c++11 -o $(TARGET) $(OBJECTS) ${CFLAGS} ${LDFLAGS}
$(OBJECTS):$(SOURCES)
	g++ -std=c++11 -c $(SOURCES) ${CFLAGS} ${LDFLAGS}
clean:
	rm *.gch $(OBJECTS) $(TARGET)


