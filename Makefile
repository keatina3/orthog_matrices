CC = gcc

CFLAGS = -g -Wall

LDFLAGS = 

OBJECTS = prog.o utils.o dense_orthogon.o

prog: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) $(OBJECTS) prog
