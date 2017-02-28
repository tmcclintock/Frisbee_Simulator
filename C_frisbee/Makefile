EXEC = test.exe

OBJS = glue.o equations_of_motion.o coefficient_model.o

CC = gcc

CFLAGS =

OFLAGS = -lm

INCL =

LIBS =

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~