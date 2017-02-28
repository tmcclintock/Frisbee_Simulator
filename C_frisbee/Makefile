
CC = gcc
ifdef ALONE
ifeq ($(ALONE),yes)
$(info Building executable)
OBJS = glue.o driver.o rk4.o equations_of_motion.o coefficient_model.o
EXEC = test.exe
CFLAGS = -g
OFLAGS = 
endif
else
$(info Building shared library)
OBJS = driver.o rk4.o equations_of_motion.o coefficient_model.o
EXEC = c_frisbee.so
CFLAGS = -fPIC
OFLAGS = -shared
endif

INCL =

LIBS = -lm

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~ *.exe *.so