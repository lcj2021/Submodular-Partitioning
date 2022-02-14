CC = g++
CFLAG = -std=c++11 -g

objdir = 

exedir = 

TARGET = $(exedir)sp
all: $(TARGET)
	@echo "Compilation ends successfully!"

$(exedir)sp: $(objdir)sp.cpp
	$(CC) $(CFLAG) sp.cpp -o $(exedir)sp
	rm -rf *.o


