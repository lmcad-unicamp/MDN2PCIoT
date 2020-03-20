CC = g++
FLAGS = -O3
DIR = src
DIROBJ = obj
SRC = $(DIR)/SourceGraph.cpp $(DIR)/TargetGraph.cpp $(DIR)/ReadInputFiles.cpp $(DIR)/GenericPartitioningAlgorithm.cpp $(DIR)/Comm.cpp $(DIR)/main.cpp
OBJ = $(SRC:.c=.o)

run: $(DIR)/main.cpp $(OBJ)
	g++ -c $(FLAGS) $(DIR)/SourceGraph.cpp    			     -o $(DIROBJ)/SourceGraph.o
	g++ -c $(FLAGS) $(DIR)/TargetGraph.cpp    			     -o $(DIROBJ)/TargetGraph.o
	g++ -c $(FLAGS) $(DIR)/ReadInputFiles.cpp 			     -o $(DIROBJ)/ReadInputFiles.o
	g++ -c $(FLAGS) $(DIR)/GenericPartitioningAlgorithm.cpp -o $(DIROBJ)/GenericPartitioningAlgorithm.o
	g++ -c $(FLAGS) $(DIR)/Comm.cpp	 			     -o $(DIROBJ)/Comm.o
	g++ -c $(FLAGS) $(DIR)/main.cpp						 -o $(DIROBJ)/main.o
	g++    $(FLAGS) $(OBJ)				             -o KLP

run-debug: 
	g++ -c -g SourceGraph.cpp    			   -o SourceGraph.o
	g++ -c -g TargetGraph.cpp    			   -o TargetGraph.o
	g++ -c -g ReadInputFiles.cpp 			   -o ReadInputFiles.o
	g++ -c -g GenericPartitioningAlgorithm.cpp -o GenericPartitioningAlgorithm.o
	g++ -c -g Comm.cpp	 			   -o Comm.o
	g++    -g $(OBJ)	   	     			   -o KLP

run-thirty: KLP
	for number in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30; do \
		./KLP $$number ; \
	done
                            
#prog.x: main.o math.o log.o
#	ld -dynamic-linker /lib64/ld-linux-x86-64.so.2 /usr/lib/x86_64-linux-gnu/crt1.o /usr/lib/x86_64-linux-gnu/crti.o -L/usr/lib64 main.o math.o log.o -lc /usr/lib/x86_64-linux-gnu/crtn.o -o prog.x
	
clean: 
	rm -f *.o KLP
