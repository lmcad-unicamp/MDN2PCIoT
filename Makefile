CC = g++
FLAGS = -O3 -fopenmp
DIR = src
DIROBJ = obj
SRC = $(DIR)/SourceGraph.cpp $(DIR)/TargetGraph.cpp $(DIR)/ReadInputFiles.cpp $(DIR)/GenericPartitioningAlgorithm.cpp $(DIR)/Comm.cpp $(DIR)/main.cpp $(DIR)/InferenceRate.cpp
OBJ = $(SRC:.c=.o)
DIRWIF = write-input-file

KLP: $(DIR)/main.cpp $(OBJ)
	g++ -c $(FLAGS) $(DIR)/SourceGraph.cpp    			     -o $(DIROBJ)/SourceGraph.o
	g++ -c $(FLAGS) $(DIR)/TargetGraph.cpp    			     -o $(DIROBJ)/TargetGraph.o
	g++ -c $(FLAGS) $(DIR)/ReadInputFiles.cpp 			     -o $(DIROBJ)/ReadInputFiles.o
	g++ -c $(FLAGS) $(DIR)/GenericPartitioningAlgorithm.cpp -o $(DIROBJ)/GenericPartitioningAlgorithm.o
	g++ -c $(FLAGS) $(DIR)/Comm.cpp	 			     -o $(DIROBJ)/Comm.o
	g++ -c $(FLAGS) $(DIR)/InferenceRate.cpp	 			     -o $(DIROBJ)/InferenceRate.o
	g++ -c $(FLAGS) $(DIR)/main.cpp						 -o $(DIROBJ)/main.o
	g++    $(FLAGS) $(OBJ)				             -o DN2PCIoT

KLP-debug: 
	g++ -c -O0 -g -fopenmp $(DIR)/SourceGraph.cpp    			   -o $(DIROBJ)/SourceGraph.o
	g++ -c -O0 -g -fopenmp $(DIR)/TargetGraph.cpp    			   -o $(DIROBJ)/TargetGraph.o
	g++ -c -O0 -g -fopenmp $(DIR)/ReadInputFiles.cpp 			   -o $(DIROBJ)/ReadInputFiles.o
	g++ -c -O0 -g -fopenmp $(DIR)/GenericPartitioningAlgorithm.cpp -o $(DIROBJ)/GenericPartitioningAlgorithm.o
	g++ -c -O0 -g -fopenmp $(DIR)/Comm.cpp	 			   -o $(DIROBJ)/Comm.o
	g++ -c -O0 -g -fopenmp $(DIR)/InferenceRate.cpp	 			     -o $(DIROBJ)/InferenceRate.o
	g++    -O0 -g -fopenmp $(OBJ)	   	     			   -o DN2PCIoT-debug

#prog.x: main.o math.o log.o
#	ld -dynamic-linker /lib64/ld-linux-x86-64.so.2 /usr/lib/x86_64-linux-O0 -gnu/crt1.o /usr/lib/x86_64-linux-O0 -gnu/crti.o -L/usr/lib64 main.o math.o log.o -lc /usr/lib/x86_64-linux-O0 -gnu/crtn.o -o prog.x
	
clean: 
	rm -f $(DIROBJ)/*.o DN2PCIoT DN2PCIoT-debug $(DIRWIF)/*.o writeInputFile*

writeInputFile_LeNet_1x1: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet_1x1.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet_1x1.c   -o $(DIRWIF)/write_input_file_LeNet_1x1.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet_1x1.o -o writeInputFile_LeNet_1x1

writeInputFile_LeNet_2x2: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet_2x2.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet_2x2.c   -o $(DIRWIF)/write_input_file_LeNet_2x2.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet_2x2.o -o writeInputFile_LeNet_2x2

writeInputFileSBAC: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.c   -o $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.o -o writeInputFile_LeNet-SBAC-PAD-2018
