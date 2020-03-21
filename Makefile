CC = g++
FLAGS = -O3 -fopenmp
DIR = src
DIROBJ = obj
SRC = $(DIR)/SourceGraph.cpp $(DIR)/TargetGraph.cpp $(DIR)/ReadInputFiles.cpp $(DIR)/GenericPartitioningAlgorithm.cpp $(DIR)/Comm.cpp $(DIR)/main.cpp $(DIR)/InferenceRate.cpp $(DIR)/CoarsenGraph.cpp
OBJ = $(SRC:.c=.o)
DIRWIF = write-input-file

MDN2PCIoT: $(DIR)/main.cpp $(OBJ)
	g++ -c $(FLAGS) $(DIR)/SourceGraph.cpp    			     -o $(DIROBJ)/SourceGraph.o
	g++ -c $(FLAGS) $(DIR)/TargetGraph.cpp    			     -o $(DIROBJ)/TargetGraph.o
	g++ -c $(FLAGS) $(DIR)/ReadInputFiles.cpp 			     -o $(DIROBJ)/ReadInputFiles.o
	g++ -c $(FLAGS) $(DIR)/CoarsenGraph.cpp 			     -o $(DIROBJ)/CoarsenGraph.o
	g++ -c $(FLAGS) $(DIR)/GenericPartitioningAlgorithm.cpp -o $(DIROBJ)/GenericPartitioningAlgorithm.o
	g++ -c $(FLAGS) $(DIR)/Comm.cpp	 			     -o $(DIROBJ)/Comm.o
	g++ -c $(FLAGS) $(DIR)/InferenceRate.cpp	 			     -o $(DIROBJ)/InferenceRate.o
	g++ -c $(FLAGS) $(DIR)/main.cpp						 -o $(DIROBJ)/main.o
	g++    $(FLAGS) $(OBJ)				             -o MDN2PCIoT

MDN2PCIoT-debug: 
	g++ -c -O0 -g3 -fopenmp $(DIR)/SourceGraph.cpp    			   -o $(DIROBJ)/SourceGraph.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/TargetGraph.cpp    			   -o $(DIROBJ)/TargetGraph.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/ReadInputFiles.cpp 			   -o $(DIROBJ)/ReadInputFiles.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/CoarsenGraph.cpp 			     -o $(DIROBJ)/CoarsenGraph.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/GenericPartitioningAlgorithm.cpp -o $(DIROBJ)/GenericPartitioningAlgorithm.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/Comm.cpp	 			   -o $(DIROBJ)/Comm.o
	g++ -c -O0 -g3 -fopenmp $(DIR)/InferenceRate.cpp	 			     -o $(DIROBJ)/InferenceRate.o
	g++    -O0 -g3 -fopenmp $(OBJ)	   	     			   -o MDN2PCIoT-debug

#prog.x: main.o math.o log.o
#	ld -dynamic-linker /lib64/ld-linux-x86-64.so.2 /usr/lib/x86_64-linux-O0 -gnu/crt1.o /usr/lib/x86_64-linux-O0 -gnu/crti.o -L/usr/lib64 main.o math.o log.o -lc /usr/lib/x86_64-linux-O0 -gnu/crtn.o -o prog.x
	
clean: 
	rm -f $(DIROBJ)/*.o MDN2PCIoT MDN2PCIoT-debug $(DIRWIF)/*.o writeInputFile*

writeInputFileCNN2D: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_CNN_2D.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_CNN_2D.c   -o $(DIRWIF)/write_input_file_CNN_2D.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_CNN_2D.o -o writeInputFileCNN2D

writeInputFileCNN2DLeNetDirected: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_CNN_2D_LeNet_directed.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_CNN_2D_LeNet_directed.c   -o $(DIRWIF)/write_input_file_CNN_2D_LeNet_directed.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_CNN_2D_LeNet_directed.o -o writeInputFileCNN2D_LeNet_directed

writeInputFileCNN2DAlexNet: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_CNN_2D_AlexNet.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_CNN_2D_AlexNet.c   -o $(DIRWIF)/write_input_file_CNN_2D_AlexNet.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_CNN_2D_AlexNet.o -o writeInputFileCNN2DAlexNet

writeInputFileCNN2DAlexNetDirected: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_CNN_2D_AlexNet_directed.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_CNN_2D_AlexNet_directed.c   -o $(DIRWIF)/write_input_file_CNN_2D_AlexNet_directed.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_CNN_2D_AlexNet_directed.o -o writeInputFileCNN2DAlexNetDirected

writeInputFile_LeNet_1x1: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet_1x1.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet_1x1.c   -o $(DIRWIF)/write_input_file_LeNet_1x1.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet_1x1.o -o writeInputFile_LeNet_1x1

writeInputFile_LeNet_2x2_directed: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet_2x2_directed.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet_2x2_directed.c   -o $(DIRWIF)/write_input_file_LeNet_2x2_directed.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet_2x2_directed.o -o writeInputFile_LeNet_2x2_directed

writeInputFile_LeNet_2x2: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet_2x2.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet_2x2.c   -o $(DIRWIF)/write_input_file_LeNet_2x2.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet_2x2.o -o writeInputFile_LeNet_2x2

writeInputFileSBAC: $(DIRWIF)/sourceGraph.c $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.c
	gcc -c $(DIRWIF)/sourceGraph.c 		  -o $(DIRWIF)/sourceGraph.o
	gcc -c $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.c   -o $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.o
	gcc    $(DIRWIF)/sourceGraph.o $(DIRWIF)/write_input_file_LeNet-SBAC-PAD-2018.o -o writeInputFile_LeNet-SBAC-PAD-2018
