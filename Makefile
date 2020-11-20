#HPP = $(shell find . -name *.hpp) ./src/
#INC = $(dir $(HPP)) 
#INC_FLAG = $(foreach d, $(INC) , -I$d) 
CPP = $(shell find ./src -type f -name *.cpp)
BASHRCPATH = alias Sam10K-solver\=\'$(shell pwd)/run\'

compile: 
	g++ -o run -O3 -I./src/ $(CPP) 
	echo $(BASHRCPATH) >> ~/.bashrc
	source ~/.bashrc
	
clean:
	#rm -rf *.o run
	
	










