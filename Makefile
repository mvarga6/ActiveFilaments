
build:
	nvcc -std=c++11 -x cu -Xptxas -O3,-v -Wno-deprecated-gpu-targets -lcurand deps/easylogging++.cc examples/simple.cpp -o examples/linux/simple

build-debug:
	nvcc -std=c++11 -x cu -lcurand -G -g -Wno-deprecated-gpu-targets deps/easylogging++.cc examples/simple.cpp  -o examples/linux/simple

run:
	./examples/linux/simple 1000 10 110 110 0 1000 10
