
build:
	nvcc -x cu -Xptxas -O3,-v -Wno-deprecated-gpu-targets -lcurand deps/easylogging++.cc examples/simple.cpp -o examples/linux/simple

build-debug:
	nvcc -x cu -lcurand -G -g -Wno-deprecated-gpu-targets deps/easylogging++.cc examples/simple.cpp  -o examples/linux/simple