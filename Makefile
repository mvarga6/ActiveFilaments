
build:
	nvcc -x cu -lcurand deps/easylogging++.cc examples/simple.cu -o examples/linux/simple

build-debug:
	nvcc -x cu -lcurand -G -g -Wno-deprecated-gpu-targets deps/easylogging++.cc examples/simple.cpp  -o examples/linux/simple