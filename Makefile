CFLAGS = -std=c++17 -O2

LDFLAGS = -lglfw -ldl -lpthread -lX11 -lXxf86vm -lXrandr -lXi -lGL -lGLEW -lGLU

ComputeShader: *.cpp *h
	g++ $(CFLAGS) -o Computer *.cpp $(LDFLAGS)
