all:
	gcc main.cpp utils.c bwt.cpp kthread.c -lstdc++ -lz -lpthread -o panlive
