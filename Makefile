build:
	gcc -o rtx.out main.c render.c aabb.c vector.c ray.c bvh.c clock.c color.c triangle.c deque.c -lpthread -lm

build_debug: 
	gcc -Wall -Wextra -fsanitize=address -g -o rtx.out main.c render.c aabb.c vector.c ray.c bvh.c clock.c color.c triangle.c deque.c -lpthread -lm

run: rtx.out
	./rtx.out ./meshes/powerplant/powerplant_2.obj

all: build run
debug: build_debug run