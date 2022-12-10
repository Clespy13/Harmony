CC = gcc
CPPFLAGS =
CFLAGS = -D__NO_INLINE__ -fsanitize=address -g -O3 -Wall -Wextra -lm `pkg-config --cflags sdl2 SDL2_image`
LDFLAGS = -lm
LDLIBS = -lm `pkg-config --libs sdl2 SDL2_image`

SOURCE_DIR := .
SRC =  ./image_transformation/image/binarization.c ./image_transformation/image/image_processing.c ./image_transformation/image/image_rotation.c ./image_transformation/image/image_resize.c ./image_transformation/image/hough_transform.c ./detection/processing.c main.c
OBJ = ${SRC:.c=.o}
DEP = ${SRC:.c=.d}
EXE = ${SRC:.c=}

all: main

main : ${OBJ}
		$(CC) -o main $(CFLAGS) $^ $(LDLFLAGS) $(LDLIBS)
clean:
	${RM} ${DEP}
	${RM} ${OBJ}
	${RM} ${EXE}
	${RM} binarized_and_rotated.png
	${RM} main

# END
