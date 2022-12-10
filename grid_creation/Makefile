# Makefile

CC = gcc
CPPFLAGS =
CFLAGS = -D__NO_INLINE__ -Wall -Wextra -O3 `pkg-config --cflags sdl2 SDL2_image SDL2_ttf`
LDFLAGS = -lm -ldl
LDLIBS = -lm `pkg-config --libs sdl2 SDL2_image SDL2_ttf`

all: grid_creation

SRC = grid_creation
OBJ = ${SRC:.c=.o}
EXE = ${SRC:.c=}

image_processing : ${OBJ}
		$(CC) -o grid_creation $(CFLAGS) $^ $(LDLFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	${RM} ${OBJ}
	${RM} ${EXE}

# END
