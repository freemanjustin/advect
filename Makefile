###################################################################
#
# particle advection 
#
# freeman.justin@gmail.com
#
##################################################################


OBJ=	./src/main.o 

# Compliler flags
INC=    -I./

CFLAGS=	-O3 -g -Wall -DNETCDF

CC=	gcc $(CFLAGS) $(INC) 

# Libraries
LFLAGS=	-lnetcdf

# Executable

EXEC=	./bin/advect

$(EXEC):$(OBJ)
	$(CC) -o $(EXEC) $(OBJ) $(LFLAGS)

clean:
	rm $(OBJ)
	rm $(EXEC)
