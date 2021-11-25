# Grandi model make
CC       = clang
CFLAGS   = -Wall -g -O0
INCLUDE  = /usr/include
MY_APP	 = CVode_Pandit_Rat
LIB	 = /usr/lib

.PHONY: all, clean
all: ${MY_APP}

${MY_APP}: CVode_Pandit_Rat.cc
	${CC} ${CFLAGS} -c -I$(INCLUDE) -L$(LIB) CVode_Pandit_Rat.cc -o CVode_Pandit_Rat.o 
	${CC} ${CFLAGS} -c -I$(INCLUDE) -L$(LIB) APD.c -o APD.o
	${CC} ${CFLAGS} CVode_Pandit_Rat.o  APD.o -I$(INCLUDE) -L$(LIB) -lsundials_cvode -lsundials_nvecserial -o ${MY_APP} -lm

run: grandi_cell.c
	rm -f ${MY_APP}
	${CC} ${CFLAGS} -c CVode_Pandit_Rat.cc -o CVode_Pandit_Rat.o -I$(INCLUDE) -L$(LIB) 
	${CC} ${CFLAGS} -c -I$(INCLUDE) -L$(LIB) APD.c -o APD.o
	${CC} ${CFLAGS} CVode_Pandit_Rat.o APD.o -I$(INCLUDE) -L$(LIB) -lsundials_cvode -lsundials_nvecserial -o ${MY_APP}
	./${MY_APP} 1000

clean:
	rm -r *.dat ${MY_APP} *.o *~ APD_measure_out.txt
