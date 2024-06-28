SHELL := /bin/bash
EXPECTED_ENV=stdpar-bh

D ?= 2
DIM_SIZE=$D
OUTFILE=./nbody_d${DIM_SIZE}
MAIN=./src/main.cpp

CONDA_FILES=-I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib

STDEXEC=./stdexec/include

VERSION=-std=c++20


# nvidia flags
FLAGS=${VERSION} -fast -DDIM_SIZE=${DIM_SIZE}
GPU_DEBUG=-gpu=lineinfo -gpu=keep
GPU_FLAGS=-stdpar=gpu

# other flags
CXX_FLAGS=${VERSION} -march=native -Ofast -ltbb -lpthread -DDIM_SIZE=${DIM_SIZE}
PEDANTIC=-Wall -Wextra -pedantic

# GET_TBB=source ~/intel/oneapi/setvars.sh


gcc: check-env
	g++ ${CXX_FLAGS} ${CONDA_FILES} ${MAIN} -o ${OUTFILE}_${@}

clang: check-env
	clang++ ${CXX_FLAGS} ${CONDA_FILES} ${MAIN} -o ${OUTFILE}_${@}

# icpx: intel
# 	${GET_TBB} && icpx ${CXX_FLAGS} ${CONDA_FILES} ${MAIN} -o ${OUTFILE}_${@}

nvcpp: check-env
	nvc++ ${FLAGS} -stdpar=multicore ${CONDA_FILES} ${MAIN} -o ${OUTFILE}_${@}

gpu: check-env
	nvc++ ${FLAGS} ${GPU_FLAGS} ${CONDA_FILES} ${MAIN} -o ${OUTFILE}_gpu

debug: check-env
	g++ ${VERSION} -DDIM_SIZE=${DIM_SIZE} -g -O0 -ltbb ${MAIN} -o ${OUTFILE}_${@}

cpu: gcc clang nvcpp

default: gcc

all: default

nothing:

check-env:
	@if [ "${CONDA_DEFAULT_ENV}" != "${EXPECTED_ENV}" ]; then \
	    echo "Error: Incorrect mamba environment. Expected ${EXPECTED_ENV}, got ${CONDA_DEFAULT_ENV}. Please run: mamba activate ${EXPECTED_ENV}"; \
	    exit 1; \
	fi

clean:
	rm -f pgcuda* .main* main.n* nbody* report*

