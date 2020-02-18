srcdir = ./src
bindir = ./bin

ana_rw_source = $(srcdir)/ana_rw/ana_rw.c
ana_rw_bin = $(bindir)/ana_rw_quiet
ana_rw_cpp = -DQUIET
ana_rw_lib = -fopenmp -lm

ana_rw: $(ana_rw_source)
	gcc $(ana_rw_cpp) $< -o $(ana_rw_bin) $(ana_rw_lib)

all: ana_rw

install: 
	@echo "No need to make install!"
