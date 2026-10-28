# gfortran is the only compiler tested; -std=legacy silences the warnings about PAUSE and tab-formatted lines
FC     = gfortran
FFLAGS = -O2 -std=legacy

surya: surya_2.f
	$(FC) $(FFLAGS) -o $@ $<

run: surya
	./run.sh

clean:
	rm -f surya

.PHONY: run clean
