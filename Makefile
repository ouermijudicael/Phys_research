# different vatiables
FC = gfortran
FFLAGS = -o

#SOURCES = matrixOperation.f90 sampling_MRT.f90 work.f 
SOURCES = work.f 
EXECUTABLE = montecarlo


all: $(SOURCES)
	$(FC) $(FFLAGS) $(EXECUTABLE) $(SOURCES)
	./$(EXECUTABLE)

clean:
	rm -rf $(EXECUTABLE) fort.*
