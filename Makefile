# different vatiables
FC = gfortran
FFLAGS = -o

#SOURCES = matrixOperation.f90 sampling_MRT.f90 work.f 
SOURCESTEST = matrixOperation.f90 im_scat.f sampling_MRT.f90 sample.f90 tests.f90 work_test.f 
SOURCES = matrixOperation.f90 im_scat.f sampling_MRT.f90 sample.f90 work.f
EXECUTABLE = montecarlo
EXEXCUTABLETEST = testmonte

run: $(SOURCES)
	$(FC) $(FFLAGS) $(EXECUTABLE) $(SOURCES)
	./$(EXECUTABLE)

test: $(SOURCESTEST)
	$(FC) $(FFLAGS) testrun $(SOURCESTEST)
	./testrun

clean:
	rm -rf $(EXECUTABLE) $(EXECUTABLETEST) testrun fort.*
