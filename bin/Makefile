ifndef PYTHIA
$(error *** define PYTHIA; must point to location of pythia8)
endif

# name of Pythia
NAME=pythiaRun
CXX=g++

bin: $(NAME)

$(NAME): 
	$(CXX) $(CXXFLAGS) -Wno-shadow -I$(PYTHIA)/include \
	$@.cc -o $@ \
	-L$(PYTHIA)/lib -lpythia8 -lLHAPDF -lHepMC

clean:
	rm -rf $(NAME) $(NAME).o
