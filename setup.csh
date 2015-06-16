setenv PATH ${PWD}/bin:${PATH}
setenv DELPHES ${HOME}/external/delphes
setenv PYTHIA ${HOME}/external
setenv DYLD_LIBRARY_PATH ${PYTHIA}/lib
setenv LD_LIBRARY_PATH ${DELPHES}/lib:${PYTHIA}/lib:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH ${PWD}/lib:${LD_LIBRARY_PATH}
setenv PYTHONPATH ${PWD}/python:${PYTHONPATH}


