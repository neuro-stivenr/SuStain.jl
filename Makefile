default: install

progdir=$(shell pwd)

install:
	julia --project=. -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
	@echo "#!/bin/bash" > SuStain.sh
	@echo "" >> SuStain.sh
	@echo 'NTHREADS=$$1; shift;' >> SuStain.sh
	@echo 'julia -t$$NTHREADS --project='"$(progdir) $(progdir)/src/SuStain.jl" '$$@' >> SuStain.sh
	chmod +x SuStain.sh
