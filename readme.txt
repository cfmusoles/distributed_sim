*** DISTRIBUTED SIMULATOR ***

Compiling instructions:

In order to run the code, you must compile it issuing the command:

make

Which will generate the distSim executable in Linux. To compile it, you must have the library METIS installed and configured in your system (ParMETIS and SCALASCA are optional and should not be a problem if you do not have them. METIS is quite simple to install and configure, see:

Download: http://glaros.dtc.umn.edu/gkhome/metis/metis/download
Manual: http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

Note that the makefile expects METIS to be installed in the following folder: $(HOME)/metis_build. If you choose to have it elsewhere, please modify the Makefile to point towards the correct folder.

Running instructions:

To run the simulation, issue the following command from the project folder:

mpirun -np XX distSim test_name comm_pattern partitioning

Where XX is the number of processes you want to use, and the parameters after distSim can be specified as follows:

* test_name: name of the filename that contains simulation results (including sim time, sync time, etc.).
* comm_pattern: option to specify the MPI communication strategy to be used. The options are: allGather (to use MPI_AllGatherv to send messages between processes), p2p (to use the point-2-point with master node coordination) and subscriber (to use the point-2-point with preloaded sender lists). See the PDF attached for more info on each of the strategies.
* partitioning: specify whether to use graph partitioning to distribute workload (using graphPartitioning) or random allocation (using random).


GENERAL NOTES

This folder contains a prototype for a distributed neuronal simulator. 

INSTRUMENTATION AND TRACING WITH SCALASCA

Instrumentation and tracing is done via Scalasca (www.scalasca.org). See installation dependencies at the bottom. 

To enable instrumentation and tracing, compile prepending the command 'scorep', eg:

scorep [--user] mpicxx -o app source.cpp # use -user if code has manual instrumentation

To run a profiling experiment prepend mpirun with 'scalasca -analyze -t', eg:


scalasca -analyze -t -e experiment_name mpirun -np 2 app # can also be filtered with -f filter_name; -t option is for enabling tracing (more expensive)


For trace capacity requirements run 'scorep-score -r [filename].cubex'
	can use a scorep filter file to exclude include functions (see Score-P User manual)

Analysis report after tracing and experimenting, scalasca -examine [experiment_name dir]
	data per function, including MPI bytes sent!


INSTRUMENTATION AND TRACING WITH EXTRAE AND PARAVER

Source code instrumentation, include "${EXTRAE HOME}/include/extrae user events.h" 

Get a copy of the following from the EXTRAE_HOME/share/example/WHICHEVER_EXAMPLE folder:
	extrae.xml --> configuration file
	trace.sh --> necessary to run the trace (sets env variables and LD_PRELOADS)
		define EXTRAE_USE_POSIX_CLOCK inside

Place them in the apps folder

Run the trace with mpirun -np X ./trace.sh ./app [app params]

It should generate and merge trace files into app.prv file, ready for Paraver



SHARC instructions

To use scorep / scalasca, need to load both MPI and GCC (version < 6.2, same as used when compiling)
	module load dev/gcc/5.4
	module load mpi/openmpi/2.0.1/gcc-4.9.4 (actually, use same version, gcc-5.4)!



SCALASCA INSTALLATION DEPENDENCIES


download www.scalasca.org
install Qt --> sudo apt-get install qt-sdk
	wget https://download.qt.io/archive/qt/4.8/4.8.4/qt-everywhere-opensource-src-4.8.4.tar.gz
	tar xvf qt-everywhere-opensource-src-4.8.4.tar.gz
	cd qt-everywhere-opensource-src-4.8.4
	module load dev/gcc/5.4 (compile only with version < 6.2)
	./configure -embedded -prefix /home/cop12c/qt_build -nomake examples -nomake tests
	gmake
	gmake install
	add to .bashrc
		export PATH="$PATH:/home/cop12c/qt_build/bin"
install scorep 2.0
	if fortran error, reinstall openmpi
	wget http://www.vi-hps.org/upload/packages/scorep/scorep-3.0.tar.gz
	tar xvf scorep-3.0.tar.gz
	cd scorep-3.0
	mkdir build
	cd build
	../configure --prefix=/home/cop12c/scorep-3.0/build
	make
	make install
	add to .bashrc
		export PATH="$PATH:/home/cop12c/scorep-3.0/build/bin"
install OTF2 library v2.0
	wget http://www.vi-hps.org/upload/packages/otf2/otf2-2.0.tar.gz
	tar xvf otf2-2.0.tar.gz
	cd otf2-2.0
	./configure --prefix=/home/cop12c/otf2-2.0/build/
	make
	make install
install cube 4.3
	download
	./configure --prefix=/home/cop12c/cube-4.3.4/build/
	make
	make install
	add to .bashrc
		export PATH="$PATH:/home/cop12c/cube-4.3.4/build/bin"
install scalasca
	wget http://apps.fz-juelich.de/scalasca/releases/scalasca/2.3/dist/scalasca-2.3.1.tar.gz
	tar xvf scalasca-2.3.1.tar.gz
	cd scalasca-2.3.1
	./configure --prefix=/home/cop12c/scalasca2_build/ --with-otf2=/home/cop12c/otf2-2.0/build --with-cube=/home/cop12c/cube-4.3.4/build
	make
	make install
add to .bashrc
	export PATH="$PATH:/home/cop12c/scalasca2_build/bin"
	export PATH="$PATH:/home/cop12c/scorep-3.0/build/bin"	

To do manual source-code instrumentation, add to the includes: scorep/SCOREP_User.h
	$HOME/[scorep build folder]

install Zoltan
	compile instructions: link to include/zoltan.h and link with -lzoltan
	wget http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v3.83.tar.gz
	tar -xvf zoltan_distrib_v3.83.tar.gz
	cd Zoltan_v3.83
	mkdir build
	cd build
	../configure -C --prefix=$HOME/Zoltan_v3.83/build --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ulong
		for support with metis, ../configure -C --prefix=$HOME/Zoltan_v3.83/build --enable-mpi --with-mpi-compilers --with-gnumake --with-id-type=ulong --with-parmetis --with-parmetis-incdir=$HOME/parmetis__ompi_build/include --with-parmetis-libdir=$HOME/parmetis_ompi_build/lib CFLAGS="-I$HOME/metis_build/include/" LDFLAGS="-L$HOME/metis_build/lib"
	make everything
	make install

install PAPI http://icl.cs.utk.edu/papi/software/
	wget http://icl.utk.edu/projects/papi/downloads/papi-5.5.1.tar.gz
	tar -xvf papi-5.5.1.tar.gz
	cd papi-5.5.1
	cd src	
	./configure --prefix=$HOME/papi_build
	make
	make install
	add to .bashrc
		export PATH="$PATH:/home/cop12c/papi_build/bin"

install binutils (required by Extrae)
	wget http://ftp.gnu.org/gnu/binutils/binutils-2.25.tar.bz2
	tar -xvf binutils-2.25.tar.bz2
	cd binutils-2.25
	./configure --prefix=$HOME/binutils-2.25_build --enable-shared --enable-install-libiberty
	make
	make install

install libxml2
	wget ftp://xmlsoft.org/libxml2/libxml2-sources-2.7.8.tar.gz
	tar -xvf libxml2-sources-2.7.8.tar.gz
	cd libxml2-2.7.8
	./configure --prefix=$HOME/libxml2_build	
	make
	make install

install Extrae https://tools.bsc.es/downloads
	wget https://ftp.tools.bsc.es/extrae/extrae-3.4.3-src.tar.bz2
	tar -xvf extrae-3.4.3-src.tar.bz2
	cd extrae-3.4.3
	./configure --prefix=$HOME/extrae-3.4.3_build --with-mpi=$MPI_HOME --without-unwind --without-dyninst --with-papi=$HOME/papi_build --with-binutils=$HOME/binutils-2.25_build --with-xml-prefix=$HOME/libxml2_build --enable-posix-clock
	make
	make install
	add to .bashrc
		export PATH="$PATH:/home/cop12c/extrae-3.4.3_build/bin"

To install Ravel MPI visualization tool, need to have Muster, Qt, OTF2.1+ and OTF1
	https://github.com/LLNL/ravel
	install Muster (requires boost!)
		BOOST: sudo apt-get install cmake libblkid-dev e2fslibs-dev libboost-all-dev libaudit-dev
		wget https://github.com/LLNL/muster/archive/master.zip
		unzip master.zip -d .
		mkdir linux-x86_64 && cd linux-x86_64
		cmake -D CMAKE_INSTALL_PREFIX=~/muster-master/linux-x86_64 -D MUSTER_USE_PMPI=TRUE ..
		make
		make install
		add to bashrc
			export LD_LIBRARY_PATH="/home/$USER/muster-master/linux-x86_64/lib/:$LD_LIBRARY_PATH"	
	install OTF1
		download from http://wwwpub.zih.tu-dresden.de/%7Emlieber/dcount/dcount.php?package=otf&get=OTF-1.12.5salmon.tar.gz
		untar
		cd OTF-1.12.5salmon
		mkdir build
		./configure --prefix=$HOME/OTF-1.12.5salmon/build
		make
		make install
	install ravel
		git clone https://github.com/scalability-llnl/ravel.git
		mkdir ravel/build
		cd ravel/build
		cmake -DCMAKE_INSTALL_PREFIX=~/ravel/build -D Muster_INCLUDE_DIRS=~/muster-master/linux-x86_64/include -D Muster_LIBRARIES=~/muster-master/linux-x86_64/lib/libmuster.so -D OTF2_INCLUDE_DIRS=~/otf2-2.0/build/include -D OTF2_LIBRARIES=~/otf2-2.0/build/lib/libotf2.a -D OTF_INCLUDE_DIRS=~/OTF-1.12.5salmon/build/include/open-trace-format -D OTF_LIBRARIES=~/OTF-1.12.5salmon/build/lib/libopen-trace-format.a ..
		make
		make install


