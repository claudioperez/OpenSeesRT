all: #-DUSE_TCL_STUBS
	# cc -shared -fPIC $$CONDA_PREFIX/lib/libtcl8.6.so HTTP.cpp  peerNGA.cpp  Socket.cpp  StringContainer.cpp -lstdc++ -o peerNGA.so
	cc -shared -fPIC $$CONDA_PREFIX/lib/libtcl8.6.so peerNGA.cpp  socket.c  StringContainer.cpp -lstdc++ -o peerNGA.so
	python -m opensees -c 'import ./peerNGA.so; puts "[searchPeerNGA -fault San -magLo 6.5 -magHi 6.75]"'
